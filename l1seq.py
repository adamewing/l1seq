#!/usr/bin/env python

import os
import re
import pysam
import argparse

import numpy as np
import align

from collections import Counter

import logging
logger = logging.getLogger(__name__)
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.INFO)


class Cluster:
    def __init__(self, bam, reads, cstrand):
        self.reads   = sorted(reads, key=lambda x: x.reference_start-x.query_alignment_start)

        self.samples = self.getRG()
        self.scount  = len(self.samples.split(','))

        self.cstrand = cstrand
        self.chrom   = bam.getrname(reads[0].reference_id)

        self.astarts = [read.reference_start for read in reads] 
        self.astops  = [read.reference_end for read in reads]

        self.minpos  = min(self.astarts)
        self.maxpos  = max(self.astops)

        # number of unique alignments
        self.uniqaln = len(list(set(zip(self.astarts, self.astops))))

        self.start_spread = max(self.astarts) - min(self.astarts)
        self.stop_spread  = max(self.astops)  - min(self.astops)

        self.mean_mq = np.mean([read.mapq for read in self.reads])

        self.mean_match = np.mean([read_matchpct(read) for read in self.reads])


    def ins_strand(self):
        if self.cstrand == '+': return '-'
        return '+'


    def getRG(self):
        ''' return read groups from RG aux tag '''
        RGs = []
        for read in self.reads:
            for tag, val in read.tags:
                if tag == 'RG':
                    RGs.append(val)

        if len(RGs) > 0:
            return ','.join(['%s|%d' % (rg,count) for rg, count in Counter(RGs).iteritems()])
        else:
            return 'NA'


    def sorted_mapped_seqs(self):
        ''' return reads in genomic order '''
        return [read.query_sequence for read in self.reads]


    def mapscore(self, maptabix):
        ''' return average mappability across chrom:start-end region; maptabix = pysam.Tabixfile'''
        scores = []

        if self.chrom in maptabix.contigs:
            for rec in maptabix.fetch(self.chrom, self.minpos, self.maxpos):
                mchrom, mstart, mend, mscore = rec.strip().split()
                mstart, mend = int(mstart), int(mend)
                mscore = float(mscore)

                while mstart < mend and mstart:
                    mstart += 1
                    if mstart >= self.minpos and mstart <= self.maxpos:
                        scores.append(mscore)

            if len(scores) > 0:
                return sum(scores) / float(len(scores))
            else:
                return 0.0
        else:
            return 0.0


    def refelt(self, reftabix, window=500):
        if self.chrom in reftabix.contigs:
            for rec in reftabix.fetch(self.chrom, self.minpos-window, self.maxpos+window):
                rchrom, rstart, rend, rfam, rmillidiv, rstr = rec.strip().split()

                if rstr == self.ins_strand() == '+' and int(rstart) < self.minpos:
                    return '|'.join(rec.strip().split())

                if rstr == self.ins_strand() == '-' and int(rend) > self.maxpos:
                    return '|'.join(rec.strip().split())
        
        return 'NA'


    def nonref(self, nrtabix, window=500):
        if self.chrom in nrtabix.contigs:
            for rec in nrtabix.fetch(self.chrom, self.minpos-window, self.maxpos+window):
                return '|'.join(rec.strip().split())

        return 'NA'

    def __str__(self):
        fields = (self.chrom, self.minpos, self.maxpos, str(self.scount), self.samples, self.maxpos-self.minpos, 
                  self.ins_strand(), len(self.reads), self.uniqaln, self.start_spread, 
                  self.stop_spread, self.mean_mq, self.mean_match) 

        fields = map(str, fields)

        return '\t'.join(fields)



def consensus(seqs, minscore=0.95):
    ''' build consensus from sorted aligned reads iteratively, expects seqs to be sorted in ref genome order '''

    S = -np.ones((256, 256)) + 2 * np.identity(256)
    S = S.astype(np.int16)

    if len(seqs) == 1: # no consensus necessary
        return seqs[0], 1.0

    uniq_seqs = [seqs[0]]
    for i, seq in enumerate(seqs[1:], start=1):
        if seq != seqs[i-1]:
            uniq_seqs.append(seq)

    if len(uniq_seqs) == 1: # all seqs were the same!
        return uniq_seqs[0], 1.0

    cons = uniq_seqs[0]
    scores = []

    if len(uniq_seqs) > 1000: uniq_seqs = np.random.choice(uniq_seqs, size=1000)

    for seq in uniq_seqs[1:]:

        s1 = align.string_to_alignment(cons)
        s2 = align.string_to_alignment(seq)

        (s, a1, a2) = align.align(s1, s2, -2, -2, S, local=True)
        a1 = align.alignment_to_string(a1)
        a2 = ''.join([b for b in list(align.alignment_to_string(a2)) if b != '-'])

        score = float(len(a1) - (len(a1)-s)) / float(len(a1))
        scores.append(score)

        if re.search(a1, cons):
            cons_start, cons_end = locate_subseq(cons, a1)

            if score >= minscore and cons_end > len(cons)-5:
                align_end = locate_subseq(seq, a2)[1]
                cons += seq[align_end:]

    return cons, np.mean(scores)


def locate_subseq(longseq, shortseq):
    ''' return (start, end) of shortseq in longseq '''
    assert len(longseq) >= len(shortseq), 'orient_subseq: %s < %s' % (longseq, shortseq)
 
    match = re.search(shortseq, longseq)
    if match is not None:
        return sorted((match.start(0), match.end(0)))
 
    return None


def read_matchpct(read):
    ''' return number of mismatches / aligned length of read '''
    nm = [value for (tag, value) in read.tags if tag == 'NM'][0]
    return 1.0 - (float(nm)/float(read.alen))


def pass_read(read, minq=1, max_distal_clip=5):
    if read.is_secondary: return False
    if read.is_unmapped:  return False
    if read.is_duplicate: return False
    if read.mapq < minq:  return False

    if read.is_reverse:
        if read.query_alignment_end - read.alen > max_distal_clip: return False
    if not read.is_reverse:
        if read.query_alignment_start > max_distal_clip: return False

    if read_matchpct(read) < 0.98: return False
    
    return True


def build_clusters(bam, window=200, minq=1):
    ''' single pass to group reads into clusters '''

    clusters = []

    pile_plus  = []
    pile_minus = []

    for read in bam.fetch():
        if pass_read(read, minq=minq):
            if read.is_reverse:
                if len(pile_minus) == 0:
                    pile_minus.append(read)

                elif bam.getrname(read.reference_id) != bam.getrname(pile_minus[-1].reference_id) or read.reference_start - window > pile_minus[-1].reference_end:
                    clusters.append(Cluster(bam, pile_minus, '-'))
                    pile_minus = []
                    pile_minus.append(read)

                else:
                    pile_minus.append(read)

            else:
                if len(pile_plus) == 0:
                    pile_plus.append(read)

                elif bam.getrname(read.reference_id) != bam.getrname(pile_plus[-1].reference_id) or read.reference_start - window > pile_plus[-1].reference_end:
                    clusters.append(Cluster(bam, pile_plus, '+'))
                    pile_plus = []
                    pile_plus.append(read)

                else:
                    pile_plus.append(read)

    if len(pile_plus) > 0:  clusters.append(Cluster(bam, pile_plus, '+'))
    if len(pile_minus) > 0: clusters.append(Cluster(bam, pile_minus, '-'))

    return clusters


def homopol_filter(cons):
    ''' return True if consensus is dominated by one base '''
    return Counter(list(cons)).most_common(1)[0][1]/float(len(cons))


def main(args):
    assert os.path.exists(args.bam + '.bai'), 'please index %s with samtools index' %  args.bam

    bam = pysam.AlignmentFile(args.bam, 'rb')

    maptbx = pysam.Tabixfile(args.map)
    reftbx = pysam.Tabixfile(args.ref)
    nrtbx  = pysam.Tabixfile(args.nonref)

    clusters = build_clusters(bam)

    fields = ['Chrom', 'Peak_Start', 'Peak_End', 'Sample_Count', 'Samples', 'Peak_Width', 'Ins_Strand', 'Total_Reads', 'Unique_Alignments', 'Start_Spread', 'End_Spread', 'Mean_MapQ', 'Mean_Matchpct']
    fields += ['Mappability', 'Ref_Ins', 'Nonref_Ins', 'Consensus_Score', 'Consensus_Homopolymer_Frac', 'Consensus_Seq']

    print '\t'.join(fields)

    for cluster in clusters:
        seqs = cluster.sorted_mapped_seqs()
        logger.info('Building consensus: %s' % str(cluster))
        cons, cons_score = consensus(seqs)
        print '\t'.join((str(cluster), str(cluster.mapscore(maptbx)), cluster.refelt(reftbx), cluster.nonref(nrtbx), str(cons_score), str(homopol_filter(cons)), cons))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyse L1-seq data (e.g. Ewing and Kazazian Genome Res. 2010)')
    parser.add_argument('-b', '--bam', required=True, help='Input BAM (for multiple samples merge with distinct readgroups)')
    parser.add_argument('-m', '--map', required=True, help='Mappability Tabix')
    parser.add_argument('--nonref', required=True, help='Nonref insertion tabix')
    parser.add_argument('--ref', required=True, help='Reference insertion tabix')

    args = parser.parse_args()
    main(args)
