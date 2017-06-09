#!/usr/bin/env python

import sys
import os
import logging
import numpy as np

from collections import defaultdict as dd
from uuid import uuid4
from bx.intervals.intersection import Intersecter, Interval # pip install bx-python

FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class InsGroup:
    def __init__(self, rec, fn, header):
        fn = os.path.basename(fn)
        self._header = header
        self._recs = [rec]
        self._fns = [fn]
        self.info = {}


    def add(self, rec, fn):
        self._recs.append(rec)
        self._fns.append(os.path.basename(fn))


    def make_info(self):
        self._fns = list(set(self._fns))

        self.info['Chrom']                      = self._recs[0]['Chrom']
        self.info['Peak_Start']                 = min([int(rec['Peak_Start']) for rec in self._recs])
        self.info['Peak_End']                   = max([int(rec['Peak_End']) for rec in self._recs])
        self.info['Sample_Count']               = sum([int(rec['Sample_Count']) for rec in self._recs])
        self.info['Samples']                    = ','.join(self._fns)
        self.info['Peak_Width']                 = self.info['Peak_End'] - self.info['Peak_Start']
        self.info['Ins_Strand']                 = self._recs[0]['Ins_Strand']
        self.info['Total_Reads']                = sum([int(rec['Total_Reads']) for rec in self._recs])
        self.info['Unique_Alignments']          = sum([int(rec['Unique_Alignments']) for rec in self._recs])
        self.info['Start_Spread']               = np.mean([float(rec['Start_Spread']) for rec in self._recs])
        self.info['End_Spread']                 = np.mean([float(rec['End_Spread']) for rec in self._recs])
        self.info['Mean_MapQ']                  = np.mean([float(rec['Mean_MapQ']) for rec in self._recs])
        self.info['Mean_Matchpct']              = np.mean([float(rec['Mean_Matchpct']) for rec in self._recs])
        self.info['Mappability']                = np.mean([float(rec['Mappability']) for rec in self._recs])
        self.info['Ref_Ins']                    = ','.join(list(set([rec['Ref_Ins'] for rec in self._recs])))
        self.info['Nonref_Ins']                 = ','.join(list(set([rec['Nonref_Ins'] for rec in self._recs])))
        self.info['Consensus_Score']            = np.mean([float(rec['Consensus_Score']) for rec in self._recs])
        self.info['Consensus_Homopolymer_Frac'] = np.mean([float(rec['Consensus_Homopolymer_Frac']) for rec in self._recs])
        self.info['Consensus_Seq']              = ','.join(list(set([rec['Consensus_Seq'] for rec in self._recs])))


    def __str__(self):
        return '\t'.join(map(str, [self.info[h] for h in self._header]))


    def __lt__(self, other):
        assert self.info, 'must run make_info() before sorting InsCall objects'
        if self.info['Chrom'] != other.info['Chrom']:
            return self.info['Chrom'] < other.info['Chrom']

        return self.info['Peak_Start'] < other.info['Peak_Start']


if len(sys.argv) > 2:
    forest = dd(Intersecter)

    header = []

    glob_ins = {}

    logger.info('comparing sites from %d files...' % (len(sys.argv)-1))

    for fn in sys.argv[1:]:

        with open(fn, 'r') as tsv:
            logger.info('processing file: %s' % fn)
            for i, line in enumerate(tsv):
                if i == 0:
                    if header:
                        assert len(header) == len(line.strip().split()), 'header mismatch on file: %s' % fn
                    else:
                        header = line.strip().split()

                    continue

                rec = {}

                cols = line.strip().split()

                if len(cols) != len(header):
                    logger.warn('WARNING, improperly formatted record on %s line %d : %s' % (fn, i, line.strip()))
                    continue

                for j, c in enumerate(cols):
                    rec[header[j]] = c

                if rec['Chrom'] not in forest:
                    glob_uuid = str(uuid4())
                    glob_ins[glob_uuid] = InsGroup(rec, fn, header)
                    forest[rec['Chrom']].add_interval(Interval(int(rec['Peak_Start']), int(rec['Peak_End']), value=glob_uuid))

                else:
                    hits = forest[rec['Chrom']].find(int(rec['Peak_Start']), int(rec['Peak_End']))

                    if len(hits) == 0:
                        glob_uuid = str(uuid4())
                        glob_ins[glob_uuid] = InsGroup(rec, fn, header)
                        forest[rec['Chrom']].add_interval(Interval(int(rec['Peak_Start']), int(rec['Peak_End']), value=glob_uuid))

                    else:
                        glob_uuid = hits[0].value
                        glob_ins[glob_uuid].add(rec, fn)

    logger.info('combining %d insertion globs...' % len(glob_ins))

    for uuid, ins in glob_ins.iteritems():
        ins.make_info()

    logger.info('outputting insertions to stdout...')

    print '\t'.join(header)

    for i in sorted(glob_ins.values()):
        print str(i)


else:
    sys.exit('usage: %s <tsv1> <tsv2> ...' % sys.argv[0])
