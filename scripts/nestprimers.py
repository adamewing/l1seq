#!/usr/bin/env python

import sys
import os

from primer3 import bindings
from pysam import Fastafile
from string import maketrans


def makeprimers(name, template, minsize, maxsize):
    try:
        binding_res = bindings.designPrimers(
            {
                'SEQUENCE_ID': name,
                'SEQUENCE_TEMPLATE': template,
            },
            {
                'PRIMER_OPT_SIZE': 24,
                'PRIMER_MIN_SIZE': 21,
                'PRIMER_MAX_SIZE': 27,
                'PRIMER_OPT_TM': 58.0,
                'PRIMER_MIN_TM': 56.0,
                'PRIMER_MAX_TM': 63.0,
                'PRIMER_GC_CLAMP': 1,
                'PRIMER_PRODUCT_SIZE_RANGE': [[minsize,maxsize]],
            }
        )

        leftprimers  = {}
        rightprimers = {}

        for k, v in binding_res.items():
            if k.endswith('SEQUENCE'):
                if k.startswith('PRIMER_LEFT'):
                    leftprimers[int(k.split('_')[2])] = v
                if k.startswith('PRIMER_RIGHT'):
                    rightprimers[int(k.split('_')[2])] = v

        return leftprimers, rightprimers
    except:
        sys.stderr.write("warning: primer design failed for " + name + "\n")
        return {}, {}

def rc(dna):
    ''' reverse complement '''
    complements = maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    return dna.translate(complements)[::-1]


if len(sys.argv) == 3:
    ref = Fastafile(sys.argv[1])
    with open(sys.argv[2], 'r') as l1seq:
        for line in l1seq:
            if not line.startswith('Chr'): # header
                c = line.strip().split('\t')
                chrom  = c[0]
                strand = c[6]
                pos = min(int(c[1]), int(c[2]))

                if strand == '+':
                    pos = max(int(c[1]), int(c[2]))

                iname   = 'c' + chrom + 'p' + str(pos) + 'IN'
                oname   = 'c' + chrom + 'p' + str(pos) + 'OUT'

                outerstart = pos - 600
                outerend   = pos + 600
                innerstart = pos - 450
                innerend   = pos + 450

                if outerstart < 0: outerstart = 0
                if innerstart < 0: innerstart = 0

                outerseq = ref.fetch(chrom, outerstart, outerend)
                innerseq = ref.fetch(chrom, innerstart, innerend)

                if strand == 'plus':
                    outerseq = rc(outerseq)
                    innerseq = rc(innerseq)

                outerprimers = makeprimers(oname, outerseq, 1000, 1200)
                innerprimers = makeprimers(iname, innerseq, 700, 900)

                leftouterprimer  = ''
                rightouterprimer = ''
                leftinnerprimer  = ''
                rightinnerprimer = ''

                primertext = ''
                primertext += oname + 'FS' + '\t'
                if len(outerprimers[0]) > 0:
                    leftouterprimer = outerprimers[0][0]
                    primertext += leftouterprimer + '\t'
                else:
                    primertext += 'None\t'

                primertext += oname + 'ES' + '\t'

                if len(outerprimers[1]) > 0:
                    rightouterprimer = outerprimers[1][0]
                    primertext += rightouterprimer + '\t'

                else:
                    primertext += 'None\t'

                primertext += iname + 'FS' + '\t'

                if len(innerprimers[0]) > 0:
                    leftinnerprimer = innerprimers[0][0]
                    key = 0
                    while leftinnerprimer == leftouterprimer:
                        key += 1
                        if key in innerprimers[1]:
                            leftinnerprimer = innerprimers[0][key]
                        else:
                            leftinnerprimer = 'None'

                    primertext += leftinnerprimer + '\t'

                else:
                    primertext += 'None\t'

                primertext += iname + 'ES' + '\t'
                if len(innerprimers[1]) > 0:
                    rightinnerprimer = innerprimers[1][0]
                    key = 0
                    while rightinnerprimer == rightouterprimer:
                        key += 1
                        if key in innerprimers[1]:
                            rightinnerprimer = innerprimers[1][key]
                        else:
                            rightinnerprimer = 'None'

                    primertext += rightinnerprimer + '\t'

                else:
                    primertext += 'None\t'

                print line.strip() + '\t' + primertext.strip()

            else:
                sys.stdout.write(line.strip() + "\touter_left_primer_name\touter_left_primer_seq\t" 
                                              + "outer_right_primer_name\touter_right_primer_seq\t"
                                              + "inner_left_primer_name\tinner_left_prime_seq\t"
                                              + "inner_right_primer_name\tinner_right_prime_seq\n")

else:
    print "usage:", sys.argv[0], "<reference> <l1seq.py output>"



