#!/usr/bin/env python

import sys
import gzip
import logging

from collections import defaultdict as dd

logger = logging.getLogger(__name__)
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)
logger.setLevel(logging.INFO)

def qualtrim(qstring):
    q = [ord(b)-41 for b in list(qstring)]
    for i in range(0,len(q)-4): # sliding window, 4bp
        if avg(q[i:i+4]) < 5:
            return qstring[:i]
    return qstring


def avg(L): return sum(L) / float(len(L))


if len(sys.argv) == 2:

    n = 0
    reads = 0
    wrote = 0

    assert sys.argv[1].endswith('.gz'), "files should be gzipped and end in .gz: %s" % sys.argv[1]

    outfn = '.'.join(sys.argv[1].split('.')[:-1]) + '.trimmed.fq.gz'

    out = gzip.open(outfn, 'wb')

    with gzip.open(sys.argv[1], 'rb') as fq:
        seqname = ''
        trimseq = ''
        tagseq  = ''

        for line in fq:
            if n == 0: # name
                seqname = line.strip()
                n += 1

            elif n == 1: # seq
                origseq = line.strip()
                trimseq = origseq[21:]
                n += 1
                
            elif n == 2: n += 1 # comment

            elif n == 3: #qual
                origqual = line.strip()

                trimqual = qualtrim(origqual[21:])
                trimseq  = trimseq[:len(trimqual)]

                if len(trimseq) >= 50:
                    out.write('\n'.join((seqname, trimseq, '+', trimqual)) + '\n')
                    wrote += 1

                n = 0
                reads += 1

                if reads > 0 and reads % 10000 == 0: logger.info('%d reads parsed, %d reads written to %s' % (reads, wrote, outfn))

    out.close()

    logger.info("done parsing %d reads, wrote %d reads." % (reads, wrote))

else:
    print "usage:", sys.argv[0], "<fq.gz>"
