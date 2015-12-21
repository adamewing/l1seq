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


def levenshtein(a,b):
    ''' from http://hetland.org/coding/python/levenshtein.py '''
    n, m = len(a), len(b)
    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        a,b = b,a
        n,m = m,n
        
    current = range(n+1)
    for i in range(1,m+1):
        previous, current = current, [i]+[0]*n
        for j in range(1,n+1):
            add, delete = previous[j]+1, current[j-1]+1
            change = previous[j-1]
            if a[j-1] != b[i-1]:
                change = change + 1
            current[j] = min(add, delete, change)
            
    return current[n]


if len(sys.argv) >= 3:

    tags = {} # write filehandles
    tagstats = dd(int) # tracker

    sys.argv.append('uncat')

    for tag in sys.argv[2:]:
        tags[tag] = gzip.open('tag.%s.fq.gz' % tag, 'wb')

    with gzip.open(sys.argv[1], 'rb') as fq:
        n = 0
        reads = 0
        tag = None

        seqname = ''
        trimseq = ''
        tagseq  = ''

        for line in fq:
            if n == 0: # name
                seqname = line.strip()
                n += 1

            elif n == 1: # seq
                origseq = line.strip()
                trimseq = origseq[16:]
                tagseq  = origseq[5:11]

                tagstats[tagseq] += 1

                n += 1
                
            elif n == 2: n += 1

            elif n == 3: #qual
                best_tag = 'uncat' 
                origqual = line.strip()

                if tagseq in tags:
                    best_tag = tagseq

                elif len(tags) == 1:
                    best_tag = tags.keys()[0]

                else:
                    t_dist = [levenshtein(tag, tagseq) for tag in tags]
                    if min(t_dist) == 1 and sorted(t_dist)[1] > 2:
                        best_tag = tags.keys()[t_dist.index(min(t_dist))]

                if best_tag is not None:
                    trimqual = qualtrim(origqual[16:])
                    trimseq  = trimseq[:len(trimqual)]
                    if len(trimseq) >= 50:
                        tags[best_tag].write('\n'.join((seqname, trimseq, '+', trimqual)) + '\n')

                n = 0
                reads += 1

                if reads > 0 and reads % 10000 == 0: logger.info('parsed %d reads.' % reads)


    for tag in tags: tags[tag].close()

    logger.info("done.")
    for tagseq, tagcount in tagstats.iteritems():
        print '%s\t%d' % (tagseq, tagcount)

else:
    print "usage:", sys.argv[0], "<fq.gz> [tag list]"
