#!/usr/bin/python

import sys

FAST = open("FASTA.out", 'r')
MAP = open('MAP.out', 'r')
seq = open('MAPSEQ.out', 'w')

for line in MAP:
    faline = FAST.readline()
    if line[0] == '~':
        seq.write('\n' + faline)
    else:
        ct = 0;
        lenline=len(faline)
        while ct < lenline:
            ma = line[ct]
            fa = faline[ct]
            if (ma == '!'):
                seq.write(fa)
            ct += 1;
FAST.close()
MAP.close()
seq.close()