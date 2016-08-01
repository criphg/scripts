#!/usr/bin/python

import sys

MAP = open('MAP.out', 'r')
bed = open("BED.out", 'w')


ctgenome = 0

#busca por linhas - ok
for line in MAP:

    #obtencao do cromossomo
    if line[0] == '~':
        pos = (line.index('|') - 1)
        chrm = line[1:pos]

    #foco na linhas das sequencias
    else:
        ct = 0
        lenline=len(line)

        #varredura na sequencia
        while ct < lenline:
            ma = line[ct]
            if (ma == '!'):
                st = ctgenome + ct
                end = -1
                while True:
                    if (ma != '!'):
                        bed.write(chrm + '\t' + str(st) + '\t' + str(st+end) + '\n')
                        break

                    ct += 1
                    end += 1
                    ma = line[ct]

            else:
                ct += 1

        ctgenome += lenline





