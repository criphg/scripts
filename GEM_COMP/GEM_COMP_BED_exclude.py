#!/usr/bin/python

import sys


MAP = open('MAP.out', 'r')
bed = open("BED_exclude.out", 'w')




#busca por linhas - ok
for line in MAP:

    ctgenome = 1

    #obtencao do cromossomo
    if line[0] == '~':
        pos = (line.index('|') - 1)
        chrm_r = line[1:pos]
        print(chrm_r)

    #foco na linhas das sequencias
    elif (len(line) > 1):
        ct = 0
        lenline=len(line)
        print(lenline)

        if lenline != 0:
            #varredura na sequencia
            while ct < lenline-2:
                ma = line[ct]

                if (ma == '!'):
                    ct += 1
                else:
                    st = ctgenome + ct
                    end = -1
                    while True:
                        if (ma == '!' or ct == (lenline-1)):
                            bed.write(chrm_r + '\t' + str(st) + '\t' + str(st+end) + '\n')
                            break

                        ct += 1
                        end += 1
                        ma = line[ct]




            ctgenome += lenline



MAP.close()
bed.close()
