#!/usr/bin/python

import sys

file_MAP = open(sys.argv[1],'r')
file_FASTA = open(sys.argv[2],'r')
file_MAP_out = open('MAP.out','w')
file_FASTA_out = open('FASTA.out','w')


for line in file_MAP:
    if line[0] == "~":
        file_MAP_out.write("\n" + line)
    else:
        file_MAP_out.write(line.replace("\n", ''))



for line in file_FASTA:
    if line[0] == ">":
        file_FASTA_out.write("\n" + line)
    else:
        file_FASTA_out.write(line.strip())


