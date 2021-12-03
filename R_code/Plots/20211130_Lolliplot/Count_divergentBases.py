#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 14:36:20 2021

@author: criph
"""
from Bio import AlignIO
import pandas as pd

align = AlignIO.read("LineageBvsB11529Consensu_Spike_transl.fa", "fasta")
print(len(align[1,:])) # sequence , possition

coord_v = 0
list_ct = list()
list_coord= list()
list_mutid=list()
for pos_v in range(0,len(align[0,:])):
    count_v = 0
    if(align[0,pos_v] != "-"):
        coord_v = coord_v + 1
    for seq_v in range(1,len(align[:,0])):
        print(align[seq_v,pos_v])
        if(align[seq_v,pos_v] != align[0,pos_v]):
            count_v = count_v + 1
            
    if(count_v != 0):
        list_coord.append(coord_v)
        list_ct.append(count_v)
        list_mutid.append(str(align[0,pos_v])+str(coord_v)+str(align[1,pos_v]) )


df = pd.DataFrame({"coord":list_coord, "category":"mut","value":list_ct,"MutID":list_mutid})
df.to_csv("DF_S_mutation_count.tab", sep="\t", index=False)


align = AlignIO.read("LineageBvsB11529Consensu_Nucleocapsid_transl.fa", "fasta")
print(len(align[1,:])) # sequence , possition

coord_v = 0
list_ct = list()
list_coord= list()
list_mutid=list()
for pos_v in range(0,len(align[0,:])):
    count_v = 0
    if(align[0,pos_v] != "-"):
        coord_v = coord_v + 1
    for seq_v in range(1,len(align[:,0])):
        print(align[seq_v,pos_v])
        if(align[seq_v,pos_v] != align[0,pos_v]):
            count_v = count_v + 1
    if(count_v != 0):
        list_coord.append(coord_v)
        list_ct.append(count_v)
        list_mutid.append(str(align[0,pos_v])+str(coord_v)+str(align[1,pos_v]) )


df = pd.DataFrame({"coord":list_coord, "category":"mut","value":list_ct,"MutID":list_mutid})
df.to_csv("DF_N_mutation_count.tab", sep="\t", index=False)