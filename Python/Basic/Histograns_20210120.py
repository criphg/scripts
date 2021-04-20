#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 08:01:03 2021

Word Histograms from a text file.

@author: criph
"""

#Functions
def histogram_table(f_name):
    fn = open(f_name, "r")

    dic = {}
    for line in fn:
        word_list = line.split()
        for word in word_list:
            dic[word] = dic.get(word, 0) + 1

    #Most prevalent key
    b_wk = None
    b_count = None
    for wk,count in dic.items():
        if b_count is None or b_count < count :
            b_count = count
            b_wk = wk

    #Sort by number of counts
    tmp = []
    for (k,v) in dic.items():
        tump_tmp=(v,k)
        tmp.append(tump_tmp)
        
        tmp = sorted(tmp, reverse=True)
        
        #new sorted dict
        dic_sorted_by_value = {}
        for v2, k2 in tmp:
            dic_sorted_by_value[k2] = v2

    #print results - ##It can be moved to the main and adapted to dict_out
    print("The most prevalent word is \""+b_wk+"\" with "+str(b_count)+" counts.")        
    
    print("\nAll words sorted:")
    print(dic_sorted_by_value)
    
    print("\nTop 10 words:")
    for val, key in tmp[:10]:
        print(key, val)
    ##    
    
    return(dic_sorted_by_value)
    

##MAIN
f_name_in = input("Enter name file =")

dict_out = histogram_table(f_name_in)


#f_name_in.close()


#print("\n\nOut:\n")
#print(dict_out)