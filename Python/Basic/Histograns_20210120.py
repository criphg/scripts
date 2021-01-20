#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 08:01:03 2021

Criação de Fluxograma

@author: criph
"""

f_name = input("Enter name file =")
fn = open(f_name, "r")

dic = {}
for line in fn:
    word_list = line.split()
    for word in word_list:
        dic[word] = dic.get(word, 0) + 1

b_wk = None
b_count = None
for wk,count in dic.items():
    if b_count is None or b_count < count :
        b_count = count
        b_wk = wk