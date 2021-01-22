#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 23:23:11 2021

@author: criph
"""

import urllib
from bs4 import BeautifulSoup

fhand = urllib.request.urlopen("http://www.ign.com")

html = fhand.read()
soup = BeautifulSoup(html, 'html.parser')

tags = soup('a')
for tag in tags:
    print(tag.get('href',None))

counts = dict()
for line in fhand:
    words = line.decode().split()
    for word in words:
        counts[word] = counts.get(word, 0) + 1
#print(counts)

