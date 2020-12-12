#!/usr/bin/env python

"""

# extract_PES.py
# usage: python %prog infile.txt
# python2.7

# infile.txt should contain scores of seq of 14 species:
'name' +'\t'+ 'chr_location' +'\t'+ 'panTro4_scores' +'\t'+ 'nomLeu3_scores' +'\t'+ ... 'loxAfr3_scores'

"""

import os,sys,re

f_path=sys.argv[1]

EVE_candidate=[]
with open(f_path, 'r') as infile:
    for line in infile:
        if re.search(r'[6-9].,.,[6-9].', line) is not None:
            ls=line.strip().split()
            n=0
            for i in range(2,15):
                if re.search(r'[0-9]+,1[0-9],[0-9]+|[0-9]+,.,[0-9]+', ls[i]) is not None:
                    n=n+1
                elif n >= 1:
                    n=-1
                    break
            if n >= 1:
                EVE_candidate.append(line)
print(''.join(EVE_candidate)).strip()

