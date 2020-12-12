#!/usr/bin/env python

"""

# batch2_remove_pseudo.py
# usage: python %prog infile.txt
# python2.7

# infile.txt should contain length of seq of 100 species:
'name' +'\t'+ 'hg38_length' +'\t'+ 'panTro4_length' +'\t'+ 'gorGor3_length' +'\t'+ ... 'petMar2_length'

"""

import os,sys,re
import numpy as np

f_path=sys.argv[1]

nonmammal_pseudo=[]
armadillo_pseudo=[]
non_pseudo=[]
with open(f_path, 'r') as infile:
    for line in infile:
        lsplit=line.strip().split()
        df=np.array(lsplit[1:], dtype='int64')
        df=df*100/df[0]
        n=0
        for i in range(62,100):
            if df[i] >= 50:
                n=n+1
        if n >= 14:
            nonmammal_pseudo.append(lsplit[0])
        elif df[61] >= 50:
            armadillo_pseudo.append(lsplit[0])
        else:
            non_pseudo.append(line)
with open('pseudo_non_mammal_name.txt', 'w') as out_non_mammal:
    out_non_mammal.write('\n'.join(nonmammal_pseudo))
with open('pseudo_armadillo_name.txt', 'w') as out_armadillo:
    out_armadillo.write('\n'.join(armadillo_pseudo))
with open('not_pseudo.txt', 'w') as out_not_pseudo:
    out_not_pseudo.write(''.join(non_pseudo))


