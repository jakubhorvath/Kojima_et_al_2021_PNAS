#!/usr/bin/env python

"""

# usage: python %prog
# python3

scipy      1.1.0

"""

import os,sys,re,shutil
from statistics import mean,stdev
import scipy.stats as st


f_path='pseudogene_perc_ident.txt'
idents=[]
with open(f_path) as infile:
    for line in infile:
        ls=line.split()
        idents.append(float(ls[1]))
        pass

mean=mean(idents)
sd=stdev(idents)

value=float('76.92')  # PVI

zscore= (value - mean) / sd
print(zscore)

pvalues=st.norm.sf(abs(zscore)) * 2  # two-sided
print(pvalues)
