#!/usr/bin/env python

"""

# usage: python %prog
# python3

pandas       0.24.1
numpy        1.15.4

"""

import os,sys,re
from statistics import mean
import pandas as pd
import numpy as np

mean_linesine={}
mean_herv={}
mean_dnate={}

f_path='calculated_ext_branch_len_Alu.txt'
with open(f_path) as infile:
    for line in infile:
        ls=line.strip().split('\t')
        vals=[float(i) for i in ls[1:11]]
        mean_linesine[ls[0]]=mean(vals)

f_path='calculated_ext_branch_len_LINE1s.txt'
with open(f_path) as infile:
    for line in infile:
        ls=line.strip().split('\t')
        vals=[float(i) for i in ls[1:11]]
        mean_linesine[ls[0]]=mean(vals)

f_path='calculated_ext_branch_len_herv.txt'
with open(f_path) as infile:
    for line in infile:
        ls=line.strip().split('\t')
        vals=[float(i) for i in ls[1:11]]
        mean_herv[ls[0]]=mean(vals)

f_path='calculated_ext_branch_len_DNATE.txt'
with open(f_path) as infile:
    for line in infile:
        ls=line.strip().split('\t')
        vals=[float(i) for i in ls[1:11]]
        mean_dnate[ls[0]]=mean(vals)



kimura_path='hg38.fa.align.gz.divergence'

kimura_d={}
with open(kimura_path) as infile:
    next(infile)
    next(infile)
    next(infile)
    next(infile)
    next(infile)
    next(infile)
    for line in infile:
        ls=line.strip().split('\t')
        if len(ls) < 5:
            break
        if not (ls[4] == '----'):
            element=ls[1] +'#'+ ls[0]
            kimura_d[element]=float(ls[4])




kimura_linesine={}
kimura_herv={}
kimura_dnate={}

mean_l=[mean_linesine, mean_herv, mean_dnate]
kimura_l=[kimura_linesine, kimura_herv, kimura_dnate]

for m,k in zip(mean_l, kimura_l):
    for name in m:
        if name in kimura_d:
            k[name]=kimura_d[name]
        pass



# pearson's r
mean=[]
kimura=[]
for m,k in zip(mean_l, kimura_l):
    tmp=[ m[i] for i in m ]
    mean.extend(tmp)
    tmp=[ k[i] for i in k ]
    kimura.extend(tmp)
    pass

s1=pd.Series(mean)
s2=pd.Series(kimura)
print(s1.corr(s2, method='pearson'))



