#!/usr/bin/env python

"""

# usage: python %prog
# python3

"""

import os,sys,re,shutil

f_path='Retrogenes_v9_pos_maf_cov.txt'

cand=[]
with open(f_path) as infile:
    for line in infile:
        ls=line.split()
        if (int(ls[10]) >= 50) and (int(ls[11]) >= 50) and (int(ls[12]) <= 5) and (int(ls[19]) <= 5) and (int(ls[38]) <= 5) and (int(ls[52]) <= 5): # calJac3 and saiBol1, otoGar3 and mm10 and felCat8 and loxAfr3
            cand.append(ls[0])
        pass

f_path='Retrogenes_v9_pos_maf_cov_tarsier.txt'

cand_tar=[]
with open(f_path) as infile:
    for line in infile:
        ls=line.split()
        if int(ls[2]) <= 5: # tarSyr2
            cand_tar.append(ls[0])
        pass

out=[]
for i in cand:
    if i in cand_tar:
        out.append(i)

print(len(out))
with open('simian_specific_pseudogenes.txt', 'w') as outfile:
    outfile.write('\n'.join(out) +'\n')
