#!/usr/bin/env python

"""

# usage: python %prog
# python3

"""

import os,sys,re,shutil

f_path='RetroGenes_v9_ucscRetroInfo9'  # downloaded from UCSC

out=[]
with open(f_path) as infile:
    next(infile)
    for line in infile:
        ls=line.strip().split('\t')
        if ('NM_' in ls[37]) or ('NR_' in ls[37]):
            if not '_' in ls[0]:
                out.append(ls[0] +'\t'+ ls[1] +'\t'+ ls[2] +'\t'+ ls[0] +':'+ ls[1] +'-'+ ls[2] +':'+ ls[37])
        pass

print(len(out))
with open('Retrogenes_v9_pos.bed', 'w') as outfile:
    outfile.write('\n'.join(out) +'\t')
