#!/usr/bin/env python

"""
# usage: python %prog in.fa
# python3.7
"""

import os,sys,re,shutil
from statistics import mean


def parse_fasta(path_to_file):
    tmp={}
    seq=''
    with open(path_to_file) as infile:
        for line in infile:
            if '>' in line and seq:
                tmp[header]=seq
                header=line.strip().replace('>', '')
                seq=''
            elif '>' in line and not seq:
                header=line.strip().replace('>', '')
            else:
                seq += line.strip()
        tmp[header]=seq
    return tmp


def count_CG_CpG_perc(seq):
    seq=seq.upper()
    l=len(seq)
    c=seq.count('C')
    g=seq.count('G')
    cpg=seq.count('CG')
    gc_perc= 100 * ((c + g) / l)
    cpg_perc= 100 * (cpg / (l - 1))
    return [gc_perc, cpg_perc]



f_path=sys.argv[1]
fa=parse_fasta(f_path)

res_GC=[]
res_CpG=[]
for id in fa:
    gc,cpg=count_CG_CpG_perc(fa[id])
    res_GC.append(gc)
    res_CpG.append(cpg)
    pass

perc_GC=mean(res_GC)
perc_CpG=mean(res_CpG)

with open('out.txt', 'w') as outfile:
    outfile.write('Percent GC = ' + str(perc_GC) + '\n' + 'Percent CpG = ' + str(perc_CpG) + '\n')
