#!/usr/bin/env python

"""

# usage: python %prog
# python3

"""

import os,sys,re,shutil

f_path='gencode.v27.annotation.gff3'
out=[]
total=0
pre_id=''
with open(f_path) as infile:
    for line in infile:
        if 'ID=exon' in line:
            lsplit=re.split('\t|:', line)
            id=lsplit[9]
            if id == pre_id:
                total += int(lsplit[4]) - int(lsplit[3]) + 1
                out.append(id +'\t'+ str(total) +'\t'+ lsplit[6])
            else:
                pre_id = id
                total = int(lsplit[4]) - int(lsplit[3]) + 1
                out.append(id +'\t'+ str(total) +'\t'+ lsplit[6])
        pass

with open('splicing_junction_pos.txt', 'w') as outfile:
    outfile.write('\n'.join(out))
