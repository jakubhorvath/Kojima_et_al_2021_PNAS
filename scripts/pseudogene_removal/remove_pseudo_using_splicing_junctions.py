#!/usr/bin/env python

"""

# usage: python %prog
# python3

"""

import os,sys,re,shutil

# list up splicing junctions
f_path='splicing_junction_pos.txt'
tmp={}
with open(f_path) as infile:
    for line in infile:
        lsplit=line.split()
        if not lsplit[0] in tmp:
            tmp[lsplit[0]]=[int(lsplit[1])]
        else:
            tmp[lsplit[0]].append(int(lsplit[1]))
junc={}
with open(f_path) as infile:
    for line in infile:
        lsplit=line.split()
        if not lsplit[0] in junc:
            junc[lsplit[0]]=[lsplit[2], tmp[lsplit[0]]]
        pass

# judge pseudo or not
f_path='blastn_result_outfmt6.txt'
pseudo=[]
pseudo_set=set()
with open(f_path) as infile:
    for line in infile:
        lsplit=line.split()
        
        if int(float(lsplit[11])) >= 100:  # blast score >= 100
            id=lsplit[1]
            if int(lsplit[8]) < int(lsplit[9]):
                start= int(lsplit[8])
                end  = int(lsplit[9])
            else:
                start= int(lsplit[9])
                end  = int(lsplit[8])
        
            p=False
            if junc[id][0] == '+':
                for i in junc[id][1]:
                    if start < i < end:
                        pseudo_set.add(lsplit[0])
                        pseudo.append(lsplit[0] +'\tputative_origin: '+ lsplit[1])
                        p=True
                        break
            else:
                transcript_len=junc[id][1][-1]
                for i in junc[id][1]:
                    i_c= transcript_len - i
                    if start < i_c < end:
                        pseudo_set.add(lsplit[0])
                        pseudo.append(lsplit[0] +'\tputative_origin: '+ lsplit[1])
                        p=True
                        break
        pass
pseudo_list=list(pseudo_set)

all=set()
f_path='query_for_blastn_from_human_transcripts.fa'
with open(f_path) as infile:
    for line in infile:
        if '>' in line:
            pos=line.strip().replace('>', '')
            all.add(pos)
        pass
all=list(all)

out=[]
for i in all:
    if not i in pseudo_list:
        out.append(i)

with open('non_pseudo.txt', 'w') as outfile:
    outfile.write('\n'.join(out) +'\n')
with open('pseudo.txt', 'w') as outfile:
    outfile.write('\n'.join(pseudo) +'\n')

