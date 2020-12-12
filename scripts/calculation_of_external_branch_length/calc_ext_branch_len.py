#!/usr/bin/env python

"""

# usage: python %prog TE_name
# python3

pandas          0.24.1
pybedtools      0.8.0

"""

import os,sys,re,subprocess,random
import pandas as pd
from pybedtools import BedTool

TEclass=sys.argv[1]
bed_path='hg38.fa.align.gz.bed'
mao_path='infer_NJ_nucleotide_kimura_partial90.mao'
hu_fasta='hg38.fa'
linsi_dir='/usr/local/bin/linsi'
megacc_dir='/usr/local/bin/megacc'

num_thread=4
min_overlap=500
random_seq_sample=100
seq_num_for_calc=10
TEs_fa=TEclass.replace('/', '_').replace('#', '_') + '.fa'
megacc_out=TEclass.replace('/', '_').replace('#', '_') + '_megacc_out'
external_branch_length_out=TEclass.replace('/', '_').replace('#', '_') + '_ext_branch_len.txt'

allTE=[]
with open(bed_path, 'r') as inbed:
    for line in inbed:
        te='\t' + TEclass + '\t'
        if te in line:
            lsplit=line.split()
            if ((int(lsplit[2]) - int(lsplit[1])) > min_overlap) and not ('_' in lsplit[0]):
                allTE.append(line.strip())

if len(allTE) < seq_num_for_calc:
    with open('error.txt', 'w') as outfile:
        outfile.write('%s: len(allTE) < seq_num_for_calc\n' % TEclass)
    sys.exit()




# calc.
exist_fa=False
        
if len(allTE) > random_seq_sample:
    TE_for_analysis=random.sample(allTE, random_seq_sample)
else:
    TE_for_analysis=allTE
bed=BedTool('\n'.join(TE_for_analysis), from_string=True)
fa=bed.sequence(fi=hu_fasta, s=True)
input_fa = str(open(fa.seqfn).read())
with open(TEs_fa, 'w') as out:
    out.write(open(fa.seqfn).read())


cmd='%s --thread %d %s' %(linsi_dir, num_thread, TEs_fa)
aligned_fa=subprocess.check_output(cmd.split())
aligned_fa=aligned_fa.decode(sys.stdout.encoding).split('\n')
fa={}
k=''
for line in aligned_fa:
    if '>' in line:
        if k: fa[k]=''.join(seq)
        k=line.replace('>', '')
        seq=[]
    else:
        seq.append(line)
fa[k]=''.join(seq)
os.remove(TEs_fa)


for i in fa.keys():
    seqlen=len(fa[i])
    break
single={}
for l in range(seqlen):
    n=0
    for header in fa:
        if fa[header][l] is not '-':
            n=n+1
    single[l]=n
max_cov=max(single.values())
max_cov_key_l=[]
for k,v in single.items():
    if v == max_cov:
        max_cov_key_l.append(k)
max_cov_pos=random.choice(max_cov_key_l)
fa_exist={}
for k,s in fa.items():
    if s[max_cov_pos] is not '-':
        fa_exist[k]=s
if max_cov > seq_num_for_calc:
    while True:
        pop_res={}
        fa_exist_keys=list(fa_exist.keys())
        for test in fa_exist_keys:
            fa_test={}
            for k,s in fa_exist.items():
                if not k == test:
                    fa_test[k]=s
            fa_test_len=len(fa_test)
            n_all=0
            for l in range(seqlen):
                n=0
                for header in fa_test:
                    if fa_test[header][l] is not '-':
                        n=n+1
                if n == fa_test_len:
                    n_all=n_all+1
            pop_res[test]=n_all
        max_pop=max(pop_res.values())
        max_pop_key_l=[]
        for k,v in pop_res.items():
            if v == max_pop:
                max_pop_key_l.append(k)
        max_pop_element=random.choice(max_pop_key_l)
        print('fasta_alignment_seq_num= ' + str(len(fa_exist)))
        print('complete_deletion_nt_num= ' + str(max(pop_res.values())) +'\n')
        fa_exist.pop(max_pop_element)
        if (len(fa_exist) <= seq_num_for_calc):
            fa_selected=fa_exist
            exist_fa=True
            break    # proceed 
elif len(allTE) <= random_seq_sample:
    with open('error.txt', 'w') as outfile:
        outfile.write('%s: (len(allTE) <= random_seq_sample) and not (max_cov > seq_num_for_calc)\n' % TEclass)

                    
if exist_fa is True:
    with open(TEs_fa, 'w') as out:
        for i in fa_selected:
            out.write('>'+ i +'\n'+ fa_selected[i] +'\n')
    cmd='%s -a %s -d %s -o %s' % (megacc_dir, mao_path, TEs_fa, megacc_out)
    success=subprocess.call(cmd.split())
    if success == 0:
        out=[]
        with open('%s.nwk' % megacc_out, 'r') as innwk:
            for line in innwk:
                lsplit=line.split(")\':")
                for i in lsplit[1:]:
                    out.append(i.split(',')[0].split(')')[0])
        with open(external_branch_length_out, 'w') as outfile:
            outfile.write(TEclass +'\t'+ '\t'.join(out) +'\n')



