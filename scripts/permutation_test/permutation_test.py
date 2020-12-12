#!/usr/bin/env python

"""
# usage: python %prog
# python3.7
"""

# permutation test; test enrichment
# linux

import os,sys,re,math,subprocess
from statistics import mean

original_pos='in.bed'  # original position; bed file
target_poss=['target_1.bed', 'target_2.bed']  # target position; bed file
chrom_size='hg38.chrom.sizes'

len_extend=0  # extend +/- from original element
n_permut=1000

len_extend=str(len_extend)


def count_intersect_num(original_pos, target_pos, chrom_size, len_extend):
    cmd='bedtools slop -i %s -g %s -b %s | bedtools intersect -a "stdin" -b %s -wa | sort | uniq | wc -l' % (original_pos, chrom_size, len_extend, target_pos)
    count=subprocess.check_output(cmd, shell=True).decode().strip()
    return int(count)

def count_randomized_intersect_num(original_pos, target_pos, chrom_size, len_extend):
    cmd='bedtools slop -i %s -g %s -b %s | bedtools shuffle -noOverlapping -i "stdin" -g %s | bedtools intersect -a "stdin" -b %s -wa | sort | uniq | wc -l'  % (original_pos, chrom_size, len_extend, chrom_size, target_pos)
    count=subprocess.check_output(cmd, shell=True).decode().strip()
    return int(count)



out=['target_bed\tintersected_num\trandomized_intersect_num\tfold_enrichment\tpval']
permut_num_out={}

for target_pos in target_poss:
    orig_count=count_intersect_num(original_pos, target_pos, chrom_size, len_extend)   # count intersected element num
    random_l=[]
    for i in range(n_permut):  # randomize
        count_random=count_randomized_intersect_num(original_pos, target_pos, chrom_size, len_extend)
        random_l.append(count_random)
    permut_num_out[target_pos]=random_l
    mean_random=mean(random_l)
    enrich= (orig_count / mean_random)  # enrich
    count_exceed_num=0
    for r in random_l:
        if r >= orig_count:
            count_exceed_num += 1
    if count_exceed_num > 0:
        pval=str(count_exceed_num / n_permut)  # pval
    else:
        pval_min= (1 / n_permut)
        pval='<' + str(pval_min)
    tmp=[target_pos, str(orig_count), str(mean_random), str(enrich), str(pval)]
    out.append('\t'.join(tmp))
    pass

with open('out_permutation_test.txt', 'w') as outfile:  # pval result
    outfile.write('\n'.join(out) +'\n')

print('\n'.join(out))

out=[]
for t in permut_num_out:
    l=[ str(i) for i in permut_num_out[t] ]
    out.append(t +'\t'+ '\t'.join(l))
with open('out_permutation_test.tsv', 'w') as outfile:  # randomized intersect num result
    outfile.write('\n'.join(out) +'\n')

