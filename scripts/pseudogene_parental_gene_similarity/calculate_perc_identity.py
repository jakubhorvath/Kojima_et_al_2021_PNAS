#!/usr/bin/env python

"""

# usage: python %prog
# python3

biopython     1.72

"""

import os,sys,re,shutil
from Bio.Blast import NCBIXML

out=[]
blast_records=NCBIXML.parse(open('blastn_pseudogenes_out.txt'))
for blast_record in blast_records:
    pseudo_id=blast_record.query.split(':')[0]
    for alignment in blast_record.alignments:
        hit_id=alignment.title.split('|')[1].split('.')[0]
        if pseudo_id == hit_id:
            no_gap_lens=[]
            idents=[]
            for hsp in alignment.hsps:
                no_gap_lens.append(hsp.align_length - hsp.gaps)
                idents.append(hsp.identities)
            percent_ident=100 * (sum(idents) / sum(no_gap_lens))
            out.append(blast_record.query +'\t'+ str(percent_ident))

print(len(out))

with open('pseudogene_perc_ident.txt', 'w') as outfile:
    outfile.write('\n'.join(out) +'\n')






