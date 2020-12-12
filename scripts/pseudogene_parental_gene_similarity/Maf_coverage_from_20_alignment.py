#!/usr/bin/env python
# coding: UTF-8

"""

# usage: python %prog infile.bed
# python2.7

# example of input file:
chr1 [tab] 36352 [tab] 36365 [tab] 314
'chr_pos' + '\t' + 'start' + '\t' + 'end' + '\t' + 'ID_num'

# if input file is 'infile.bed', output file name will be 'infile_out.txt'
# output example
'ID_num' + '\t' + 'hg38_coverage' + '\t' + 'panTro4_coverage' + ....

# prerequisite
pybedtools (0.8.0)
pandas (0.23.4)
biopython (1.71)
bx-python (0.7.4)
"""

import os,sys,re
import pandas as pd
from pybedtools import BedTool
from Bio import AlignIO
from Bio.AlignIO import MafIO
from Bio.AlignIO.MafIO import MafIndex
from bx.align import maf

f_path=sys.argv[1]
cwd=os.getcwd()
outfile_dir=cwd
temporary_file_dir=cwd
basename=os.path.splitext(os.path.basename(f_path))[0]
outfile_name=basename + '_tarsier_out.txt'

maf_dir='/path/to/maf_file/dir'

species=['hg38','tarSyr2'] # tarsier
df_res=pd.DataFrame(index=[species], columns=[])

try:
    with open(f_path, 'r') as infile:
        for line in infile:
            LMR=[]
            emt=re.split(r'\t', line.strip())
            emt[1]=int(emt[1])
            emt[2]=int(emt[2])
            LMR.append(emt[0]+'\t'+str(emt[1])+'\t'+str(emt[2])+'\t'+emt[3]+'_M')

            LMR_pos=[]
            for i in LMR:
                LMR_emt=i.split()
                base=LMR_emt[3]
                LMR_pos.append(base+'_1'+'\t'+LMR_emt[1]+'\t'+LMR_emt[2])
            
            df2=pd.DataFrame(index=[species], columns=[])
            idx = MafIO.MafIndex('%s/%s.mafindex' % (maf_dir,emt[0]), '%s/%s_hu_tarsier.maf' % (maf_dir,emt[0]), 'hg38.%s' % emt[0])
            for LMR_pos_line in LMR_pos:
                emt_LMR_pos=re.split(r'\t', LMR_pos_line)
                start=int(emt_LMR_pos[1])
                end=int(emt_LMR_pos[2])
                multiple_alignment=idx.search([start], [end])
                name=LMR_pos_line.replace('\t', '.')
                AlignIO.write(multiple_alignment, '%s/%s.maf' % (temporary_file_dir, name), 'maf')
                if sum(1 for line_num in open('%s/%s.maf' % (temporary_file_dir, name), 'r')) > 4:
                    texts={}
                    for s in species: texts[s]=[]
                    with open('%s/%s.maf' % (temporary_file_dir, name), 'r') as maf_file:
                        maf_reader=maf.Reader(maf_file)
                        for m in maf_reader:
                            for s in species:
                                c=m.get_component_by_src_start(s)
                                if c: texts[s].append(c.text)
                                else: texts[s].append('-' * m.text_size)
                
                    hg38_before=''.join(texts['hg38'])
                    if hg38_before.find('-') == -1:
                        texts2=texts
                    else:
                        iter=re.finditer(r'[ATGCatgc]+', hg38_before)
                        texts2={}
                        for s in species:
                            texts2[s]=[]
                        for match in iter:
                            for s in species:
                                whole_seq=''.join(texts[s])
                                texts2[s].append(whole_seq[match.start():match.end()])

                    with open('%s/%s.maf' % (temporary_file_dir, name), 'r') as maf_file:
                        maf_file_list=[maf_file_l for maf_file_l in maf_file]
                    for i in maf_file_list:
                        if 'hg38' in i:
                            maf_start=int(i.split()[2])
                            break
                    for i in maf_file_list:
                        if 'hg38' in i:
                            k=i
                    emt_hg38=re.split(r' +', k)
                    maf_end=(int(emt_hg38[2]) + int(emt_hg38[3]))
                    fa_start=(start - maf_start)
                    seq_len=len(''.join(texts2['hg38']))
                    fa_end=(seq_len + end - maf_end)
                    count1=[]
                    for s in species:
                        seq=''.join(texts2[s])[fa_start:fa_end]
                        count1.append(len(seq.replace('-', '')))
                    count1[0]=(end - start)
                    df1=pd.Series(count1, index=[species])
                else:
                    hu_len=(end - start)
                    count1=[hu_len]
                    for s in species[1:]:
                        count1.append(0)
                    df1=pd.Series(count1, index=[species])
                df2[emt_LMR_pos[0]]=df1
                os.remove('%s/%s.maf' % (temporary_file_dir, name))

            emt_LMR_pos=re.split(r'_', LMR_pos[0])
            df3=pd.DataFrame(index=[species], columns=[])
            df3['M']=df2.filter(regex=('_M_')).sum(axis=1)
            df_res[emt_LMR_pos[0]+'_'+emt_LMR_pos[1]]=df3.apply(lambda x: 100*x/df3.iloc[0], axis=1).astype(int).astype(str)

finally:
    df_res.T.to_csv('%s/%s' %(outfile_dir, outfile_name), sep='\t')


