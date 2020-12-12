#!/usr/bin/env python

"""

# Maf_coverage.py
# usage: python %prog infile.bed
# python2.7

# example of input file:
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
outfile_name=basename + '_out.txt'

hu_fasta='/path/to/hg38_RepeatMasked_fasta_file'
maf_dir='/path/to/maf_file'

species=['hg38','panTro4','nomLeu3','macFas5','chlSab2','calJac3','saiBol1','otoGar3','tupChi1','mm10','ochPri3','oviAri3','canFam3','loxAfr3']
df_res=pd.DataFrame(index=[species], columns=[])

with open(f_path, 'r') as infile:
    for line in infile:
        LMR=[]
        emt=re.split(r'\t', line.strip())
        emt[1]=int(emt[1])
        emt[2]=int(emt[2])
        LMR.append(emt[0]+'\t'+str(emt[1]-1000)+'\t'+str(emt[1])+'\t'+emt[3]+'_L')
        LMR.append(emt[0]+'\t'+str(emt[1])+'\t'+str(emt[2])+'\t'+emt[3]+'_M')
        LMR.append(emt[0]+'\t'+str(emt[2])+'\t'+str(emt[2]+1000)+'\t'+emt[3]+'_R')

        inbed=BedTool('\n'.join(LMR), from_string=True)
        fa=inbed.sequence(fi=hu_fasta, name=True)

        LMR_pos=[]
        with open(fa.seqfn) as seq:
            for line in seq:
                num=1
                if line.find('>') == 0:
                    emt=re.split(r':|-', line.replace('>', '').strip())
                else:
                    iter=re.finditer(r'[ATGC]+', line)
                    if iter is not None:
                        for match in iter:
                            if ( match.span()[1] - match.span()[0] ) > 100:
                                LMR_pos.append(emt[0]+'_'+str(num)+'\t'+str(int(emt[3])+match.span()[0])+'\t'+str(int(emt[3])+match.span()[1]))
                                num=num+1

        joined_LMR_pos='\n'.join(LMR_pos)
        if not ( '_L' in joined_LMR_pos and '_M' in joined_LMR_pos and '_R' in joined_LMR_pos ):
            del LMR_pos[:]
            for i in LMR:
                LMR_emt=re.split(r'\t', i)
                LMR_pos.append(LMR_emt[3]+'_1'+'\t'+LMR_emt[1]+'\t'+LMR_emt[2])
                    
        df2=pd.DataFrame(index=[species], columns=[])
        idx = MafIO.MafIndex('%s/%s.mafindex' % (maf_dir,emt[2]), '%s/%s.maf' % (maf_dir,emt[2]), 'hg38.%s' % emt[2])
        for LMR_pos_line in LMR_pos:
            emt_LMR_pos=re.split(r'\t', LMR_pos_line)
            start=int(emt_LMR_pos[1])
            end=int(emt_LMR_pos[2])
            multiple_alignment=idx.search([start], [end])
            name=LMR_pos_line.replace('\t', '.')
            AlignIO.write(multiple_alignment, '%s/%s.maf' % (temporary_file_dir, name), 'maf')
            if sum(1 for line_num in open('%s/%s.maf' % (temporary_file_dir, name), 'r')) > 3:
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
                df1=pd.DataFrame([hu_len,0,0,0,0,0,0,0,0,0,0,0,0,0], index=[species])
            df2[emt_LMR_pos[0]]=df1
            os.remove('%s/%s.maf' % (temporary_file_dir, name))

        emt_LMR_pos=re.split(r'_', LMR_pos[0])
        df3=pd.DataFrame(index=[species], columns=[])
        df3['L']=df2.filter(regex=('_L_')).sum(axis=1)
        df3['M']=df2.filter(regex=('_M_')).sum(axis=1)
        df3['R']=df2.filter(regex=('_R_')).sum(axis=1)
        df4=df3.apply(lambda x: 100*x/df3.iloc[0], axis=1).astype(int).astype(str)
        df_res[emt_LMR_pos[0]]=df4.apply(lambda x : ','.join(x), axis=1)

df_res.T.to_csv('%s/%s' %(outfile_dir, outfile_name), sep='\t')

