#!/usr/bin/env python

"""

# maf_retrieve_to_fa.py
# usage: python %prog infile.bed
# python2.7

# infile.bed requires 4 rows; chr* +'\t'+ start +'\t'+ end +'\t'+ name

"""

import os,sys,re
from Bio import AlignIO
from Bio.AlignIO import MafIO
from Bio.AlignIO.MafIO import MafIndex
from bx.align import maf

f_path=sys.argv[1]
maf_path='/path/to/directory_where_contains_maf_files'
species=['hg38','panTro4','gorGor3','ponAbe2','nomLeu3','rheMac3','macFas5','papAnu2','chlSab2','calJac3','saiBol1','otoGar3','tupChi1','speTri2','jacJac1','micOch1','criGri1','mesAur1','mm10','rn6','hetGla2','cavPor3','chiLan1','octDeg1','oryCun2','ochPri3','susScr3','vicPac2','camFer1','turTru2','orcOrc1','panHod1','bosTau8','oviAri3','capHir1','equCab2','cerSim1','felCat8','canFam3','musFur1','ailMel1','odoRosDiv1','lepWed1','pteAle1','pteVam1','myoDav1','myoLuc2','eptFus1','eriEur2','sorAra2','conCri1','loxAfr3','eleEdw1','triMan1','chrAsi1','echTel2','oryAfe1','dasNov3','monDom5','sarHar1','macEug2','ornAna1','falChe1','falPer1','ficAlb2','zonAlb1','geoFor1','taeGut2','pseHum1','melUnd1','amaVit1','araMac1','colLiv1','anaPla1','galGal4','melGal1','allMis1','cheMyd1','chrPic2','pelSin1','apaSpi1','anoCar2','xenTro7','latCha1','tetNig2','fr3','takFla1','oreNil2','neoBri1','hapBur1','mayZeb1','punNye1','oryLat2','xipMac1','gasAcu1','gadMor1','danRer10','astMex1','lepOcu1','petMar2']

# retrieve maf from intervals; then convert maf to fasta
with open(f_path, 'r') as infile:
    for line in infile:
        lsplit=line.split()
        start=int(lsplit[1])
        end=int(lsplit[2])
        idx=MafIO.MafIndex('%s/%s.mafindex' % (maf_path, lsplit[0]), '%s/%s.maf' % (maf_path, lsplit[0]), 'hg38.%s' % lsplit[0])
        multiple_alignment=idx.search([start], [end])
        AlignIO.write(multiple_alignment, './retrieved_maf/%s.maf' % lsplit[3], 'maf')

        texts={}
        for s in species: texts[s]=[]
        with open('./retrieved_maf/%s.maf' % lsplit[3], 'r') as inmaf:
            maf_reader=maf.Reader(inmaf)
            for m in maf_reader:
                for s in species:
                    c=m.get_component_by_src_start(s)
                    if c: texts[s].append(c.text)
                    else: texts[s].append('-' * m.text_size)
        
        outfa=[]
        for s in species:
            outfa.append('>' + s)
            outfa.append(''.join(texts[s]))
        with open('./retrieved_fa/%s.fa' % lsplit[3], 'w') as outfilefa:
            outfilefa.write('\n'.join(outfa) + '\n')

