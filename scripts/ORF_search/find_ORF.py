#!/usr/bin/env python

"""

# usage: python %prog
# python3

matplotlib    3.0.3

"""

# map start codon, stop codons and ORF

import os,sys,re,shutil
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

seq=''  # DNA seq of interest
seq=seq.upper()

def parse_to_three_frames(pos_list):
    frame1_list=[]
    frame2_list=[]
    frame3_list=[]
    for pos in pos_list:
        if   pos[0] % 3 == 0: frame1_list.append(pos)
        elif pos[0] % 3 == 1: frame2_list.append(pos)
        elif pos[0] % 3 == 2: frame3_list.append(pos)
    return [frame1_list, frame2_list, frame3_list]

seq_len=len(seq)
full_pos=[[0, seq_len]]

atg_pos=[]
iter=re.finditer('ATG', seq)
if iter is not None:
    for match in iter:
        pos=[match.span()[0], match.span()[1]]
        atg_pos.append(pos)

term_pos=[]
for codon in ['TAG', 'TGA', 'TAA']:
    iter=re.finditer(codon, seq)
    if iter is not None:
        for match in iter:
            pos=[match.span()[0], match.span()[1]]
            term_pos.append(pos)

orf_pos=[]
for start in atg_pos:
    tmp_orf_pos=[]
    orf_len=seq_len
    for stop in term_pos:
        if stop[0] > start[0]:
            new_orf_len= (stop[0] - start[0])
            if (new_orf_len < orf_len) and (new_orf_len % 3 == 0):
                tmp_orf_pos=[start[0], stop[0]]
                orf_len=new_orf_len
    if tmp_orf_pos:
        orf_pos.append(tmp_orf_pos)
    pass

frame1atg,frame2atg,frame3atg = parse_to_three_frames(atg_pos)
frame1term,frame2term,frame3term = parse_to_three_frames(term_pos)
frame1orf,frame2orf,frame3orf = parse_to_three_frames(orf_pos)


# plot
plt.figure(figsize=(5,1.5))
gs=gridspec.GridSpec(3, 1, height_ratios=[1,1,1])
gs.update(hspace=0.1)

pos_list=[[frame1orf, frame1term, frame1atg],
          [frame2orf, frame2term, frame2atg],
          [frame3orf, frame3term, frame3atg]]
colors=['lightskyblue', 'coral', 'limegreen']

n=0
for f in pos_list:
    ax=plt.subplot(gs[n])
    for l,c in zip(f, colors):
        for p in l:
            rect=matplotlib.patches.Rectangle((p[0],0), p[1]-p[0], 1, color=c, ec=None)
            ax.add_patch(rect)
        ax.set_xlim(full_pos[0])
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(0.5)
        ax.yaxis.set_ticks([])
        if not n == (len(pos_list) - 1):
            ax.xaxis.set_ticks([])
    n+=1

plt.savefig('plot_out.pdf')


