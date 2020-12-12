#!/usr/bin/env python

"""
# usage: python %prog aligned.fa
# python3.7
"""

# make alignment from aligned.fa

import os,sys,re
from statistics import mean
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

cwd=os.getcwd()
outfile_dir=cwd   # outputs results in this dir
f_path=sys.argv[1]

split_num=10  # bin for alignment identity




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




# parse fasta
fa_for_plot=parse_fasta(f_path)
label=list(fa_for_plot.keys())


# prepare for plot
rect={}
for s in label:
    rect[s]={}
    iter=re.finditer(r'[atgcnATGCN]+', fa_for_plot[s])
    if iter is not None:
        rect[s]['start']=[]
        rect[s]['end']=[]
        for match in iter:
            rect[s]['start'].append(match.span()[0])
            rect[s]['end'].append(match.span()[1])
        
count=[]
atgc=['A','G','T','C','N']
seqlen=len(fa_for_plot[label[0]])
for l in range(seqlen):
    chars=''
    for s in label:
        if s in fa_for_plot:
            chars=chars+fa_for_plot[s][l]
    ntcount=[]
    for nt in atgc:
        ntcount.append(chars.upper().count(nt))
    count.append(max(ntcount))
pos=[]
ave=[]
seqpos=range(len(count))
for i in range(0, len(count), split_num):
    pos.append(mean(seqpos[i:i+split_num]))
    ave.append(mean(count[i:i+split_num]))



# plt
full_length=len(count)

grid_num=len(label) + 1
ratio=[]
for i in range(grid_num-1):
    ratio.append(5)
ratio.append(8)

plt.figure(figsize=(5,grid_num*0.2))
gs=gridspec.GridSpec(grid_num, 1, height_ratios=ratio)
gs.update(hspace=0.2)
plt.rcParams['font.size']=5


# plot alignment
n=0
for s in rect:
    l=label[n]
    ax=plt.subplot(gs[n])
    if rect[s]:
        for s,e in zip(rect[s]['start'], rect[s]['end']):
            r=matplotlib.patches.Rectangle((s,0), e-s, 1, color='black', ec=None)
            ax.add_patch(r)
    ax.set_xlim([0,full_length])
    ax.text(0, 1, l, ha='left', va='bottom', fontsize=2)
    ax.axis('off')
    n += 1


# plot identity
ax=plt.subplot(gs[n])
ax.fill_between(pos, ave, facecolor='#6495ED80', edgecolor='black', linewidth=0.5)
ax.text(0, grid_num-1, 'Identity', ha='left', va='bottom', fontsize=2)  # if figure, comment out
for axis in ['top','right']:
    ax.spines[axis].set_linewidth(0)
ax.spines['bottom'].set_linewidth(0.5)
ax.spines['left'].set_linewidth(0.5)
ax.set_xlim([0,full_length])
ax.set_ylim([0,grid_num-2])
ax.tick_params('both', width=0.5)
ax.yaxis.set_ticks([0,grid_num-2])
ax.set_yticklabels([0, 100])

plt.subplots_adjust(left=0.05, right=0.95, bottom=0.1, top=0.95)
plt.savefig('plot_out.pdf')



