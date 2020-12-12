#!/bin/sh

# merge_bed.sh
# usage: sh prog file1.bed file2.bed

# prerequisite
# bedtools version: v2.26.0

# file1 should be better condition than file2

basename1=`basename $1 .bed`
basename2=`basename $2 .bed`

bedtools subtract -A -a $2 -b $1 -f 0.95 -F 0.95 > ${basename2}_clean.bed

cat $1 ${basename2}_clean.bed | sort -n > ${basename1}plus${basename2}.bed

