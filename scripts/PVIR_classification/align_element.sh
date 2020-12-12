#!/bin/sh

# align_element.sh
# usage: sh prog

# prerequisites
# mafft: v7.407 (2018/Jul/23)
# MEGA: Version: 10.0.5, Build: 10180924-x86_64

# 'test_element.fa' is a fasta file which contains a sequence in your interest

species=species_in_interest
element_name=element_name_in_interest
mkdir ${element_name}
cd ./${element_name}

prealigned_fasta_1=/path/to/prealigned_fasta_1.fa
prealigned_fasta_2=/path/to/prealigned_fasta_2.fa
prealigned_fasta_3=/path/to/prealigned_fasta_3.fa
failed_file=/path/to/failed_file.txt
mao_file=/path/to/mao_file_generated_by_megaproto.mao

linsi --thread 2 --add test_element.fa --reorder ${prealigned_fasta_1} > ${element_name}_linsi_1.fa
megacc -a ${mao_file} -d ${element_name}_linsi_1.fa -o ${element_name}_linsi_1
num=`ls | grep .nwk | wc -l | tr -d ' '`
if [ "${num}" == 0 ]; then
linsi --thread 2 --add test_element.fa --reorder ${prealigned_fasta_2} > ${element_name}_linsi_2.fa
megacc -a ${mao_file} -d ${element_name}_linsi_2.fa -o ${element_name}_linsi_2
num=`ls | grep .nwk | wc -l | tr -d ' '`
if [ "${num}" == 0 ]; then
linsi --thread 2 --add test_element.fa --reorder ${prealigned_fasta_3} > ${element_name}_linsi_3.fa
megacc -a ${mao_file} -d ${element_name}_linsi_3.fa -o ${element_name}_linsi_3
num=`ls | grep .nwk | wc -l | tr -d ' '`
if [ "${num}" == 0 ]; then
cat test_element_name.txt >> ${failed_file}
fi
fi
fi

cd ../
