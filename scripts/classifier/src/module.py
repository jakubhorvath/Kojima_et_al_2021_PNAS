# coding: UTF-8

import numpy as np
import re
from Bio import SeqIO


def load(data_list_name):
    """
    Load .fasta data and store gene sequences to output gene list

        Parameters
        ----------
        data_list_name : list
            list of target data set names

        Returns
        -------
        genes_list : list
            Parsed dataset
        genes_idlist : List
            Parsed id of the dataset
    """
    print('# Loading Dataset...')
    # initialize Genes list dictionary
    genes_list = {}
    genes_idlist = {}

    for i in range(len(data_list_name)):
        print('# Loading "{}.fasta"'.format(data_list_name[i]))
        filename = '../data/%s.fasta' % data_list_name[i]
        gene_list = []
        gene_idlist = []
        with open(filename, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                id = record.id
                seq = record.seq
                gene_idlist.append(">"+str(id))
                gene_list.append(str(seq))

        genes_list[data_list_name[i]] = list(gene_list)
        genes_idlist[data_list_name[i]] = list(gene_idlist)

    return genes_list, genes_idlist


def get_length(input_data):
    """
        Obtain the number of elements in each dataset

        Parameters
        ----------
        input_data : list
            Parsed dataset of genes

        Returns
        -------
        length : list
            the length of each dataset
    """
    length = []
    for data in input_data.values():
        length.append(len(data))

    return length


def k_mer(input_data, k):
    """
        Generate k-mer features

        Parameters
        ----------
        input_data : list
            Parsed dataset of genes
        k : int
            k the numebr of k of k-mer.

        Returns
        -------
        features_data : ndarray
            k-mer dataset
    """
    print('# Generating Feature Matrix... : K = {}'.format(k))
    all_data = []
    for key in input_data:
        all_data = all_data + input_data[key]

    # A->0, T->1, G->2, C->3, others->4
    trans_dict = str.maketrans('ATGC', '0123')
    itr_list = list((map((lambda x: x.translate(trans_dict)), all_data)))
    itr_list = list((map((lambda x: re.sub(r'\D', '4', x)), itr_list)))

    size = [len(itr_list)]
    size.extend([5]*k)
    data = np.zeros(tuple(size))
    for itr in range(len(itr_list)):
        sample = itr_list[itr]
        if k == 1:
            for i in range(len(sample)):
                data[itr][int(sample[i])] += 1
        if k == 2:
            for i in range(len(sample)-1):
                data[itr][int(sample[i])][int(sample[i+1])] += 1
        if k == 3:
            for i in range(len(sample) - 2):
                data[itr][int(sample[i])][int(sample[i + 1])
                                          ][int(sample[i + 2])] += 1
        if k == 4:
            for i in range(len(sample) - 3):
                data[itr][int(sample[i])][int(sample[i + 1])
                                          ][int(sample[i + 2])][int(sample[i + 3])] += 1
        if k == 5:
            for i in range(len(sample) - 4):
                data[itr][int(sample[i])][int(
                    sample[i + 1])][int(sample[i + 2])][int(sample[i + 3])][int(sample[i + 4])] += 1
        data[itr] = data[itr] / np.sum(data[itr])
    for i in range(k):
        data = np.delete(data, [4], axis=(i+1))
    features_data = data.reshape((len(itr_list), 4**k))

    return features_data
