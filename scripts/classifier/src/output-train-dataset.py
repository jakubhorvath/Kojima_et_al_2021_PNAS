from module import *
import pandas as pd
import sys


def main():
    argv = sys.argv
    k = int(argv[1])

    # train datalist
    # 7 groups
    data_list = ['EVE_training',
                 'coding_exon_training',
                 'InterGenic_training',
                 'Intron_training',
                 'noncoding_exon_training',
                 'pseudogene_exon_training',
                 'TSS_training',
                 ]

    # Load dataset
    genes_List, _ = load(data_list)

    # Obtain the number of elements in each dataset
    genes_length = get_length(genes_List)

    # Generate k-mer features
    X = k_mer(genes_List, k)

    # Convert to data.frame type with pandas
    dfX = pd.DataFrame(X)
    dfLength = pd.DataFrame(genes_length)
    dfDatalist = pd.DataFrame(data_list)

    # Output data.frame in csv format
    dfX.to_csv("tmpdata/dataset.csv", header=False, index=False)
    dfLength.to_csv("tmpdata/length.csv", header=False, index=False)
    dfDatalist.to_csv("tmpdata/datalist.csv", header=False, index=False)


if __name__ == "__main__":
    main()
