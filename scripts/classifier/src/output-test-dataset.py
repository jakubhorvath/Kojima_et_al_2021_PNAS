from module import *
import pandas as pd
import sys


def main():
    argv = sys.argv
    k = int(argv[1])

    # test datalist, test.fasta
    data_list = ['test']

    # Load dataset
    genes_List, genes_IdList = load(data_list)

    # Obtain the number of elements in each dataset
    genes_length = get_length(genes_List)

    # Generate k-mer features
    X = k_mer(genes_List, k)

    # Convert to data.frame type with pandas
    dfX = pd.DataFrame(X)
    dfLength = pd.DataFrame(genes_length)
    dfIdList = pd.DataFrame(genes_IdList)

    # Output data.frame in csv format
    dfX.to_csv("tmpdata/dataset.csv", header=False, index=False)
    dfLength.to_csv("tmpdata/length.csv", header=False, index=False)
    dfIdList.to_csv("tmpdata/geneid.csv", header=False, index=False)


if __name__ == "__main__":
    main()
