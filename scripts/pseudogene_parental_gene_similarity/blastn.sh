
# blastn to calculate pseudogene-original_gene nucleotide identity
# GRCh38_latest_rna.fna was downloaded from NCBI

cat simian_specific_pseudogenes.bed | bedtools getfasta -fi hg38.fa -bed 'stdin' -name -fo simian_specific_pseudogenes.fa
makeblastdb -in GRCh38_latest_rna.fna -out /GRCh38_latest_rna.fna.blastdb -dbtype nucl -parse_seqids
blastn -db GRCh38_latest_rna.fna.blastdb -query simian_specific_pseudogenes.fa -evalue 1e-10 -num_alignments 10 -word_size 7 -num_threads 6 -out blastn_pseudogenes_out.txt -outfmt 5


