#!/usr/bin/env bash

#Launch Module-1: 1) Rules-based nucleotide translation and peptide fasta/tsv generation
#Rules: 1) For riboseq GENCODE sequences specifically: 1-frame translation
#       2) Strand known + ORF unknown: 3-frame translation
#	3) Neither strand nor ORF known: 6-frame translation
#These rules apply to ANY sequence


#2) Output diagnostic proteogenomic visualizations of the translated petpide (2 randomly chosen peptides)

fasta_path="/data/neopipe/path_c_ddbj_dataset/Module_1/module1"

#nextflow run Module_1.nf --fasta $fasta_path/test5.fasta --sample demo --tissue tumor

nextflow run Module_1.nf --fasta $fasta_path/ddbj_out.fasta --sample mel04 --tissue tumor
