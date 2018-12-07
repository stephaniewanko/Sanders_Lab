#!/bin/sh
#author: Stephanie Wankowicz	
#created:12/05/2018
#purpose: from VEP annotated variant list, 1) create and output bedfile of these variants, 2)create fasta file for each variant
#toosl used: python, bedtools,human genome fasta (from:), 

##########inputs#############
flanking_size='100' #size of the region you want around your mutation (ie if you want a window of 100bp, put in 100) 
variant_list='/Users/stephaniewankowicz/Downloads/Sanders_Lab/list_gene.hqDNV_sscwgs.20171115.P231_WGS519_256.hg38.vep_gene_20171229.txt'
output_name=181205_test_5kb
#############################


python 181024_Jasper_PreProcessing.py $flanking_size $variant_list $output_name


#output_fasta_name=output+'.fa'
output_fasta_name=$output_name'_fasta.fa'
input_bed=$output_name'.bed'
bedtools getfasta -fi /Users/stephaniewankowicz/Downloads/Sanders_Lab/hg38/hg38_combined.fasta -bed $input_bed -fo $output_fasta_name


#for file in /Users/stephaniewankowicz/Downloads/Sanders_Lab/hg38/* #where your hg38 fasta files are located
#do
#    output_fasta_name=$output_name'_fasta_'${file##/*/}
#    printf $output_fasta_name
#    bedtools getfasta -fi $file -bed 181205_test_5kb.bed -fo $output_fasta_name
#done



##de-duplicate fasta list
python remove_duplicates.py $output_fasta_name

#read in fasta file, if duplicate, write down the name of each header, write only one version of fasta file with both headers and delete other copy
#add up the number of duplicates
