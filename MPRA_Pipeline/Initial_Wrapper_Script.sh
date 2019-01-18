#!/bin/bash
#Created by: Stephanie Wankowicz
#Last edited: 1/14/2019

#Please note, this script can take up to 15 hours.
#___________________________________________________________________________#
#This script takes in 3 associated fastqs and spit out a pickle with barcodes-variant. It is currently configured to look at exact matches of both the variant expected fasta and the barcodes.
#The hard coded portions are for the Sanders' lab MPRA experiment from Fall 2018.
#Outputs: Variant Sequences aligned to expected sequence (bam), pickle with variant-barcode, feather files with statistics on the variant-barcode pickle

###INPUTS####
R1= #Fastq R1
R2= #Fastq R2
R3= #Fastq R3
Date=          #for the output; example 181205
dir= #directory where your files are located. Coded as tools (BOWTIE and bamtools are in seperate folders)

#________________________________________________________________________#
###BOWTIE###
printf('Creating Bowtie Library File')
$dir/bowtie2-2.3.4.3-linux-x86_64/bowtie2-build sequence.Sanders_renamed_20181128.fasta Sanders_expected_sequences
printf('Aligning FASTQs. Please note this can take up to 8 hours')
$dir/bowtie2-2.3.4.3-linux-x86_64/bowtie2 -q --very-sensitive -t -x Sanders_expected_sequences -1 $R1 -2 $R3 >$Date+_Sanders_Asso.sam\

#________________________________________________________________________#
###Bamtools###
printf('Converting SAM file to BAM file')
samtools view -bS $Date+_Sanders_Asso.sam  > $Date+_Sanders_Asso.bam
printf('Selecting only exact matches from Associated Bam')

$dir/bamtools/build/bin/bamtools filter -tag XM:0 -in $Date+_Sanders_Asso.bam  -out Perfect_$Date+_Sanders_Asso.bam

#________________________________________________________________________#
###Mapping Barcodes###
printf('Mapping Barcodes to Variants.')
python map_barcodes.py #BAM, #R2 #re-write script to take in arguements
printf('DONE! Move onto wrapper script #2.')
