#!/bin/bash
#Replication Wrapper Script
#Created by: Stephanie Wankowicz
#Last edited: 1/18/2018

#___________________________________________________________________________#
#This script should be run for each replicate.
#The first part will extract the barcodes from the DNA and RNA.
#The second script will output some QC graphics
#The third script will find significant or outlier variants
#The forth script will output expected sequences for the outlier variants (you can change it to output expected sequences for all)
  #This can be feed into MEMO.
#The hard coded portions are for the Sanders' lab MPRA experiment from Fall 2018.
#Outputs: QC graphics, barcode counts for DNA and RNA, outlier variant files, expected sequences of outlier variants

##REQUIREMENTS##
#R version 3.x
#Python version 3 

###INPUTS####
RNA_FASTQ= #Fastq R1 RNA 
DNA_FASTQ= #Fastq R1 DNA
Pickle= #pickle output from initial MPRA wrapper script
Date=          #for the output; example 181205
dir= #directory where your files are located. Coded as tools (BOWTIE and bamtools are in seperate folders)

#________________________________________________________________________#
###Extracting Barcodes###
printf('Extracting Barcode Counts')
python barcode_counter2.py #need to edit to be able to accept args
#this will spit out raw counts for RNA and DNA. The next script will calculate an 'MRPA scrore' of RNA/DNA. 
printf('Calculating relative barcodes from DNA/RNA')
python barcode_pipeline.py #need to edit to be able to accept args
#________________________________________________________________________#

###QC Graphics###
printf('Creating some QC metrics and graphics')
R MRPA_QC_graphics.R #need to edit to be able to accept args (input=CSV file)

#________________________________________________________________________#
###Extracting Outliers###
printf('Extacting Outliers.')
R Extract_Outliers.R #need to edit to be able to accept args (input=CSV file) Also, combine with script above.

#________________________________________________________________________#
###Expected Sequences###
./Jaspar_preprorcessing_wrapper_script.sh #need to take inputs
#This also calls 181024_Jasper_PreProcessing.py

printf('DONE! You can put expected sequences into MEMO.')
printf('Also! Make sure you re-run this script for every ')
