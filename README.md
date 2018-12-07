# Sanders Lab @ UCSF https://sanderslab.github.io/
This repositiory includes scripts that were created during my rotation in the Sander's lab during the fall of 2018. 

Graphic Scripts:
181112_Candlestick_Plot.R: This script was created to make a candlestick plot of mutations in cases versus controls over a moving window in the genome. While the script is currently hard coded to look at upstream/downstream conserved variants, it can easily be changed to subset down to other variables within the VEP annotation file and create the plot (ie only variants that are in ChmmState15_E1_Brain). 

For your input, you will need a VEP annotated mutation file. Additionally, you will need to decide what annotation you would like to examine, how large around a window you would like to look, and how far away from the transcription start site (TSS) you would like to examine. 


JASPAR/Motif PreProcessing Scripts:
The goal of this script is to take a variant file input (from VEP), create (1) a bed file of variants, (2) a fasta file of a certain number of basepairs around each variant.
The main script is the wrapper script. That then calls the other two python scripts.

