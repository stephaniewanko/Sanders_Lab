# Sanders Lab @ UCSF https://sanderslab.github.io/
This repositiory includes scripts that were created during my rotation in the Sander's lab during the fall of 2018. 

181112_Candlestick_Plot.R: This script was created to make a candlestick plot of mutations in cases versus controls over a moving window in the genome. While the script is currently hard coded to look at upstream/downstream conserved variants, it can easily be changed to subset down to other variables within the VEP annotation file and create the plot (ie only variants that are in ChmmState15_E1_Brain). 

For your input, you will need a VEP annotated mutation file. Additionally, you will need to decide what annotation you would like to examine, how large around a window you would like to look, and how far away from the transcription start site (TSS) you would like to examine. 
