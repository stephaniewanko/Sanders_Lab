The goal of this script is to take a variant file input (from VEP), create (1) a bed file of variants, (2) a fasta file of a certain number of basepairs around each variant.
The main script is the wrapper script. That then calls the other two python scripts. For your input, you will need a VEP annotated mutation file. Additionally, you will need to decide what annotation you would like to examine, how large around a window you would like to look, 
and how far away from the transcription start site (TSS) you would like to examine. 
