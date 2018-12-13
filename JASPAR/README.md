The goal of these scripts is to take a variant file input (from VEP), create (1) a bed file of variants, (2) a fasta file of a certain number of basepairs around each variant. This script was used to put variants into transcription factor binding site software, such as the MEME Suite of tools (http://meme-suite.org/index.html). 


Prerequisites:
(1)Python-3.0
(2)Bedtools
(3)Human FASTA



Wrapper Script:
    This will call the other python scripts and BEDTOOLs scripts.
    Inputs: 
      (1) VEP annotated mutation file
      (2) Size of window you would like to examine around the variant (ie. if you want to look at 100bp on either side of the       SNV, size=100.
