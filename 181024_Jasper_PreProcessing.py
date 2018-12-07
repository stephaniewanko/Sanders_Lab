
##JASPER Pre Processing
#author: Stephanie Wankowicz
#date: 10/24/2-18

import numpy as np
import pandas as pd
import sys
import os
import Bio
from Bio import SeqIO
import sys

size=int(sys.argv[1])

#import variant list
print('Loading in Variant List.')
variant_list=pd.read_table(sys.argv[2])

print('Preprocessing Input File')
variant_list['Chromosome'] = variant_list.ID.str[0:5]
variant_list['Chromosome']=variant_list['Chromosome'].map(lambda x:x.rstrip(':')) #this is removing all trailing ':' characters from the function
variant_list['Position']=variant_list.ID.str[5:]


print('Formatting Variant List.')
dd_chrom=['chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']
for index, row in variant_list.iterrows():
    if variant_list.loc[index,'Chromosome'] in dd_chrom:
        variant_list.loc[index,'Position']=variant_list.loc[index,'Position'].split(':', 2)[1]
    else:
        variant_list.loc[index,'Position']=variant_list.loc[index,'Position'].split(':',1)[0]
variant_list=variant_list.reset_index(drop=True)



print('Creating Bedfile.')
variant_list_bed=pd.DataFrame(columns=['chrom','chromStart','chromEnd'], index=range(0,len(variant_list.index)))
for index, row in variant_list.iterrows():
    variant_list_bed.loc[index]['chromStart']=int(variant_list.loc[index]['Position'])
    variant_list_bed.loc[index]['chrom']=variant_list.loc[index]['Chromosome']


variant_list_bed['chromEnd']=variant_list_bed['chromStart']+size
variant_list_bed['chromStart']=variant_list_bed['chromStart']-size

variant_list_bed.to_csv(sys.argv[3]+'.bed',index=False, sep='\t', header=False)

print('Moving onto generate FASTA file.')



##I am keeping this in here, but bedtools is 100x faster
'''
#print('Loading in HG38 FASTA Sequence')
#human_fasta=list(SeqIO.parse("Homo_sapiens.GRCh38.dna.primary_assembly.fa", 'fasta'))
#seqDict=SeqIO.to_dict(human_fasta)

#print('Creating new fasta!')
#new_fasta = []
#print(variant_list_bed.head())
#for index, row in variant_list_bed.iterrows():
#    chr=str(variant_list_bed.loc[index,"chrom"])
    start=int(variant_list_bed.loc[index,"chromStart"])
    end=int(variant_list_bed.loc[index,"chromEnd"])
    ID=chr+str(variant_list_bed.loc[index,"chromStart"])
    print(ID)
    chr_num=chr[3:]
    print(chr_num)
    long_seq_record=seqDict[chr_num]
    #print(long_seq_record)
    long_seq = long_seq_record.seq
    alphabet = long_seq.alphabet
    short_seq = str(long_seq)[start-1:end]
    variant_list_bed.loc[index,"Sequence"]=short_seq
    #print(short_seq)
    new_fasta.append('>%s\n%s' % (ID, short_seq))
#output_file.close()
with open('Variant_200_bp.fasta', 'w') as f:
    f.write('\n'.join(new_fasta))
'''
