#written by: Stephanie Wankowicz
#10/12/2018 
#take 2

#packages
import pandas as pd
import numpy as np
import pickle
import collections
from collections import defaultdict
import sys
import math
from collections import Counter
import pysam
from tqdm import tqdm

#set up
seq_list=[]
#fastq=sys.argv[1]

f=open("12_Variant_coords_to_barcodes.pickle", "rb")
barcode_pickle = pickle.load(f)

print(len(barcode_pickle.values()))
print(barcode_pickle.values()[1])
print(len(barcode_pickle.keys()))
print(barcode_pickle.keys()[1])
fastq = pysam.FastxFile("Steph-DNA-rep1_S18_R1_001.fastq")

n_fastq_records = 44633047 #wc -l Ste /4
barcode_seqs=set(barcode_pickle.keys())
coords_to_barcodes = defaultdict(list)
append = seq_list.append

for i,barcode in tqdm(enumerate(fastq), 'barcodes',total=n_fastq_records):
	if barcode.sequence in barcode_seqs:
		append(barcode.sequence)
print(len(seq_list))
barcode_count=Counter(seq_list)
df_barcodes=pd.DataFrame.from_dict(barcode_count, orient='index').reset_index()
df_barcodes.columns=['Barcode','count']
print(df_barcodes.head())
df_barcodes.to_csv('Variant_barcodes.csv')


def barcode_to_variant(df):
    df["Variant"] = ''
    for index,row in df.iterrows():
	variant=barcode_pickle.get(df.loc[index,'Barcode'])
        df.loc[index,"Variant"]=str(variant)
    return (df)

barcode_variant=barcode_to_variant(df_barcodes)
output=fastq+'R1_output'
barcode_variant.to_csv(output+'.csv')
