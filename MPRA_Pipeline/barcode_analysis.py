#!/usr/bin/env python
#last edited: 12-08-2018
#last edited by: Stephanie Wankowicz

import os
import pandas as pd
import pickle


RNA = pd.read_csv("RNA_Rep3_combined_pickle.csv",sep=',')
DNA = pd.read_csv("DNA_Rep3_combined_pickle.csv",sep=',')
print(DNA.head())
DNA.columns=['X','Barcode','DNA_count', 'Variant']
RNA.columns=['X','Barcode','RNA_count', 'Variant']
print(RNA.head())


DNA_RNA=pd.merge(DNA, RNA,how='outer',left_on="Barcode", right_on="Barcode")
DNA_RNA=DNA_RNA.drop(['X_x', 'X_y'], axis=1)
DNA_RNA['RNA/DNA']=DNA_RNA['RNA_count']/DNA_RNA['DNA_count']
print(DNA_RNA.head())
output=pd.DataFrame(columns=['Number_of_Barcodes', 'Variant','Mean_DNA_RNA'])
print(output.head())
print(DNA_RNA.groupby('Variant_x').sum())
summary_table=DNA_RNA.groupby('Variant_x').mean() #rename
summary_table.to_csv('Rep3_SummaryTable.csv')
