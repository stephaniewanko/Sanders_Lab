
##Remove Duplicate FASTA files
#author: Stephanie Wankowicz
#date: 11/20/2018

import numpy as np
import pandas as pd
import sys
import os
import Bio
from Bio import SeqIO
import sys
from collections import defaultdict

dedup_records = defaultdict(list)

for record in SeqIO.parse(sys.argv[1],'fasta'):
    dedup_records[str(record.seq)].append(record.id)


with open("de_duplicated"+sys.argv[1], 'w') as output:
    for seq, ids in dedup_records.items():
        # Join the ids and write them out as the fasta
        output.write(">{}\n".format('|'.join(ids)))
        output.write(seq + "\n")


