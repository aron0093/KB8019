import sys
import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

orf_fasta = sys.argv[1]

orf_data = pd.read_csv(orf_fasta, header = None)

for i in range(len(orf_data)):

    if orf_data.iloc[i][0][0] == '>':
    
        continue
    
    else:
    
        tran = Seq(orf_data.iloc[i][0]).translate(table="Standard", to_stop=True)
        
        orf_data.iloc[i] = str(tran)
        
orf_data.to_csv(orf_fasta[:-4]+'_translated.fasta', header = None, index = None)
