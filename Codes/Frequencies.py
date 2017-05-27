import sys
import pandas as pd

genome_fil = open (sys.argv [1])

genome = pd.read_csv(genome_fil)


chars = list(genome[genome.columns[0]].iloc[0])

st_nuc = ['A', 'T', 'C', 'G']

# Nucleotide frequencies

st_count = 0

nuc_freq = {}

for char in chars:
    
    if char in st_nuc:
    
        st_count += 1
    
        try:
    
            nuc_freq[char] += 1

        except KeyError:
    
            nuc_freq[char] = 1

nuc_freq.update(((k, v / float(st_count))) for k,v in nuc_freq.items())
    
print('Nucleotide frequencies are')
print('\n')
print(nuc_freq)
print('\n')

for k, v in nuc_freq.items():

    open('nuc_freq_out.txt', 'a').write(genome.columns[0][3:5]+'\t'+k+'\t'+str(v)+'\n')

# Dinucleotide frequencies

dinuc_freq = {}

for i in range(len(chars)-1):

    if char in st_nuc:

        try:
    
            dinuc_freq[chars[i]+chars[i+1]] += 1
    
        except KeyError:
    
            dinuc_freq[chars[i]+chars[i+1]] = 1
        
dinuc_freq.update((k, v / float((st_count)-1)) for k,v in dinuc_freq.items())

print('Dinucleotide frequencies are')
print('\n')
print(dinuc_freq)
print('\n')

for k, v in dinuc_freq.items():

    open('dinuc_freq_out.txt', 'a').write(genome.columns[0][3:5]+'\t'+k+'\t'+str(v)+'\n')
    
# GC content

print('GC content is ', (nuc_freq['G']+nuc_freq['C']))

open('GC_con_out.txt', 'a').write(genome.columns[0][3:5]+'\t'+str(nuc_freq['G']+nuc_freq['C'])+'\n')

# Amino acid freq

from Bio import SeqIO
from Bio import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

sequence = SeqIO.read(open(sys.argv[1]), "fasta")

translation = [Seq.translate(sequence.seq[i:], table='Standard', stop_symbol='*', to_stop=False, cds=False) for i in range(3)]

print('Amino acid frequencies in all 3 reading frames are : ')
print('\n')

f_count = 0

for frames in translation:

    f_count += 1    
    
    analysed_seq = ProteinAnalysis(str(frames))

    anal = analysed_seq.count_amino_acids()
    
    anal.update((k, v / float(len(str(frames)))) for k,v in anal.items())
    
    print('Amino frequencies are')
    print('\n')
    print(anal)  
    print('\n')
    
    for k, v in anal.items():

        open('anal_out.txt', 'a').write(genome.columns[0][3:5]+'\t'+str(f_count)+'\t'+k+'\t'+str(v)+'\n')
    
    diamino_freq = {}
    
    for i in range(len(str(frames))-1):

        try:
    
            diamino_freq[str(frames)[i]+str(frames)[i+1]] += 1
    
        except KeyError:
    
            diamino_freq[str(frames)[i]+str(frames)[i+1]] = 1
        
    diamino_freq.update((k, v / float((len(str(frames))-1))) for k,v in diamino_freq.items())

    print('Diamino frequencies are')
    print('\n')
    print(diamino_freq)
    print('\n')

    for k, v in diamino_freq.items():

        open('diamino_freq_out.txt', 'a').write(genome.columns[0][3:5]+'\t'+str(f_count)+'\t'+k+'\t'+str(v)+'\n')
























