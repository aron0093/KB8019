# ORF Finder:

# Comparative Genomics

# Authors: Javier Lanillos,Revant Gupta, Siddharth Thomar

# This script is a complete package for detecting and analysing potential ORFs and generates a fasta file with predictions.

# The ORFs are labelled in the following format 

'''> Species name_reading frame_ORF number_prediction score_reverse/forward'''

# Run: python ORFinder.py <input genome in FASTA file> 

##############################

import sys
import re
import time
import pandas as pd
from Bio import SeqIO
from Bio import Seq
from Bio import pairwise2 as pw2
from Bio.pairwise2 import format_alignment

genome_fil = sys.argv[1]

genome = SeqIO.read(open(fname), "fasta")
genome_name = genome_fil.partition('.')[0]

		####### FORWARD STRAND #######

start_time = time.time()

print('Forward direction...')

blah = []
prot=''

# List to store the ORFs found in the forward direction no matter at which Reading Frame that is provided by the list below (fRF)

forwardORFseq=[]  

## List to store the accepted ORFs shorter than a given minimum length
#shortforwardORFseq = [] 
## This is just to store them separatedly from the longest ones
#shortfRF= []

# Reading Frame Positions Series --> y[RF+(j = 1,2,3)] = y[0]RFj + 3
# Reading Frame +1: 0, 3, 6, 9... so, for y being the base position, y[0] = 0 --> RF+1 when y % 3 == 0
# Reading Frame +2: 1, 4, 7, 10...so, for y being the base position, y[0] = 1 --> RF+1 when y % 3 == 1
# Reading Frame +3: 2, 5, 8, 11...so, for y being the base position, y[0] = 2 --> RF+1 when y % 3 == 2

ORFcoord_frame_forward = list() #Store ORF coordinates (begin, end) and reading frame

for y in range(0,len(sequence)-2):

	begin = sequence[y]+sequence[y+1]+sequence[y+2]
	if begin !='ATG':

		ORF_begin = y
		
		if y % 3 == 0:
			actualRF = 1 # actualRF is the Reading Frame corresponding to the ORF we have found
		elif y % 3 == 1:
			actualRF = 2
		elif y % 3 == 2:
			actualRF = 3
		continue
	else:
		for x in range(y,(len(sequence)-2),3):
			codon=sequence[x]+sequence[x+1]+sequence[x+2]

			if codon in ['TAA', 'TAG', 'TGA']:
				ORF_end = x
				break

			blah.append(codon)

		prot=''.join(blah)
	if prot=='':
		print ("No start codon found")
	else:
		
		if ORF_end-ORF_begin  >= 100: #len(prot) >= 100: 
			ORFcoord_frame_forward.append([ORF_begin+1,ORF_end,actualRF])
			forwardORFseq.append(prot)

		prot=''
		blah=[]

print ("Total number of possible ORFs is ", len(forwardORFseq))

# Save forwards ORFs coordinates into a file

orf_count = 0
with open(genome_name + '_ORFs_forward.txt', 'w') as output: 

	for coord_Frame in ORFcoord_frame_forward: #Writting forward ORFs by coordinates

		orf_count = orf_count + 1

		ORFname = [genome_name,'RF'+str(coord_Frame[2]),'ORF'+str(orf_count)]
		ORFname =  '_'.join(ORFname)

		output.write("{}\t{}\t{}\n".format(ORFname,coord_Frame[0],coord_Frame[1]))


print("--- %s seconds ---" % (time.time() - start_time))

       	####### REVERSED STRAND #######

print('Reversed direction...')
start_time = time.time()
sequence = sequence[::-1]
ORFcoord_frame_reversed = list()
#shortORF_coord_reversed= list()
# Repeat the same for the reversed strand
blah = []
prot=''

reversedORFseq=[]  


## List to store the accepted ORFs shorter than a given minimum length
#shortreversedORFseq = [] 
## This is just to store them separatedly from the longest ones
#shortrRF= []

for y in range(0,len(sequence)-2):


	begin = sequence[y]+sequence[y+1]+sequence[y+2]
	if begin !='ATG':

		ORF_begin = y
		if y % 3 == 0:
			actualRF = 1 # actualRF is the Reading Frame corresponding to the ORF we have found
		elif y % 3 == 1:
			actualRF = 2
		elif y % 3 == 2:
			actualRF = 3
		continue
	else:
		for x in range(y,(len(sequence)-2),3):
			codon=sequence[x]+sequence[x+1]+sequence[x+2]

			if codon in ['TAA', 'TAG', 'TGA']:
				ORF_end = x
				break

			blah.append(codon)

		prot=''.join(blah)
	if prot=='':
		print ("No start codon found")
	else:
		
		if ORF_end-ORF_begin  >= 100: #len(prot) >= 100: 
			ORFcoord_frame_reversed.append([ORF_begin+1,ORF_end,actualRF])
			reversedORFseq.append(prot)

		prot=''
		blah=[]
print ("Total number of possible ORFs is ", len(reversedORFseq))
#print ("Total number of short possible ORFs (50 to 100 bp.) ", len(shortreversedORFseq))

# Save reversed ORFs coordinates into a file

orf_count = 0
with open(genome_name + '_ORFs_reversed.txt', 'w') as output: 

	for coord_Frame in ORFcoord_frame_reversed: #Writting reversed ORFs by coordinates

		orf_count = orf_count + 1

		ORFname = [genome_name,'RF'+str(coord_Frame[2]),'ORF'+str(orf_count),'rev']
		ORFname =  '_'.join(ORFname)

		output.write("{}\t{}\t{}\n".format(ORFname,coord_Frame[0],coord_Frame[1]))

print("--- %s seconds ---" % (time.time() - start_time))

##################################################################

# This part of the code is in case we need to translate the sequences

"""
#comp_dic = {'A':'T', 'T':'A', 'G':'C', 'C':'G' }
#rna_dic = {'A':'A', 'T':'U', 'G':'G', 'C':'C' }
#prot_dic = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
	}




# Make the three reading frames (forward and reversed) and translate them into peptide sequence.
# Resources used: http://biopython.org/DIST/docs/api/Bio.Seq.Seq-class.html#translate

forwardRFs_prot = [Seq.translate(sequence.seq[i:], table='Standard', stop_symbol='*', to_stop=False, cds=False) for i in range(3)]

#forwardRFs_prot = [Seq.translate(sequence.seq[i:], table='Standard', stop_symbol='*', to_stop=False, cds=False) for i in range(3)]



reversedRFs = [Seq.translate(sequence[::-1].seq[i:], table='Standard', stop_symbol='*', to_stop=False, cds = False) for i in range(3)]

forwardpept = [[] for x in range(len(forwardRFs))] # List of lists to save peptide sequences from ORF in forward direction. Each nested list contains all the ORFs for one of the three reading frames

for framenumber in range(len(forwardRFs)):

	frame = forwardRFs[framenumber]

	for peptide in frame.split('*'): #Split the translated sequences when there is a stopcodon
		if len(peptide) > 30:
			forwardpept[framenumber].append(peptide)




reversedpept = [] # List to save peptide sequences from ORF in reversed direction


reversedpept = [[] for x in range(len(reversedRFs))] # List of lists to save peptide sequences from ORF in forward direction. Each nested list contains all the ORFs for one of the three reading frames

for framenumber in range(len(reversedRFs)):

	frame = reversedRFs[framenumber]

	for peptide in frame.split('*'): #Split the translated sequences when there is a stopcodon
		if len(peptide) > 30:
			reversedpept[framenumber].append(peptide)

"""


