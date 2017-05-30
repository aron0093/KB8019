# Script to generate database with the coordinates of Promoters and RBSs

# https://en.wikipedia.org/wiki/Ribosome-binding_site

# Eukaryotic : Kozak Sequence https://en.wikipedia.org/wiki/Kozak_consensus_sequence (gcc)gccRccAUGG \ R is purine (A, G)

# Prokaryotic : Shine-Dalgarno sequence https://en.wikipedia.org/wiki/Shine-Dalgarno_sequence AGGAGG

import sys
import re
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2 as pw2
from Bio.pairwise2 import format_alignment

genome_file = sys.argv[1]

genome_type = input('Select genome type :   1 - Eukaryotic   0 - Prokaryotic ')

genome_data = pd.read_csv(genome_file)

genome = genome_data[genome_data.columns[0]].iloc[0]

species_dic = {'05': 'Chlamydia trachomatis', '08': 'Dictyoglomus turgidum', '14': 'Mycobacterium tuberculosis', '16': 'Rhodopirellula baltica', '25': 'Saccharomyces cerevisiae'}

species = species_dic[genome_data.columns[0][3:5]]

out_file = species+'_trans_eles.txt'

if genome_type == '1':

    out_data = pd.DataFrame(columns = ['Translation linked elements', 'Alignment', 'Score', 'Start Pos', 'End Pos', 'Frame' ], index = list(range(2*len(genome)-20)))

    for i in range(0, len(genome)-10): # For first Kozak
    
        rf = (i%3)+1
    
        frame = genome[i:i+10]
        
        if frame[:-4] == 'ATGG' and frame[3] == 'A':
    
            alignment = pw2.align.globalms("GCCACCATGG", frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 20

            align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))
            
            out_data['Translation linked elements'][i] = 'GCCACCATGG'
            
            out_data['Alignment'][i] = align_data[10:20]
            
            out_data['Score'][i] = align_data[25:]
            
            out_data['Start Pos'][i] = str(i+1)
            
            out_data['End Pos'][i] = str(i+11)
            
            out_data['Frame'][i] = rf
            
    for i in range(0, len(genome)-10): # For second Kozak
    
        rf = (i%3)+1
    
        frame = genome[i:i+10]        
        
        if frame[:-4] == 'ATGG' and frame[3] == 'G':
    
            alignment = pw2.align.globalms("GCCGCCATGG", frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 20

            align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))
            
            out_data['Translation linked elements'][i+len(genome)-10] = 'GCCGCCATGG'
            
            out_data['Alignment'][i+len(genome)-10] = align_data[10:20]
            
            out_data['Score'][i+len(genome)-10] = align_data[25:]
            
            out_data['Start Pos'][i+len(genome)-10] = str(i+1)
            
            out_data['End Pos'][i+len(genome)-10] = str(i+11)
            
            out_data['Frame'][i+len(genome)-10] = rf
                   
    out_data.dropna(axis = 0, how = 'any', inplace=True)
    
    out_data.reset_index(drop = True, inplace = True)
    
    out_data.to_csv(out_file, sep = '\t', index = None)
    
elif genome_type == '0':

    out_data = pd.DataFrame(columns = ['Translation linked elements', 'Alignment', 'Score', 'Start Pos', 'End Pos', 'Directionality', 'Frame' ], index = list(range(len(genome)-6)))

    for i in range(0, len(genome)-6): # For both Pibnow boxes
    
        rf = (i%3)+1
        
        frame = genome[i:i+6]
        
        # First
        
        alignment = pw2.align.globalms("AGGAGG", frame, 2,0,-1,-1,one_alignment_only = True)     # Max Score 12

        align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))
        
        out_data['Translation linked elements'][i] = 'AGGAGG'
        
        out_data['Alignment'][i] = align_data[6:12]
        
        out_data['Score'][i] = align_data[17:]
        
        out_data['Start Pos'][i] = str(i+1)
        
        out_data['End Pos'][i] = str(i+7)
    
        out_data['Frame'][i] = rf

        
    out_data.dropna(axis = 0, how = 'any', inplace=True) 
    
    out_data.reset_index(drop = True, inplace = True)         
        
    out_data.to_csv(out_file, sep = '\t', index = None)
            
else:

    print('Genome type not recognised! Please try again')            
            
