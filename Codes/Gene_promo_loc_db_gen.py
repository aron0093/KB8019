# Script to generate database with the coordinates of Promoters and RBSs

# https://en.wikipedia.org/wiki/Promoter_(genetics)#Identification_of_relative_location

# Eukaryotic :
 
#   TATA Box - 5' TATAAA 3' 25-35   https://en.wikipedia.org/wiki/TATA_box
#   CAAT Box - 5' GGCCAATCT 3'  60-100   Bidirectional  https://en.wikipedia.org/wiki/CAAT_box

# Prokaryotic : 

#   Pribnow Box - 5' TATAAT 3'  10  https://en.wikipedia.org/wiki/Pribnow_box
#   Probnow Box - 5' TTGACA 3'  35  https://en.wikipedia.org/wiki/Promoter_(genetics)#Bacterial

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

out_file = species+'_promoters.txt'

if genome_type == '1':

    out_data = pd.DataFrame(columns = ['Promoter', 'Alignment', 'Score', 'Start Pos', 'End Pos', 'Directionality', 'Frame' ], index = list(range(2*len(genome)-15)))

    for i in range(0, len(genome)-6): # For TATA Box
    
        rf = (i%3)+1
    
        frame = genome[i:i+6]
        
        if frame[:4] == 'TATA':
    
            alignment = pw2.align.globalms("TATAAA", frame, 2,0,-1,-1)   # Max score 12  

            align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))
            
            out_data['Promoter'][i] = 'TATAAA'
            
            out_data['Alignment'][i] = align_data[6:12]
            
            out_data['Score'][i] = align_data[17:]
            
            out_data['Start Pos'][i] = str(i+1)
            
            out_data['End Pos'][i] = str(i+7)
            
            out_data['Directionality'][i] = 'forward'
            
            out_data['Frame'][i] = rf
        
         
    for i in range(0, len(genome)-9): # For CAAT Box
    
        rf = (i%3)+1
        
        frame = genome[i:i+9]
        
        if frame[0] == 'G' or frame[0] == 'A' and frame[3:7] == 'CCAA':
        
            alignment = pw2.align.globalms("GGCCAATCT", frame, 2,0,-1,-1)     # Max score 18

            align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))
            
            out_data['Promoter'][i+len(genome)-6] = 'GGCCAATCT'
            
            out_data['Alignment'][i+len(genome)-6] = align_data[9:18]
            
            out_data['Score'][i+len(genome)-6] = align_data[23:]
            
            out_data['Start Pos'][i+len(genome)-6] = str(i+1)
            
            out_data['End Pos'][i+len(genome)-6] = str(i+10)
            
            out_data['Directionality'][i+len(genome)-6] = 'bi'
            
            out_data['Frame'][i+len(genome)-6] = rf    
    
    out_data.dropna(axis = 0, how = 'any', inplace=True)
    
    out_data.reset_index(drop = True, inplace = True)
    
    out_data.to_csv(out_file, sep = '\t', index = None)
    
elif genome_type == '0':

    out_data = pd.DataFrame(columns = ['Promoter', 'Alignment', 'Score', 'Start Pos', 'End Pos', 'Directionality', 'Frame' ], index = list(range(2*len(genome)-12)))

    for i in range(0, len(genome)-6): # For both Pibnow boxes
    
        rf = (i%3)+1
        
        frame = genome[i:i+6]
        
        # First
        
        alignment = pw2.align.globalms("TATAAT", frame, 2,0,-1,-1)     # Max Score 12

        align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))
        
        out_data['Promoter'][i] = 'TATAAT'
        
        out_data['Alignment'][i] = align_data[6:12]
        
        out_data['Score'][i] = align_data[17:]
        
        out_data['Start Pos'][i] = str(i+1)
        
        out_data['End Pos'][i] = str(i+7)
        
        out_data['Directionality'][i] = 'forward'   
    
        out_data['Frame'][i] = rf
        
        # Second
    
        alignment = pw2.align.globalms("TTGACA", frame, 2,0,-1,-1)     # Max Score 12

        align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))
        
        out_data['Promoter'][i+len(genome)-6] = 'TTGACA'
        
        out_data['Alignment'][i+len(genome)-6] = align_data[6:12]
        
        out_data['Score'][i+len(genome)-6] = align_data[17:]
        
        out_data['Start Pos'][i+len(genome)-6] = str(i+1)
        
        out_data['End Pos'][i+len(genome)-6] = str(i+7)
        
        out_data['Directionality'][i+len(genome)-6] = 'forward'
        
        out_data['Frame'][i+len(genome)-6] = rf
        
    out_data.dropna(axis = 0, how = 'any', inplace=True) 
    
    out_data.reset_index(drop = True, inplace = True)         
        
    out_data.to_csv(out_file, sep = '\t', index = None)
            
else:

    print('Genome type not recognised! Please try again')            
            
