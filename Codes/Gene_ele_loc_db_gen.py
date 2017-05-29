# Script to generate database with the coordinates of Promoters and RBSs

# https://en.wikipedia.org/wiki/Promoter_(genetics)#Identification_of_relative_location

# Eukaryotic :
 
#   TATA Box - 5' TATAAA 3' 25-35   https://en.wikipedia.org/wiki/TATA_box
#   CAAT Box - 5' GGCCAATCT 3'  60-100   Bidirectional  https://en.wikipedia.org/wiki/CAAT_box

# Prokaryotic : 

#   Pribnow Box - 5' TATAAT 3'  10  https://en.wikipedia.org/wiki/Pribnow_box
#   Probnow Box - 5' TTGACA 3'  35  https://en.wikipedia.org/wiki/Promoter_(genetics)#Bacterial

import sys
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

out_file = open(species+'_promoters.txt', 'a')

out_data = pd.DataFrame(columns = ['Promoter', 'Alignment', 'Score', 'Start Pos', 'End Pos', 'Directionality' ])

if genome_type == '1':

    for i in range(0, len(genome)-6): # For TATA Box
    
        rf = (i % 3) + 1
    
        frame = genome[i:i+6]
        
        if frame[:4] == 'TATA':
    
            alignment = pw2.align.globalms("TATAAA", frame, 2,0,-1,-1)     

            align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0])))
            
            out_data['Promoter'][i] = 'TATAAA'
            
            out_data['Alignment'][i] = align_data[6:12]
            
            out_data['Score'][i] = align_data[17:]
            
            out_data['Start Pos'][i] = str(i+1)
            
            out_data['Stop Pos'][i] = str(i+7)
            
            out_data['Directionality'][i] = 'forward'
            
         
    for i in range(0, len(genoem)-         
            
            
            
            
            
