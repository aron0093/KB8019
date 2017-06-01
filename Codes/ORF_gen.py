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

genome_type = input('Select genome type :   1 - Eukaryotic   0 - Prokaryotic ')

genome_record = SeqIO.read(open(genome_fil), "fasta")

genome = genome_record.seq

reverse_genome  = genome.reverse_complement()

species = genome_fil.partition('.')[0]

species_dic = {'05': 'Chlamydia trachomatis', '08': 'Dictyoglomus turgidum', '14': 'Mycobacterium tuberculosis', '16': 'Rhodopirellula baltica', '25': 'Saccharomyces cerevisiae'}



out_file = open(species+'_orfs.fasta', 'a')

start_time = time.time()

if genome_type == '1':

        ####### FORWARD STRAND #######

    print('Forward direction...')

    # List to store the ORFs found in the forward direction no matter at which Reading Frame that is provided by the list below (fRF)

    # Reading Frame Positions Series --> y[RF+(j = 1,2,3)] = y[0]RFj + 3
    
    # Reading Frame +1: 0, 3, 6, 9... so, for y being the base position, y[0] = 0 --> RF+1 when y % 3 == 0
    # Reading Frame +2: 1, 4, 7, 10...so, for y being the base position, y[0] = 1 --> RF+1 when y % 3 == 1
    # Reading Frame +3: 2, 5, 8, 11...so, for y being the base position, y[0] = 2 --> RF+1 when y % 3 == 2

    ORF_name_counter = 0
    tata_begin  = 0
    caat_begin  = 0
    kozak_a_begin = 0
    kozak_b_begin = 0
    
    
    for y in range(0,len(genome)-2):
    
    # Proto ORF positions

        begin = genome[y]+genome[y+1]+genome[y+2]
        
        if begin =='ATG':

            ORF_begin = y
   
            actualRF = (y%3)+1 # actualRF is the Reading Frame corresponding to the ORF we have found
                      
            for x in range(y,(len(genome)-2),3):
            
                codon=genome[x]+genome[x+1]+genome[x+2]

                if codon in ['TAA', 'TAG', 'TGA']:
                
                    ORF_end = x
                    
                    break
                    
    # Length check and search space definitions            
                                
            if ORF_end-ORF_begin  >= 100: #len(prot) >= 100:
            
                if ORF_begin >= 99: 
            
                    promo_search_space = range((ORF_begin-99), ORF_begin-34)
                                
                else:
                
                    promo_search_space = range(0, ORF_begin-34)
                    
                if ORF_begin >= 19: 
            
                    rbs_search_space = range((ORF_begin-19), ORF_begin+19)
                                
                else:
                
                    rbs_search_space = range(0, ORF_end+19)  
                    
    # Check for promoters                
                                    
                for i in promo_search_space: 
                
                # Check TATA

                    tata_frame = str(genome[i:i+6])
        
                    if tata_frame[:4] == 'TATA':

                        alignment = pw2.align.globalms("TATAAA", tata_frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 12  

                        align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))
                        
                        if int(align_data[17:]) >= 10:
                        
                            tata_begin = i
                            
                # Check CAAT
                    
                    caat_frame = str(genome[i:i+9])
        
                    if caat_frame[0] == 'G' or caat_frame[0] == 'A' and caat_frame[3:7] == 'CCAA':            
                  
                        alignment = pw2.align.globalms("GGCCAATCT", caat_frame, 2,0,-1,-1,one_alignment_only = True)     # Max score 18

                        align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))
                
                        if int(align_data[23:]) >= 10:
                
                            caat_begin  = i 
                              
    # Check RBS
                
                for j in rbs_search_space:
                
                # First Kozak
                
                    rbs_frame = str(genome[j:j+10])
                            
                    if rbs_frame[:-4] == 'ATGG' and rbs_frame[3] == 'A':

                        alignment = pw2.align.globalms("GCCACCATGG", rbs_frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 20

                        align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))        
                            
                        if int(align_data[25:]) >= 12:    
                            
                            kozak_a_begin  = 10
                            
                    if rbs_frame[:-4] == 'ATGG' and rbs_frame[3] == 'G':

                        alignment = pw2.align.globalms("GCCGCCATGG", rbs_frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 20

                        align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))        
                            
                        if int(align_data[25:]) >= 12:    
                            
                            kozak_b_begin  = 10        
                   
   # Output Orfs
   
                score = 0            
                                 
                if tata_begin > 0:
                                
                    score += 1
                                 
                if caat_begin > 0:
                
                    score += 1
                 
                if kozak_a_begin > 0 or kozak_b_begin > 0:
                
                    score += 1
                    
                if score > 0 and tata_begin > 0:
                            
                    out_file.write('> '+species+'_RF'+str(actualRF)+'_Sco_'+str(score)+'_ORF'+str(ORF_name_counter)+'\n')
                    out_file.write(str(genome[ORF_begin:ORF_end])+'\n')
                    
                    ORF_name_counter += 1
                    
                elif score > 0 and caat_begin > 0: 
                
                    out_file.write('> '+species+'_RF'+str(actualRF)+'_Sco_'+str(score)+'_ORF'+str(ORF_name_counter)+'\n')
                    out_file.write(str(genome[ORF_begin:ORF_end])+'\n')
                    
                    ORF_name_counter += 1
                    
                elif score > 0:
                
                    out_file.write('> '+species+'_RF'+str(actualRF)+'_Sco_'+str(score)+'_ORF'+str(ORF_name_counter)+'\n')
                    out_file.write(str(genome[ORF_begin:ORF_end])+'\n')
                    
                    ORF_name_counter += 1
        

    print("--- %s seconds ---" % (time.time() - start_time))

               ####### REVERSED STRAND #######

    print('Reversed direction...')
    

    tata_begin  = 0
    caat_begin  = 0
    kozak_a_begin = 0
    kozak_b_begin = 0
    

    for y in range(0,len(reverse_genome)-2):
    
    # Proto ORF positions

        begin = reverse_genome[y]+reverse_genome[y+1]+reverse_genome[y+2]
        
        if begin =='ATG':

            ORF_begin = y
   
            actualRF = (y%3)+1 # actualRF is the Reading Frame corresponding to the ORF we have found
                    
            for x in range(y,(len(reverse_genome)-2),3):
            
                codon=reverse_genome[x]+reverse_genome[x+1]+reverse_genome[x+2]

                if codon in ['TAA', 'TAG', 'TGA']:
                
                    ORF_end = x
                    
                    break
                    
        # Length check and search space definitions            
                                    
            if ORF_end-ORF_begin  >= 100: #len(prot) >= 100:
            
                if ORF_begin >= 99: 
            
                    promo_search_space = range((ORF_begin-99), ORF_begin-34)
                                
                else:
                
                    promo_search_space = range(0, ORF_begin-34)
                    
                if ORF_begin >= 19: 
            
                    rbs_search_space = range((ORF_begin-19), ORF_begin+19)
                                
                else:
                
                    rbs_search_space = range(0, ORF_end+19)  
                    
    # Check for promoters                
                                    
                for i in promo_search_space: 
                
                # Check TATA

                    tata_frame = str(genome[i:i+6])
        
                    if tata_frame[:4] == 'TATA':

                        alignment = pw2.align.globalms("TATAAA", tata_frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 12  

                        align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))
                        
                        if int(align_data[17:]) >= 10:
                        
                            tata_begin = i
                            
                # Check CAAT
                    
                    caat_frame = str(genome[i:i+9])
        
                    if caat_frame[0] == 'G' or caat_frame[0] == 'A' and caat_frame[3:7] == 'CCAA':            
                  
                        alignment = pw2.align.globalms("GGCCAATCT", caat_frame, 2,0,-1,-1,one_alignment_only = True)     # Max score 18

                        align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))
                
                        if int(align_data[23:]) >= 10:
                
                            caat_begin  = i 
                              
    # Check RBS
                
                for j in rbs_search_space:
                
                # First Kozak
                
                    rbs_frame = str(genome[j:j+10])
                            
                    if rbs_frame[:-4] == 'ATGG' and rbs_frame[3] == 'A':

                        alignment = pw2.align.globalms("GCCACCATGG", rbs_frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 20

                        align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))        
                            
                        if int(align_data[25:]) >= 12:    
                            
                            kozak_a_begin  = 10
                            
                    if rbs_frame[:-4] == 'ATGG' and rbs_frame[3] == 'G':

                        alignment = pw2.align.globalms("GCCGCCATGG", rbs_frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 20

                        align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))        
                            
                        if int(align_data[25:]) >= 12:    
                            
                            kozak_b_begin  = 10        
                           
       # Output Orfs
       
       
                score = 0            
                                 
                if tata_begin > 0:
                                
                    score += 1
                                 
                if caat_begin > 0:
                
                    score += 1
                 
                if kozak_a_begin > 0 or kozak_b_begin > 0:
                
                    score += 1
                    
                if score > 0 and tata_begin > 0:
                            
                    out_file.write('> '+species+'_RF'+str(actualRF)+'_Sco_'+str(score)+'_ORF'+str(ORF_name_counter)+'_rev'+'\n')
                    out_file.write(str(reverse_genome[ORF_begin:ORF_end])+'\n')
                    
                    ORF_name_counter += 1
                    
                elif score > 0 and caat_begin > 0: 
                
                    out_file.write('> '+species+'_RF'+str(actualRF)+'_Sco_'+str(score)+'_ORF'+str(ORF_name_counter)+'_rev'+'\n')
                    out_file.write(str(reverse_genome[ORF_begin:ORF_end])+'\n')
                    
                    ORF_name_counter += 1
                    
                elif score > 0:
                
                    out_file.write('> '+species+'_RF'+str(actualRF)+'_Sco_'+str(score)+'_ORF'+str(ORF_name_counter)+'_rev'+'\n')
                    out_file.write(str(reverse_genome[ORF_begin:ORF_end])+'\n')
                    
                    ORF_name_counter += 1
                

    print("--- %s seconds ---" % (time.time() - start_time))

elif genome_type == '0':

    print('Forward direction...')

    # List to store the ORFs found in the forward direction no matter at which Reading Frame that is provided by the list below (fRF)

    # Reading Frame Positions Series --> y[RF+(j = 1,2,3)] = y[0]RFj + 3
    
    # Reading Frame +1: 0, 3, 6, 9... so, for y being the base position, y[0] = 0 --> RF+1 when y % 3 == 0
    # Reading Frame +2: 1, 4, 7, 10...so, for y being the base position, y[0] = 1 --> RF+1 when y % 3 == 1
    # Reading Frame +3: 2, 5, 8, 11...so, for y being the base position, y[0] = 2 --> RF+1 when y % 3 == 2

    ORF_name_counter = 0
    promo_a_begin  = 0
    promo_b_begin = 0
    tle_begin = 0
    

    for y in range(0,len(genome)-2):
    
    # Proto ORF positions

        begin = genome[y]+genome[y+1]+genome[y+2]
        
        if begin =='ATG':

            ORF_begin = y
   
            actualRF = (y%3)+1 # actualRF is the Reading Frame corresponding to the ORF we have found
                     
            for x in range(y,(len(genome)-2),3):
            
                codon=genome[x]+genome[x+1]+genome[x+2]

                if codon in ['TAA', 'TAG', 'TGA']:
                
                    ORF_end = x
                    
                    break
                    
        # Length check and search space definitions            
                                    
            if ORF_end-ORF_begin  >= 100: #len(prot) >= 100:
            
                if ORF_begin >= 99: 
            
                    promo_search_space = range((ORF_begin-99), ORF_begin-34)
                                
                else:
                
                    promo_search_space = range(0, ORF_begin-34)
                    
                if ORF_begin >= 7: 
            
                    rbs_search_space = range((ORF_begin-7), ORF_begin+7)
                                
                else:
                
                    rbs_search_space = range(0, ORF_end+7)  
                    
    # Check for promoters                
                                    
                for i in promo_search_space: 
                
                # Check Pribnow

                    promo_frame = str(genome[i:i+6])
        
                    alignment_a = pw2.align.globalms("TATAAT", promo_frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 12  

                    align_data_a = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment_a[0]))
                    
                    alignment_b = pw2.align.globalms("TTGACA", promo_frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 12  

                    align_data_b = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment_b[0]))
                        
                    if int(align_data_a[17:]) >= 10:
                    
                        promo_a_begin = i 
                    
                    if int(align_data_b[17:]) >= 10:
                        
                        promo_b_begin = i                          
                              
    # Check TLE
                
                for j in rbs_search_space:
                               
                    tle_frame = str(genome[j:j+6])
                              
                    alignment = pw2.align.globalms("AGGAGG", tle_frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 20

                    align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))        
                            
                    if int(align_data[17:]) >= 10:    
                            
                        tle_begin  = 10
                                
                           
       # Output Orfs
       
                score = 0            
                                 
                if promo_a_begin > 0:
                                
                    score += 1
                
                if promo_b_begin > 0:
                                
                    score += 1
                                            
                if tle_begin > 0:
                
                    score += 1
                    
                if score > 0 and promo_a_begin > 0:
                            
                    out_file.write('> '+species+'_RF'+str(actualRF)+'_Sco_'+str(score)+'_ORF'+str(ORF_name_counter)+'\n')
                    out_file.write(str(genome[ORF_begin:ORF_end])+'\n')
                    
                    ORF_name_counter += 1
                    
                elif score > 0 and promo_b_begin > 0:
                            
                    out_file.write('> '+species+'_RF'+str(actualRF)+'_Sco_'+str(score)+'_ORF'+str(ORF_name_counter)+'\n')
                    out_file.write(str(genome[ORF_begin:ORF_end])+'\n')
                    
                    ORF_name_counter += 1
                    
                elif score > 0 and tle_begin > 0: 
                
                    out_file.write('> '+species+'_RF'+str(actualRF)+'_Sco_'+str(score)+'_ORF'+str(ORF_name_counter)+'\n')
                    out_file.write(str(genome[ORF_begin:ORF_end])+'\n')
                    
                    ORF_name_counter += 1
                    

    print("--- %s seconds ---" % (time.time() - start_time))


               ####### REVERSED STRAND #######

    print('Reversed direction...')

    # List to store the ORFs found in the forward direction no matter at which Reading Frame that is provided by the list below (fRF)

    # Reading Frame Positions Series --> y[RF+(j = 1,2,3)] = y[0]RFj + 3
    
    # Reading Frame +1: 0, 3, 6, 9... so, for y being the base position, y[0] = 0 --> RF+1 when y % 3 == 0
    # Reading Frame +2: 1, 4, 7, 10...so, for y being the base position, y[0] = 1 --> RF+1 when y % 3 == 1
    # Reading Frame +3: 2, 5, 8, 11...so, for y being the base position, y[0] = 2 --> RF+1 when y % 3 == 2


    promo_a_begin  = 0
    promo_b_begin = 0
    tle_begin = 0
    

    for y in range(0,len(reverse_genome)-2):
    
    # Proto ORF positions

        begin = reverse_genome[y]+reverse_genome[y+1]+reverse_genome[y+2]
        
        if begin =='ATG':

            ORF_begin = y
   
            actualRF = (y%3)+1 # actualRF is the Reading Frame corresponding to the ORF we have found
                       
            for x in range(y,(len(reverse_genome)-2),3):
            
                codon=reverse_genome[x]+reverse_genome[x+1]+reverse_genome[x+2]

                if codon in ['TAA', 'TAG', 'TGA']:
                
                    ORF_end = x
                    
                    break
                    
        # Length check and search space definitions            
                                
            if ORF_end-ORF_begin  >= 100: #len(prot) >= 100:
            
                if ORF_begin >= 99: 
            
                    promo_search_space = range((ORF_begin-99), ORF_begin-34)
                                
                else:
                
                    promo_search_space = range(0, ORF_begin-34)
                    
                if ORF_begin >= 7: 
            
                    rbs_search_space = range((ORF_begin-7), ORF_begin+7)
                                
                else:
                
                    rbs_search_space = range(0, ORF_end+7)  
                    
    # Check for promoters                
                                    
                for i in promo_search_space: 
                
                # Check Pribnow

                    promo_frame = genome[i:i+6]
        
                    alignment_a = pw2.align.globalms("TATAAT", promo_frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 12  

                    align_data_a = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))
                    
                    alignment_b = pw2.align.globalms("TTGACA", promo_frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 12  

                    align_data_b = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))
                        
                    if int(align_data_a[17:]) >= 10:
                    
                        promo_a_begin = i 
                    
                    if int(align_data_b[17:]) >= 10:
                        
                        promo_b_begin = i                          
                              
    # Check TLE
                
                for j in rbs_search_space:
                               
                    rbs_frame = genome[j:j+6]
                              
                    alignment = pw2.align.globalms("AGGAGG", rbs_frame, 2,0,-1,-1,one_alignment_only = True)   # Max score 20

                    align_data = re.sub('[^A-Za-z0-9]+', '', format_alignment(*alignment[0]))        
                            
                    if int(align_data[17:]) >= 10:    
                            
                        tle_begin  = 10
                                
                       
       # Output Orfs
       
                score = 0            
                                 
                if promo_a_begin > 0:
                                
                    score += 1
                
                if promo_b_begin > 0:
                                
                    score += 1
                                            
                if tle_begin > 0:
                
                    score += 1
                    
                if score > 0 and promo_a_begin > 0:
                            
                    out_file.write('> '+species+'_RF'+str(actualRF)+'_Sco_'+str(score)+'_ORF'+str(ORF_name_counter)+'_rev'+'\n')
                    out_file.write(str(reverse_genome[ORF_begin:ORF_end])+'\n')
                    
                    ORF_name_counter += 1
                    
                elif score > 0 and promo_b_begin > 0:
                            
                    out_file.write('> '+species+'_RF'+str(actualRF)+'_Sco_'+str(score)+'_ORF'+str(ORF_name_counter)+'_rev'+'\n')
                    out_file.write(str(reverse_genome[ORF_begin:ORF_end])+'\n')
                    
                    ORF_name_counter += 1
                    
                elif score > 0 and tle_begin > 0: 
                
                    out_file.write('> '+species+'_RF'+str(actualRF)+'_Sco_'+str(score)+'_ORF'+str(ORF_name_counter)+'_rev'+'\n')
                    out_file.write(str(reverse_genome[ORF_begin:ORF_end])+'\n')
                    
                    ORF_name_counter += 1
                        

    print("--- %s seconds ---" % (time.time() - start_time))
    
else:

    print('Genome type not recognised! Please try again')


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


