# ORF Finder: discard proteins shorter than 30 aa

# python ORFinder.py <input genome in FASTA file> 



import sys
from Bio import SeqIO
from Bio import Seq

##############################3333
    

        
comp_dic = {'A':'T', 'T':'A', 'G':'C', 'C':'G' }
rna_dic = {'A':'A', 'T':'U', 'G':'G', 'C':'C' }
prot_dic = {
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




fname = sys.argv[1]

sequence = SeqIO.read(open(fname), "fasta")
genome_name = fname.replace('.fa','')

blah = []
prot=''
comp=''
rna=''
protein=[]

for y in range(0,len(sequence)-2):

	begin = sequence[y]+sequence[y+1]+sequence[y+2]
	if begin !='ATG':
			continue
	else:

		for x in range(y,(len(sequence)-2),3):
			codon=sequence[x]+sequence[x+1]+sequence[x+2]

			if codon in ['TAA', 'TAG', 'TGA']:
				break

			blah.append(codon)

		prot=''.join(blah)

	if prot=='':
		print ("No start codon found")
	else:
		if len(prot) >= 50: 
			protein.append(prot)
		prot=''
		blah=[]
print ("Total number of possible proteins is ", len(protein))
#print (protein[:1], protein[-1:])

##################################################################

# Make the three reading frames (forward and reversed) and translate them into peptide sequence.
# Resources used: http://biopython.org/DIST/docs/api/Bio.Seq.Seq-class.html#translate

forwardRFs = [Seq.translate(sequence.seq[i:], table='Standard', stop_symbol='*', to_stop=False, cds=False) for i in range(3)]



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




#Use PotentialORFs.txt as output           

with open('genome_'+ genome_name + '_PotentialORFs.txt', 'w') as output: 

	for frame in range(len(forwardpept)): #Writting forward ORFs
		frame_name = str(frame)
		orf_count = 0
		for peptide in forwardpept[frame]:
			ofr_count = orf_count + 1
			name_aux = ['>'+genome_name,'ORF'+str(orf_count),'RF'+str(frame)]
			ORFarbitraryname =  '_'.join(name_aux)
			output.write("{}\n".format(ORFarbitraryname))
			output.write("{}\n".format(peptide))



	for frame in range(len(reversedpept)): # Writting reversed ORFs
		frame_name = str(frame)
		orf_count = 0
		for peptide in reversedpept[frame]:
			ofr_count = orf_count + 1
			name_aux = ['>'+genome_name,'ORF'+str(orf_count),'RF'+str(frame),'rev'] #reversed '_rev'
			ORFarbitraryname =  '_'.join(name_aux)
			output.write("{}\n".format(ORFarbitraryname))
			output.write("{}\n".format(peptide))

