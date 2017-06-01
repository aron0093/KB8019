import sys
import re
from Bio.Blast import NCBIXML

blastOutputXMLFile = sys.argv [1]
#ref_tag = sys.argv[2]
#target_tag = sys.argv[3]

blastOutputXMLHandle = open (blastOutputXMLFile)

listOfBlastRecords = NCBIXML.parse (blastOutputXMLHandle) #Generator to iterate through each genome
title_list = list()

for aSingleBlastRecord in listOfBlastRecords:

    for i in range (len (aSingleBlastRecord.alignments)):   # Iterate through each query seq
    
        if i == 0: # Confine results to first, i.e. best hit

            description = aSingleBlastRecord.descriptions [i]

            title = re.compile ("gnl\|BL_ORD_ID\|\d* ").sub ("", description.title)
            title_list.append(title)
            #ref = aSingleBlastRecord.query
            
            #print (ref_tag+' '+ref+' '+target_tag+' '+title+' ') # Redirect to get file

print(len(title_list))
print(len(set(title_list)))
