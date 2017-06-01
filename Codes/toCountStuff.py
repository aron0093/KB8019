# To count stuff in a FASTA file

import sys

f = open(sys.argv[1],'r') 

#count = 0

count_1 = 0
count_2 = 0
count_3 = 0

for i in f.readlines():

    if i [0] == '>':
#        count += 1
        
        words = i.partition('_')
        
        score = int(words[2][8])
        
        if score == 1:
        
            count_1 += 1
        
        elif score == 2:
        
            count_2 += 1
            
        elif score == 3:
        
            count_3 += 1
            
        else:
        
            print('Error')
f.close()
#print(count)

print(count_1)
print(count_2)
print(count_3)
