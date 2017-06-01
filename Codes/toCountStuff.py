# To count stuff in a FASTA file
import sys

f = open(sys.argv[1],'r') 
t = list()

for i in f.readlines():
	if i [0] == '>':
		t.append(i)

f.close()
print(len(t))

