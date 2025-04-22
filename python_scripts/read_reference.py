import sys

from pprint import pprint

File = open(sys.argv[1],'r')

Genome = File.readlines()

File.close()

startPos = int(sys.argv[2]) ### start position of the gene on the reference genome

endPos = int(sys.argv[3]) ### end position of the gene on the reference genome

length = endPos - startPos

outFile = open(sys.argv[4], 'w')

gene = {} ### create a dictionary for the refenence sequence of the gene. Key is base position, value is the nucleotide

for i in range(startPos,endPos+1): ## populate the dict the base calls for each position
	gene[i] = Genome[1][i]

seq= "".join(gene.values()) ## join the values of the dictionary to create a sequence

outFile.write(Genome[0] + seq + "\n")
outFile.close()