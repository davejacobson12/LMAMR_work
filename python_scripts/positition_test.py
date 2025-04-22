import sys
import re

vcfFile = open(sys.argv[1], "r")

fileHandle = sys.argv[1].split('.')[0] ### get the name of the file so you can have the output in fasta format at the end

vcfLines = vcfFile.readlines()

vcfFile.close()

fastaOut = open(sys.argv[5], 'w') ### name of output file

refFile = open(sys.argv[4],'r') ### get the length of the reference genome

refGenome = refFile.readlines() 

refFile.close()

genomeLength = len(refGenome[1]) ### the second index is the genome, so get the length

startPos = int(sys.argv[2]) ### start position of the gene on the reference genome

endPos = int(sys.argv[3]) ### end position of the gene on the reference genome

vcfSplit = [None] * genomeLength ### create empty list. I'm not sure that this length needs to be this long
vcfSeq = ['N'] * genomeLength ### populate the final sequence list with Ns at each position for the length of the reference genome. The length will change depending on which genome you mapped to
#posVCF = ['-'] * len(vcfLines)

for i in range(1,len(vcfLines)):
	vcfSplit[i] = vcfLines[i].split() #### populate the empty list with strings seperated by space, so each column can be indexed
	if vcfSplit[i][0] == 'NC_015385.1': ### make sure that I'm matching just the succinifaciens genome
		position = int(vcfSplit[i][1])
		infoColumn = vcfSplit[i][7].split(';') ## create object that holds the column info with depth and forward/reverse info
		dpString = infoColumn[0].split('=') ## depth of coverage is the first index within that object and the numerical values are after the =, so split by =
		depth = int(dpString[1])
		if depth >=5: ### if depth of coverage is greater than 10, check to see if there is an alternate base
			altLen = len(vcfSplit[i][4]) ## check to see if there is an alternate nucleotide
			if altLen == 3:  ### if the alternate length column is equal to 3, then the base call is the same as the reference
				vcfSeq[position] = vcfSplit[i][3]
				#fastaOut.write(vcfSeq[position])
			elif altLen >= 6: ### to get rid of any amibigous alternative base calls
				vcfSeq[position] = "-"
				#fastaOut.write(vcfSeq[position])
			else: ### if the alternative length is not three or greater than 5
				I16 = infoColumn[1].split('=')[1] ### pull out the string of the I16 column
				ref_alt = I16.split(',')[0:4] ### this string is seperated by commas, so split. We are concerned with the first 4 indices, 1 = reference forward, 2 = rerference reverse, 3 = alt forward, 4 = alt reverse
				heterozygosity= (int(ref_alt[2]) + int(ref_alt[3]))/(int(ref_alt[0]) + int(ref_alt[1])+ int(ref_alt[2]) + int(ref_alt[3])) ### find the proportion of reads that are from the alternate
				if heterozygosity >=0.8: ### for calls that are definitely the alternative nucleotide
					heteroColumn = vcfSplit[i][4] ### pull out the column with the alt bases and <*> syntax
					bases = heteroColumn.split('<')[0] ### split by < to isolate the individual base calls
					altCall = bases.replace(',','') ### because there can be multiple alternative calls, but we only want sites with one base call, so the sites with multiple alts I will put a hyphen
					if len(altCall) >= 2:
						vcfSeq[position] = "-"
						#fastaOut.write(vcfSeq[position])
					elif altCall == 'A':
						vcfSeq[position] = altCall
						#fastaOut.write(vcfSeq[position]) 
					elif altCall == 'T':
						vcfSeq[position] = altCall
						#fastaOut.write(vcfSeq[position]) 
					elif altCall == 'G':
						vcfSeq[position] = altCall
						#fastaOut.write(vcfSeq[position]) 
					elif altCall == 'C':
						vcfSeq[position] = altCall
						#fastaOut.write(vcfSeq[position]) 
					#elif altCall == 'CT' or altCall == "TC":
					#	vcfSeq[position] = '-'
					#	fastaOut.write(vcfSeq[position]) 
					#elif altCall == 'AG' or altCall== "GA":
					#	vcfSeq[position] = '-'
					#	fastaOut.write(vcfSeq[position])
					#elif altCall == 'GC' or  altCall == "CG":
					#	vcfSeq[position] = '-'
					#	fastaOut.write(vcfSeq[position])
					#elif altCall == 'AT' or  altCall == "TA":
					#	vcfSeq[position] ='-'
					#	fastaOut.write(vcfSeq[position])
					#elif altCall == 'GT' or altCall== "TG":
					#	vcfSeq[position] ='-'
					#	fastaOut.write(vcfSeq[position])
					#elif altCall == 'AC' or altCall == "CA":
					#	vcfSeq[position] ='-'
					#	fastaOut.write(vcfSeq[position])
				#	elif altCall == 'CGT' or altCall== "CTG" or altCall == "TGC" or altCall == "TCG" or altCall == "GTC" or altCall =="GCT":
				#		vcfSeq[position] = "B"
				#		fastaOut.write(vcfSeq[position]) 
				#	elif altCall == 'AGT' or altCall == "ATG" or  altCall == "GAT" or altCall == "GTA" or altCall== "TAG" or altCall =="TGA":
				#		vcfSeq[position] = "D"
				#		fastaOut.write(vcfSeq[position]) 
				#	elif altCall == 'ACT' or altCall == "ATC" or altCall == "TAC" or altCall== "TCA" or  altCall == "CAT" or altCall == "CTA":
				#		vcfSeq[position] = "H"
				#		fastaOut.write(vcfSeq[position]) 
				#	elif altCall == 'ACG' or altCall == "AGC" or  altCall == "GAC" or altCall == "GCA" or altCall == "CAG" or altCall == "CGA":
				#		vcfSeq[position] = "V"
				#		fastaOut.write(vcfSeq[position]) 
				elif heterozygosity < 0.2: ### few alternative nucleotide reads
					vcfSeq[position] = vcfSplit[i][3]
					#fastaOut.write(vcfSplit[i][3]) ### because there are few alternative reads, just write the reference base call
				else: ### now to write the code for the calls that true heterozygotes
					heteroColumn = vcfSplit[i][4] ### pull out the column with the alt bases and <*> syntax
					bases = heteroColumn.split('<')[0] ### split by < to isolate the individual base calls
					altCall = bases.replace(',','')
					heterozygotes = vcfSplit[i][3] + altCall
					if len(heterozygotes) >= 3:
						vcfSeq[position] = "-"
						#fastaOut.write(vcfSeq[position])
					elif heterozygotes == 'CT' or heterozygotes== "TC":
						vcfSeq[position] = 'Y'
						#fastaOut.write(vcfSeq[position]) 
				#	elif heterozygotes == 'CGT' or heterozygotes=="CTG" or heterozygotes=="TGC" or heterozygotes=="TCG" or heterozygotes=="GTC" or heterozygotes=="GCT":
				#		vcfSeq[position] = "B"
				#		fastaOut.write(vcfSeq[position]) 
				#	elif heterozygotes == 'AGT' or heterozygotes=="ATG" or heterozygotes=="GAT" or heterozygotes=="GTA" or heterozygotes=="TAG" or heterozygotes=="TGA":
				#		vcfSeq[position] = "D"
				#		fastaOut.write(vcfSeq[position]) 
				#	elif heterozygotes == 'ACT' or heterozygotes=="ATC" or heterozygotes=="TAC" or heterozygotes=="TCA" or heterozygotes=="CAT" or heterozygotes=="CTA":
				#		vcfSeq[position] = "H"
				#		fastaOut.write(vcfSeq[position]) 
				#	elif heterozygotes == 'ACG' or heterozygotes=="AGC" or heterozygotes== "GAC" or heterozygotes=="GCA" or heterozygotes=="CAG" or heterozygotes=="CGA":
				#		vcfSeq[position] = "V"
				#		fastaOut.write(vcfSeq[position]) 
					elif heterozygotes == 'AG' or heterozygotes=="GA":
						vcfSeq[position] = 'R'
						#fastaOut.write(vcfSeq[position])
					elif heterozygotes == 'GC' or heterozygotes=="CG":
						vcfSeq[position] = 'S'
						#fastaOut.write(vcfSeq[position])
					elif heterozygotes == 'AT' or heterozygotes=="TA":
						vcfSeq[position] ='W'
						#fastaOut.write(vcfSeq[position])
					elif heterozygotes == 'GT' or heterozygotes== "TG":
						vcfSeq[position] ='K'
						#fastaOut.write(vcfSeq[position])
					elif heterozygotes == 'AC' or heterozygotes=="CA":
						vcfSeq[position] ='M'
						#fastaOut.write(vcfSeq[position])
		else: ### if depth is less than 10
			vcfSeq[position] = 'N'
			#fastaOut.write(vcfSeq[position])	

fastaOut.write(">" + fileHandle +"\n" + "".join(vcfSeq) + "\n") ### format for fasta. Need to join each element of the list
temp = "".join(vcfSeq)
print(fileHandle + "\n" + str((temp.count('N')/genomeLength)*100)) #print out the proportion of bases that are Ns, just to see depth of mapping
fastaOut.close()