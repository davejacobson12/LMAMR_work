import sys

vcfFile = open(sys.argv[1], "r")

fileHandle = sys.argv[1].split('.')[0] ### get the name of the file so you can have the output in fasta format at the end

vcfLines = vcfFile.readlines()

vcfFile.close()

fastaOut = open(sys.argv[3], 'w') ### name of output file

refFile = open(sys.argv[2],'r') ### the reference genome must be in single linge fasta and have just the genome with no plasmid

refGenome = refFile.readlines() 

refHandle = refGenome[0].split()[0].split(">")[1] ### get the genome name of the reference, need to split the first line a couple of times

refFile.close()

genomeLength = len(refGenome[1]) ### the second index is the genome, so get the length

vcfSplit = [None] * genomeLength ### create empty list. I'm not sure that this length needs to be this long
vcfSeq = ['N'] * genomeLength ### populate the final sequence list with Ns at each position for the length of the reference genome. The length will change depending on which genome you mapped to

for i in range(1,len(vcfLines)):
	vcfSplit[i] = vcfLines[i].split() #### populate the empty list with strings seperated by space, so each column can be indexed
	if vcfSplit[i][0] == refHandle: ### make sure that I'm matching just the succinifaciens genome
		position = int(vcfSplit[i][1])
		infoColumn = vcfSplit[i][7].split(';') ## create object that holds the column info with depth and forward/reverse info
		dpString = infoColumn[0].split('=') ## depth of coverage is the first index within that object and the numerical values are after the =, so split by =
		depth = int(dpString[1])
		if depth >=5: ### if depth of coverage is greater than 10, check to see if there is an alternate base
			altLen = len(vcfSplit[i][4]) ## check to see if there is an alternate nucleotide
			if altLen == 3:  ### if the alternate length column is equal to 3, then the base call is the same as the reference
				vcfSeq[position] = vcfSplit[i][3]
			elif altLen >= 6: ### to get rid of any amibigous alternative base calls
				vcfSeq[position] = "-"
			else: ### if the alternative length is not three or greater than 6. This means there is one, non-ambiguous alternative base call
				I16 = infoColumn[1].split('=')[1] ### pull out the string of the I16 column
				ref_alt = I16.split(',')[0:4] ### this string is seperated by commas, so split. We are concerned with the first 4 indices, 1 = reference forward, 2 = rerference reverse, 3 = alt forward, 4 = alt reverse
				heterozygosity= (int(ref_alt[2]) + int(ref_alt[3]))/(int(ref_alt[0]) + int(ref_alt[1])+ int(ref_alt[2]) + int(ref_alt[3])) ### find the proportion of reads that are from the alternate
				#print(int(heterozygosity))
				if heterozygosity >=0.8: ### for calls that are definitely the alternative nucleotide
					heteroColumn = vcfSplit[i][4] ### pull out the column with the alt bases and <*> syntax
					bases = heteroColumn.split('<')[0] ### split by < to isolate the individual base calls
					altCall = bases.replace(',','') ### because there can be multiple alternative calls, but we only want sites with one base call, so the sites with multiple alts I will put a hyphen
					if len(altCall) >= 2: ### If there are any sites with multiple reference alleles, write a hyphen
						vcfSeq[position] = "-"
					elif altCall == 'A' or altCall == 'T' or altCall == 'G' or altCall == 'C': ### if the alternative call is a single nucleotide, write that nucleotide
						vcfSeq[position] = altCall
				elif heterozygosity < 0.2: ### few alternative nucleotide reads
					vcfSeq[position] = vcfSplit[i][3]### because there are few alternative reads, just write the reference base call
				else: ### now to write the code for the calls that true heterozygotes
					heteroColumn = vcfSplit[i][4] ### pull out the column with the alt bases and <*> syntax
					bases = heteroColumn.split('<')[0] ### split by < to isolate the individual base calls
					altCall = bases.replace(',','')
					heterozygotes = vcfSplit[i][3] + altCall
					if len(heterozygotes) >= 3: ### if there are multiple alternative alleles, write a hyphen
						vcfSeq[position] = "-"
					elif heterozygotes == 'CT' or heterozygotes== "TC": ### Go through all of the combinations of heterozygotes and write the IUPAC code for each combination
						vcfSeq[position] = 'Y' 
					elif heterozygotes == 'AG' or heterozygotes=="GA":
						vcfSeq[position] = 'R'
					elif heterozygotes == 'GC' or heterozygotes=="CG":
						vcfSeq[position] = 'S'
					elif heterozygotes == 'AT' or heterozygotes=="TA":
						vcfSeq[position] ='W'
					elif heterozygotes == 'GT' or heterozygotes== "TG":
						vcfSeq[position] ='K'
					elif heterozygotes == 'AC' or heterozygotes=="CA":
						vcfSeq[position] ='M'
		else: ### if depth is less than 5, write N
			vcfSeq[position] = 'N'	

fastaOut.write(">" + fileHandle +"\n" + "".join(vcfSeq) + "\n") ### format for fasta. Need to join each element of the list
temp = "".join(vcfSeq)
print(fileHandle + "\n" + str((temp.count('N')/genomeLength)*100)) #print out the proportion of bases that are Ns, just to see depth of mapping
fastaOut.close()
