originalFile  = open('/Users/dave/Desktop/nisha_lmamr15_sheet.txt')  ###path to file with Barcodes that are 12bp

lines = originalFile.readlines()  ### read the lines of that file

header = lines[0].split()  ### pull out the first row, which should be the headers with one of the columns labeled Barcode

num_lines = len(lines)  ### calculate the number of lines in the file


barcode_index = header.index('Barcode')  ### find the index of the column with Barcode as the header

revComp = [0] * 12  ### make an empty array that is 12 characters long, because we are using 12bp barcodes

barcodeFile = open("/Users/dave/Desktop/nisha_lmamr15_revComp.txt", 'w')  ### create a file on the desktop, can change path to any directory you want
for i in range(1,num_lines):  ### for each line, starting with the first non-header row
	each_line = lines[i].split()  ###  split each line by tabs into character/words
	barcodes = each_line[barcode_index] ### using the barcode index, find the barcode in each line
	revSeq = barcodes[::-1] ### reverse the sequence of the barcode
	start = 0 ### initialize a loop
	while start <= 12:  ### run the loop 12 times (for 12bp barcode)
		start = start + 1  ### add one to start each loop, so it stops after 12
		for i, j in enumerate(revSeq):  ### basically print each letter of the reverse sequence (until the 12th letter)
			if j == "G":
				revComp[i] = "C"
			if j == "C":
				revComp[i] = "G"
			if j == "A":
				revComp[i] = "T"
			if j == "T":
				revComp[i] = "A"
	finalRevCompSeq= ''.join(revComp) ### join the sequences together, so they are not separated by commas
	barcodeFile.write(finalRevCompSeq + '\n') ### write to the new file, with a new line after each barcode


barcodeFile.close() ### close the file
	
	


