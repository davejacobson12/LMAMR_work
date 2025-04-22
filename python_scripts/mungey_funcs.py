import pandas as pd
from pandas import DataFrame
import argparse
parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(description='Mungy description here')

parser.add_argument('-t', '--table', help="The table you wish to format", required=True, metavar="") # Goes to Allie's function
parser.add_argument('-m', '--metadata', help="The metadata file you wish to merge", required=True) # Goes to Dave's functionparser.add_argument('-x', '--transpose', help="Use this if you want to transpose your table") # Richard's function y/n
parser.add_argument('-c', '--collapsedOut', help='Desired collapsed OTU table name, default is collapsed_otu_table.txt', metavar='', default='collapsed_otu_table.txt')

args = parser.parse_args()

class colors:
	COMPLETE = '\033[92m'

def tableCheck():
	global cleanedTable
	print('Checking OTU table format...')
	with open(args.table, 'r') as f:
		first = f.readline() #pull out only first line of file
		if '# ' in first: #check to see if leading header is there
			inTable = pd.read_csv(args.table, sep='\t', skiprows=1) #if so, skip
		else:
			inTable = pd.read_csv(args.table, sep='\t') #else read in normally

	inTable.replace(to_replace=' *', value='', regex=True, inplace=True) #remove spaces from taxa strings
	inTable.columns = [x.strip().replace(' ', '_') for x in inTable.columns] #replace any spaces from column names and change to underscores

	cleanedTable = inTable.groupby('#OTU_ID').sum() #group by taxonomy string, get sum across all rows
	print("Number of duplicate taxa strings to be collapsed: %i" % (inTable.shape[0] - cleanedTable.shape[0]-1))	
	
	with open(args.collapsedOut, 'w') as outfile:
		cleanedTable.to_csv(outfile, sep='\t')
	outfile.close()

	print(colors.COMPLETE + "Table checking complete. Collapsed OTU table written to: %s " % args.collapsedOut)

def transpose(cleanedTable): #the input to this function should be Allie's cleaned table
	global transposedTable
	transposedTable = cleanedTable.transpose()
	transposedTable.to_csv('transposedTable.tsv', sep='\t',)

tableCheck()

transpose(cleanedTable)

def merging():
	metaDF = pd.read_csv(args.metadata, sep='\t')
	taxaDF = pd.read_csv('transposedTable.tsv', sep = '\t')
	mergedDF = pd.merge(metaDF, taxaDF, left_index = True, right_index = True, how = 'outer') ## merge the metadata and transposed files using the first header in each file as the criteria
	### should we do merge inner or outer? 
	del mergedDF['Unnamed: 0'] ### delete the column that has sample names from taxa table
	mergedDF.to_csv('mergedTable.txt', sep = '\t', index = False)

merging()


