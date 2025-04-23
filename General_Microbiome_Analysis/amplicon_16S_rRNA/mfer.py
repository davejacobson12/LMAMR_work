#!/usr/local/bin/python3

'''This script formats and transposes OTU tables, merges them with metadata files, and performs categorical analyses. Usage: mfer.py -t otu_table -m metadata_file -n categorical_variable'''
__author__ = "Richard Hagan, David Jacobson, Allison Mann, and Krithivasan Sankaranarayanan"
__license__ = "GPL"
__version__ = "1.0"

import pandas as pd
from pandas import DataFrame
import argparse
import subprocess
import glob
import os
import shutil

parser = argparse.ArgumentParser(description='This script formats and transposes OTU tables, merges them with metadata files, and performs categorical analyses. Usage: mfer.py -t otu_table -m metadata_file -n categorical_variable')

#Required arguments
parser.add_argument('-t', '--table', help="The table you wish to format.", required=True)
parser.add_argument('-m', '--metadata', help="The metadata file you wish to merge", required=True) 

#Optional arguments
parser.add_argument('-n', '--categoryName', help='Metadata category used to compare taxa significance. Default is all categories')
parser.add_argument('-c', '--collapsedOut', help='Desired collapsed OTU table name, default is collapsed_otu_table.txt', default='collapsed_otu_table.txt')
parser.add_argument('-o', '--compareOut', help='Output file name for cateogry compare step, default is compareOut.txt', default='compareOut.txt')
parser.add_argument('-v', '--mergedOut', help='Output file name for merged table step, default is mergedOut.txt', default='mergedOut.txt')
parser.add_argument('-r', '--transpose', help='Output file name for transpose table step, default is tranposedTable.txt', default='transposedTable.txt')
parser.add_argument('--fdr', help='Optional. Choose fdr value for significance, default is 0.1',default=0.1)
parser.add_argument('--pval', help='Optional. Choose p-value for significant, default is 0.05',default=0.05)
parser.add_argument('-s', '--significanceOutput', help = 'the output from the R script, filtered for significant pval and fdr' , default = 'significanceOutput.txt')
parser.add_argument('--taxSum', help = 'Optional. If you wish to filter out any taxa which have a total number of reads less than a number, please use this flag and enter your chosen integer. Default is 0', default = 0)
parser.add_argument('--taxMax', help = 'Optional. Similar to taxSum. If you wish to filter out any taxa who have who have a maximum count of less than a defined number, please use this flag and enter your chosen integer. Default is 0', default = 0)
parser.add_argument('--summedOut', help = 'Do not use this flag, it simply produces an intermediate file with the total number of reads for each taxa', default = 'summed.txt')
parser.add_argument('--maxOut', help = 'Do not use this flag, it simply produces an intermediate file with the max count for each taxa', default = 'maxOut.txt')
parser.add_argument('--intermediate', help = 'The name of the directory in which to store the intermediate files, default is "intermediateFiles"', default = 'intermediateFiles')
args = parser.parse_args()

class colors:
	RUN = '\033[94m'
	COMPLETE = '\033[92m'
	ENDC = '\033[0m'
	FAIL = '\033[91m'

def prelimCheck():
	assert os.path.exists(args.table), colors.FAIL + 'ERROR: File does not exist: %s. Is the path correct?' % args.table
	assert os.path.exists(args.metadata), colors.FAIL + 'ERROR: File does not exist: %s. Is the path correct?' % args.metadata
	if os.path.isfile('Compare.R'):
		print('Running...')
	else:
		print(colors.FAIL + 'ERROR: R script not found, is it in your current working directory?' + colors.ENDC)


def tableCheck():
	global cleanedTable
	print(colors.RUN + 'Checking OTU table format...' + colors.ENDC)
	with open(args.table, 'r') as f:
		first = f.readline() #pull out only first line of file
		if '# ' in first: #check to see if leading header is there
			inTable = pd.read_csv(args.table, sep='\t', skiprows=1) #if so, skip
		else:
			inTable = pd.read_csv(args.table, sep='\t') #else read in normally

	inTable.replace(to_replace=' *', value='', regex=True, inplace=True) #remove spaces from taxa strings
	inTable.columns = [x.strip().replace(' ', '_') for x in inTable.columns] #replace any spaces from column names and change to underscores

	cleanedTable = inTable.groupby('#OTU_ID').sum() #group by taxonomy string, get sum across all rows
	print(colors.COMPLETE + "Number of duplicate taxa strings to be collapsed: %i" % (inTable.shape[0] - cleanedTable.shape[0]+1) + colors.ENDC) #number of taxa collapsed, added one to account for differences in shape	
	
	with open(args.collapsedOut, 'w') as outfile:
		cleanedTable.to_csv(outfile, sep='\t')
	outfile.close()

	print(colors.RUN + "Table checking complete. Collapsed OTU table written to: %s..." % args.collapsedOut + colors.ENDC)

def addTaxa():
	userSum = int(args.taxSum)
	df = cleanedTable 
	df['summedUp'] = df.sum(axis=1)
	summedDF = df[df.summedUp> userSum ]
	summedDF.to_csv(args.summedOut, sep="\t")

def maxTaxa():
	userMax = int(args.taxMax)
	df = pd.read_csv(args.summedOut, sep = "\t")
	df = df.drop('summedUp', 1)
	df['taxaMax'] = df.max(axis=1)
	maxedDF = df[df.taxaMax > userMax]
	maxedDF.to_csv(args.maxOut, sep ="\t", index = False)



def transpose(): #the input to this function should be Allie's cleaned table
	global transposedTable
	table = pd.read_csv(args.maxOut, sep ="\t")
	table = table.drop('taxaMax', 1)
	table.set_index('#OTU_ID', inplace= True)
	transposedTable = table.transpose()
	transposedTable.to_csv(args.transpose, sep='\t')
	print(colors.RUN + "Table transposed. Written to %s..." % args.transpose + colors.ENDC)


def merging():
	global metadataHeaders
	print(colors.RUN + 'Merging %s and %s...' % (args.metadata, args.transpose) + colors.ENDC)
	with open(args.metadata, 'r') as m:
		metadataHeaders = list(m.readline().split())[:]
	x = len(metadataHeaders)
	metaDF = pd.read_csv(args.metadata, dtype = {metadataHeaders[0]: object} , sep='\t')
	taxaDF = pd.read_csv(args.transpose, sep = '\t')
	new_columns = taxaDF.columns.values; new_columns[0] = 'SampID'
	taxaDF.columns = new_columns
	mergedDF = pd.merge(metaDF, taxaDF, left_on = metadataHeaders[0], right_on = 'SampID', how = 'inner') ## merge the metadata and transposed files using the first header in each file as the criteria
	### should we do merge inner or outer? 
	mergedDF= mergedDF.drop(mergedDF.columns[x], axis =1) ### delete the column that has sample names from taxa table
	mergedDF.columns = [x.strip().replace('#', '') for x in mergedDF.columns] #remove hash symbol	
	mergedDF.to_csv(args.mergedOut, sep = '\t', index = False)


	
def Rcomp():
	if args.categoryName:
		print(colors.RUN + 'Comparing groups, split by %s..' % args.categoryName + colors.ENDC)
		subprocess.call(args=['Rscript Compare.R %s %s %s' %(args.mergedOut, args.categoryName, args.compareOut)], shell=True)
		print(colors.COMPLETE + 'Complete! Output written to %s' % args.compareOut + colors.ENDC)
	else:
		for category in metadataHeaders:
			print(colors.RUN + 'Comparing groups, split by %s..' % category + colors.ENDC)
			subprocess.call(args=['Rscript Compare.R %s %s %s' %(args.mergedOut, category, category + "Compared.txt")], shell=True)
			print(colors.COMPLETE + 'Complete! Output written to %s' % (category + "Compared.txt") + colors.ENDC)	
	

def filterR():
	if args.categoryName:
		df =  pd.read_csv(args.compareOut, sep = '\t')
		userPval = float(args.pval)
		userFDR = float(args.fdr)
		newDF = df.ix[(df['pval']<=userPval) & (df['fdr']<=userFDR)]
		newDF.to_csv(args.categoryName + "_" + args.significanceOutput, sep = '\t', index = False)
	else:
		fileList = glob.glob('*Compared.txt')
		for file in fileList:
			print(file)
			df = pd.read_csv(file, sep = "\t")
			userPval = float(args.pval)
			userFDR = float(args.fdr)
			newDF = df.ix[(df['pval']<=userPval) & (df['fdr']<=userFDR)]
			newDF.to_csv("significance_" + file , sep = "\t", index = False)	

def directoryMaker():
	workDir = os.getcwd()
	intermidatePath = workDir + '/' + args.intermediate
	outputFiles = [args.transpose, args.collapsedOut, args.mergedOut, args.summedOut, args.maxOut]
	while True:
		if os.path.exists(intermidatePath) == True:
			print('You currently have a directory named ' + args.intermediate + '.' )
			print('Due to this, your intermediate files are currently in your working directory ' + '(' + workDir + ').')
			print('Please choose another name for the intermeidate files directory to ensure that no files are overwritten')
			break
		else:
			os.makedirs(workDir + "/" + args.intermediate)
			break
	for files in outputFiles:
		if os.path.exists(intermidatePath + "/" + str(outputFiles)) == True:
			break
		else:
			shutil.move(workDir + "/" + files, intermidatePath + "/" + files)
		print(colors.COMPLETE + 'Intermediate files stored in '+ args.intermediate)


	
	


prelimCheck()
tableCheck()
addTaxa()
maxTaxa()
transpose()
merging()
Rcomp()
filterR()
#directoryMaker()



