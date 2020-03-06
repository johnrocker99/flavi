#Author: Thulek@gmail.com
#Estimating the model using Iqtree
#Input: a set of alignment
#Output: a substitution model that represent the evolutionary process of input
 
from os import listdir
from os.path import isfile, join
import os
import sys
import platform
import subprocess
import time
import array as arr
import numpy

import argparse

text = "===================\nEstimate the model by IQTREE\n===================="

parser = argparse.ArgumentParser(description = text)  
#parser.parse_args()  

jobIds = arr.array('i')

parser.add_argument("-d", "--directory", help="The Phylip aligments folder (* required)")
parser.add_argument("-t", "--thread", help="The number of threads (sub-jobs) ")
parser.add_argument("-i", "--iter", help="iter ")
parser.add_argument("-m", "--model", help="The model that be used to train (optional, default: GTR20)")
parser.add_argument("-mset", "--mset", help="Set of models, default: FLU,HIVb,HIVw,JTT,WAG,LG")
parser.add_argument("-th", "--threshold", help = "Correlation threshold, range: 0 - 1")
# read arguments from the command line
args = parser.parse_args()


if args.iter:
	iter = int(args.iter)
else:
	iter = 3

maxJob = 4

iqtree_path = "$IQTREE/" ##### enter path to IQTREE folder here #####
mset = "FLU,HIVb,HIVw,JTT,WAG,LG"
if args.mset:
	mset = args.mset
THRESHOLD = 0.99
if args.threshold:
	if (args.threshold>0) and (args.threshold<=1):
		THRESHOLD = args.threshold

#########CLONE A FILE FROM ANOTHER FILE WITH A RANGE##############

def writeFileFrom(fromFile, toFile, startText, endText):
	with open(fromFile) as infile, open(toFile, 'w+') as outfile:
		copy = False
		for line in infile:
			if startText in line.strip():
#				outfile.write(line) # add this
				copy = True
			elif endText in line.strip():
#				outfile.write(line) # add this
				copy = False
			elif copy:
				outfile.write(line)

def correlMatrix(matrix1, matrix2):
	list1 = []
	list2 = []
	rmatrix = open(matrix1,"r")
	for l in rmatrix:
		if len(l.strip()) > 0:
			elet = l.strip().split(" ")
			for ele in elet:
				list1.append(float(ele))
	rmatrix.close()

	cmatrix = open(matrix2,"r")
	for l in cmatrix:
		if len(l.strip()) > 0:
			elet = l.strip().split(" ")
			for ele in elet:
				list2.append(ele)
	cmatrix.close()

	return numpy.corrcoef(list1, list2)[0, 1]

####GET PARAMETERS AND CREATE TEMP FOLDERS FOR IQTREE#####
alnFolder = args.directory
nThread = 1
if args.thread:
	nThread = int(args.thread)

files = [f for f in listdir(alnFolder) if isfile(join(alnFolder, f))]
nFilesPerJob = -(-len(files) // nThread)

########CREATE FOLDER FUNCTION############
def createFolder(folder):
	if not os.path.isdir(folder):
		cmd = 'mkdir '+folder
		os.system(cmd)
	else:
		print(folder+" is exists")


##########LIST FILES IN FOLDER ORDER BY SIZE##################
def get_files_by_file_size(dirname, reverse=False):
	""" Return list of file paths in directory sorted by file size """
	# Get list of files
	filepaths = []
	for basename in os.listdir(dirname):
		filename = os.path.join(dirname, basename)
		if os.path.isfile(filename):
			filepaths.append(filename)
	# Re-populate list with filename, size tuples
	for i in xrange(len(filepaths)):
		filepaths[i] = (filepaths[i], os.path.getsize(filepaths[i]))

	# Sort list by file size
	# If reverse=True sort from largest to smallest
	# If reverse=False sort from smallest to largest
	filepaths.sort(key=lambda filename: filename[1], reverse=reverse)
	# Re-populate list with just filenames
	for i in xrange(len(filepaths)):
		filepaths[i] = filepaths[i][0]
	return filepaths			

#Model-joint
model = ""
if args.model:
	model = str(args.model)
else:
	model = "GTR20"

NEWREV = ""

numJob = 0
it = 1

#In case of using several jobs, folders are created to handle alignments for each job
if nThread > 1:
	for x in range(1, nThread+1):
		if not os.path.isdir(alnFolder+"."+str(x)):
			cmd = 'mkdir ' + alnFolder+"."+str(x)
			os.system(cmd)
		else:
			cmd = 'rm -r ' + alnFolder+"."+str(x)
			os.system(cmd)
			cmd = 'mkdir ' + alnFolder+"."+str(x)
			os.system(cmd)

	createFolder("logs")
	i=0
	y=1
	z=1

	##########Divide The ALIGMENTS to Temp Folders##################
	for f in get_files_by_file_size(alnFolder):
		z = (i % nThread) + 1
		cmd = "cp " + f + " " + alnFolder+"."+str(z)+"/"
		os.system(cmd)
		i+=1
####Divided alignments into folders######

#Get folder's name
alnName = ""
if "/" in alnFolder:
	alnTmpName = alnFolder.split("/")
	alnName = alnTmpName[len(alnTmpName)-1].strip()
else:
	alnName = alnFolder
		
###Loop - iter times
while (it <= iter):
	####MULTIPLE SUB-JOBS	
	if nThread > 1: #if split multiple sub-jobs
		for x in range(1,nThread+1):
			if it == 1:
				command = iqtree_path+"iqtree2 --seed 1 -S "+alnFolder+"."+str(x)+" --no-seq-comp -mset "+mset+"  -cmax 16 -nt 36 -safe -pre "+alnFolder+"."+str(it)+"."+str(x)+""
			else:
				command = iqtree_path+"iqtree2 --seed 1 -S "+alnFolder+"."+str(x)+" -t "+alnFolder+"."+str(it-1)+"."+str(x)+".treefile --no-seq-comp -mset "+mset+"  -cmax 16 -nt 36 -safe -pre "+alnFolder+"."+str(it)+"."+str(x)+""
			os.system(command)
			
		cmd = "python joinNex_os.py -d "+alnFolder+"."+str(it)+" -t "+str(nThread) #MERGE TREEFILE & NEXFILE
		os.system(cmd)
		print("Finish Step 1")

		##################ESTIMATE MODEL##################################################
		if(os.path.isfile(alnFolder+"."+str(it)+".trees") and os.path.isfile(alnFolder+"."+str(it)+".nex")):
			if it > 1:
				init_model="--init-model "+NEWREV
			else:
				init_model = ""
			command = iqtree_path+"iqtree2 --seed 1 -S "+alnFolder+"."+str(it)+".nex -te "+alnFolder+"."+str(it)+".trees --no-seq-comp -nt 36 -safe --model-joint "+model+" "+init_model+" -pre "+alnFolder+"."+str(it)+"."+model
			os.system(command)
				
		#read new paml
		if(os.path.isfile(alnFolder+"."+str(it)+"."+model+".iqtree")):
			writeFileFrom(alnFolder+"."+str(it)+"."+model+".iqtree", alnFolder+".M"+str(it)+".PAML", "Substitution parameters (lower-diagonal) and state frequencies in PAML format (can be used as input for IQ-TREE)", "State frequencies: (empirical counts from alignment)")
			if(os.path.isfile(alnFolder+".M"+str(it)+".PAML")):
				NEWREV = alnFolder+".M"+str(it)+".PAML"
				print("NEWREV: "+alnFolder+".M"+str(it)+".PAML")
				if it >= 2:
					correlValue = float(correlMatrix(alnFolder+".M"+str(it-1)+".PAML",alnFolder+".M"+str(it)+".PAML"))
					if correlValue > THRESHOLD:
						it = iter + 1 # quit if correlation > THRESHOLD
						print ("Correlation of NEWREV and CURREV is "+ str(correlValue) )
	
	#NO SUB-JOB
	else: #if only one job
		if it == 1:
			command = iqtree_path+"iqtree2 --seed 1 -S "+alnFolder+" --no-seq-comp -mset "+mset+"  -cmax 16 -nt 36 -safe -pre "+alnFolder+"."+str(it)+""
		else:
			command = iqtree_path+"iqtree2 --seed 1 -S "+alnFolder+" -t "+alnFolder+"."+str(it-1)+".treefile --no-seq-comp -mset "+mset+"  -cmax 16 -nt 36 -safe -pre "+alnFolder+"."+str(it)+""
		os.system(command)
		
		print("Finish Step 1")

		##################ESTIMATE MODEL##################################################
		if(os.path.isfile(alnFolder+"."+str(it)+".treefile") and os.path.isfile(alnFolder+"."+str(it)+".best_model.nex")):
			if it > 1:
				init_model="--init-model "+NEWREV
			else:
				init_model = ""
			command = iqtree_path+"iqtree2 --seed 1 -S "+alnFolder+"."+str(it)+".best_model.nex -te "+alnFolder+"."+str(it)+".treefile --no-seq-comp -nt 36 -safe --model-joint "+model+" "+init_model+" -pre "+alnFolder+"."+str(it)+"."+model
			os.system(command)
		
		#Read new paml		
		if(os.path.isfile(alnFolder+"."+str(it)+"."+model+".iqtree")):
			writeFileFrom(alnFolder+"."+str(it)+"."+model+".iqtree", alnFolder+".M"+str(it)+".PAML", "Substitution parameters (lower-diagonal) and state frequencies in PAML format (can be used as input for IQ-TREE)", "State frequencies: (empirical counts from alignment)")
			if(os.path.isfile(alnFolder+".M"+str(it)+".PAML")):
				NEWREV = alnFolder+".M"+str(it)+".PAML"
				print("NEWREV: "+alnFolder+".M"+str(it)+".PAML")
				if it >= 2:
					correlValue = float(correlMatrix(alnFolder+".M"+str(it-1)+".PAML",alnFolder+".M"+str(it)+".PAML"))
					if correlValue > THRESHOLD:
						it = iter + 1 # quit if correlation > THRESHOLD
						print ("Correlation of NEWREV and CURREV is "+ str(correlValue))
	it = it + 1	

###CLEAN TEMPORARY DATA
nJob = 1
if nThread > 1:
	while nJob <= nThread:
		os.system("rm -rf "+alnFolder+"."+str(nJob))
		nJob = nJob + 1

nIt = 1
while nIt <= iter:
	os.system("rm "+alnFolder+"."+str(nIt)+".*")
	nIt += 1

if(os.path.isfile("tempPart")):
	os.system("rm tempPart")

if(os.path.isfile("tempBestMF")):
	os.system("rm tempBestMF")