#Author: Thulek@gmail.com
#Join the nex files and treefiles

from os import listdir
from os.path import isfile, join
import os
import sys
import platform
import subprocess
import time
import array as arr

import argparse

text = "joinNex"

parser = argparse.ArgumentParser(description = text)  

jobIds = arr.array('i')

parser.add_argument("-d", "--directory", help="The Phylip alignment folder")
parser.add_argument("-t", "--thread", help="The number of threads")

args = parser.parse_args()

def writeFileFrom(fromFile, toFile, startText, endText):
	with open(fromFile) as infile, open(toFile, 'a+') as outfile:
		copy = False
		for line in infile:
			if startText in line.strip():
				#outfile.write(line) # add this
				copy = True
			elif endText in line.strip():
				#outfile.write(line) # add this
				copy = False
			elif copy:
				outfile.write(line)

####GET PARAMETERS AND CREATE TEMP FOLDERS FOR IQTREE#####
alnFolder = args.directory
nThread = int(args.thread)
##########CREATE TREE FILE#################
command = "cat "
for x in range(1,nThread+1):
	#if(x==nThread):
	#	print(str(x))
	command += ""+alnFolder+"."+str(x)+".treefile "
command += " > "+alnFolder+".trees"
os.system(command)

print(command)

##########CREATE NEX FILE###################
#fullPart 
parFile = "tempPart"
bestMFile = "tempBestMF"
open(parFile,"w+")
open(bestMFile,"w+")
if(os.path.isfile(parFile)):
	os.system("rm "+parFile)
if(os.path.isfile(bestMFile)):
	os.system("rm "+bestMFile)
for x in range(1,nThread+1):
	writeFileFrom(alnFolder+"."+str(x)+".best_model.nex", parFile, "begin sets", "charpartition mymodels")
	writeFileFrom(alnFolder+"."+str(x)+".best_model.nex", bestMFile, "charpartition mymodels", "end;")
	cmd = "sed -i 's/"+alnFolder+"."+str(x)+"/"+alnFolder+"/g' "+parFile
	os.system(cmd)
	if(x < nThread):
		cmd = "sed -i 's/;/,/g' "+bestMFile
		os.system(cmd)

nexFile=open(alnFolder+".nex","w+");
nexFile.write("#nexus\n")
nexFile.write("begin sets;\n")
parF = open(parFile,"r")
nexFile.write(parF.read())
parF.close()
nexFile.write("  charpartition mymodels =\n")
parF = open(bestMFile,"r")
nexFile.write(parF.read())
parF.close()
nexFile.write("end;")
nexFile.close()


#print("There are " + str(i) + " aligments and " + str(nFolder) + " threads." + str(nFilesPerJob))