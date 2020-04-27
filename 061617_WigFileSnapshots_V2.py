#############################################################
######### 061617_WigFileSnapshots_V2.py                #########
######### Have a wig file in four columns and       #########
######### then pull out wig signal for a given range#########
#############################################################

#!/usr/bin/env python


import os
from math import *
from optparse import OptionParser
import sys
import re

################# Argument Import : Uncomment for actual running script  #################

def parse_options():
	parser = OptionParser()
	parser.add_option("-W", "--wigfile", dest="wigfile",
					  help="4 columned TXT file of normalized wig signal: chr, start, end, wigvalue")
	parser.add_option("-R", "--regfile", dest="regfile",
					  help="5 columned TXT file of wanted genomic regions: Number, chr, start, end, Value")
	parser.add_option("-l", "--length", dest="length",
					  help="total genomic length returned")
	(options, args) = parser.parse_args()
	return options
            
options = parse_options()
parser = OptionParser()
if not options.wigfile:
    print("Wigfile option is missing\n")
    parser.print_help()
    exit(-1)
if not options.regfile:
    print("Region file option is missing\n")
    parser.print_help()
    exit(-1)
if not options.length:
    print("Length option is missing\n")
    parser.print_help()
    exit(-1)

################# Data Import   #################
# for jupyter Notebook pass these arguments
# wigfile = '061417_multiBigWigCompare_bins_041417_Bowtie2_ATACPools_merged.txt'
# regfile = 'MergedALL2.txt'
# length = 2000

####### Read in 4 columned Wig file
#'chr'	'start'	'end'	'041417_Bowtie2_ATACPools_merged_10Whole_10MDmel.bw'	'041417_Bowtie2_ATACPools_merged_20Anterior_10MDmel.bw'	'041417_Bowtie2_ATACPools_merged_20Posterior_10MDmel.bw'	'041417_Bowtie2_ATACPools_merged_AplusP_10MDmel.bw'	'101416_10MNorm_DNase-I_dmels5r1.bw'
#chr2R	0	10	nan	nan	nan	nan	nan
#chr2R	10	20	nan	nan	nan	nan	nan

WigFile_file = open(str(options.wigfile), 'r')
WigList = WigFile_file.readlines()
WigFile_file.close()
print('Wig File Imported')
                    
####### Read in File with regions I want to get
                    
#"Sample"	"chr"	"start"	"end"	"Merged_Whole"	"Merged_20Ant"	"Merged_20Post"	"Merged_AplusP"	"DnaseI"	"Type"	"Location"	"Name"	"Dotsize"
#"1"	"chr2R"	17598962	17599122	"118.7752312"	"135.1978857"	"71.82917128"	"164.9578282"	"67.41856544"	"Promoter"	"Posteriorly Primarily"	"CG30403"	1.5
#"2"	"chr3R"	7702946	7702963	"718.3436746"	"436.3480655"	"775.9198213"	"575.2080976"	"696.9214442"	"Promoter"	"OnlyA"	"sad"	1.5
#"3"	"chr2R"	8146899	8146969	"116.2616414"	"98.17407118"	"160.4658203"	"104.0412238"	"130.1133728"	"Promoter"	"Mostly A"	"Cam"	1.5

RegionsFile_file = open(str(options.regfile), 'r')
RegionsList = RegionsFile_file.readlines()
RegionsFile_file.close()
print('Regions File Imported')

WigFirstLine = WigList.pop(0)
RegionFirstLine = RegionsList.pop(0)

####### Create a list data frame where each list is a list containing information for one enhancer location  #################
RegionsListElement = []

# For each enhancer, strip the line ending, and split up line by \t and add that new list to a list of lists
WigDict = {}
for line in WigList:
    line = line.strip('\n')
    #WigListElement.append(line.split('\t'))
    wigregion = line.split('\t')
    #print(wigregion)
    wigChrom = wigregion[0]
    wigStart = int(wigregion[1])
    wigEnd = int(wigregion[2])
    WigKey = '_'.join([str(wigChrom),str(wigStart),str(wigEnd)])
    WigDict[WigKey] = line
#print(WigDict)
print("wig dictionary made")
print(WigDict["chrX_9583490_9583500"])

for line in RegionsList:
    line = line.strip('\n')
    RegionsListElement.append(line.split('\t'))
print(RegionsListElement[0:10])
print("Files made into lists")

################# Make a dictionary of all the keys between each region 

# Then for each region make a key and then compare it to the wigDict. If it's there, print the key value to a file. 
for region in RegionsListElement:
    #make variables 
    chrom = region[0]
    Start = int(region[1])
    End = int(region[2])
    TotalLength = End - Start
    #adjust the start and end to reflect the user input length
    if TotalLength < float(options.length):
        AddedLength = (float(options.length) - TotalLength)/2
        ExtendedStart = Start - AddedLength
        ExtendedEnd = End + AddedLength
    else:
        ExtendedStart = Start
        ExtendedEnd = End
    #round start and end so they match wig
    RoundedStart = int(int(ExtendedStart/10) * 10)
    #print(Start)
    #print(RoundedStart)
    RoundedEnd = int(int(ExtendedEnd/10) * 10) + 10
#     print(RoundedEnd)
#     print(chrom)
#     print(RoundedStart)
#     print(RoundedEnd)
    #make output files
    OutputFileNameString = str(region[5]) + "_WigFileSnapshot.txt"
    print(OutputFileNameString)
    outputFile = open(OutputFileNameString, 'w')
    
    #make a key for each 10bp wig interval in the extended region and compare that to wig dict
    while RoundedStart < RoundedEnd: # if start plus 10 is less than the end,
        RoundedStart += 10
        key = '_'.join([str(chrom),str(RoundedStart - 10),str(RoundedStart)]) #add this key
        outputFile.write(str(WigDict[key]))
        outputFile.write('\n')
    outputFile.close()
print("All done! :)")