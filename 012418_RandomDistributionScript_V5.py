####################################################################
######### 012418_RandomDistributionScript_V5.py             #########
######### Take a bed file and find random regions           #########
######### that have same overall wig signal                 #########
######### output : Random Region overlap file, and histograms#########
#####################################################################

#!/usr/bin/env python


import os
from math import *
from optparse import OptionParser
import sys
import re
from pylab import figure, title, xlabel, ylabel, hist, axis, grid, savefig
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from random import shuffle

################# Argument Import : Uncomment for actual running script  #################

def parse_options():
	parser = OptionParser()
	parser.add_option("-W", "--wigfile", dest="wigfile",
					  help="Random region file - 4 columned TXT file (unix line endings) with columns : #'chr'	'start'	'end'	'041417_Bowtie2_ATACPools_merged_20Anterior_10MDmel.bw'	'062817_LinregNorm_041417_Bowtie2_ATACPools_merged_20Posterior_10MDmel.bw'")
	parser.add_option("-R", "--regfile", dest="regfile",
					  help="5 columned TXT file of wanted genomic regions columns: chr	start	end	Type	New Location Assignment	Name	Blind use?	Old Location Assignment	Insitu?	Blind Notes	Length	FlybaseID	Source	#'chr'	'start'	'end'	'041417_Bowtie2_ATACPools_merged_10Whole_10MDmel.bw'	'081217_LinregNormWhole_041417_Bowtie2_ATACPools_merged_20Anterior_10MDmel.bw'	'081217_LinregNormWhole_041417_Bowtie2_ATACPools_merged_20Posterior_10MDmel.bw'	'081217_LinregNormWhole_041417_Bowtie2_ATACPools_merged_AplusP_10MDmel.bw'	'101416_10MNorm_DNase-I_dmels5r1.bw'")
	parser.add_option("-o", "--output", dest="output",
					  help="name of output file")
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
if not options.output:
    print("output file is missing\n")
    parser.print_help()
    exit(-1)
    
################# Data Import   #################
enhancerList_file = open(str(options.regfile), 'r')
enhancerList = enhancerList_file.readlines()
enhancerList_file.close()

#   B. Read in Random region ATAC file (zeros included)
#Header : #'chr'	'start'	'end'	'041417_Bowtie2_ATACPools_merged_20Anterior_10MDmel.bw'	'062817_LinregNorm_041417_Bowtie2_ATACPools_merged_20Posterior_10MDmel.bw'
randomList_file = open(str(options.wigfile), 'r')
randomList = randomList_file.readlines()
randomList_file.close()

#   D. Take out the header of both random and enhancer ATAC files. it is a list
enhFirstLine = enhancerList.pop(0)
enhFirstLine = enhFirstLine.strip("/n")
randFirstLine = randomList.pop(0)
randFirstLine = randFirstLine.strip("/n")

print("files imported")

################# Create Dictionaries for region and random files   #################
enhListElement = []
randomListElement = []
enhScoreDict = {}
enhInfoDict = {}
enhLocDict = {}
enhTypeDict = {}
enhAntDict = {}
enhPostRDict = {}
enhNameList = []
Insitu_use ={}


# For each enhancer, strip the line ending, and split up line by \t 
# make a dictionary for score and dictionary for the line by name

for line in enhancerList:
    line2 = line.strip('\n')
    enhListElement = line2.split('\t')
    enhNameList.append(enhListElement[5])
    EnhancerTotalScore = float(enhListElement[25]) + float(enhListElement[27])
    enhScoreDict[enhListElement[5]] = EnhancerTotalScore
    enhInfoDict[enhListElement[5]] = line2
    enhLocDict[enhListElement[5]] = enhListElement[4]
    enhTypeDict[enhListElement[5]] = enhListElement[3]
    enhAntDict[enhListElement[5]] = enhListElement[25]
    enhPostRDict[enhListElement[5]] = enhListElement[27]
    Insitu_use[enhListElement[5]] = enhListElement[6]


# Calculate the ATACskewScore and add it to the output file
EnhancerATACSkewScoreList =[]
AnteriorATACSkewScoreList = []
PosteriorATACSkewScoreList =[]
PosteriorPromoterATACSkewScoreList =[]
AnteriorPromoterATACSkewScoreList =[]
PromoterATACSkewScoreList =[]
Insitu_use_List =[]
ATACSkewScoreDict = {}



for enhancer in enhNameList:
    if Insitu_use[enhancer] == "yes": #if used in the final list
        if enhTypeDict[enhancer] == "Enhancer":
            if enhLocDict[enhancer] == "Anterior" or enhLocDict[enhancer] == "Mostly Ant" :
                ATACSkewScore = (float(enhAntDict[enhancer]) - float(enhPostRDict[enhancer])) / (float(enhAntDict[enhancer]) + float(enhPostRDict[enhancer]))
                print("Anterior: ", enhancer)
                EnhancerATACSkewScoreList.append(ATACSkewScore) #add it to the general enhancers list
                AnteriorATACSkewScoreList.append(ATACSkewScore) #add it to the anterior enhancers only list
            if enhLocDict[enhancer] == "Posterior" or enhLocDict[enhancer] == "Mostly Post" :
                ATACSkewScore = (float(enhPostRDict[enhancer]) - float(enhAntDict[enhancer])) / (float(enhAntDict[enhancer]) + float(enhPostRDict[enhancer]))
                print("Posterior: ", enhancer)
                EnhancerATACSkewScoreList.append(ATACSkewScore) #add it to the general enhancers list
                PosteriorATACSkewScoreList.append(ATACSkewScore) #add it to the posterior enhancers list
            else: #DV enhancers
                ATACSkewScore = (float(enhAntDict[enhancer]) - float(enhPostRDict[enhancer])) / (float(enhAntDict[enhancer]) + float(enhPostRDict[enhancer]))
                EnhancerATACSkewScoreList.append(ATACSkewScore)
        else: #promoter
            if enhLocDict[enhancer] == "Anterior" or enhLocDict[enhancer] == "Mostly Ant" :
                ATACSkewScore = (float(enhAntDict[enhancer]) - float(enhPostRDict[enhancer])) / (float(enhAntDict[enhancer]) + float(enhPostRDict[enhancer]))
                PromoterATACSkewScoreList.append(ATACSkewScore) #add it to the promoters list
                AnteriorPromoterATACSkewScoreList.append(ATACSkewScore) #add it to the anterior promoters list
            if enhLocDict[enhancer] == "Posterior" or enhLocDict[enhancer] == "Mostly Post" :
                ATACSkewScore = (float(enhPostRDict[enhancer]) - float(enhAntDict[enhancer])) / (float(enhAntDict[enhancer]) + float(enhPostRDict[enhancer]))
                PromoterATACSkewScoreList.append(ATACSkewScore) #add it to the promoters list
                PosteriorPromoterATACSkewScoreList.append(ATACSkewScore) #add it to the posterior promoters list
            else: #DV promoters
                ATACSkewScore = (float(enhAntDict[enhancer]) - float(enhPostRDict[enhancer])) / (float(enhAntDict[enhancer]) + float(enhPostRDict[enhancer]))
                PromoterATACSkewScoreList.append(ATACSkewScore)
        ATACSkewScoreDict[enhancer] = ATACSkewScore
print("ATACSkewScores Calculated!")


###### Need to shuffle the random region list and peform the following standard deviation calculation
#Import Random regions, calculate the Random skew score for each region
RandScoreDict = {}
RandInfoDict = {}
EnhrandATACSkewScoreList = []
RandomATACSkewScoreDict = {}
PromrandATACSkewScoreList =[]
randNamelist =[]
perm_counter = 0

##### Mu and Std Lists
mu_enh_list =[]
mu_prom_list = []
std_enh_list = []
std_prom_list = []

for line in randomList:
    line2 = line.strip('\n')
    randomListElement = line2.split('\t')
    RandomTotalScore = float(randomListElement[3]) + float(randomListElement[4])
    RandName = "-".join([randomListElement[0], randomListElement[1], randomListElement[2]])
    randNamelist.append(RandName)
    RandScoreDict[RandName] = RandomTotalScore
    RandInfoDict[RandName] = line2
    if RandomTotalScore != 0:
        randATACSkewScore = (float(randomListElement[3]) - float(randomListElement[4])) / RandomTotalScore
        RandomATACSkewScoreDict[RandName] = randATACSkewScore
            

RandUniqueDict = {}
enhRandregionMatchDict ={}
MatchedRandomATACSkewScoreDict ={}

shuffle(randNamelist)
for enhancer in enhNameList:
    if Insitu_use[enhancer] == "yes": #if used in the final list
        enhancerCounter = 0
        TotalScore = enhScoreDict[enhancer]
        for randomkey in randNamelist:
            if enhancerCounter < 1:
                if randomkey not in RandUniqueDict:
                    randValue = RandScoreDict[randomkey]
                    if TotalScore * 0.8 < randValue < TotalScore * 1.1:
                        enhancerCounter += 1
                        EnhrandATACSkewScoreList.append(RandomATACSkewScoreDict[randomkey])
                        RandUniqueDict[randomkey] = True
                        enhRandregionMatchDict[enhancer] = RandInfoDict[randomkey]
                        MatchedRandomATACSkewScoreDict[enhancer] = RandomATACSkewScoreDict[randomkey]
                        print("Match!")
                        break
# Calculate Mu and Std from RandATAC Score List
mu_enh, std_enh = norm.fit(EnhrandATACSkewScoreList)
print("mu_enh: ", mu_enh, "std_enh :",  std_enh)

# Making the output file
print("Making output file")
outputfileName = str(options.output)
outputFile = open(str(outputfileName), 'w')
outputheader = "\t".join([enhFirstLine.strip("\n"), 'ATACSkewScore', randFirstLine.strip("\n"), 'RandSkewScore' , 'ZScore', 'PValue', "\n"])
outputFile.write(outputheader)

for enhancer in enhNameList:
    if Insitu_use[enhancer] == "yes":
        outputFile.write(enhInfoDict[enhancer])
        outputFile.write("\t")
        outputFile.write(str(ATACSkewScoreDict[enhancer]))
        outputFile.write("\t")
        EnhancerZscore = (ATACSkewScoreDict[enhancer] - float(mu_enh))/float(std_enh)
        enhpValue = norm.sf(abs((EnhancerZscore)))*2
        outputFile.write(enhRandregionMatchDict[enhancer])
        outputFile.write("\t")
        outputFile.write(str(MatchedRandomATACSkewScoreDict[enhancer]))
        outputFile.write("\t")
        outputFile.write(str(EnhancerZscore))
        outputFile.write('\t')
        outputFile.write(str(enhpValue))
        outputFile.write('\n')
outputFile.close() 

#Making the graphs

#random region histogram
plt.hist(EnhrandATACSkewScoreList, bins=20, normed=True, alpha=0.6, color='r')
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu_enh, std_enh)
plt.plot(x, p, 'k', linewidth=2)
title = "Fit results: mu = %.2f,  std = %.2f" % (mu_enh, std_enh)
plt.title(title)
savefig('RandomRegionHistogram.png')
plt.close()
print("Done calculating Enhancer_RandomRegion Curve")