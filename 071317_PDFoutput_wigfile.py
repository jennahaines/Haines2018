#############################################################
######### 071317_PDFoutput_wigfile.py               #########
######### Take a directory of wig snapshots,        #########
#########  insitus, and data and compile it into    #########
#########  one pdf per gene                         #########
#############################################################

#!/usr/bin/env python


import os
from math import *
from optparse import OptionParser
import sys
import re
from pylab import figure, title, xlabel, ylabel, hist, axis, grid, savefig
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.stats import norm
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch

def parse_options():
	parser = OptionParser()
	parser.add_option("-W", "--wigfileDest", dest="wigfileDest",
					  help="Directory where wig_snapshot.png are")
	parser.add_option("-R", "--regfile", dest="regfile",
					  help="14 columned TXT file of wanted genomic regions columns: Sample	Chr	Start	End	Whole	Ant	LinregPost	AplusP	Dnase1	Type	Location	Name	Zscore	Pvalue")
	parser.add_option("-S", "--insitu", dest="insitu",
					  help="Directory where insitu.png are")
	(options, args) = parser.parse_args()
	return options
            
options = parse_options()
parser = OptionParser()
if not options.wigfileDest:
    print("Wigfile option is missing\n")
    parser.print_help()
    exit(-1)
if not options.regfile:
    print("Region file option is missing\n")
    parser.print_help()
    exit(-1)
if not options.insitu:
    print("output file is missing\n")
    parser.print_help()
    exit(-1)

################# Data Import   #################

RegionList_file = open(str(options.regfile), 'r')
WigfileDirectory = str(options.wigfileDest)
InsituDirectory = str(options.insitu)


RegionList = RegionList_file.readlines()
RegionList_file.close()

#print(RegionList)

#Make a test regionList
regionList_test = RegionList[0:10]

enhInfoDict={}
enhLocDict={}
enhTypeDict={}
enhZScoreDict={}
enhPValueDict={}
enhNameList =[]



#for reportlab: Justifying the text for some reason. I probably don't need to do this
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))


with PdfPages('multipage_pdf.pdf') as pdf: # making a pdf called multipage_pdf.pdf
    for enhancer in RegionList:
        line2 = enhancer.strip('\n')
        enhListElement = line2.split('\t')
        enhNameList.append(enhListElement[10])
        enhInfoDict[enhListElement[10]] = line2 #Dictionary that returns the line without line ending
        enhLocDict[enhListElement[10]] = enhListElement[9] #Dictionary that returns Location
        enhTypeDict[enhListElement[10]] = enhListElement[8] #Dictionary that returns Type
        enhZScoreDict[enhListElement[10]] = enhListElement[11] #Dictionary that returns ZScore
        enhPValueDict[enhListElement[10]] = enhListElement[12] #Dictionary that returns Pvalue
        Wigname = "".join([str(enhListElement[10]), "_WigFileSnapshot.png"]) #making the wigname
        PDFTitle = str(enhListElement[10]) # Making the PDF title the name of the enhancer
        PDFName = "".join([PDFTitle, "_Report.pdf"]) #Making the PDFName for the PDF File name the name of the enhancer + .pdf
        #for reportlab: first specify the document format
        doc = SimpleDocTemplate(PDFName,pagesize=letter,
                            rightMargin=72,leftMargin=72,
                            topMargin=72,bottomMargin=18)
        #for reportlab: make a story variable that is a list that contains all images + text
        Story=[]
        ptext = '<font size=20>%s</font>' % PDFTitle #creating the title from the enhancer name
        Story.append(Paragraph(ptext, styles["Normal"])) #Adding the title to the story
        Story.append(Spacer(2, 12)) #adding a space of one line below the title
        for file in os.listdir(WigfileDirectory): #list the contents of the wigfile directory
            #print(file) #TEST
            if str(file) == Wigname:  #if the file is the same as the enhancer 
                WigSnapshotimage = str(file) #make a string of the file name
                WigSnapshotimage_Path = "".join([str(WigfileDirectory) , "/", WigSnapshotimage]) #add the path
                WigSnapshotimage_object = Image(WigSnapshotimage_Path, 5*inch, 4*inch) #make an image object
                #print(WigSnapshotimage_object) # TEST
                Story.append(WigSnapshotimage_object) #add it to the story
                Story.append(Spacer(1, 12)) # add a space 
        for insitu in os.listdir(InsituDirectory): #For each insitu in the insitu directory
            #print(insitu) #TEST
            insituName = insitu.split("_insit")
            insituName = insituName[0] #Obtain the insitu name 
            #print(insituName) #TEST
            if str(insituName) == str(enhListElement[10]): #if the insitu name is the same as the enhancer name
                Insituimage = str(insitu) #make the insitu file a string 
                Insituimage_Path = "".join([str(InsituDirectory) , "/", Insituimage]) #add it's path
                Insituimage_object = Image(Insituimage_Path, 5*inch, 4*inch) #make it into an image object
                Story.append(Insituimage_object) #append it to the PDF
        #Creating the string from the pvalue and Zscore
        DataString = "".join(["Location: ", str(enhListElement[9]), " Type: ", str(enhListElement[8]), " ZScore: ", str(enhListElement[11]), " PValue: ", str(enhListElement[12])])
        DataText = '<font size=12>%s</font>' % DataString #create the text with the enhancer data 
        Story.append(Paragraph(DataText, styles["Normal"])) #Adding the data text to the story
        doc.build(Story) #build the story list into the document
        print(Story) #TEST