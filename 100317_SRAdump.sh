######################################################################
##########   100317_DnaseISRAdump.sh                              ##########
######################################################################


#! /bin/bash

# already did all of stage 5

# stage 9 Rep 1
inputname[1]=SRR060800

# stage 9 Rep 2
inputname[2]=SRR060801  
inputname[3]=SRR060802  
inputname[4]=SRR060803  
inputname[5]=SRR060804  
inputname[6]=SRR060805  

# stage 10 Rep 1
inputname[7]=SRR060769  
inputname[8]=SRR060770  
inputname[9]=SRR060771  

# stage 10 Rep 2
inputname[10]=SRR060772 
inputname[11]=SRR060773 
inputname[12]=SRR060774 

#Stage 11 rep1
inputname[13]=SRR060775 
inputname[14]=SRR060776 

#Stage 11 rep2
inputname[15]=SRR060778 
inputname[16]=SRR060779 

#stage 14 rep1
inputname[17]=SRR060780 
inputname[18]=SRR060781 

#stage 14 rep2
inputname[19]=SRR060782 
inputname[20]=SRR060783 

for k in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
    input=${inputname[${k}]}
    fastq-dump ${input}
    echo "done"
done