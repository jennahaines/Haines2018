################################################################################################
##########   011718_10Mnormalization_mergedfiles.sh                                   ##########
##########   Input:  merged bed files  + merged wig files                             ##########
#########    Output: merged wig files normalized to 10M reads                         ##########
################################################################################################

#bash 011718_10Mnormalization_mergedfiles.sh 2> 011718_10Mnormalization_mergedfiles_STERR.txt | cat > 011718_10Mnormalization_mergedfiles_STDOUT.txt &

#! /bin/bash
### First I concatenated individual bedfiles into a merged file. 
cat 011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-02-20p_shifted_lessthan130.bed 011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-02-20Post_shifted_lessthan130.bed > 011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130.bed &
cat 011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-01-20Ant_shifted_lessthan130.bed 011617_101117_RemDUP_PE_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-01-20a_shifted_lessthan130.bed > 011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130.bed &
cat 011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-7-10whole_shifted_lessthan130.bed 011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-04-10whole_shifted_lessthan130.bed > 011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_10whole_shifted_lessthan130.bed &
### make a new directory and move merged files into that directory 
mkdir 011718_mergedfiles
mv 011718_MERGED 011718_mergedfiles
perl ~/scripts/xl-bed2wig-dirproc-HOA-2014-new-chrsz 011718_mergedfiles

### Then I calculate the normalization factor for 10 million reads by dividing 10,000,000 by the number of lines in the merged bed file. I multiply the wig signal by that calculated number to normalize the signal between wig files.
filesetname[1]=011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_10whole
filesetname[2]=011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant
filesetname[3]=011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post

for k in 1 2 3 
do
  cd 011718_mergedfiles
  filesetname=${filesetname[${k}]}
  NumberofReads=$(wc -l < "${filesetname}_shifted_lessthan130.bed")
  echo "$NumberofReads"
  export NumberofReads
  echo "Normalizing to 10M reads"
  perl -n -e ' $Scale = 10000000/$ENV{NumberofReads} ; if (/^\d+/) {@A = split (/\s+/, $_); $norm = $A[1] * $Scale ; print "$A[0]\t $norm \n"} else { print "$_"}' ${filesetname}_shifted_lessthan130.wig > ${filesetname}_shifted_lessthan130_10MNorm.wig
  echo "${filesetname} Normalized wig File made"

done