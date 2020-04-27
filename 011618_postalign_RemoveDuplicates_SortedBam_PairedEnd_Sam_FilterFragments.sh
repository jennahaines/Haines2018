###################################################################################################
#########  011618_postalign_RemoveDuplicates_SortedBam_PairedEnd_Sam_FilterFragments.sh  ##########
#########    Input:  sam file                                                            ##########
#########    Output: DupRemovedBam file                                                  ##########
#########    ATAC shifted, 130 and below bed file, wig file, 10Mnormalized wig file      ##########
#########    Dependecy: samtools,
################################################################################################

#bash 011618_postalign_RemoveDuplicates_SortedBam_PairedEnd_Sam_FilterFragments.sh 2> 011618_postalign_RemoveDuplicates_SortedBam_PairedEnd_Sam_FilterFragments_STERR.txt | cat > 011618_postalign_RemoveDuplicates_SortedBam_PairedEnd_Sam_FilterFragments_STDOUT.txt &

#! /bin/bash


# original sam file
Bowtie2SamFile[1]=041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-01-20a.sam
Bowtie2SamFile[2]=041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-02-20p.sam
Bowtie2SamFile[3]=041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-5-1A.sam
Bowtie2SamFile[4]=041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-6-1p.sam
Bowtie2SamFile[5]=041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-7-10whole.sam
Bowtie2SamFile[6]=041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-01-20Ant.sam
Bowtie2SamFile[7]=041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-02-20Post.sam
Bowtie2SamFile[8]=041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-03-1plus2.sam
Bowtie2SamFile[9]=041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-04-10whole.sam

#Core name that will be carried with the sample
filesetname[1]=101117_RemDUP_PE_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-01-20a 
filesetname[2]=081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-02-20p
filesetname[3]=081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-5-1A   
filesetname[4]=081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-6-1p
filesetname[5]=081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-7-10whole         
filesetname[6]=081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-01-20Ant
filesetname[7]=081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-02-20Post 
filesetname[8]=081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-03-1plus2
filesetname[9]=081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-04-10whole

for k in 1 2 3 4 5 6 7 8 9
do
    samfile=${Bowtie2SamFile[${k}]}
    echo 'sam file imported'
    filesetname=${filesetname[${k}]}
    samtools view -bS -F 1804 -f 2 -q 30 ${samfile} | samtools sort - ${filesetname}
    java -jar ~/scripts/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE \
      I=${filesetname}.bam \
      O=${filesetname}_nodups.bam \
      M=${filesetname}_marked_dup_metrics.txt
    echo 'Duplicates removed'
    samtools sort -n ${filesetname}_nodups.bam ${filesetname}.srt.temp 
    bedtools bamtobed -bedpe -i ${filesetname}.srt.temp.bam | perl -n -e '@A = split (/\t/, $_); $start = $A[1] +4 ; $end = $A[5] - 6 ; print "$A[0]\t$start\t$end\t$A[6]\t$A[7]\t$A[8]\t$[9]\n"' > 011617_${filesetname}_shifted.bed
    echo 'shifted bed PE file made'
    perl -n -e '@A = split (/\t/, $_); $readsize = abs($A[2] - $A[1]); if ($readsize <= 130) {print $_}' 011617_${filesetname}_shifted.bed > 011617_${filesetname}_shifted_lessthan130.bed
    echo 'filtered based on 130 bp threshold'
   
done
echo 'making wig files'
mkdir Bedfiles_011617
mv *.bed Bedfiles_011617
perl xl-bed2wig-dirproc-HOA-2014-new-chrsz Bedfiles

#This script will take the filtered bed files.. count the number of reads and then use that to noramlize the wig files

for k in 1 2 3 4 5 6 7 8 9
do
  cd Bedfiles
  filesetname=${filesetname[${k}]}
  NumberofReads=$(wc -l < "011617_${filesetname}_shifted_lessthan130.bed")
  echo "$NumberofReads"
  export NumberofReads
  echo "Normalizing to 1M reads"
  perl -n -e ' $Scale = 1000000/$ENV{NumberofReads} ; if (/^\d+/) {@A = split (/\s+/, $_); $norm = $A[1] * $Scale ; print "$A[0]\t $norm \n"} else { print "$_"}' 011617_${filesetname}_shifted_lessthan130.wig > 011617_${filesetname}_shifted_lessthan130_1MNorm.wig
  echo "${filesetname} Normalized wig File made"

done