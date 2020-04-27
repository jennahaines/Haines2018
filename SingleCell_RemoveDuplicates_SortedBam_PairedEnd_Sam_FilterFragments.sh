####################################################################################################
##########   SingleCell_RemoveDuplicates_SortedBam_PairedEnd_Sam_FilterFragments.sh       ##########
##########   Input:  bam file from Cusanovich 2017                                        ##########
#########    Output: shifted bed file, merged bed files, wig file, 10Mnormalized wig file ##########
####################################################################################################

#bash 040418_postalign_RemoveDuplicates_SortedBam_PairedEnd_Sam_FilterFragments.sh 2> 040418_postalign_RemoveDuplicates_SortedBam_PairedEnd_Sam_FilterFragments_STERR.txt | cat > 040418_postalign_RemoveDuplicates_SortedBam_PairedEnd_Sam_FilterFragments_STOUT.txt &

#! /bin/bash

# bam file from Cusanovich 2018
Bowtie2SamFile[1]=SCatac_tha.bowtie.2to4.fixed.nodups.15.bam
Bowtie2SamFile[2]=SCatac_tha.bowtie.2to4.fixed.nodups.16.bam
Bowtie2SamFile[3]=SCatac_tha.bowtie.2to4.fixed.nodups.4.bam
Bowtie2SamFile[4]=SCatac_tha.bowtie.2to4.fixed.nodups.6.bam
Bowtie2SamFile[5]=SCatac_tha.bowtie.2to4.fixed.nodups.7.bam


#Core name that will be carried with the sample

filesetname[1]=SCatac_tha.bowtie.2to4.fixed.nodups.15
filesetname[2]=SCatac_tha.bowtie.2to4.fixed.nodups.16
filesetname[3]=SCatac_tha.bowtie.2to4.fixed.nodups.4
filesetname[4]=SCatac_tha.bowtie.2to4.fixed.nodups.6
filesetname[5]=SCatac_tha.bowtie.2to4.fixed.nodups.7
                         

for k in 1 2 3 4 5
do
    samfile=${Bowtie2SamFile[${k}]}
    echo 'sam file imported'
    filesetname=${filesetname[${k}]}
    samtools sort -n ${filesetname}.bam ${filesetname}.srt.temp 
    bedtools bamtobed -bedpe -i ${filesetname}.srt.temp.bam | perl -n -e '@A = split (/\t/, $_); $start = $A[1] +4 ; $end = $A[5] - 6 ; print "$A[0]\t$start\t$end\t$A[6]\t$A[7]\t$A[8]\t$[9]\n"' > 040418_${filesetname}_shifted.bed
    echo 'shifted bed PE file made'
    # did not filter based on the 130bp threshold..since they did not do that in their pipeline
    #perl -n -e '@A = split (/\t/, $_); $readsize = abs($A[2] - $A[1]); print $_' 011617_${filesetname}_shifted.bed > 011617_${filesetname}_shifted.bed
    #echo 'filtered based on 130 bp threshold'
   
done

#merge anterior (clusters 6 and 15) and posterior (4,7,16) beds together

cat 040418_SCatac_tha.bowtie.2to4.fixed.nodups.6_shifted.bed 040418_SCatac_tha.bowtie.2to4.fixed.nodups.15_shifted.bed > 040418_MERGED_ANTERIOR_SCatac_tha.bowtie.2to4.fixed.nodups_shifted.bed
cat 040418_SCatac_tha.bowtie.2to4.fixed.nodups.4_shifted.bed 040418_SCatac_tha.bowtie.2to4.fixed.nodups.7_shifted.bed 040418_SCatac_tha.bowtie.2to4.fixed.nodups.16_shifted.bed > 040418_MERGED_POSTERIOR_SCatac_tha.bowtie.2to4.fixed.nodups_shifted.bed
echo 'Merged Bed files'


mkdir Bedfiles_040418
mv 040418_MERGED_*.bed Bedfiles_040418
perl ~/scripts/xl-bed2wig-dirproc-HOA-2014-new-chrsz Bedfiles_040418
echo 'made wig files'

#This script will take the filtered bed files.. count the number of reads and then use that to noramlize the wig files to 1M filtered reads

MERGEDfilesetname[1]=040418_MERGED_ANTERIOR_SCatac_tha.bowtie.2to4.fixed.nodups_shifted
MERGEDfilesetname[2]=040418_MERGED_POSTERIOR_SCatac_tha.bowtie.2to4.fixed.nodups_shifted

for k in 1 2
do
  cd Bedfiles_040418
  filesetname=${MERGEDfilesetname[${k}]}
  NumberofReads=$(wc -l < "${filesetname}.bed")
  echo "$filesetname" > outputfile.txt
  export NumberofReads
  echo "$NumberofReads" > outputfile.txt
  echo "Normalizing to 1M reads"
  perl -n -e ' $Scale = 1000000/$ENV{NumberofReads} ; if (/^\d+/) {@A = split (/\s+/, $_); $norm = $A[1] * $Scale ; print "$A[0]\t $norm \n"} else { print "$_"}' ${filesetname}.wig > ${filesetname}1MNorm.wig
  echo "${filesetname} Normalized wig File made"

done