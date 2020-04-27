#########################################################################################################
##########  DnaseI_postalign.sh                                                                ##########
########## 011818_postalign_RemoveDuplicates_SortedBam_PairedEnd_Sam_FilterFragments_DNASEI.sh ##########
##########   Input:  sam file                                                                  ##########
#########    Output: DupRemovedBam file                                                        ##########
#########    ATAC shifted, 130 and below bed file, wig file, 10Mnormalized wig file            ##########
#########################################################################################################

#bash DnaseI_postalign.sh 2> DnaseI_postalign_STERR.txt | cat > DnaseI_postalign_STDOUT.txt &

#! /bin/bash

####### First I aligned SRA Dnase 1 data downloaded in 100317_SRAdump.sh
#100417_BDTNP_DNaseI_stage5
bowtie2 -p 10 -5 5 -3 5 -N 1 -X 2000 --local --very-sensitive-local --un 100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_UNMAPPED.fastq -x ~/Indexed_Genomes/dmel-r5.57-index -U SRR060796.fastq,SRR060797.fastq,SRR060798.fastq,SRR060799.fastq -S 100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED.sam &


#100417_BDTNP_DNaseI_stage9
bowtie2 -p 10 -5 5 -3 5 -N 1 -X 2000 --local --very-sensitive-local --un 100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED_UNMAPPED.fastq -x ~/Indexed_Genomes/dmel-r5.57-index \
-U SRR060800.fastq,SRR060801.fastq,SRR060802.fastq,SRR060803.fastq,SRR060804.fastq,SRR060805.fastq \
-S 100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED.sam &

#100417_BDTNP_DNaseI_stage11
bowtie2 -p 10 -5 5 -3 5 -N 1 -X 2000 --local --very-sensitive-local --un 100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED_UNMAPPED.fastq -x ~/Indexed_Genomes/dmel-r5.57-index \
-U SRR060775.fastq,SRR060776.fastq,SRR060778.fastq,SRR060779.fastq \
-S 100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED.sam &

#100417_BDTNP_DNaseI_stage14
bowtie2 -p 10 -5 5 -3 5 -N 1 -X 2000 --local --very-sensitive-local --un 100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED_UNMAPPED.fastq -x ~/Indexed_Genomes/dmel-r5.57-index \
-U SRR060780.fastq,SRR060781.fastq,SRR060782.fastq,SRR060783.fastq \
-S 100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED.sam &

#######  made a new working directory where everything will be in
#Directory='/netdata/jhaines/BDTNP_s5_DNase-seq/100317-SRADownload/011818_AnalysisforReviews_PE'


####### Loop to process the sam files into normalized wig files.
# original sam file
Bowtie2SamFile[1]=100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED.sam
Bowtie2SamFile[2]=100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED.sam
Bowtie2SamFile[3]=100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED.sam
Bowtie2SamFile[4]=100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED.sam


#Core name that will be carried with the sample
filesetname[1]=011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED
filesetname[2]=011818_100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED
filesetname[3]=011818_100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED
filesetname[4]=011818_100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED

for k in 1 2 3 4
do
    samfile=${Bowtie2SamFile[${k}]}
    echo 'sam file imported'
    filesetname=${filesetname[${k}]}
    samtools view -bS -F 1804 -q 30 ${samfile} | samtools sort - ${filesetname}
    java -jar ~/scripts/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE \
      I=${filesetname}.bam \
      O=${filesetname}_nodups.bam \
      M=${filesetname}_marked_dup_metrics.txt
    echo 'remove duplicates from bam complete'
    bedtools bamtobed -i ${filesetname}_nodups.bam > ${filesetname}_shifted.bed
    echo 'bed file created'
done
echo 'all samples done'
mkdir 011818_Bedfiles
mv *.bed 011818_Bedfiles
perl ~/scripts/xl-bed2wig-dirproc-HOA-2014-new-chrsz 011818_Bedfiles

#This script will take the filtered bed files.. count the number of reads and then use that to noramlize the wig files

for k in 1 2 3 4
do
  cd 011818_Bedfiles
  filesetname=${filesetname[${k}]}
  NumberofReads=$(wc -l < "${filesetname}_shifted.bed")
  echo "$NumberofReads"
  export NumberofReads
  echo "Normalizing to 1M reads"
  perl -n -e ' $Scale = 10000000/$ENV{NumberofReads} ; if (/^\d+/) {@A = split (/\s+/, $_); $norm = $A[1] * $Scale ; print "$A[0]\t $norm \n"} else { print "$_"}' ${filesetname}_shifted.wig > ${filesetname}_shifted_10MNorm.wig
  echo "${filesetname} Normalized wig File made"

done