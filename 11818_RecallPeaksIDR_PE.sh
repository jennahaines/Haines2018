#########################################################################################
##########   011818_RecallPeaksIDR_PE.sh                                      ###########
##########   Calls peaks with Macs2 with bedfiles                             ###########
##########   Overlaps peaks from both replicates to find a set of peaks       ###########
##########   found in both replicates                                         ###########
#########################################################################################
#bash 011818_RecallPeaksIDR_PE.sh 2> 011818_RecallPeaksIDR_PE_stderr.txt | cat > 011818_RecallPeaksIDR_PE_stdout.txt &
#! /bin/bash

#Import files
bedfile[1]=011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-02-20p_shifted
bedfile[2]=011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-02-20p_shifted_lessthan130
bedfile[3]=011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-5-1A_shifted_lessthan130
bedfile[4]=011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-6-1p_shifted_lessthan130
bedfile[5]=011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-7-10whole_shifted_lessthan130
bedfile[6]=011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-01-20Ant_shifted_lessthan130
bedfile[7]=011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-02-20Post_shifted_lessthan130
bedfile[8]=011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-03-1plus2_shifted_lessthan130
bedfile[9]=011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-04-10whole_shifted_lessthan130
bedfile[10]=011617_101117_RemDUP_PE_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-01-20a_shifted_lessthan130
bedfile[11]=011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_10whole_shifted_lessthan130
bedfile[12]=011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130
bedfile[13]=011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130

for k in {2..13}
do
    bedfile=${bedfile[${k}]}
    perl -pe '$_ =~ tr/chr//d' ${bedfile}.bed > ${bedfile}_nochr.bed
    echo "Took out chr"
    macs2 callpeak --nomodel -f BEDPE -g dm -p 1e-3 --call-summits --bdg -t ${bedfile}_nochr.bed -n ${bedfile}_nochr
    echo "Called peaks"
    Peakfile=${bedfile[${k}]}_nochr_peaks.narrowPeak
    sort -k 8gr,8gr $Peakfile | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -c > Sorted_$Peakfile.gz
done

#Make a directory for all of the peak files
mkdir 011818_RecallPeaksIDR_PE_Dir_PEAKS
mv *_nochr* 011818_RecallPeaksIDR_PE_Dir_PEAKS
cd 011818_RecallPeaksIDR_PE_Dir_PEAKS

#Anterior High Confidence replicate Overlap peaks
Antrep1=Sorted_011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-01-20Ant_shifted_lessthan130_nochr_peaks.narrowPeak.gz
Antrep2=Sorted_011617_101117_RemDUP_PE_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-01-20a_shifted_lessthan130_nochr_peaks.narrowPeak.gz
Antmerged=Sorted_011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_nochr_peaks.narrowPeak.gz

#Intersect merged with rep 1 then intersect that with rep2
intersectBed -wo -a ${Antmerged} \
-b ${Antrep1} \
| awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5))
{print $0}}' | cut -f 1-10 | sort | uniq | \
intersectBed -wo -a stdin -b ${Antrep2} \
| awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5))
{print $0}}' | cut -f 1-10 | sort | uniq > 011818_ANTERIORPooledInRep1AndRep2.narrowPeak.gz

echo "Anterior High confidence peak file made"

# Posterior High Confidence replicate Overlap peaks
Postrep1=Sorted_011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-02-20p_shifted_lessthan130_nochr_peaks.narrowPeak.gz
Postrep2=Sorted_011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-02-20Post_shifted_lessthan130_nochr_peaks.narrowPeak.gz
Postmerged=Sorted_011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_nochr_peaks.narrowPeak.gz

#Intersect merged with rep 1 
intersectBed -wo -a ${Postmerged} \
-b ${Postrep1} \
| awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5))
{print $0}}' | cut -f 1-10 | sort | uniq | \
intersectBed -wo -a stdin -b ${Postrep2} \
| awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5))
{print $0}}' | cut -f 1-10 | sort | uniq > 011818_POSTERIORPooledInRep1AndRep2.narrowPeak.gz

echo "Posterior High confidence peak file made"

# Whole High Confidence replicate Overlap peaks
Wholerep1=Sorted_011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105-ATACSlice2-7-10whole_shifted_lessthan130_nochr_peaks.narrowPeak.gz
Wholerep2=Sorted_011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109-ATACSlice03-04-10whole_shifted_lessthan130_nochr_peaks.narrowPeak.gz
Wholemerged=Sorted_011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_10whole_shifted_lessthan130_nochr_peaks.narrowPeak.gz

#Whole merged with rep 1 
intersectBed -wo -a ${Wholemerged} \
-b ${Wholerep1} \
| awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5))
{print $0}}' | cut -f 1-10 | sort | uniq | \
intersectBed -wo -a stdin -b ${Wholerep2} \
| awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5))
{print $0}}' | cut -f 1-10 | sort | uniq > 011818_WholePooledInRep1AndRep2.narrowPeak.gz

echo "Whole High confidence peak file made"

echo "All Done! :)"