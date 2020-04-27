####################################################################################################################
##########   FigS8_PeakAnalysis_pipeline.sh                                                              ###########
##########   Uses wig_sig_around_bedfile_013018.pl to summarize wig signal around all 10bp regions       ###########
##########   Runs 020618_RandomRegions_AroundPeaks.py                                                    ###########
##########   Overlap and intersect peaks with bedtools                                                   ###########
##########   Run  020718_Pivotlistconvert.py to change datatable format                                  ###########
####################################################################################################################


#Wig Signal from all TF and ATAC data around these whole peaks
perl /Users/jennahaines/'Box Sync'/Eisen_Lab/Scripts/wig_sig_around_bedfile_013018.pl \
020818_WholePooledInRep1AndRep2.narrowPeak_JustRegions.bed \
~/020218_TFwigfiles \
021618_wig_sig_around_Wholepeaks_TFbindingDirectory.txt

#random regions to get ATAC skew/ Positional scores and significance
python 020618_RandomRegions_AroundPeaks.py \
-W 013118_wig_sig_around_bedfile_randomregions.txt \
-R 021618_wig_sig_around_Wholepeaks_TFbindingDirectory.txt \
-o 020618_030628WholePeaks_Randomregv2_output.txt


# ############################# Peak overlap #############################
perl -pe '$_ =~ tr/chr//d' 020818_WholePooledInRep1AndRep2.narrowPeak_JustRegions.bed > 020818_WholePooledInRep1AndRep2.narrowPeak_JustRegions_nochr.bed
#took out Uextra, and Many het chromosomes to do this analysis

bedtools intersect -wao -names Zelda1hr Zelda2hr Zelda3hr Hb Kr Gt Bcd Cad Kni Hoskins Kvon Prom H3K27acc14a H3K27acc14c H3K4me1c14a H3K4me1c14c \
-a 020718_WholePeaks_plusATACScore_PValue_nochr.bed \
-b \
100817s6_1hr_ZLD-Inpc_peaks_nochr.bed \
100817s7_2hr_ZLD-Inpc_peaks_nochr.bed \
100817s8_3hr_ZLD-Inpc-chr_peaks_nochr.bed \
DvS2TG_090826s6_HB_peaks_nochr.bed \
DvS2TG_090826s7_KR_peaks_nochr.bed \
DvS2TG_090826s8_GT_peaks_nochr.bed \
GSM511083_Dmel-BCD_peaks_nochr.bed \
GSM511087_Dmel-CADpeaks_nochr.bed \
GSM511088_Dmel-KNIpeaks_nochr.bed \
011618_ReviewsRevised_RevisedAPDV_enhPromRegions_nochr.bed \
Kvon_AllCRMS_nochr.bed \
Supplementary_data_file3_IntegratedPromoters_nochr_shortened.bed \
GSM1424902_Dmel-H3K27ac-c14a-peaks_nochr.bed \
GSM1424903_Dmel-H3K27ac-c14c-peaks_nochr.bed \
GSM1424906_Dmel-H3K4me1-c14a-peaks_nochr.bed \
GSM1424907_Dmel-H3K4me1-c14c-peaks_nochr.bed > 020718_ReplicatedwholePeaks_BedtoolsOverlap_Kitchensink.txt


# First I want to sort them to make everything go faster
sort -k1,1 -k2,2n 

#### Run 020718_Pivotlistconvert.py to change the bedtools output into a file format that makes sense
python 020718_Pivotlistconvert.py -W 020718_ReplicatedwholePeaks_BedtoolsOverlap_Kitchensink.txt -o 020718_WholePeaks_Pivotlistconvert_output.txt -n Names.txt