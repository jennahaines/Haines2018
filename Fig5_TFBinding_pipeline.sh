###############################################################################################################
##########   Fig5_TFBinding_pipeline.sh                                                             ###########
##########   Uses wig_sig_around_bedfile_013018.pl to summarize wig signal around all 10bp regions  ###########
##########   Runs 061617_WigFileSnapshots_V2.py                                                     ###########
##########   Pulls out the wig signal around a given range                                          ###########
###############################################################################################################
# ############################# Tf binding figure #############################

# Generate text file of average wig signal around 10 bp regions
perl wig_sig_around_bedfile_013018_for10bpwindows.pl \
/Volumes/'JHHD 1'/020118_dm3_10bpwindows.bed \
~/020218_TFwigfiles \
020218_TF_wig_sig_10bp_windows_XY.txt

## Make a text file containing only wig signal for supplied bed intervals
mkdir 020218_TF_Wigfilesplitter
cd 020218_TF_Wigfilesplitter

python /Users/jennahaines/'Box Sync'/Eisen_Lab/Scripts/061617_WigFileSnapshots_V2.py \
-W ~/020218_TF_wig_sig_10bp_windows_XY.txt \
-R ~/012318_ReviewsRevised_RevisedAPDV_enhPromRegions_FINAL_YESONLY.txt \
-l 3000