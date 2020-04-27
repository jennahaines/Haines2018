############################################################################################################
##########   Random_regions_pipeline.sh                                                          ###########
##########   uses bedtools to make a list of random regions                                      ###########
##########   Uses wig_sig_around_bedfile_013018.pl to summarize wig signal around random regions ###########
##########   Runs 063017_RandomDistributionScript_V5.py                                          ###########
############################################################################################################
# Make a list of random regions excluding genes and patterning regions
bedtools intersect -v -a 081517_Dmel5_1E6_randomregions.bed -b 081517_newgenelist_stripped_plus_FlybaseGenes_coordinates.bed > 081517_RandomRegions_Excl_flybasegenes_Enhancers.bed

# Random region match up 

perl ~/Scripts/wig_sig_around_bedfile_013018.pl \
013118_081517_RandomRegions_Excl_flybasegenes_Enhancers.bed \
~/013118_RandomWigDirectory \
013118_wig_sig_around_bedfile_randomregions.txt


#063017_RandomDistributionScript_V5.py
python 012418_RandomDistributionScript_V5.py \
-W 013118_wig_sig_around_bedfile_randomregions.txt \
-R 013118_wig_sig_around_bedfile_012318GENELIST_matched.txt \
-o 013118_012318GENELIST_Randomregv2_output.txt

echo "randomization script complete"