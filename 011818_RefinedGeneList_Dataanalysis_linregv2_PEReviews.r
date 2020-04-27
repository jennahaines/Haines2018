############### ################ ################ ################ ################
################ 011818_RefinedGeneList_Dataanalysis_linregv2_PEReviews.r     ################ 
################ make scatterplots with ATAC Score and position               ################ 
################calculated in 051117_randomdistribution_script                ################
################ ################ ################ ################ ################


library("plyr")
library("dplyr")
library("ggplot2")


############################ Linear Regression Normalization ############################ 

###### ### ###  Working directory
setwd("~/Box Sync/Eisen_Lab/Experiments/ATAC-seq/Halves/ATAC-seq_Pools/040617_Analysis/061417_BrowserTraces/062317_RegressionNormalization/011618_PairedEndAnalysisforReviews/")

###### ### ###  Input files
# Opens the entire merged wig file made in bedtools that binned every 1Kb
# first input file is 011918_multiBigWigCompare_1KB_041417_Bowtie2_ATACPools_merged.txt
# multiBigwigSummary bins -bs 1000 -b \
# /Users/jennahaines/'Box Sync'/Eisen_Lab/Experiments/ATAC-seq/Halves/ATAC-seq_Pools/040617_Analysis/041717_NormalizedWigs/bigwig/011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_10whole_shifted_lessthan130_10MNorm.bw \
# /Users/jennahaines/'Box Sync'/Eisen_Lab/Experiments/ATAC-seq/Halves/ATAC-seq_Pools/040617_Analysis/041717_NormalizedWigs/bigwig/011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm.bw \
# /Users/jennahaines/'Box Sync'/Eisen_Lab/Experiments/ATAC-seq/Halves/ATAC-seq_Pools/040617_Analysis/041717_NormalizedWigs/bigwig/011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm.bw \
# /Users/jennahaines/'Box Sync'/Eisen_Lab/Experiments/ATAC-seq/Halves/ATAC-seq_Pools/040617_Analysis/041717_NormalizedWigs/bigwig/011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm.bw \
# -out 011918_multiBigWigCompare_1KBbins_041417_Bowtie2_ATACPools_merged.npz --outRawCounts 011918_multiBigWigCompare_1KB_041417_Bowtie2_ATACPools_merged.txt

merged_norm_1kb = read.delim2("011918_multiBigWigCompare_1KB_041417_Bowtie2_ATACPools_merged.txt", sep = "\t", stringsAsFactors = FALSE) 

colnames(merged_norm_1kb)  = c("chr", "start", "end", "Whole", "Ant", "Post", "DnaseI")

allSamples_NoNAs = merged_norm_1kb %>%
  select(chr, start, end, Whole, Ant, Post, DnaseI) %>%
  filter(Whole != 'nan') %>%
  filter(Ant != 'nan') %>%
  filter(Post != 'nan') %>%
  filter(DnaseI != 'nan') %>%
  mutate(Type = "NA") %>%
  mutate(Location = "NA") %>%
  mutate(Name = "NA") %>%
  mutate(Dotsize = 1)
### make linear regressions

merged_norm_1kb.lm = lm(as.numeric(allSamples_NoNAs$Post) ~ as.numeric(allSamples_NoNAs$Ant))
### make linear regressions
merged_norm_1kb_AntvsWhole.lm = lm(as.numeric(allSamples_NoNAs$Ant) ~ as.numeric(allSamples_NoNAs$Whole))
merged_norm_1kb_PostvsWhole.lm = lm(as.numeric(allSamples_NoNAs$Post) ~ as.numeric(allSamples_NoNAs$Whole))

AP_NoNormCor = cor(as.numeric(allSamples_NoNAs$Ant),as.numeric(allSamples_NoNAs$Post), method = "spearman")
###### ### ###  Ant vs Whole
# Call:
#   lm(formula = as.numeric(allSamples_NoNAs$Ant) ~ as.numeric(allSamples_NoNAs$Whole))
# 
# Coefficients:
#   (Intercept)  
# 36.0167  
# as.numeric(allSamples_NoNAs$Whole)  
# 0.5787

###### ### ###  Post vs Whole
# Call:
#   lm(formula = as.numeric(allSamples_NoNAs$Post) ~ as.numeric(allSamples_NoNAs$Whole))
# 
# Coefficients:
#   (Intercept)  
# 50.0868  
# as.numeric(allSamples_NoNAs$Whole)  
# 0.6063 


# Call:
#   lm(formula = as.numeric(allSamples_NoNAs$Post) ~ as.numeric(allSamples_NoNAs$Ant))
# 
# Coefficients:
#   (Intercept)  
# 16.939  
# as.numeric(allSamples_NoNAs$Ant)  
# 0.996  

############## No Normalization ############## 

### Ant vs Post - No Lin reg normalization
png('011918_AllRegions1KB_Scatter_Merged_PE_10M_NoNormalization_AntVsPost.png', width = 2000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Ant),y= as.numeric(Post))) +
  #scale_color_manual(name = 'Legand', values=c("darkorange1" , "dodgerblue3")) +
  geom_point(colour="grey35", alpha = 0.2) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('011918_AllRegions1KB_Scatter_Merged_PE_10M_NoNormalization_AntVsPost') +
  #scale_x_log10() +
  #scale_y_log10() +
  xlab('Anterior Halves merged') +
  ylab('Posterior Halves merged') +
  geom_smooth(method='lm', formula = y~x) +
  ylim(0,NA) + 
  xlim(0,NA) +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=100, y=500, label= "slope = 0.996 ", colour = 'darkorange1') +
  #annotate("text", x=500, y=1500, label= Post_Prom_Label, colour = 'dodgerblue3') +
  theme(panel.border = element_rect(fill = NA),
        panel.background = element_rect(fill= NA),
        axis.title.x = element_text(vjust = 0, size = 12),
        axis.title.y = element_text(vjust = 1, size = 12),
        axis.text.x = element_text(size=8),
        axis.text.y  = element_text(size=8),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()

############ Whole Vs Ant ############ 

png('011918_AllRegions1KB_Scatter_Merged_PE_10M_NoNormalization_AntVsWhole.png', width = 2000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Whole),y= as.numeric(Ant))) +
  geom_point(colour="grey35", alpha = 0.2) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('011918_AllRegions1KB_Scatter_Merged_PE_10M_NoNormalization_WholevsAnt') +
  xlab('Whole Halves merged') +
  ylab('Ant Halves merged') +
  geom_smooth(method='lm', formula = y~x) +
  ylim(0,NA) + 
  xlim(0,NA) +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  theme(panel.border = element_rect(fill = NA),
        panel.background = element_rect(fill= NA),
        axis.title.x = element_text(vjust = 0, size = 12),
        axis.title.y = element_text(vjust = 1, size = 12),
        axis.text.x = element_text(size=8),
        axis.text.y  = element_text(size=8),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()

############ Whole Vs Post ############ 

png('011918_AllRegions1KB_Scatter_Merged_PE_10M_NoNormalization_WholevsPost.png', width = 2000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Whole),y= as.numeric(Post))) +
  geom_point(colour="grey35", alpha = 0.2) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('011918_AllRegions1KB_Scatter_Merged_PE_10M_NoNormalization_WholevsPost') +
  xlab('Whole Halves merged') +
  ylab('Post Halves merged') +
  geom_smooth(method='lm', formula = y~x) +
  ylim(0,NA) + 
  xlim(0,NA) +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  theme(panel.border = element_rect(fill = NA),
        panel.background = element_rect(fill= NA),
        axis.title.x = element_text(vjust = 0, size = 12),
        axis.title.y = element_text(vjust = 1, size = 12),
        axis.text.x = element_text(size=8),
        axis.text.y  = element_text(size=8),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()

############## Lin Reg Normalization ############## 
png('011918_AllRegions1KB_Scatter_Merged_PE_10M_Linregwhole_WholevsAnt.png', width = 2000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Whole),y= (as.numeric(Ant)/ as.numeric(merged_norm_1kb_AntvsWhole.lm$coefficients[2])))) +
  geom_point(colour="grey35", alpha = 0.2) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('011918_AllRegions1KB_Scatter_Merged_PE_10M_Linregwhole_WholevsAnt') +
  # scale_x_log10() +
  # scale_y_log10() +
  ylim(0,NA) + 
  xlim(0,NA) +
  xlab('Whole Halves merged') +
  ylab('Ant LinRegWhole merged') +
  geom_smooth(method='lm', formula = y~x) +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  theme(panel.border = element_rect(fill = NA),
        panel.background = element_rect(fill= NA),
        axis.title.x = element_text(vjust = 0, size = 12),
        axis.title.y = element_text(vjust = 1, size = 12),
        axis.text.x = element_text(size=8),
        axis.text.y  = element_text(size=8),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()

png('011918_AllRegions1KB_Scatter_Merged_PE_10M_Linregwhole_WholevsPost.png', width = 2000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Whole),y= (as.numeric(Post)/ as.numeric(merged_norm_1kb_PostvsWhole.lm$coefficients[2])))) +
  geom_point(colour="grey35", alpha = 0.2) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('011918_AllRegions1KB_Scatter_Merged_PE_10M_Linregwhole_WholevsPost') +
  xlab('Whole Halves merged') +
  ylab('Post linregWhole') +
  geom_smooth(method='lm', formula = y~x) +
  ylim(0,3000) + 
  xlim(0,3000) +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  theme(panel.border = element_rect(fill = NA),
        panel.background = element_rect(fill= NA),
        axis.title.x = element_text(vjust = 0, size = 12),
        axis.title.y = element_text(vjust = 1, size = 12),
        axis.text.x = element_text(size=8),
        axis.text.y  = element_text(size=8),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()


######### ######### ######### Input ######### ######### ######### ######### ######### ######### 

setwd("~/Box Sync/Eisen_Lab/Experiments/ATAC-seq/Halves/ATAC-seq_Pools/040617_Analysis/061417_BrowserTraces/062317_RegressionNormalization/011618_PairedEndAnalysisforReviews/")

allSamples <-read.delim2("012318_2reps_nodups_PE_130Cutoff__NewgeneList_Randomregv2_output.txt", sep = "\t", stringsAsFactors = FALSE)
allSamples_1kb_bins <- read.delim2("012318_2reps_nodups_PE_130Cutoff__multiBigWigCompare_1kbbins.txt", sep = "\t", stringsAsFactors = FALSE)
allPeaksFile <- read.delim2("012418_2reps_REppeaks_overlap_REVISEDGENELIST.txt", sep ="\t", stringsAsFactors = FALSE, header = FALSE)

####### Peak matching #######
# Parse all the peaks out by sample, then only keep the highest fold change per sample, then match it back to the original dataframe
allPeaks = allPeaksFile %>%
  select(V1, V2, V3, V6, V11, V12, V13, V15, V21)
colnames(allPeaks) = c("chr", "start", "end", "name", "chrPeak",	"startPeak",	"endPeak", "fold_enrichment", "Peakname")

allPeaksWhole <- allPeaks %>%
  dplyr::filter(grepl(pattern = "whole", Peakname)) %>%
  group_by(name) %>%
  arrange(as.numeric(fold_enrichment)) %>%
  dplyr::top_n(1, as.numeric(fold_enrichment))
allPeaksAnt <- allPeaks %>%
  dplyr::filter(grepl(pattern = "Anterior", Peakname)) %>%
  group_by(name) %>%
  arrange(as.numeric(fold_enrichment)) %>%
  dplyr::top_n(1, as.numeric(fold_enrichment))
allPeaksPost <- allPeaks %>%
  dplyr::filter(grepl(pattern = "Posterior", Peakname)) %>%
  group_by(name) %>%
  arrange(as.numeric(fold_enrichment)) %>%
  dplyr::top_n(1, as.numeric(fold_enrichment))
allPeaksNone <- allPeaks %>%
  dplyr::filter(Peakname == ".")

#Add back to allSamples
allSamples$WholePeaks = as.numeric(allPeaksWhole$fold_enrichment)[match(allSamples$Name, allPeaksWhole$name)]
allSamples$AntPeaks = as.numeric(allPeaksAnt$fold_enrichment)[match(allSamples$Name, allPeaksAnt$name)]
allSamples$PostPeaks = as.numeric(allPeaksPost$fold_enrichment)[match(allSamples$Name, allPeaksPost$name)]
allSamples$NoPeak = as.character(allPeaksNone$fold_enrichment)[match(allSamples$Name, allPeaksNone$name)]

################# Filter out regions that were yes and had peaks (that is kind of redundant but whatevs - there are only two positions that were yes but didn't have a peak)
allSamples_yes_peaks = allSamples %>%
  filter(X011617.use == "yes") %>%
  filter(is.na(NoPeak) == TRUE) %>%
    select(X..chr., X.start., X.end., X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_10whole_shifted_lessthan130_10MNorm.bw., X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole.bw., X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole.bw., X.011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.03.1plus2_shifted_lessthan130_1MNorm.bw., X.011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm.bw., Name, New.Location.Assignment, Type) %>%
  mutate(Dotsize = 1.5)
colnames(allSamples_yes_peaks) = c("Chr", "Start", "End", "Whole", "Ant", "LinregPost", "AplusP", "Dnase1", "Name", "Location", "Type", "Dotsize")


AP_Enhancers = allSamples_yes_peaks %>%
  filter(Type == "Enhancer") %>%
  filter(Location %in% c('Anterior', "Posterior", "Mostly Ant", "Mostly Post"))

AP_Promoters = allSamples_yes_peaks %>%
  filter(Type == "Promoter") %>%
  filter(Location %in% c('Anterior', "Posterior", "Mostly Ant", "Mostly Post"))

DV_Enhancers = allSamples_yes_peaks %>%
  filter(Type == "Enhancer") %>%
  filter(Location %in% c('Dorsal', "Ventral"))

DV_Promoters = allSamples_yes_peaks %>%
  filter(Type == "Promoter") %>%
  filter(Location %in% c('Dorsal', "Ventral"))

################# filter our NANs in the 1kb bin file and make standard column names
allSamples_NoNAs = allSamples_1kb_bins %>%
  select(X..chr., X.start., X.end., X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_10whole_shifted_lessthan130_10MNorm.bw., X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole.bw., X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole.bw.,  X.011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.03.1plus2_shifted_lessthan130_1MNorm.bw., X.011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm.bw.) %>%
  filter(X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_10whole_shifted_lessthan130_10MNorm.bw.!= 'nan') %>%
  filter(X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole.bw.!= 'nan') %>%
  filter(X.011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.03.1plus2_shifted_lessthan130_1MNorm.bw.!= 'nan') %>%
  filter(X.011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm.bw.!= 'nan') %>%
  filter(X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole.bw.!= 'nan') %>%
  mutate(Type = "NA") %>%
  mutate(Location = "NA") %>%
  mutate(Name = "NA") %>%
  mutate(Dotsize = 1)
colnames(allSamples_NoNAs) = c("Chr", "Start", "End", "Whole", "Ant", "LinregPost", "AplusP", "Dnase1", "Name", "Location", "Type", "Dotsize")


################# merge the regions with the 1kb windows
MergedALL = rbind(allSamples_NoNAs, allSamples_yes_peaks)


######### ######### ######### AP ENHANCER SCATTERPLOT
AP_Enhancers_windows = MergedALL %>%
  filter(Type %in% c("Enhancer", "NA")) %>%
  filter(Location %in% c("NA", 'Anterior', "Posterior", "Mostly Post"))
AP_Enhancers_windows$Location <- factor(AP_Enhancers_windows$Location,
                                        levels = c("NA", 'Anterior', "Mostly Post", 'Posterior') ,ordered = TRUE)

png('012518_PE_2reps_nodups_AllRegions1KB_RegNormWhole_APEnhancers_Scatter_0118Newgenelist.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(AP_Enhancers_windows,
       aes(x = as.numeric(Ant),y= as.numeric(LinregPost), colour = as.factor(Location))) +
  scale_color_manual(name = 'Legand', values=c("grey70","darkorange1",  "dodgerblue3", "dodgerblue3")) +
  geom_point(alpha = 0.8, aes(size = Dotsize)) +
  scale_size_continuous(range = c(1, 2)) +
  #geom_density2d(colour="lightblue1") + 
  ggtitle('012518_PE_2reps_nodups_AllRegions1KB_RegNormWhole_APEnhancers_Scatter_0118Newgenelist') +
  #scale_x_log10() +
  #scale_y_log10() +
  xlab('Reg Norm Whole Anterior') +
  ylab('Reg Norm Whole Posterior') +
  #geom_smooth(method='lm', formula = y~x) +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  #annotate("text", x=100, y=500, label= Ant_Label_prom, colour = 'darkorange1') +
  #annotate("text", x=500, y=1500, label= Post_Prom_Label, colour = 'dodgerblue3') +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()

######### ######### ######### DV ENHANCER SCATTERPLOT
DV_Enhancers_windows = MergedALL %>%
  filter(Type %in% c("Enhancer", "NA")) %>%
  filter(Location %in% c("NA", 'Dorsal', "Ventral"))
DV_Enhancers_windows$Location <- factor(DV_Enhancers_windows$Location,
                                        levels = c("NA", 'Dorsal', 'Ventral') ,ordered = TRUE)

png('012518_PE_2reps_nodups_AllRegions1KB_DVEnhancers_LinRegWhole_Scatter_0118newgeneList.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(DV_Enhancers_windows,
       aes(x = as.numeric(Ant),y= as.numeric(LinregPost), colour = as.factor(Location))) +
  scale_color_manual(name = 'Legand', values=c( "grey70", "orchid4", "seagreen4")) +
  geom_point(alpha = 0.9, aes(size = Dotsize)) +
  scale_size_continuous(range = c(1, 2)) +
  #geom_density2d(colour="lightblue1") + 
  ggtitle('012518_PE_2reps_nodups_AllRegions1KB_DVEnhancers_LinRegWhole_Scatter_0118newgeneList') +
  #scale_x_log10() +
  #scale_y_log10() +
  xlab('LinregWhole Anterior') +
  ylab('LinregWhole  Posterior') +
  #geom_smooth(method='lm', formula = y~x) +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  #annotate("text", x=100, y=500, label= Ant_Label_prom, colour = 'darkorange1') +
  #annotate("text", x=500, y=1500, label= Post_Prom_Label, colour = 'dodgerblue3') +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20))
dev.off()


######### ######### ######### AP Promoters SCATTERPLOT

AP_Promoters_windows = MergedALL %>%
  filter(Type %in% c("Promoter", "NA")) %>%
  filter(Location %in% c("NA", 'Anterior', "Posterior", "Mostly Ant", "Mostly Post"))
AP_Promoters_windows$Location <- factor(AP_Promoters_windows$Location,
                                        levels = c("NA", 'Anterior', "Mostly Ant", "Mostly Post", 'Posterior') ,ordered = TRUE)


png('012518_PE_2reps_nodups_AllRegions1KB_APPromoters_LinRegWhole_Scatter_0118newgeneList.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(AP_Promoters_windows,
       aes(x = as.numeric(Ant),y= as.numeric(LinregPost), colour = as.factor(Location))) +
  scale_color_manual(name = 'Legand', values=c("grey70","darkorange1" , "darkorange1",  "dodgerblue3", "dodgerblue3")) +
  geom_point(alpha = 1, aes(size = Dotsize)) +
  scale_size_continuous(range = c(1, 2)) +
  #geom_density2d(colour="lightblue1") + 
  ggtitle('012518_PE_2reps_nodups_AllRegions1KB_APPromoters_LinRegWhole_Scatter_0118newgeneList') +
  #scale_x_log10() +
  #scale_y_log10() +
  xlab('LinRegWhole Anterior') +
  ylab('LinRegWhole Posterior') +
  #geom_smooth(method='lm', formula = y~x) +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  #annotate("text", x=100, y=500, label= Ant_Label_prom, colour = 'darkorange1') +
  #annotate("text", x=500, y=1500, label= Post_Prom_Label, colour = 'dodgerblue3') +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()

######### ######### ######### DV Promoters SCATTERPLOT
DV_Promoter_windows = MergedALL %>%
  filter(Type %in% c("Promoter", "NA")) %>%
  filter(Location %in% c("NA", 'Dorsal', "Ventral"))
DV_Promoter_windows$Location <- factor(DV_Promoter_windows$Location,
                                       levels = c("NA", 'Dorsal', 'Ventral') ,ordered = TRUE)

png('012518_PE_2reps_nodups_AllRegions1KB_DVPromoters_linRegWhole_Scatter_0118newGeneList.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(DV_Promoter_windows,
       aes(x = as.numeric(Ant),y= as.numeric(LinregPost), colour = as.factor(Location))) +
  scale_color_manual(name = 'Legand', values=c( "grey70", "orchid4", "seagreen4")) +
  geom_point(alpha = 0.9, aes(size = Dotsize)) +
  scale_size_continuous(range = c(1, 2)) +
  #geom_density2d(colour="lightblue1") + 
  ggtitle('012518_PE_2reps_nodups_AllRegions1KB_DVPromoters_linRegWhole_Scatter_0118newGeneList') +
  #scale_x_log10() +
  #scale_y_log10() +
  xlab('LinregWhole Anterior') +
  ylab('LinRegWhole Posterior') +
  #geom_smooth(method='lm', formula = y~x) +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  #annotate("text", x=100, y=500, label= Ant_Label_prom, colour = 'darkorange1') +
  #annotate("text", x=500, y=1500, label= Post_Prom_Label, colour = 'dodgerblue3') +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20))
dev.off()


######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### #########
######### ######### ######### Bargraphs!!!!


######### ######### ######### AP Enhancer Bargraphs 
#########  Filter out A vs P enhancers to calculate the positional score
allSamples_yes_peaks_ATACSkew = allSamples %>%
  filter(X011617.use == "yes") %>%
  filter(is.na(NoPeak) == TRUE) %>%
  select(X..chr., X.start., X.end., X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_10whole_shifted_lessthan130_10MNorm.bw., X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole.bw., X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole.bw., X.011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.03.1plus2_shifted_lessthan130_1MNorm.bw., X.011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm.bw., Name, New.Location.Assignment, Type, ATACSkewScore) %>%
  mutate(Dotsize = 1.5)
colnames(allSamples_yes_peaks_ATACSkew) = c("Chr", "Start", "End", "Whole", "Ant", "LinregPost", "AplusP", "Dnase1", "Name", "Location", "Type", "ATACSkewScore", "Dotsize")

JustAEnhancers = allSamples_yes_peaks_ATACSkew %>%
  select(Location, Type, Name, Ant, LinregPost, ATACSkewScore) %>%
  filter(Location %in%  c('Anterior')) %>%
  filter(Type == 'Enhancer') %>%
  mutate(PositionalScore = (as.numeric(Ant) - as.numeric(LinregPost)) / (as.numeric(Ant) + as.numeric(LinregPost))) %>%
  distinct()

JustPEnhancers = allSamples_yes_peaks_ATACSkew %>%
  select(Location, Type, Name, Ant, LinregPost, ATACSkewScore) %>%
  filter(Location %in%  c('Posterior', 'Mostly Post')) %>%
  filter(Type == 'Enhancer') %>%
  mutate(PositionalScore = (as.numeric(Ant) - as.numeric(LinregPost)) / (as.numeric(Ant) + as.numeric(LinregPost))) %>%
  distinct()

JustAPEnhancers = bind_rows(JustAEnhancers, JustPEnhancers)

JustAPEnhancers$Location <- factor(JustAPEnhancers$Location,
                                   levels = c('Anterior', 'Posterior', 'Mostly Post'),ordered = TRUE)
ggplot(JustAPEnhancers,
       aes(x = reorder(as.factor(Name),- as.numeric(PositionalScore)),
           y = as.numeric(PositionalScore),
           fill = Location)) +
  #guides(fill= FALSE) +
  geom_col(position = "stack", colour = "black") +
  ylim(-0.6, 0.6) +
  xlab("A-P Enhancers") +
  ylab("PositionalScore") +
  ggtitle("0125_PE_2reps_nodups_LinRegWhole_AP_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList") +
  scale_fill_manual(name = "Legand", values=c("darkorange1", "dodgerblue3", "dodgerblue3")) +
  #geom_abline(slope = 0, intercept=0.2, linetype = 'dashed', size = 1) +
  #geom_abline(slope = 0, intercept=-0.2, linetype = 'dashed', size = 1) +
  #annotate("rect", xmin = 0 , xmax = Inf,   ymin = -0.2, ymax = 0.2,   alpha = 0.2, fill = "slategray4") +
  #scale_fill_gradient(high = "darkorange1", low = "dodgerblue3") +
  #annotate("text", x=500, y=0, label= Ant_Label, colour = 'darkorange1') +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 15),
        axis.title.y = element_text(vjust = 1, size = 15),
        axis.text.x = element_text(size = 7, angle = -90),
        axis.text.y  = element_text(size= 15),
        plot.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12)) +
  png('0125_PE_2reps_nodups_LinRegWhole_AP_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList.png', width = 5000, height = 2000, units = "px",  res=300)
dev.off()


######### ######### ######### AP Promoters Bargraphs 
#########  Filter out A vs P enhancers to calculate the positional score

JustAPromoters = allSamples_yes_peaks_ATACSkew %>%
  select(Location, Type, Name, Ant, LinregPost, ATACSkewScore) %>%
  filter(Location %in%  c('Anterior', 'Mostly Ant')) %>%
  filter(Type == 'Promoter') %>%
  mutate(PositionalScore = (as.numeric(Ant) - as.numeric(LinregPost)) / (as.numeric(Ant) + as.numeric(LinregPost))) %>%
  distinct()

JustPPromoters = allSamples_yes_peaks_ATACSkew %>%
  select(Location, Type, Name, Ant, LinregPost, ATACSkewScore) %>%
  filter(Location %in%  c('Posterior', 'Mostly Post')) %>%
  filter(Type == 'Promoter') %>%
  mutate(PositionalScore = (as.numeric(Ant) - as.numeric(LinregPost)) / (as.numeric(Ant) + as.numeric(LinregPost))) %>%
  distinct()

JustAPPromoters = bind_rows(JustAPromoters, JustPPromoters)

JustAPPromoters$Location <- factor(JustAPPromoters$Location,
                                   levels = c('Anterior', 'Mostly Ant', 'Posterior', 'Mostly Post'),ordered = TRUE)

png('012518_2reps_nodups_AP_Promoters_PositionalScore_Bargraph_Position_0118newgeneList.png', width = 5000, height = 2000, units = "px",  res=300)
ggplot(JustAPPromoters,
       aes(x = reorder(as.factor(Name),- as.numeric(PositionalScore)),
           y = as.numeric(PositionalScore),
           fill = Location)) +
  #guides(fill= FALSE) +
  geom_col(position = "stack", colour = "black") +
  ylim(-0.6, 0.6) +
  xlab("A-P Promoters") +
  ylab("PositionalScore") +
  ggtitle("012518_2reps_nodups_AP_Promoters_PositionalScore_Bargraph_Position_0118newgeneList") +
  scale_fill_manual(name = "Legand", values=c("darkorange1", "darkorange1", "dodgerblue3", "dodgerblue3")) +
  #geom_abline(slope = 0, intercept=0.2, linetype = 'dashed', size = 1) +
  #geom_abline(slope = 0, intercept=-0.2, linetype = 'dashed', size = 1) +
  #annotate("rect", xmin = 0 , xmax = Inf,   ymin = -0.2, ymax = 0.2,   alpha = 0.2, fill = "slategray4") +
  #scale_fill_gradient(high = "darkorange1", low = "dodgerblue3") +
  #annotate("text", x=500, y=0, label= Ant_Label, colour = 'darkorange1') +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 15),
        axis.title.y = element_text(vjust = 1, size = 15),
        axis.text.x = element_text(size = 7, angle = -90),
        axis.text.y  = element_text(size= 15),
        plot.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))

dev.off()

######### ######### ######### DV Enhancer Bargraphs 
JustDorsalEnhancers = allSamples_yes_peaks_ATACSkew %>%
  select(Location, Type, Name, Ant, LinregPost, ATACSkewScore) %>%
  filter(Location == "Dorsal") %>%
  filter(Type == 'Enhancer') %>%
  mutate(PositionalScore = (as.numeric(Ant) - as.numeric(LinregPost)) / (as.numeric(Ant) + as.numeric(LinregPost))) %>%
  distinct()

JustVentralEnhancers = allSamples_yes_peaks_ATACSkew %>%
  select(Location, Type, Name, Ant, LinregPost, ATACSkewScore) %>%
  filter(Location == "Ventral") %>%
  filter(Type == 'Enhancer') %>%
  mutate(PositionalScore = (as.numeric(Ant) - as.numeric(LinregPost)) / (as.numeric(Ant) + as.numeric(LinregPost))) %>%
  distinct()

JustDVEnhancers = bind_rows(JustDorsalEnhancers, JustVentralEnhancers)

JustDVEnhancers$Location <- factor(JustDVEnhancers$Location,
                                   levels = c("Dorsal", "Ventral"),ordered = TRUE)

ggplot(JustDVEnhancers,
       aes(x = reorder(as.factor(Name),- as.numeric(PositionalScore)),
           y = as.numeric(PositionalScore),
           fill = Location)) +
  #guides(fill= FALSE) +
  geom_col(position = "stack", colour = "black") +
  ylim(-0.6, 0.6) +
  xlab("D-V Enhancers") +
  ylab("PositionalScore") +
  ggtitle("012518_PE_2reps_nodups_LinRegWhole_DV_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList") +
  scale_fill_manual(name = "Legand", values=c("orchid4", "seagreen4")) +
  #geom_abline(slope = 0, intercept=0.2, linetype = 'dashed', size = 1) +
  #geom_abline(slope = 0, intercept=-0.2, linetype = 'dashed', size = 1) +
  #annotate("rect", xmin = 0 , xmax = Inf,   ymin = -0.2, ymax = 0.2,   alpha = 0.2, fill = "slategray4") +
  #scale_fill_gradient(high = "darkorange1", low = "dodgerblue3") +
  #annotate("text", x=500, y=0, label= Ant_Label, colour = 'darkorange1') +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 15),
        axis.title.y = element_text(vjust = 1, size = 15),
        axis.text.x = element_text(size = 7, angle = -90),
        axis.text.y  = element_text(size= 15),
        plot.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12)) +
  png('012518_PE_2reps_nodups_LinRegWhole_DV_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList.png', width = 5000, height = 2000, units = "px",  res=300)
dev.off()
######### ######### ######### DV Promoter Bargraphs 

JustDorsalPromoter = allSamples_yes_peaks_ATACSkew %>%
  select(Location, Type, Name, Ant, LinregPost, ATACSkewScore) %>%
  filter(Location == "Dorsal") %>%
  filter(Type == 'Promoter') %>%
  mutate(PositionalScore = (as.numeric(Ant) - as.numeric(LinregPost)) / (as.numeric(Ant) + as.numeric(LinregPost))) %>%
  distinct()

JustVentralPromoter = allSamples_yes_peaks_ATACSkew %>%
  select(Location, Type, Name, Ant, LinregPost, ATACSkewScore) %>%
  filter(Location == "Ventral") %>%
  filter(Type == 'Promoter') %>%
  mutate(PositionalScore = (as.numeric(Ant) - as.numeric(LinregPost)) / (as.numeric(Ant) + as.numeric(LinregPost))) %>%
  distinct()

JustDVPromoter = bind_rows(JustDorsalPromoter, JustVentralPromoter)

JustDVPromoter$Location <- factor(JustDVPromoter$Location,
                                  levels = c("Dorsal", "Ventral"),ordered = TRUE)

ggplot(JustDVPromoter,
       aes(x = reorder(as.factor(Name),- as.numeric(PositionalScore)),
           y = as.numeric(PositionalScore),
           fill = Location)) +
  #guides(fill= FALSE) +
  geom_col(position = "stack", colour = "black") +
  ylim(-0.6, 0.6) +
  xlab("D-V Promoters") +
  ylab("PositionalScore") +
  ggtitle("012518_2reps_nodups_LinRegWhole_DV_Promoters_PositionalScore_Bargraph_Position_0118newgeneList") +
  scale_fill_manual(name = "Legand", values=c("orchid4", "seagreen4")) +
  #geom_abline(slope = 0, intercept=0.2, linetype = 'dashed', size = 1) +
  #geom_abline(slope = 0, intercept=-0.2, linetype = 'dashed', size = 1) +
  #annotate("rect", xmin = 0 , xmax = Inf,   ymin = -0.2, ymax = 0.2,   alpha = 0.2, fill = "slategray4") +
  #scale_fill_gradient(high = "darkorange1", low = "dodgerblue3") +
  #annotate("text", x=500, y=0, label= Ant_Label, colour = 'darkorange1') +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 15),
        axis.title.y = element_text(vjust = 1, size = 15),
        axis.text.x = element_text(size = 7, angle = -90),
        axis.text.y  = element_text(size= 15),
        plot.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12)) +
  png('012518_2reps_nodups_LinRegWhole_DV_Promoters_PositionalScore_Bargraph_Position_0118newgeneList.png', width = 5000, height = 2000, units = "px",  res=300)
dev.off()

######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### #########
######### ######### ######### Boxplots!!!!


######### ######### Enhancers

allSamples$Sample <- seq.int(nrow(allSamples)) # Add an sampleID row at the end

RandomRegions = allSamples %>%
  select(Sample, X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole.bw..1, X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole.bw..1, New.Location.Assignment, Type, RandSkewScore) %>%
  filter(Type =="Enhancer") %>%
  mutate(Location.real = "Random") %>%
  mutate(Type.real = "Random") %>%
  mutate(TotalScore = as.numeric(X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole.bw..1) + as.numeric(X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole.bw..1)) %>%
  select(Sample, Location.real, Type.real, RandSkewScore, TotalScore)
colnames(RandomRegions) = c("Name", "Location", "Type", "ATACSkewScore", "TotalScore")
RandomRegions$Name = as.character(RandomRegions$Name)

JustAEnhancers_merged = allSamples_yes_peaks_ATACSkew %>%
  filter(Location %in%  c('Anterior')) %>%
  filter(Type == 'Enhancer') %>%
  mutate(TotalScore = as.numeric(Ant) + as.numeric(LinregPost)) %>%
  mutate(merged_location = "Anterior") %>%
  select(Name, merged_location, Type, ATACSkewScore, TotalScore) %>%
  distinct()

JustPEnhancers_merged = allSamples_yes_peaks_ATACSkew %>%
  filter(Location %in%  c('Posterior', 'Mostly Post')) %>%
  filter(Type == 'Enhancer') %>%
  mutate(TotalScore = as.numeric(Ant) + as.numeric(LinregPost)) %>%
  mutate(merged_location = "Posterior") %>%
  select(Name, merged_location, Type, ATACSkewScore, TotalScore) %>%
  distinct()

JustDorsalEnhancers_merged = allSamples_yes_peaks_ATACSkew %>%
  filter(Location == "Dorsal") %>%
  filter(Type == 'Enhancer') %>%
  mutate(TotalScore = as.numeric(Ant) + as.numeric(LinregPost)) %>%
  mutate(merged_location = "Dorsal") %>%
  select(Name, merged_location, Type, ATACSkewScore, TotalScore) %>%
  distinct()

JustVentralEnhancers_merged = allSamples_yes_peaks_ATACSkew %>%
  filter(Location == "Ventral") %>%
  filter(Type == 'Enhancer') %>%
  mutate(TotalScore = as.numeric(Ant) + as.numeric(LinregPost)) %>%
  mutate(merged_location = "Ventral") %>%
  select(Name, merged_location, Type, ATACSkewScore, TotalScore) %>%
  distinct()

colnames(JustAEnhancers_merged) = c("Name", "Location", "Type", "ATACSkewScore", "TotalScore")
colnames(JustPEnhancers_merged  ) = c("Name", "Location", "Type", "ATACSkewScore", "TotalScore")
colnames(JustDorsalEnhancers_merged) = c("Name", "Location", "Type", "ATACSkewScore", "TotalScore")
colnames(JustVentralEnhancers_merged) = c("Name", "Location", "Type", "ATACSkewScore", "TotalScore")

rand_and_Enhancers = bind_rows(RandomRegions, JustAEnhancers_merged, JustPEnhancers_merged, JustDorsalEnhancers_merged, JustVentralEnhancers_merged)
rand_and_Enhancers$Location <- factor(rand_and_Enhancers$Location,
                                      levels = c('Anterior',  'Posterior', 'Dorsal', 'Ventral', 'Random') ,ordered = TRUE)


ggplot(rand_and_Enhancers,
       aes(x = as.factor(Location),
           y = as.numeric(ATACSkewScore),
           fill = Location)) +
  #guides(fill= FALSE) +
  geom_boxplot(outlier.size = 0) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.6 , alpha = 0.3, binpositions="all") +
  xlab("Location") +
  ylab("ATAC Skew Score") +
  coord_cartesian(ylim = c(-0.5, 0.7)) +
  geom_abline(slope = 0, intercept=0, linetype = 'dotted') +
  ggtitle("012518_PE_2reps_nodups_LinRegWhole_AllLocations_Enhancers_ATACSkewScore_WithRandom_0118newgeneList") +
  scale_fill_manual(values=c("darkorange1", "dodgerblue3", "orchid4", "seagreen4", "grey50")) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.position="none") +
  png('012518_PE_2reps_nodups_LinRegWhole_AllLocations_Enhancers_ATACSkewScore_WithRandom_0118newgeneList.png', width = 2000, height = 2000, units = "px",  res=300)
dev.off()

################# ANOVAS for Enhancers  ################# 
Anova_Enhancers_dataframe = rand_and_Enhancers %>%
  select(Location, ATACSkewScore)

Anova_Enhancers_lm = lm(as.numeric(Anova_Enhancers_dataframe$ATACSkewScore) ~ Anova_Enhancers_dataframe$Location)
Anova_Enhancer = aov(formula = as.numeric(Anova_Enhancers_dataframe$ATACSkewScore) ~ Anova_Enhancers_dataframe$Location)
Tukey_Enhancers = TukeyHSD(Anova_Enhancer)

write.table(Tukey_Enhancers$`Anova_Enhancers_dataframe$Location`, "012518_Anova_ALLLocations_Enhancers.txt")
# Promoters

######### ######### Promoters
write.csv(allSamples, file = "012518_PE_2reps_nodups_randomRegv2_ouput_0118newgenelist_withpeaks.csv")


############### ################ ################ ################ ################
################ RNA-seq MATCHUP     ################ 
######### ######### ######### Input ######### ######### ######### ######### ######### ######### 


Susan_RNAseq_data <- read.delim2("~/Box Sync/Eisen_Lab/Experiments/ATAC-seq/Halves/ATAC-seq_Pools/040617_Analysis/061417_BrowserTraces/062317_RegressionNormalization/081617_2reps/Lott2011_RNA-seq_GeneExpression.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(Susan_RNAseq_data)[1] = "Name"
Joined_ATACPeaks_plus_RNAseq = full_join(allSamples, Susan_RNAseq_data, by = "Name")

JustPromoters = Joined_ATACPeaks_plus_RNAseq %>%
  filter(Type == "Promoter") %>%
  select(Name, Type, New.Location.Assignment,X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole.bw. ,X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole.bw. , ATACSkewScore, PValue, WholePeaks, AntPeaks, PostPeaks, NoPeak, CLASS, F10, F11, F12, U13, F14A, U14A, F14B, F14B_r2, F14C, F14C_r2, F14D, U14D)

write.csv(JustPromoters, file = "012518_Promoters_peaks_RNAseq.csv")

#allSamples <-read.delim2("100617_Promoters_peaks_RNAseq.csv", sep = ",", stringsAsFactors = FALSE)

#boxplot to show whether zygotic genes are more differentially accessible than maternal genes

JustAProm_merged = JustPromoters %>%
  filter(New.Location.Assignment %in%  c('Anterior', 'Mostly Ant')) %>%
  filter(Type == 'Promoter') %>%
  #mutate(TotalScore = as.numeric(Ant) + as.numeric(LinregPost)) %>%
  mutate(merged_location = "Anterior") %>%
  select(Name, merged_location, Type, ATACSkewScore, CLASS) %>%
  distinct()

JustP_prom_merged = JustPromoters %>%
  filter(New.Location.Assignment %in%  c('Posterior', 'Mostly Post')) %>%
  filter(Type == 'Promoter') %>%
  #mutate(TotalScore = as.numeric(Ant) + as.numeric(LinregPost)) %>%
  mutate(merged_location = "Posterior") %>%
  select(Name, merged_location, Type, ATACSkewScore, CLASS) %>%
  distinct()

JustDorsal_Prom_merged = JustPromoters %>%
  filter(New.Location.Assignment == "Dorsal") %>%
  filter(Type == 'Promoter') %>%
  #mutate(TotalScore = as.numeric(Ant) + as.numeric(LinregPost)) %>%
  mutate(merged_location = "Dorsal") %>%
  select(Name, merged_location, Type, ATACSkewScore, CLASS) %>%
  distinct()

JustVentral_Prom_merged = JustPromoters %>%
  filter(New.Location.Assignment == "Ventral") %>%
  filter(Type == 'Promoter') %>%
  #mutate(TotalScore = as.numeric(Ant) + as.numeric(LinregPost)) %>%
  mutate(merged_location = "Ventral") %>%
  select(Name, merged_location, Type, ATACSkewScore, CLASS) %>%
  distinct()

colnames(JustAProm_merged) = c("Name", "Location", "Type", "ATACSkewScore", "CLASS")
colnames(JustP_prom_merged  ) = c("Name", "Location", "Type", "ATACSkewScore", "CLASS")
colnames(JustDorsal_Prom_merged) = c("Name", "Location", "Type", "ATACSkewScore", "CLASS")
colnames(JustVentral_Prom_merged) = c("Name", "Location", "Type", "ATACSkewScore", "CLASS")

rand_and_Promoters = bind_rows(JustAProm_merged, JustP_prom_merged, JustDorsal_Prom_merged, JustVentral_Prom_merged)
rand_and_Promoters$Location <- factor(rand_and_Promoters$Location,
                                      levels = c('Anterior', 'Posterior', 'Dorsal', 'Ventral') ,ordered = TRUE)

ggplot(rand_and_Promoters,
       aes(x = as.factor(CLASS),
           y = as.numeric(ATACSkewScore),
           fill = Location)) +
  #guides(fill= FALSE) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = as.numeric(rand_and_Enhancers$TotalScore)/3000, alpha = 0.3, binpositions="all") +
  geom_boxplot(outlier.size = 0) +
  facet_wrap(~Location) +
  xlab("Class") +
  ylab("ATAC Skew Score") +
  coord_cartesian(ylim = c(-0.5, 0.7)) +
  geom_abline(slope = 0, intercept=0, linetype = 'dotted') +
  ggtitle("012518_MatvZygotic_ATACSkewScore_Promoters") +
  scale_fill_manual(values=c("darkorange1", "dodgerblue3", "orchid4", "seagreen4", "grey50")) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=15),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.position="none") +
  png('012518_MatvZygotic_ATACSkewScore_Promoters.png', width = 2000, height = 2000, units = "px",  res=300)
dev.off()

####### ANOVAS to determine if maternal and zygotic promoters have different skew scores.
Anova_AntPromoters_dataframe = JustAProm_merged %>%
  select(ATACSkewScore, CLASS)
Anova_AntPromoters = aov(formula = as.numeric(Anova_AntPromoters_dataframe$ATACSkewScore) ~ Anova_AntPromoters_dataframe$CLASS)
Tukey_AntPromoters = TukeyHSD(Anova_AntPromoters)

Anova_PostPromoters_dataframe = JustP_prom_merged %>%
  select(ATACSkewScore, CLASS)
Anova_PostPromoters = aov(formula = as.numeric(Anova_PostPromoters_dataframe$ATACSkewScore) ~ Anova_PostPromoters_dataframe$CLASS)
Tukey_PostPromoters = TukeyHSD(Anova_PostPromoters)

Anova_DorPromoters_dataframe = JustDorsal_Prom_merged %>%
  select(ATACSkewScore, CLASS)
Anova_DorPromoters = aov(formula = as.numeric(JustDorsal_Prom_merged$ATACSkewScore) ~ JustDorsal_Prom_merged$CLASS)
Tukey_DorPromoters = TukeyHSD(Anova_DorPromoters)

Anova_VentralPromoters_dataframe = JustVentral_Prom_merged %>%
  select(ATACSkewScore, CLASS)
Anova_VenPromoters = aov(formula = as.numeric(JustVentral_Prom_merged$ATACSkewScore) ~ JustVentral_Prom_merged$CLASS)
Tukey_VenPromoters = TukeyHSD(Anova_VenPromoters)
AntpromANOVA = as.data.frame(Tukey_AntPromoters$`Anova_AntPromoters_dataframe$CLASS`)
Joined_Anovas = rbind(as.data.frame(Tukey_AntPromoters$`Anova_AntPromoters_dataframe$CLASS`, ), as.data.frame(Tukey_PostPromoters$`Anova_PostPromoters_dataframe$CLASS`), as.data.frame(Tukey_DorPromoters$`JustDorsal_Prom_merged$CLASS`), as.data.frame(Tukey_VenPromoters$`JustVentral_Prom_merged$CLASS`))
write.table(Joined_Anovas, "012518_Anova_RNAseq_Promoters.txt")

############## Posterior Promoter Classes and Random Regions ############## 

RandomRegions = allSamples %>%
  select(Sample, X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole.bw..1,X.011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole.bw..1, New.Location.Assignment, Type, RandSkewScore) %>%
  filter(Type =="Promoter") %>%
  mutate(Location.real = "Random") %>%
  mutate(Type.real = "Random") %>%
  mutate(CLASS = "Random") %>%
  select(Sample, Location.real, Type.real, RandSkewScore, CLASS)
colnames(RandomRegions) = c("Name", "Location", "Type", "ATACSkewScore", "CLASS")
RandomRegions$Name = as.character(RandomRegions$Name)

rand_and_PostPromoters = rbind(JustP_prom_merged, RandomRegions)
rand_and_PostPromoters$CLASS <- factor(rand_and_PostPromoters$CLASS,
                                       levels = c('mat', 'matzyg', 'zyg', 'Random') ,ordered = TRUE)

ggplot(rand_and_PostPromoters,
       aes(x = as.factor(CLASS),
           y = as.numeric(ATACSkewScore),
           fill = Location)) +
  #guides(fill= FALSE) +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize = as.numeric(rand_and_Enhancers$TotalScore)/3000, alpha = 0.3, binpositions="all") +
  geom_boxplot(outlier.size = 0) +
  #facet_wrap(~Location) +
  xlab("Class") +
  ylab("ATAC Skew Score") +
  coord_cartesian(ylim = c(-0.5, 0.7)) +
  geom_abline(slope = 0, intercept=0, linetype = 'dotted') +
  ggtitle("012518_PE_MatvZygotic_ATACSkewScore_Promoters_POSTVsRANDOM") +
  scale_fill_manual(values=c("dodgerblue3", "grey50")) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=15),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.position="none") +
  png('012518_PE_MatvZygotic_ATACSkewScore_Promoters_POSTVsRANDOM.png', width = 2000, height = 2000, units = "px",  res=300)
dev.off()


############ ANOVAS for Posterior and Random ############ 
Anova_PostPromotersRand_dataframe = rand_and_PostPromoters %>%
  select(ATACSkewScore, CLASS)
Anova_PostvRandPromoters = aov(formula = as.numeric(Anova_PostPromotersRand_dataframe$ATACSkewScore) ~ Anova_PostPromotersRand_dataframe$CLASS)
Tukey_PostvRandPromoters = TukeyHSD(Anova_PostvRandPromoters)

write.table(Tukey_PostvRandPromoters$`Anova_PostPromotersRand_dataframe$CLASS`, "012518_Anova_RNAseq_PromotersPostvRand.txt")




##############################################################################################################################
#New promoter boxplot figure


JustAProm_merged_ZYG = JustAProm_merged %>%
  filter(CLASS == "zyg")

JustP_prom_merged_ZYG = JustP_prom_merged %>%
  filter(CLASS == "zyg")

JustDorsal_Prom_merged_ZYG = JustDorsal_Prom_merged %>%
  filter(CLASS == "zyg")
  
JustVentral_Prom_merged_ZYG = JustVentral_Prom_merged %>%
  filter(CLASS == "zyg")


rand_and_PromotersZYG = bind_rows(RandomRegions, JustAProm_merged_ZYG, JustP_prom_merged_ZYG, JustDorsal_Prom_merged_ZYG, JustVentral_Prom_merged_ZYG)
rand_and_PromotersZYG$Location <- factor(rand_and_PromotersZYG$Location,
                                      levels = c('Anterior', 'Posterior', 'Dorsal', 'Ventral', 'Random') ,ordered = TRUE)

ggplot(rand_and_PromotersZYG,
       aes(x = as.factor(Location),
           y = as.numeric(ATACSkewScore),
           fill = Location)) +
  #guides(fill= FALSE) +
  geom_boxplot(outlier.size = 0) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.6 , alpha = 0.3, binpositions="all") +
  xlab("Location") +
  ylab("ATAC Skew Score") +
  coord_cartesian(ylim = c(-0.5, 0.7)) +
  geom_abline(slope = 0, intercept=0, linetype = 'dotted') +
  ggtitle("012518_PE_2reps_nodups_LinRegWhole_ALLLocations_Promoters_ATACSkewScore_WithRandom_0118newgeneList_Zyggenesonly") +
  scale_fill_manual(values=c("darkorange1", "dodgerblue3", "orchid4", "seagreen4", "grey50")) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.position="none") +
  png('012518_PE_2reps_nodups_LinRegWhole_ALLLocations_Promoters_ATACSkewScore_WithRandom_0118newgeneList_Zyggenesonly.png', width = 2000, height = 2000, units = "px",  res=300)
dev.off()

################# ANOVAS for BoxPlot Zyg only promoters  ################# 
Anova_Promoters_dataframe = rand_and_PromotersZYG %>%
  select(Location, ATACSkewScore)

Anova_Promoters_lm = lm(as.numeric(Anova_Promoters_dataframe$ATACSkewScore) ~ Anova_Promoters_dataframe$Location)
Anova_Promoter = aov(formula = as.numeric(Anova_Promoters_dataframe$ATACSkewScore) ~ Anova_Promoters_dataframe$Location)
Tukey_Promoters = TukeyHSD(Anova_Promoter)
write.table(Tukey_Promoters$`Anova_Promoters_dataframe$Location`, "012518_Anova_ALLLocations_Promoters_Zygonly.txt")

######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### #########
######### ######### ######### Correlation coeefficients for Figure 3AB and Figure 4 AB


######### ######### Figure 3AB
AenhCor = cor(as.numeric(JustAEnhancers$Ant), as.numeric(JustAEnhancers$LinregPost), method = "spearman")
PenhCor = cor(as.numeric(JustPEnhancers$Ant), as.numeric(JustPEnhancers$LinregPost), method = "spearman")
VenhCor = cor(as.numeric(JustVentralEnhancers$Ant), as.numeric(JustVentralEnhancers$LinregPost), method = "spearman")
DenhCor = cor(as.numeric(JustDorsalEnhancers$Ant), as.numeric(JustDorsalEnhancers$LinregPost), method = "spearman")

######### ######### Figure 4AB
APromCor = cor(as.numeric(JustAPromoters$Ant), as.numeric(JustAPromoters$LinregPost), method = "spearman")
PPromCor = cor(as.numeric(JustPPromoters$Ant), as.numeric(JustPPromoters$LinregPost), method = "spearman")
VPromCor = cor(as.numeric(JustVentralPromoter$Ant), as.numeric(JustVentralPromoter$LinregPost), method = "spearman")
DPromCor = cor(as.numeric(JustDorsalPromoter$Ant), as.numeric(JustDorsalPromoter$LinregPost), method = "spearman")

######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### #########
######### ######### ######### Histograms of random regions and ATACskewscore for each group to check the random distribution script

justallEnhancers = rbind(JustAPEnhancers,JustDVEnhancers )

rand_and_AP_Enhancers = rand_and_Enhancers %>%
  filter(Location %in%  c('Anterior', 'Posterior', 'Random'))

rand_and_DV_Enhancers = rand_and_Enhancers %>%
  filter(Location %in%  c('Dorsal', 'Ventral', 'Random'))

png('012518_PE_histogram_density_enhancers.png', width = 2000, height = 2000, units = "px",  res=300)
ggplot(rand_and_AP_Enhancers,
       aes(x = as.numeric(ATACSkewScore), fill = as.factor(Location), colour = as.factor(Location))) +
  #geom_histogram(aes(y=..density..), alpha = 0.3 ) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values=c("darkorange1", "dodgerblue3", "grey50")) +
  scale_colour_manual(values=c("darkorange1", "dodgerblue3","grey50")) +
  ggtitle('012518_PE_histogram_density_enhancers') +
  xlab('ATACSkewScore') +
  #coord_equal() +
  xlim(-1,1) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.position="none")
dev.off()

png('012518_PE_histogram_density_DVenhancers.png', width = 2000, height = 2000, units = "px",  res=300)
ggplot(rand_and_DV_Enhancers,
       aes(x = as.numeric(ATACSkewScore), fill = as.factor(Location), colour = as.factor(Location))) +
  #geom_histogram(aes(y=..density..), alpha = 0.3 ) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values=c("orchid4", "seagreen4", "grey50")) +
  scale_colour_manual(values=c("orchid4", "seagreen4", "grey50")) +
  ggtitle('012518_PE_histogram_density_DVenhancers') +
  xlab('ATACSkewScore') +
  #coord_equal() +
  xlim(-1,1) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.position="none")
dev.off()

rand_and_AP_Promoters = rand_and_Promoters %>%
  filter(Location %in%  c('Anterior', 'Posterior', 'Random'))

rand_and_DV_Promoters = rand_and_Promoters %>%
  filter(Location %in%  c('Dorsal', 'Ventral', 'Random'))

png('012518_PE_histogram_density_APPromoters.png', width = 2000, height = 2000, units = "px",  res=300)
ggplot(rand_and_AP_Promoters,
       aes(x = as.numeric(ATACSkewScore), fill = as.factor(Location), colour = as.factor(Location))) +
  #geom_histogram(aes(y=..density..), alpha = 0.3 ) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values=c("darkorange1", "dodgerblue3", "grey50")) +
  scale_colour_manual(values=c("darkorange1", "dodgerblue3","grey50")) +
  ggtitle('012518_PE_histogram_density_APPromoters') +
  xlab('ATACSkewScore') +
  #coord_equal() +
  xlim(-1,1) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.position="none")
dev.off()

png('012518_PE_histogram_density_DVPromoters.png', width = 2000, height = 2000, units = "px",  res=300)
ggplot(rand_and_DV_Promoters,
       aes(x = as.numeric(ATACSkewScore), fill = as.factor(Location), colour = as.factor(Location))) +
  #geom_histogram(aes(y=..density..), alpha = 0.3 ) +
  geom_density(alpha = 0.3) +
  scale_fill_manual(values=c("orchid4", "seagreen4", "grey50")) +
  scale_colour_manual(values=c("orchid4", "seagreen4", "grey50")) +
  ggtitle('012518_PE_histogram_density_DVPromoters') +
  xlab('ATACSkewScore') +
  #coord_equal() +
  xlim(-1,1) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.position="none")
dev.off()

######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### #########
######### ######### ######### Figure 1 Scatterplots of 1kb regions


######### ######### ######### Ant

### Correction new graph for figure 1 that is just AvsP
AvsP_cor = cor(as.numeric(allSamples_NoNAs$Ant),y= as.numeric(allSamples_NoNAs$LinregPost), method = "spearman")
png('012518_PE2reps_nodups_LinRegWhole_AllRegions1KB_AvsP_nolog.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(LinregPost),y= as.numeric(Ant), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  #scale_size_continuous(range = c(1, 2)) +
  geom_point(alpha = 0.05, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightskyblue2", alpha = 0.5, size = 0.4) + 
  ggtitle('012518_PE_2reps_nodups_LinRegWhole_AllRegions1KB_AvsP') +
  #scale_x_log10(limits = c(10,1500)) +
  #scale_y_log10(limits = c(10,1500)) +
  xlab('Reg Norm Whole Posterior') +
  ylab('Reg Norm Whole Anterior') +
  #geom_smooth(method='lm', formula = y~x) +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=1000, y=20, label= AvsP_cor) +
  #annotate("text", x=500, y=1500, label= Post_Prom_Label, colour = 'dodgerblue3') +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()


AvsWhole_cor = cor(as.numeric(allSamples_NoNAs$Ant),y= as.numeric(allSamples_NoNAs$Whole), method = "spearman")
png('012518_PE_2reps_nodups_LinRegWhole_AvsWhole_nolog.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Ant),y= as.numeric(Whole), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.05, size = 0.3, show.legend = FALSE) +
  #scale_size_continuous(range = c(1, 2)) +
  geom_density2d(colour="lightblue1", alpha = 0.5, size = 0.4) + 
  ggtitle('012518_PE_2reps_nodups_LinRegWhole_AvsWhole') +
  #scale_x_log10(limit = c(10,1500)) +
  #scale_y_log10(limit = c(10,1500)) +
  xlab('LinRegWhole Anterior') +
  ylab('LinRegWhole Whole') +
  #geom_smooth(method='lm', formula = y~x) +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  #annotate("text", x=20, y=1000, label= AvsWhole_cor) +
  #annotate("text", x=500, y=1500, label= Post_Prom_Label, colour = 'dodgerblue3') +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()


AvsAP_cor = cor(as.numeric(allSamples_NoNAs$Ant),y= as.numeric(allSamples_NoNAs$AplusP), method = "spearman")
png('012518_PE_2reps_nodups_LinRegWhole_AvsAP.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Ant),y= as.numeric(AplusP), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('012518_PE_2reps_nodups_LinRegWhole_AvsAP') +
  xlab('LinregWhole Anterior') +
  ylab('LinregWhole AplusP') +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=1500, y=500, label= AvsAP_cor) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()

AvsDnaseI_cor = cor(as.numeric(allSamples_NoNAs$Ant),y= as.numeric(allSamples_NoNAs$Dnase1), method = "spearman")
png('012518_PE_2reps_nodups_LinRegWhole_AvsDnaseI.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Dnase1),y= as.numeric(Ant), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('012518_PE_2reps_nodups_LinRegWhole_AvsDnaseI') +
  xlab('RegNormWhole DnaseI') +
  ylab('RegNormWhole Anterior') +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  #annotate("text", x=1500, y=500, label= AvsDnaseI_cor) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()


######### ######### ######### Posterior

PvsWhole_cor = cor(as.numeric(allSamples_NoNAs$LinregPost),y= as.numeric(allSamples_NoNAs$Whole), method = "spearman")
png('012518_PE_2reps_nodups_LinRegWhole_PostvsWhole.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(LinregPost),y= as.numeric(Whole), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('012518_PE_2reps_nodups_LinRegWhole_PostvsWhole') +
  xlab('LinregWhole Post') +
  ylab('Whole') +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=1500, y=500, label= PvsWhole_cor) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()

PvsAP_cor = cor(as.numeric(allSamples_NoNAs$LinregPost),y= as.numeric(allSamples_NoNAs$AplusP), method = "spearman")
png('012518_PE_2reps_nodups_LinRegWhole_PostvsAplusP.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(LinregPost),y= as.numeric(AplusP), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('012518_PE_2reps_nodups_LinRegWhole_PostvsAplusP') +
  xlab('LinregPost') +
  ylab('LinregWhole AplusP') +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=1500, y=500, label= PvsAP_cor) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()

PvsDnaseI_cor = cor(as.numeric(allSamples_NoNAs$LinregPost),y= as.numeric(allSamples_NoNAs$Dnase1), method = "spearman")
png('012518_PE_2reps_nodups_LinRegWhole_PostvsDnaseI.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Dnase1),y= as.numeric(LinregPost), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('012518_PE_2reps_nodups_LinRegWhole_PostvsDnaseI') +
  xlab('DnaseI') +
  ylab('LinregWhole LinregPost') +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  #annotate("text", x=1500, y=500, label= PvsDnaseI_cor) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()

######### ######### ######### Whole


WholevsDnaseI_cor = cor(as.numeric(allSamples_NoNAs$Whole),y= as.numeric(allSamples_NoNAs$Dnase1), method = "spearman")
png('012518_PE_2reps_nodups_LinRegWhole_WholevsDnaseI.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Dnase1),y= as.numeric(Whole), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('012518_PE_2reps_nodups_LinRegWhole_WholevsDnaseI') +
  xlab('Dnase1') +
  ylab('Whole') +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  #annotate("text", x=1500, y=500, label= WholevsDnaseI_cor) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()

WholevsAplusP_cor = cor(as.numeric(allSamples_NoNAs$Whole),y= as.numeric(allSamples_NoNAs$AplusP), method = "spearman")
png('012518_PE_2reps_nodups_LinRegWhole_WholevsAplusP.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Whole),y= as.numeric(AplusP), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('012518_PE_2reps_nodups_LinRegWhole_WholevsAplusP') +
  xlab('Whole') +
  ylab('RegNormWhole AplusP') +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=1500, y=500, label= WholevsAplusP_cor) +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill= NA),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(vjust = 0, size = 20),
        axis.title.y = element_text(vjust = 1, size = 20),
        axis.text.x = element_text(size=20),
        axis.text.y  = element_text(size=20),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12))
dev.off()






