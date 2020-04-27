##############################################################################################
################ LinearregNormalization.r                                     ################ 
################ 011818_RefinedGeneList_Dataanalysis_linregv2_PEReviews.r     ################ 
################ make scatterplots with ATAC Score and position               ################ 
################ calculated in 051117_randomdistribution_script                ################
##############################################################################################


library("plyr")
library("dplyr")
library("ggplot2")


############################ Linear Regression Normalization ############################ 

###### ### ###  Working directory
setwd("~/011618_PairedEndAnalysisforReviews/")

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

### calculate linear regression model
merged_norm_1kb_AntvsWhole.lm = lm(as.numeric(allSamples_NoNAs$Ant) ~ as.numeric(allSamples_NoNAs$Whole))
merged_norm_1kb_PostvsWhole.lm = lm(as.numeric(allSamples_NoNAs$Post) ~ as.numeric(allSamples_NoNAs$Whole))

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

#### Correlations 
AP_NoNormCor = cor(as.numeric(allSamples_NoNAs$Ant),as.numeric(allSamples_NoNAs$Post), method = "spearman")
AP_NoNormCor_P = cor(as.numeric(allSamples_NoNAs$Ant),as.numeric(allSamples_NoNAs$Post), method = "pearson")
AP_NoNormCor_P_sq = (AP_NoNormCor_P)^2
AP_NoNormCor_sq = (AP_NoNormCor)^2

AW_NoNormCor = cor(as.numeric(allSamples_NoNAs$Ant),as.numeric(allSamples_NoNAs$Whole), method = "spearman")
AW_NoNormCor_P = cor(as.numeric(allSamples_NoNAs$Ant),as.numeric(allSamples_NoNAs$Whole), method = "pearson")
AW_NoNormCor_P_sq = (AW_NoNormCor_P)^2
AW_NoNormCor_sq = (AW_NoNormCor)^2

PW_NoNormCor = cor(as.numeric(allSamples_NoNAs$Post),as.numeric(allSamples_NoNAs$Whole), method = "spearman")
PW_NoNormCor_P = cor(as.numeric(allSamples_NoNAs$Post),as.numeric(allSamples_NoNAs$Whole), method = "pearson")
PW_NoNormCor_P_sq = (PW_NoNormCor_P)^2
PW_NoNormCor_sq = (PW_NoNormCor)^2

AW_NormCor = cor((as.numeric(allSamples_NoNAs$Ant)/ as.numeric(merged_norm_1kb_AntvsWhole.lm$coefficients[2])),as.numeric(allSamples_NoNAs$Whole), method = "spearman")
AW_NormCor_P = cor((as.numeric(allSamples_NoNAs$Ant)/ as.numeric(merged_norm_1kb_AntvsWhole.lm$coefficients[2])),as.numeric(allSamples_NoNAs$Whole), method = "pearson")
AW_NormCor_P_sq = (AW_NormCor_P)^2
AW_NormCor_sq = (AW_NormCor)^2

PW_NormCor = cor((as.numeric(allSamples_NoNAs$Post)/ as.numeric(merged_norm_1kb_PostvsWhole.lm$coefficients[2])),as.numeric(allSamples_NoNAs$Whole), method = "spearman")
PW_NormCor_P = cor((as.numeric(allSamples_NoNAs$Post)/ as.numeric(merged_norm_1kb_PostvsWhole.lm$coefficients[2])),as.numeric(allSamples_NoNAs$Whole), method = "pearson")
PW_NormCor_P_sq = (PW_NormCor_P)^2
PW_NormCor_sq = (PW_NormCor)^2

############################ S1 Figs. ############################

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
  coord_equal() +
  ylim(0,3000) + 
  xlim(0,3000) +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  #annotate("text", x=100, y=500, label= "slope = 0.996 ", colour = 'darkorange1') +
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

############ S1Fig-A Whole Vs Ant ############ 

png('011918_AllRegions1KB_Scatter_Merged_PE_10M_NoNormalization_AntVsWhole.png', width = 2000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Whole),y= as.numeric(Ant))) +
  geom_point(colour="grey35", alpha = 0.2) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('011918_AllRegions1KB_Scatter_Merged_PE_10M_NoNormalization_WholevsAnt') +
  xlab('Whole Halves merged') +
  ylab('Ant Halves merged') +
  geom_smooth(method='lm', formula = y~x) +
  coord_equal() +
  ylim(0,3000) + 
  xlim(0,3000) +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
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

############ S1Fig-C Whole Vs Post ############ 

png('011918_AllRegions1KB_Scatter_Merged_PE_10M_NoNormalization_WholevsPost.png', width = 2000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Whole),y= as.numeric(Post))) +
  geom_point(colour="grey35", alpha = 0.2) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('011918_AllRegions1KB_Scatter_Merged_PE_10M_NoNormalization_WholevsPost') +
  xlab('Whole Halves merged') +
  ylab('Post Halves merged') +
  geom_smooth(method='lm', formula = y~x) +
  coord_equal() +
  ylim(0,3000) + 
  xlim(0,3000) +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
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

############## Lin Reg Normalization ##############
############## S1Fig-B
png('011918_AllRegions1KB_Scatter_Merged_PE_10M_Linregwhole_WholevsAnt.png', width = 2000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Whole),y= (as.numeric(Ant)/ as.numeric(merged_norm_1kb_AntvsWhole.lm$coefficients[2])))) +
  geom_point(colour="grey35", alpha = 0.2) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('011918_AllRegions1KB_Scatter_Merged_PE_10M_Linregwhole_WholevsAnt') +
  # scale_x_log10() +
  # scale_y_log10() +
  coord_equal() +
  ylim(0,3000) + 
  xlim(0,3000) +
  xlab('Whole Halves merged') +
  ylab('Ant LinRegWhole merged') +
  geom_smooth(method='lm', formula = y~x) +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
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

png('011918_AllRegions1KB_Scatter_Merged_PE_10M_Linregwhole_WholevsPost.png', width = 2000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Whole),y= (as.numeric(Post)/ as.numeric(merged_norm_1kb_PostvsWhole.lm$coefficients[2])))) +
  geom_point(colour="grey35", alpha = 0.2) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('011918_AllRegions1KB_Scatter_Merged_PE_10M_Linregwhole_WholevsPost') +
  xlab('Whole Halves merged') +
  ylab('Post linregWhole') +
  geom_smooth(method='lm', formula = y~x) +
  coord_equal() +
  ylim(0,3000) + 
  xlim(0,3000) +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
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
dev.off()`