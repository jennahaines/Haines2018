############################################################################################################
##########   Figures_1_3_4_S1_S2_S3_S6_S7.r                                                                  ###########
##########   011818_RefinedGeneList_Dataanalysis_linregv2_PEReviews.r                            ###########
##########   Inputs: summarized text file, randomregion output, spatially patterned gene list    ###########
##########   Outputs: graphs and analysis                                                        ###########
############################################################################################################


library("plyr")
library("dplyr")
library("ggplot2")

#######################################################################
########################### IMPORT  FILES #############################
#######################################################################

allSamples <-read.delim2("021218_012318GENELIST_Randomregv2_output.txt", sep = "\t", stringsAsFactors = FALSE)
allSamples_1kb_bins <- read.delim2("021218_wig_sig_around_bedfile_1kbwindows_all.txt", sep = "\t", stringsAsFactors = FALSE)
allPeaksFile <- read.delim2("021518_2reps_REppeaks_overlap_REVISEDGENELIST.txt", sep ="\t", stringsAsFactors = FALSE, header = FALSE)

#######################################################################
########################### DATA FILTERING  ###########################
#######################################################################

####### Peak matching --> Match the filtered peak file back to the sample list #######
# Parse all the peaks out by sample, then only keep the highest fold change per sample, then match it back to the original dataframe
allPeaks = allPeaksFile %>%
  select(V1, V2, V3, V6, V11, V12, V13, V17, V21)
colnames(allPeaks) = c("chr", "start", "end", "name", "chrPeak",	"startPeak",	"endPeak", "fold_enrichment", "Peakname")

allPeaksWhole <- allPeaks %>%
  dplyr::filter(grepl(pattern = "Whole", Peakname)) %>%
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

################# Filter out regions that were confirmed spatially resolved in situ expression and overlap peaks
allSamples_yes_peaks = allSamples %>%
  filter(X011617.use == "yes") %>%
  filter(is.na(NoPeak) == TRUE) %>%
  select(chr, start, end, X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_10whole_shifted_lessthan130_10MNorm, X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole, X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole, X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.03.1plus2_shifted_lessthan130_1MNorm, X101416_10MNorm_DNase.I_dmels5r1, Name,New.Location.Assignment, Type) %>%
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

#######################################################################
########################### FIGURE 1  #################################
#######################################################################

######### ######### ######### Figure 1 Scatterplots of 1kb regions

################# filter our NANs in the 1kb bin file and make standard column names
allSamples_NoNAs = allSamples_1kb_bins %>%
  select(X.chr., X.start., X.end., X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_10whole_shifted_lessthan130_10MNorm, X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole, X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole, X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.03.1plus2_shifted_lessthan130_1MNorm, X101416_10MNorm_DNase.I_dmels5r1) %>%
  filter(X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_10whole_shifted_lessthan130_10MNorm != 'nan') %>%
  filter(X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole!= 'nan') %>%
  filter(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.03.1plus2_shifted_lessthan130_1MNorm!= 'nan') %>%
  filter(X101416_10MNorm_DNase.I_dmels5r1!= 'nan') %>%
  filter(X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole!= 'nan') %>%
  mutate(Name = "NA") %>%
  mutate(Location = "NA") %>%
  mutate(Type = "NA") %>%
  mutate(Dotsize = 1)
colnames(allSamples_NoNAs) = c("Chr", "Start", "End", "Whole", "Ant", "LinregPost", "AplusP", "Dnase1", "Name", "Location", "Type", "Dotsize")



######### ######### ######### Ant

### Correction new graph for figure 1 that is just AvsP
AvsP_cor = cor(as.numeric(allSamples_NoNAs$Ant),y= as.numeric(allSamples_NoNAs$LinregPost), method = "spearman")
AvsP_cor_Pearson = cor(as.numeric(allSamples_NoNAs$Ant),y= as.numeric(allSamples_NoNAs$LinregPost), method = "pearson")
AvsP_cor_rsq = (AvsP_cor)^2
AvsP_cor_Pearson_rsq = (AvsP_cor_Pearson)^2


png('020118_PE2reps_nodups_LinRegWhole_AllRegions1KB_AvsP_log.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(LinregPost),y= as.numeric(Ant), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.05, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightskyblue2", alpha = 0.5, size = 0.4) + 
  ggtitle('020118_PE_2reps_nodups_LinRegWhole_AllRegions1KB_AvsP') +
  scale_x_log10(limits = c(10,1500)) +
  scale_y_log10(limits = c(10,1500)) +
  xlab('Reg Norm Whole Posterior') +
  ylab('Reg Norm Whole Anterior') +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=1000, y=20, label= AvsP_cor) +
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


png('020118_PE2reps_nodups_LinRegWhole_AllRegions1KB_AvsP_nolog.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(LinregPost),y= as.numeric(Ant), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.05, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightskyblue2", alpha = 0.5, size = 0.4) + 
  ggtitle('020118_PE_2reps_nodups_LinRegWhole_AllRegions1KB_AvsP') +
  xlab('Reg Norm Whole Posterior') +
  ylab('Reg Norm Whole Anterior') +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=1000, y=20, label= AvsP_cor) +
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

####################################################################################
########################### FIGURE 3A, 3B, 4A, 4B  #################################
####################################################################################

################# merge the regions with the 1kb windows
MergedALL = rbind(allSamples_NoNAs, allSamples_yes_peaks)


######### ######### ######### AP ENHANCER SCATTERPLOT
AP_Enhancers_windows = MergedALL %>%
  filter(Type %in% c("Enhancer", "NA")) %>%
  filter(Location %in% c("NA", 'Anterior', "Posterior", "Mostly Post"))
AP_Enhancers_windows$Location <- factor(AP_Enhancers_windows$Location,
                                        levels = c("NA", 'Anterior', "Mostly Post", 'Posterior') ,ordered = TRUE)

png('020118_PE_2reps_nodups_AllRegions1KB_RegNormWhole_APEnhancers_Scatter_0118Newgenelist.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(AP_Enhancers_windows,
       aes(x = as.numeric(Ant),y= as.numeric(LinregPost), colour = as.factor(Location))) +
  scale_color_manual(name = 'Legand', values=c("grey70","darkorange1",  "dodgerblue3", "dodgerblue3")) +
  geom_point(alpha = 0.8, aes(size = Dotsize)) +
  scale_size_continuous(range = c(1, 2)) +
  ggtitle('020118_PE_2reps_nodups_AllRegions1KB_RegNormWhole_APEnhancers_Scatter_0118Newgenelist') +
  xlab('Reg Norm Whole Anterior') +
  ylab('Reg Norm Whole Posterior') +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
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

######### ######### ######### DV ENHANCER SCATTERPLOT
DV_Enhancers_windows = MergedALL %>%
  filter(Type %in% c("Enhancer", "NA")) %>%
  filter(Location %in% c("NA", 'Dorsal', "Ventral"))
DV_Enhancers_windows$Location <- factor(DV_Enhancers_windows$Location,
                                        levels = c("NA", 'Dorsal', 'Ventral') ,ordered = TRUE)

png('020118_PE_2reps_nodups_AllRegions1KB_DVEnhancers_LinRegWhole_Scatter_0118newgeneList.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(DV_Enhancers_windows,
       aes(x = as.numeric(Ant),y= as.numeric(LinregPost), colour = as.factor(Location))) +
  scale_color_manual(name = 'Legend', values=c( "grey70", "orchid4", "seagreen4")) +
  geom_point(alpha = 0.9, aes(size = Dotsize)) +
  scale_size_continuous(range = c(1, 2)) +
  ggtitle('020118_PE_2reps_nodups_AllRegions1KB_DVEnhancers_LinRegWhole_Scatter_0118newgeneList') +
  xlab('LinregWhole Anterior') +
  ylab('LinregWhole  Posterior') +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
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


png('020118_PE_2reps_nodups_AllRegions1KB_APPromoters_LinRegWhole_Scatter_0118newgeneList.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(AP_Promoters_windows,
       aes(x = as.numeric(Ant),y= as.numeric(LinregPost), colour = as.factor(Location))) +
  scale_color_manual(name = 'Legand', values=c("grey70","darkorange1" , "darkorange1",  "dodgerblue3", "dodgerblue3")) +
  geom_point(alpha = 1, aes(size = Dotsize)) +
  scale_size_continuous(range = c(1, 2)) +
  ggtitle('020118_PE_2reps_nodups_AllRegions1KB_APPromoters_LinRegWhole_Scatter_0118newgeneList') +
  xlab('LinRegWhole Anterior') +
  ylab('LinRegWhole Posterior') +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
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

######### ######### ######### DV Promoters SCATTERPLOT
DV_Promoter_windows = MergedALL %>%
  filter(Type %in% c("Promoter", "NA")) %>%
  filter(Location %in% c("NA", 'Dorsal', "Ventral"))
DV_Promoter_windows$Location <- factor(DV_Promoter_windows$Location,
                                       levels = c("NA", 'Dorsal', 'Ventral') ,ordered = TRUE)

png('020118_PE_2reps_nodups_AllRegions1KB_DVPromoters_linRegWhole_Scatter_0118newGeneList.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(DV_Promoter_windows,
       aes(x = as.numeric(Ant),y= as.numeric(LinregPost), colour = as.factor(Location))) +
  scale_color_manual(name = 'Legand', values=c( "grey70", "orchid4", "seagreen4")) +
  geom_point(alpha = 0.9, aes(size = Dotsize)) +
  scale_size_continuous(range = c(1, 2)) +
  ggtitle('020118_PE_2reps_nodups_AllRegions1KB_DVPromoters_linRegWhole_Scatter_0118newGeneList') +
  xlab('LinregWhole Anterior') +
  ylab('LinRegWhole Posterior') +
  ylim(0,1500) + 
  xlim(0,1500) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
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

####################################################################################
########################### FIGURE 3D, 4D, S4  #####################################
####################################################################################

######### ######### ######### AP Enhancer Bargraphs 
#########  Filter out A vs P enhancers to calculate the positional score

allSamples_yes_peaks_ATACSkew = allSamples %>%
  filter(X011617.use == "yes") %>%
  filter(is.na(NoPeak) == TRUE) %>%
  select(chr, start, end, X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_10whole_shifted_lessthan130_10MNorm, X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole, X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole, X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.03.1plus2_shifted_lessthan130_1MNorm, X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm, Name,New.Location.Assignment, Type, ATACSkewScore) %>%
  mutate(Dotsize = 1.5)
colnames(allSamples_yes_peaks_ATACSkew) = c("Chr", "Start", "End", "Whole", "Ant", "LinregPost", "AplusP", "Dnase1", "Name", "Location", "Type", "ATACSkewScore", "Dotsize")



JustAEnhancers = allSamples_yes_peaks_ATACSkew %>%
  select(Chr, Start, End, Location, Type, Name, Ant, LinregPost, ATACSkewScore) %>%
  filter(Location %in%  c('Anterior')) %>%
  filter(Type == 'Enhancer') %>%
  mutate(PositionalScore = (as.numeric(Ant) - as.numeric(LinregPost)) / (as.numeric(Ant) + as.numeric(LinregPost))) %>%
  distinct()

JustPEnhancers = allSamples_yes_peaks_ATACSkew %>%
  select(Chr, Start, End, Location, Type, Name, Ant, LinregPost, ATACSkewScore) %>%
  filter(Location %in%  c('Posterior', 'Mostly Post')) %>%
  filter(Type == 'Enhancer') %>%
  mutate(PositionalScore = (as.numeric(Ant) - as.numeric(LinregPost)) / (as.numeric(Ant) + as.numeric(LinregPost))) %>%
  distinct()

JustAPEnhancers = bind_rows(JustAEnhancers, JustPEnhancers)

write.csv(JustAPEnhancers, "020318_2Reps_NewGenelist_PE_Analysis_APEnhancers_PositionalScore.csv")

JustAPEnhancers$Location <- factor(JustAPEnhancers$Location,
                                 levels = c('Anterior', 'Posterior', 'Mostly Post'),ordered = TRUE)
ggplot(JustAPEnhancers,
       aes(x = reorder(as.factor(Name),- as.numeric(PositionalScore)),
           y = as.numeric(PositionalScore),
           fill = Location)) +
  geom_col(position = "stack", colour = "black") +
  ylim(-0.6, 0.6) +
  xlab("A-P Enhancers") +
  ylab("PositionalScore") +
  ggtitle("020118_PE_2reps_nodups_LinRegWhole_AP_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList") +
  scale_fill_manual(name = "Legand", values=c("darkorange1", "dodgerblue3", "dodgerblue3")) +
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
  png('020118_PE_2reps_nodups_LinRegWhole_AP_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList.png', width = 5000, height = 2000, units = "px",  res=300)
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

png('020118_2reps_nodups_AP_Promoters_PositionalScore_Bargraph_Position_0118newgeneList.png', width = 5000, height = 2000, units = "px",  res=300)
ggplot(JustAPPromoters,
       aes(x = reorder(as.factor(Name),- as.numeric(PositionalScore)),
           y = as.numeric(PositionalScore),
           fill = Location)) +
  geom_col(position = "stack", colour = "black") +
  ylim(-0.6, 0.6) +
  xlab("A-P Promoters") +
  ylab("PositionalScore") +
  ggtitle("020118_2reps_nodups_AP_Promoters_PositionalScore_Bargraph_Position_0118newgeneList") +
  scale_fill_manual(name = "Legand", values=c("darkorange1", "darkorange1", "dodgerblue3", "dodgerblue3")) +
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
  geom_col(position = "stack", colour = "black") +
  ylim(-0.6, 0.6) +
  xlab("D-V Enhancers") +
  ylab("PositionalScore") +
  ggtitle("020118_PE_2reps_nodups_LinRegWhole_DV_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList") +
  scale_fill_manual(name = "Legand", values=c("orchid4", "seagreen4")) +
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
  png('020118_PE_2reps_nodups_LinRegWhole_DV_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList.png', width = 5000, height = 2000, units = "px",  res=300)
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
  ggtitle("020118_2reps_nodups_LinRegWhole_DV_Promoters_PositionalScore_Bargraph_Position_0118newgeneList") +
  scale_fill_manual(name = "Legand", values=c("orchid4", "seagreen4")) +
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
  png('020118_2reps_nodups_LinRegWhole_DV_Promoters_PositionalScore_Bargraph_Position_0118newgeneList.png', width = 5000, height = 2000, units = "px",  res=300)
dev.off()

####################################################################################
########################### FIGURE 3C  #############################################
####################################################################################

allSamples$Sample <- seq.int(nrow(allSamples)) # Add an sampleID row at the end

RandomRegions = allSamples %>%
  select(Sample, X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole.1, X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole.1, New.Location.Assignment, Type, RandSkewScore) %>%
  filter(Type =="Enhancer") %>%
  mutate(Location.real = "Random") %>%
  mutate(Type.real = "Random") %>%
  mutate(TotalScore = as.numeric(X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole.1) + as.numeric(X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole.1)) %>%
  select(Sample, Location.real, Type.real, RandSkewScore, TotalScore)
colnames(RandomRegions) = c("Name", "Location", "Type", "ATACSkewScore", "TotalScore")
RandomRegions$Name = as.character(RandomRegions$Name)

RandomRegions$ATACSkewScore = as.numeric(RandomRegions$ATACSkewScore)

JustAEnhancers_merged = JustAEnhancers %>%
  filter(Location %in%  c('Anterior')) %>%
  filter(Type == 'Enhancer') %>%
  mutate(TotalScore = as.numeric(Ant) + as.numeric(LinregPost)) %>%
  mutate(merged_location = "Anterior") %>%
  select(Name, merged_location, Type, PositionalScore, TotalScore) %>%
  distinct()


JustPEnhancers_merged = JustPEnhancers %>%
  filter(Location %in%  c('Posterior', 'Mostly Post')) %>%
  filter(Type == 'Enhancer') %>%
  mutate(TotalScore = as.numeric(Ant) + as.numeric(LinregPost)) %>%
  mutate(merged_location = "Posterior") %>%
  select(Name, merged_location, Type, PositionalScore, TotalScore) %>%
  distinct()

JustDorsalEnhancers_merged = allSamples_yes_peaks_ATACSkew %>%
  filter(Location == "Dorsal") %>%
  filter(Type == 'Enhancer') %>%
  mutate(TotalScore = as.numeric(Ant) + as.numeric(LinregPost)) %>%
  mutate(merged_location = "Dorsal") %>%
  select(Name, merged_location, Type, ATACSkewScore, TotalScore) %>%
  distinct()
JustDorsalEnhancers_merged$ATACSkewScore = as.numeric(JustDorsalEnhancers_merged$ATACSkewScore)

JustVentralEnhancers_merged = allSamples_yes_peaks_ATACSkew %>%
  filter(Location == "Ventral") %>%
  filter(Type == 'Enhancer') %>%
  mutate(TotalScore = as.numeric(Ant) + as.numeric(LinregPost)) %>%
  mutate(merged_location = "Ventral") %>%
  select(Name, merged_location, Type, ATACSkewScore, TotalScore) %>%
  distinct()
JustVentralEnhancers_merged$ATACSkewScore = as.numeric(JustVentralEnhancers_merged$ATACSkewScore)

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
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.2 , alpha = 0.3, binpositions="all") +
  xlab("Location") +
  ylab("Positional Skew Score") +
  coord_cartesian(ylim = c(-0.5, 0.7)) +
  geom_abline(slope = 0, intercept=0, linetype = 'dotted') +
  ggtitle("020118_PE_2reps_nodups_LinRegWhole_AllLocations_Enhancers_ATACSkewScore_WithRandom_0118newgeneList") +
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
  png('020118_PE_2reps_nodups_LinRegWhole_AllLocations_Enhancers_ATACSkewScore_WithRandom_0118newgeneList.png', width = 2000, height = 2000, units = "px",  res=300)
dev.off()

################# ANOVAS for Enhancers  ################# 
Anova_Enhancers_dataframe = rand_and_Enhancers %>%
  select(Location, ATACSkewScore)

Anova_Enhancers_lm = lm(as.numeric(Anova_Enhancers_dataframe$ATACSkewScore) ~ Anova_Enhancers_dataframe$Location)
Anova_Enhancer = aov(formula = as.numeric(Anova_Enhancers_dataframe$ATACSkewScore) ~ Anova_Enhancers_dataframe$Location)
Tukey_Enhancers = TukeyHSD(Anova_Enhancer)

write.table(Tukey_Enhancers$`Anova_Enhancers_dataframe$Location`, "020118_Anova_ALLLocations_Enhancers.txt")

####################################################################################
########################### FIGURE 4C  #############################################
####################################################################################

Susan_RNAseq_data <- read.delim2("~/Box Sync/Eisen_Lab/Experiments/ATAC-seq/Halves/ATAC-seq_Pools/040617_Analysis/061417_BrowserTraces/062317_RegressionNormalization/081617_2reps/Lott2011_RNA-seq_GeneExpression.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(Susan_RNAseq_data)[1] = "Name"
Joined_ATACPeaks_plus_RNAseq = full_join(allSamples, Susan_RNAseq_data, by = "Name")

JustPromoters = Joined_ATACPeaks_plus_RNAseq %>%
  filter(Type == "Promoter") %>%
  select(Name, Type, New.Location.Assignment,X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole ,X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole , ATACSkewScore, PValue, WholePeaks, AntPeaks, PostPeaks, NoPeak, CLASS, F10, F11, F12, U13, F14A, U14A, F14B, F14B_r2, F14C, F14C_r2, F14D, U14D)

write.csv(JustPromoters, file = "020118_Promoters_peaks_RNAseq.csv")

JustAProm_merged_ZYG = JustAProm_merged %>%
  filter(CLASS == "zyg")

JustP_prom_merged_ZYG = JustPProm_merged %>%
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
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.2 , alpha = 0.3, binpositions="all") +
  xlab("Location") +
  ylab("Positional Skew Score") +
  coord_cartesian(ylim = c(-0.5, 0.7)) +
  geom_abline(slope = 0, intercept=0, linetype = 'dotted') +
  ggtitle("020118_PE_2reps_nodups_LinRegWhole_ALLLocations_Promoters_ATACSkewScore_WithRandom_0118newgeneList_Zyggenesonly") +
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
  png('020118_PE_2reps_nodups_LinRegWhole_ALLLocations_Promoters_ATACSkewScore_WithRandom_0118newgeneList_Zyggenesonly.png', width = 2000, height = 2000, units = "px",  res=300)
dev.off()

################# ANOVAS for BoxPlot Zyg only promoters  ################# 
Anova_Promoters_dataframe = rand_and_PromotersZYG %>%
  select(Location, ATACSkewScore)

Anova_Promoters_lm = lm(as.numeric(Anova_Promoters_dataframe$ATACSkewScore) ~ Anova_Promoters_dataframe$Location)
Anova_Promoter = aov(formula = as.numeric(Anova_Promoters_dataframe$ATACSkewScore) ~ Anova_Promoters_dataframe$Location)
Tukey_Promoters = TukeyHSD(Anova_Promoter)
write.table(Tukey_Promoters$`Anova_Promoters_dataframe$Location`, "020118_Anova_ALLLocations_Promoters_Zygonly.txt")

####################################################################################
########################### FIGURE S6  #############################################
####################################################################################

###### Positional score replicates                                                


#first I need to calculate Positional score for all the replicates

RepSamples_NoNAs = allSamples %>%
  select(chr, start, end, Type, New.Location.Assignment, Name, X011617.use, ATACSkewScore, NoPeak, X011617_101117_RemDUP_PE_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.01.20a_shifted_lessthan130_1MNorm, X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.01.20Ant_shifted_lessthan130_1MNorm, X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.02.20p_shifted_lessthan130_1MNorm,  X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.02.20Post_shifted_lessthan130_1MNorm, X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.7.10whole_shifted_lessthan130_1MNorm, X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.04.10whole_shifted_lessthan130_1MNorm) %>%
  filter(X011617.use == "yes") %>%
  filter(is.na(NoPeak) == TRUE) %>%
  mutate(PositionalScore_1105 = (as.numeric(X011617_101117_RemDUP_PE_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.01.20a_shifted_lessthan130_1MNorm) - as.numeric(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.02.20p_shifted_lessthan130_1MNorm)) / (as.numeric(X011617_101117_RemDUP_PE_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.01.20a_shifted_lessthan130_1MNorm) + as.numeric(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.02.20p_shifted_lessthan130_1MNorm))) %>%
  mutate(PositionalScore_1109 = (as.numeric(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.01.20Ant_shifted_lessthan130_1MNorm) - as.numeric(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.02.20Post_shifted_lessthan130_1MNorm)) / (as.numeric(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.01.20Ant_shifted_lessthan130_1MNorm) + as.numeric(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.02.20Post_shifted_lessthan130_1MNorm))) %>%
  select(Name, PositionalScore_1105, PositionalScore_1109)

# Then add it back to the positional score data frame used in the bargraph
JustAEnhancers_Reps = inner_join(JustAEnhancers, RepSamples_NoNAs, by = "Name")
JustPEnhancers_Reps = inner_join(JustPEnhancers, RepSamples_NoNAs, by = "Name")

JustAPEnhancers_Reps = bind_rows(JustAEnhancers_Reps, JustPEnhancers_Reps)

JustAPEnhancers_Reps$Location <- factor(JustAPEnhancers_Reps$Location,
                                   levels = c('Anterior', 'Posterior', 'Mostly Post'),ordered = TRUE)

ggplot(JustAPEnhancers_Reps,
       aes(x = reorder(as.factor(Name),- as.numeric(PositionalScore)),
           y = as.numeric(PositionalScore),
           fill = Location)) +
  geom_col(position = "stack", alpha = 0.2) +
  geom_boxplot(outlier.size = 0) +
  ylim(-0.8, 0.8) +
  xlab("A-P Enhancers") +
  ylab("PositionalScore") +
  geom_dotplot(aes(x = reorder(as.factor(Name),- as.numeric(PositionalScore)), y = as.numeric(PositionalScore_1105)) ,binaxis='y', stackdir='center', dotsize = 0.7 , alpha = 1, binpositions="all") +
  geom_dotplot(aes(x = reorder(as.factor(Name),- as.numeric(PositionalScore)), y = as.numeric(PositionalScore_1109)) ,binaxis='y', stackdir='center', dotsize = 0.7 , alpha = 1, binpositions="all") +
  geom_hline(yintercept = 0) +
  ggtitle("020118_PE_2reps_nodups_LinRegWhole_AP_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList_REPLICATES") +
  scale_fill_manual(name = "Legand", values=c("darkorange1", "dodgerblue3", "dodgerblue3")) +
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
  png('020118_PE_2reps_nodups_LinRegWhole_AP_Enhancers_PositionalScore_Bargraph_Position_0118newgeneListREPLICATES.png', width = 5000, height = 2000, units = "px",  res=300)
dev.off()

####################################################################################
########################### FIGURE S1  #############################################
####################################################################################

RepSamples_NoNAs = allSamples_1kb_bins %>%
  select(X.chr., X.start., X.end., X011617_101117_RemDUP_PE_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.01.20a_shifted_lessthan130_1MNorm, X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.01.20Ant_shifted_lessthan130_1MNorm, X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.02.20p_shifted_lessthan130_1MNorm,  X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.02.20Post_shifted_lessthan130_1MNorm, X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.7.10whole_shifted_lessthan130_1MNorm, X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.04.10whole_shifted_lessthan130_1MNorm, X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm) %>%
  filter(X011617_101117_RemDUP_PE_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.01.20a_shifted_lessthan130_1MNorm != 'nan') %>%
  filter(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.01.20Ant_shifted_lessthan130_1MNorm != 'nan') %>%
  filter(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.02.20p_shifted_lessthan130_1MNorm != 'nan') %>%
  filter(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.02.20Post_shifted_lessthan130_1MNorm != 'nan') %>%
  filter(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.7.10whole_shifted_lessthan130_1MNorm != 'nan') %>%
  filter(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1109.ATACSlice03.04.10whole_shifted_lessthan130_1MNorm != 'nan') %>%
  filter(X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm != 'nan')
colnames(RepSamples_NoNAs) = c("Chr", "Start", "End", "1105-20a", "1109-20a", "1105-20p", "1109-20p", "1105-10w", "1109-10w", "DnaseI")


#### Anterior replicate graph
AntReps_cor = cor(as.numeric(RepSamples_NoNAs$`1105-20a`),y= as.numeric(RepSamples_NoNAs$`1109-20a`), method = "spearman")
AntReps_cor_pearson = cor(as.numeric(RepSamples_NoNAs$`1105-20a`),y= as.numeric(RepSamples_NoNAs$`1109-20a`), method = "pearson")
AntReps_cor_r_squared = (AntReps_cor)^2
AntReps_cor_pearson_Squared = (AntReps_cor_pearson)^2
png('020118_2reps_nodups_Antreps_nolog.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(RepSamples_NoNAs,
       aes(x = as.numeric(`1105-20a`),y= as.numeric(`1109-20a`))) +
  geom_point(alpha = 0.05, size = 0.3, show.legend = FALSE, colour = "lightsteelblue4") +
  #scale_size_continuous(range = c(1, 2)) +
  geom_density2d(colour="lightblue1", alpha = 0.5, size = 0.4) + 
  ggtitle('020118_2reps_nodups_Antreps_nolog') +
  #scale_x_log10(limit = c(10,1500)) +
  #scale_y_log10(limit = c(10,1500)) +
  scale_x_log10() +
  scale_y_log10() +
  xlab('nodups 1105-20ant') +
  ylab('nodups 1109-20ant') +
  #geom_smooth(method='lm', formula = y~x) +
  #ylim(0, 300) + 
  #xlim(0, 300) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=500, y=500, label= AntReps_cor) +
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

#### Posterior replicate graph
PostReps_cor = cor(as.numeric(RepSamples_NoNAs$`1105-20p`),y= as.numeric(RepSamples_NoNAs$`1109-20p`), method = "spearman")
PostReps_cor_pearson = cor(as.numeric(RepSamples_NoNAs$`1105-20p`),y= as.numeric(RepSamples_NoNAs$`1109-20p`), method = "pearson")
PostReps_cor_rsq = (PostReps_cor)^2
PostReps_cor_pearson__rsq = (PostReps_cor_pearson)^2
png('020118_2reps_nodups_Postreps_nolog.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(RepSamples_NoNAs,
       aes(x = as.numeric(`1105-20p`),y= as.numeric(`1109-20p`))) +
  geom_point(alpha = 0.05, size = 0.3, show.legend = FALSE, colour = "lightsteelblue4") +
  #scale_size_continuous(range = c(1, 2)) +
  geom_density2d(colour="lightblue1", alpha = 0.5, size = 0.4) + 
  ggtitle('020118_2reps_nodups_Postreps_nolog') +
  #scale_x_log10(limit = c(10,1500)) +
  #scale_y_log10(limit = c(10,1500)) +
  xlab('nodups 1105-20post') +
  ylab('nodups 1109-20post') +
  #geom_smooth(method='lm', formula = y~x) +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=500, y=500, label= PostReps_cor) +
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

#### Whole replicate graph
WholeReps_cor = cor(as.numeric(RepSamples_NoNAs$`1105-10w`),y= as.numeric(RepSamples_NoNAs$`1109-10w`), method = "spearman")
WholeReps_cor_pearson = cor(as.numeric(RepSamples_NoNAs$`1105-10w`),y= as.numeric(RepSamples_NoNAs$`1109-10w`), method = "pearson")
WholeReps_cor_rsq = (WholeReps_cor)^2
WholeReps_cor_pearson_rsq = (WholeReps_cor_pearson)^2
png('020118_2reps_nodups_Wholereps_nolog.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(RepSamples_NoNAs,
       aes(x = as.numeric(`1105-10w`),y= as.numeric(`1109-10w`))) +
  geom_point(alpha = 0.05, size = 0.3, show.legend = FALSE, colour = "lightsteelblue4") +
  #scale_size_continuous(range = c(1, 2)) +
  geom_density2d(colour="lightblue1", alpha = 0.5, size = 0.4) + 
  ggtitle('020118_2reps_nodups_wholereps_nolog') +
  #scale_x_log10(limit = c(10,1500)) +
  #scale_y_log10(limit = c(10,1500)) +
  xlab('nodups 1105-10whole') +
  ylab('nodups 1109-10whole') +
  #geom_smooth(method='lm', formula = y~x) +
  scale_x_log10() +
  scale_y_log10() +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=500, y=500, label= WholeReps_cor) +
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


####################################################################################
########################### FIGURE S2  #############################################
####################################################################################

######### ######### ######### Ant

### Correction new graph for figure 1 that is just AvsP
AvsP_cor = cor(as.numeric(allSamples_NoNAs$Ant),y= as.numeric(allSamples_NoNAs$LinregPost), method = "spearman")
AvsP_cor_Pearson = cor(as.numeric(allSamples_NoNAs$Ant),y= as.numeric(allSamples_NoNAs$LinregPost), method = "pearson")
AvsP_cor_rsq = (AvsP_cor)^2
AvsP_cor_Pearson_rsq = (AvsP_cor_Pearson)^2


png('020118_PE2reps_nodups_LinRegWhole_AllRegions1KB_AvsP_log.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(LinregPost),y= as.numeric(Ant), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  #scale_size_continuous(range = c(1, 2)) +
  geom_point(alpha = 0.05, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightskyblue2", alpha = 0.5, size = 0.4) + 
  ggtitle('020118_PE_2reps_nodups_LinRegWhole_AllRegions1KB_AvsP') +
  scale_x_log10(limits = c(10,1500)) +
  scale_y_log10(limits = c(10,1500)) +
  xlab('Reg Norm Whole Posterior') +
  ylab('Reg Norm Whole Anterior') +
  #geom_smooth(method='lm', formula = y~x) +
  #ylim(0,1500) + 
  #xlim(0,1500) +
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


png('020118_PE2reps_nodups_LinRegWhole_AllRegions1KB_AvsP_nolog.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(LinregPost),y= as.numeric(Ant), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  #scale_size_continuous(range = c(1, 2)) +
  geom_point(alpha = 0.05, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightskyblue2", alpha = 0.5, size = 0.4) + 
  ggtitle('020118_PE_2reps_nodups_LinRegWhole_AllRegions1KB_AvsP') +
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
AvsWhole_cor_Pearson = cor(as.numeric(allSamples_NoNAs$Ant),y= as.numeric(allSamples_NoNAs$Whole), method = "pearson")
AvsWhole_cor_RSQ = (AvsWhole_cor)^2
AvsWhole_cor_Pearson_RSQ = (AvsWhole_cor_Pearson)^2


png('020118_PE_2reps_nodups_LinRegWhole_AvsWhole_nolog.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Ant),y= as.numeric(Whole), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.05, size = 0.3, show.legend = FALSE) +
  #scale_size_continuous(range = c(1, 2)) +
  geom_density2d(colour="lightblue1", alpha = 0.5, size = 0.4) + 
  ggtitle('020118_PE_2reps_nodups_LinRegWhole_AvsWhole') +
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
AvsAP_cor_Pearson = cor(as.numeric(allSamples_NoNAs$Ant),y= as.numeric(allSamples_NoNAs$AplusP), method = "pearson")
AvsAP_cor_rq = (AvsAP_cor)^2
AvsAP_cor_Pearson_rq = (AvsAP_cor_Pearson)^2

png('020118_PE_2reps_nodups_LinRegWhole_AvsAP.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Ant),y= as.numeric(AplusP), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020118_PE_2reps_nodups_LinRegWhole_AvsAP') +
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
AvsDnaseI_cor_Pearson = cor(as.numeric(allSamples_NoNAs$Ant),y= as.numeric(allSamples_NoNAs$Dnase1), method = "pearson")
AvsDnaseI_cor_rsq = (AvsDnaseI_cor)^2
AvsDnaseI_cor_Pearson_rq = (AvsDnaseI_cor_Pearson)^2
png('021218_PE_2reps_nodups_LinRegWhole_AvsDnaseI.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Dnase1),y= as.numeric(Ant), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020118_PE_2reps_nodups_LinRegWhole_AvsDnaseI') +
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
PvsWhole_cor_Pearson = cor(as.numeric(allSamples_NoNAs$LinregPost),y= as.numeric(allSamples_NoNAs$Whole), method = "pearson")
PvsWhole_cor_rsq = (PvsWhole_cor)^2
PvsWhole_cor_Pearson_rsq = (PvsWhole_cor_Pearson)^2
png('020118_PE_2reps_nodups_LinRegWhole_PostvsWhole.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(LinregPost),y= as.numeric(Whole), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020118_PE_2reps_nodups_LinRegWhole_PostvsWhole') +
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
PvsAP_cor_Pearson = cor(as.numeric(allSamples_NoNAs$LinregPost),y= as.numeric(allSamples_NoNAs$AplusP), method = "pearson")
png('020118_PE_2reps_nodups_LinRegWhole_PostvsAplusP.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(LinregPost),y= as.numeric(AplusP), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020118_PE_2reps_nodups_LinRegWhole_PostvsAplusP') +
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
PvsDnaseI_cor_Pearson = cor(as.numeric(allSamples_NoNAs$LinregPost),y= as.numeric(allSamples_NoNAs$Dnase1), method = "pearson")
PvsDnaseI_cor_rsq = (PvsDnaseI_cor)^2
PvsDnaseI_cor_Pearson_rsq = (PvsDnaseI_cor_Pearson)^2

png('021218_PE_2reps_nodups_LinRegWhole_PostvsDnaseI.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Dnase1),y= as.numeric(LinregPost), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020118_PE_2reps_nodups_LinRegWhole_PostvsDnaseI') +
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
WholevsDnaseI_cor_Pearson = cor(as.numeric(allSamples_NoNAs$Whole),y= as.numeric(allSamples_NoNAs$Dnase1), method = "pearson")
WholevsDnaseI_cor_rsq = (WholevsDnaseI_cor)^2
WholevsDnaseI_cor_Pearson_rsq = (WholevsDnaseI_cor_Pearson)^2

png('021218_PE_2reps_nodups_LinRegWhole_WholevsDnaseI.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Dnase1),y= as.numeric(Whole), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020118_PE_2reps_nodups_LinRegWhole_WholevsDnaseI') +
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
WholevsAplusP_cor_Pearson = cor(as.numeric(allSamples_NoNAs$Whole),y= as.numeric(allSamples_NoNAs$AplusP), method = "pearson")
WholevsAplusP_cor_Pearson_rsq = (WholevsAplusP_cor_Pearson)^2
WholevsAplusP_cor_rsq = (WholevsAplusP_cor)^2
png('020118_PE_2reps_nodups_LinRegWhole_WholevsAplusP.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_NoNAs,
       aes(x = as.numeric(Whole),y= as.numeric(AplusP), colour = as.factor(Location))) +
  scale_color_manual(values=c("lightsteelblue4")) +
  geom_point(alpha = 0.4, size = 0.3, show.legend = FALSE) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020118_PE_2reps_nodups_LinRegWhole_WholevsAplusP') +
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

####################################################################################
########################### FIGURE S3  #############################################
####################################################################################
#### Import New DnaseI 1kb region graph 
allSamples_DNASEI_NoNAs = allSamples_1kb_bins %>%
  select(X.chr., X.start., X.end., X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm, X011818_100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED_shifted_10MNorm, X011818_100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED_shifted_10MNorm,  X011818_100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED_shifted_10MNorm) %>%
  filter(X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm != 'nan') %>%
  filter(X011818_100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED_shifted_10MNorm != 'nan') %>%
  filter(X011818_100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED_shifted_10MNorm != 'nan') %>%
  filter(X011818_100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED_shifted_10MNorm != 'nan') %>%
  mutate(Name = "NA") %>%
  mutate(Location = "NA") %>%
  mutate(Type = "NA") %>%
  mutate(Dotsize = 1) %>%
  select(Name, Location, Type, X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm, X011818_100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED_shifted_10MNorm, X011818_100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED_shifted_10MNorm, X011818_100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED_shifted_10MNorm, Dotsize)
#c


st5_st14_Cor_pearson = round(cor(as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm), as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED_shifted_10MNorm), method = "pearson"),2)
st5_st14_Cor = round(cor(as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm), as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED_shifted_10MNorm), method = "spearman"),2)
st5_st14_Cor_pearson_sq = (st5_st14_Cor_pearson)^2
st5_st14_Cor_sq = (st5_st14_Cor)^2

png('020318_DnaseI_Stg5-Stg14_AllRegions1KB_RegNormWhole_Scatter_Newgenelist_LOG.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_DNASEI_NoNAs,
       aes(x = as.numeric(BDTNP_DnaseAcce_S5r1),y= as.numeric(BDTNP_DnaseAcce_S14r1))) +
  scale_color_manual(name = 'Legand', values=c("grey70")) +
  geom_point(alpha = 0.3, colour = "grey44") +
  scale_size_continuous(range = c(1, 2)) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020318_DnaseI_Stg5-Stg14_AllRegions1KB_RegNormWhole_Scatter_Newgenelist_LOG') +
  scale_x_log10() +
  scale_y_log10() +
  xlab('BDTNP_DNaseI_stage5') +
  ylab('BDTNP_DNaseI_stage14') +
  #geom_smooth(method='lm', formula = y~x) +
  #ylim(0,300) + 
  #xlim(0,300) +
  #coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=100, y=300, label= st5_st14_Cor, colour = 'black') +
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


st5_st9_Cor = round(cor(as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm), as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED_shifted_10MNorm), method = "spearman"),2)
st5_st9_Cor_pearson = round(cor(as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm), as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED_shifted_10MNorm), method = "pearson"),2)

st5_st9_Cor_pearson_sq = (st5_st9_Cor_pearson)^2
st5_st9_Cor_sq = (st5_st9_Cor)^2


png('020318_DnaseI_Stg5-Stg9_AllRegions1KB_RegNormWhole_Scatter_Newgenelist_LOG.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_DNASEI_NoNAs,
       aes(x = as.numeric(X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm),y= as.numeric(X011818_100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED_shifted_10MNorm))) +
  scale_color_manual(name = 'Legand', values=c("grey70")) +
  geom_point(alpha = 0.3, colour = "grey44") +
  scale_size_continuous(range = c(1, 2)) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020318_DnaseI_Stg5-Stg9_AllRegions1KB_RegNormWhole_Scatter_Newgenelist_LOG') +
  scale_x_log10() +
  scale_y_log10() +
  xlab('BDTNP_DNaseI_stage5') +
  ylab('BDTNP_DNaseI_stage9') +
  #geom_smooth(method='lm', formula = y~x) +
  #ylim(0,300) + 
  #xlim(0,300) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=100, y=300, label= st5_st9_Cor, colour = 'black') +
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


st5_st11_Cor = round(cor(as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm), as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED_shifted_10MNorm), method = "spearman"),2)
st5_st11_Cor_pearson = round(cor(as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm), as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED_shifted_10MNorm), method = "pearson"),2)

st5_st11_Cor_pearson_sq = (st5_st11_Cor_pearson)^2
st5_st11_Cor_sq = (st5_st11_Cor)^2

png('020318_DnaseI_Stg5-Stg11_AllRegions1KB_RegNormWhole_Scatter_Newgenelist_LOG.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_DNASEI_NoNAs,
       aes(x = as.numeric(X011818_100417_BDTNP_DNaseI_stage5_ALLREPSCOMBINED_shifted_10MNorm),y= as.numeric(X011818_100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED_shifted_10MNorm))) +
  scale_color_manual(name = 'Legand', values=c("grey70")) +
  geom_point(alpha = 0.3, colour = "grey44") +
  scale_size_continuous(range = c(1, 2)) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020318_DnaseI_Stg5-Stg11_AllRegions1KB_RegNormWhole_Scatter_Newgenelist_LOG') +
  scale_x_log10() +
  scale_y_log10() +
  xlab('BDTNP_DNaseI_stage5') +
  ylab('BDTNP_DNaseI_stage11') +
  #geom_smooth(method='lm', formula = y~x) +
  #ylim(0,300) + 
  #xlim(0,300) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=100, y=300, label= st5_st11_Cor, colour = 'black') +
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


st9_st11_Cor = round(cor(as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED_shifted_10MNorm), as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED_shifted_10MNorm), method = "spearman"),2)
st9_st11_Cor_pearson = round(cor(as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED_shifted_10MNorm), as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED_shifted_10MNorm), method = "pearson"),2)

st9_st11_Cor_pearson_sq = (st9_st11_Cor_pearson)^2
st9_st11_Cor_sq = (st9_st11_Cor)^2

png('020318_DnaseI_Stg9-Stg11_AllRegions1KB_RegNormWhole_Scatter_Newgenelist_LOG.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_DNASEI_NoNAs,
       aes(x = as.numeric(X011818_100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED_shifted_10MNorm),y= as.numeric(X011818_100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED_shifted_10MNorm))) +
  scale_color_manual(name = 'Legand', values=c("grey70")) +
  geom_point(alpha = 0.3, colour = "grey44") +
  scale_size_continuous(range = c(1, 2)) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020318_DnaseI_Stg9-Stg11_AllRegions1KB_RegNormWhole_Scatter_Newgenelist_LOG') +
  scale_x_log10() +
  scale_y_log10() +
  xlab('BDTNP_DNaseI_stage9') +
  ylab('BDTNP_DNaseI_stage11') +
  #geom_smooth(method='lm', formula = y~x) +
  #ylim(0,300) + 
  #xlim(0,300) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=100, y=300, label= st9_st11_Cor, colour = 'black') +
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

st9_st14_Cor = round(cor(as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED_shifted_10MNorm), as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED_shifted_10MNorm), method = "spearman"),2)
st9_st14_Cor_pearson = round(cor(as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED_shifted_10MNorm), as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED_shifted_10MNorm), method = "pearson"),2)

st9_st14_Cor_pearson_sq = (st9_st14_Cor_pearson)^2
st9_st14_Cor_sq = (st9_st14_Cor)^2
png('020318_DnaseI_Stg9-Stg14_AllRegions1KB_RegNormWhole_Scatter_Newgenelist_LOG.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_DNASEI_NoNAs,
       aes(x = as.numeric(X011818_100417_BDTNP_DNaseI_stage9_ALLREPSCOMBINED_shifted_10MNorm),y= as.numeric(X011818_100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED_shifted_10MNorm))) +
  scale_color_manual(name = 'Legand', values=c("grey70")) +
  geom_point(alpha = 0.3, colour = "grey44") +
  scale_size_continuous(range = c(1, 2)) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020318_DnaseI_Stg9-Stg14_AllRegions1KB_RegNormWhole_Scatter_Newgenelist_LOG') +
  scale_x_log10() +
  scale_y_log10() +
  xlab('BDTNP_DNaseI_stage9') +
  ylab('BDTNP_DNaseI_stage14') +
  #geom_smooth(method='lm', formula = y~x) +
  #ylim(0,300) + 
  #xlim(0,300) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=100, y=300, label= st9_st14_Cor, colour = 'black') +
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

st11_st14_Cor = round(cor(as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED_shifted_10MNorm), as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED_shifted_10MNorm), method = "spearman"),2)
st11_st14_Cor_pearson = round(cor(as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED_shifted_10MNorm), as.numeric(allSamples_DNASEI_NoNAs$X011818_100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED_shifted_10MNorm), method = "pearson"),2)
st11_st14_Cor_pearson_sq = (st11_st14_Cor_pearson)^2
st11_st14_Cor_sq = (st11_st14_Cor)^2


png('020318_DnaseI_Stg11-Stg14_AllRegions1KB_RegNormWhole_Scatter_Newgenelist_LOG.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_DNASEI_NoNAs,
       aes(x = as.numeric(X011818_100417_BDTNP_DNaseI_stage11_ALLREPSCOMBINED_shifted_10MNorm),y= as.numeric(X011818_100417_BDTNP_DNaseI_stage14_ALLREPSCOMBINED_shifted_10MNorm))) +
  scale_color_manual(name = 'Legand', values=c("grey70")) +
  geom_point(alpha = 0.3, colour = "grey44") +
  scale_size_continuous(range = c(1, 2)) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020318_DnaseI_Stg11-Stg14_AllRegions1KB_RegNormWhole_Scatter_Newgenelist_LOG') +
  scale_x_log10() +
  scale_y_log10() +
  xlab('BDTNP_DNaseI_stage11') +
  ylab('BDTNP_DNaseI_stage14') +
  #geom_smooth(method='lm', formula = y~x) +
  #ylim(0,300) + 
  #xlim(0,300) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  annotate("text", x=100, y=300, label= st11_st14_Cor, colour = 'black') +
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
        strip.text


####################################################################################
########################### FIGURE S7  #############################################
####################################################################################

allSamples_yes_peaks_ATACSkew_1Avs1P = allSamples %>%
  filter(X011617.use == "yes") %>%
  filter(is.na(NoPeak) == TRUE) %>%
  select(chr, start, end, X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.5.1A_shifted_lessthan130_1MNorm, X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.6.1p_shifted_lessthan130_1MNorm, Name,New.Location.Assignment, Type, ATACSkewScore) %>%
  mutate(Dotsize = 1.5)
colnames(allSamples_yes_peaks_ATACSkew_1Avs1P) = c("Chr", "Start", "End", "OneA", "OneP", "Name", "Location", "Type", "ATACSkewScore", "Dotsize")



JustAEnhancers_one = allSamples_yes_peaks_ATACSkew_1Avs1P %>%
  select(Chr, Start, End, Location, Type, Name, OneA, OneP, ATACSkewScore) %>%
  filter(Location %in%  c('Anterior')) %>%
  filter(Type == 'Enhancer') %>%
  mutate(PositionalScore = (as.numeric(OneA) - as.numeric(OneP)) / (as.numeric(OneA) + as.numeric(OneP))) %>%
  distinct()

JustPEnhancers_one = allSamples_yes_peaks_ATACSkew_1Avs1P %>%
  select(Chr, Start, End, Location, Type, Name, OneA, OneP, ATACSkewScore) %>%
  filter(Location %in%  c('Posterior', 'Mostly Post')) %>%
  filter(Type == 'Enhancer') %>%
  mutate(PositionalScore = (as.numeric(OneA) - as.numeric(OneP)) / (as.numeric(OneA) + as.numeric(OneP))) %>%
  distinct()

JustAPEnhancers_one = bind_rows(JustAEnhancers_one, JustPEnhancers_one)


JustAPEnhancers_one$Location <- factor(JustAPEnhancers_one$Location,
                                       levels = c('Anterior', 'Posterior', 'Mostly Post'),ordered = TRUE)
ggplot(JustAPEnhancers_one,
       aes(x = reorder(as.factor(Name),- as.numeric(PositionalScore)),
           y = as.numeric(PositionalScore),
           fill = Location)) +
  geom_col(position = "stack", colour = "black") +
  xlab("A-P Enhancers") +
  ylab("PositionalScore") +
  ggtitle("020918_PE_2reps_nodups_SingleHalves_AP_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList") +
  scale_fill_manual(name = "Legand", values=c("darkorange1", "dodgerblue3", "dodgerblue3")) +
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
  png('020918_PE_2reps_nodups_SingleHalves_AP_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList.png', width = 5000, height = 2000, units = "px",  res=300)
dev.off()

oneAvsoneP_Cor_pearson = round(cor(as.numeric(allSamples_1kb_bins$X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.5.1A_shifted_lessthan130_1MNorm), as.numeric(allSamples_1kb_bins$X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.6.1p_shifted_lessthan130_1MNorm), method = "pearson"),2)
oneAvsoneP_Cor = round(cor(as.numeric(allSamples_1kb_bins$X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.5.1A_shifted_lessthan130_1MNorm), as.numeric(allSamples_1kb_bins$X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.6.1p_shifted_lessthan130_1MNorm), method = "spearman"),2)
oneAvsoneP_Cor_rsq = (oneAvsoneP_Cor)^2
oneAvsoneP_Cor_pearson_rsq = (oneAvsoneP_Cor_pearson)^2
oneAvsoneP_lm = lm(as.numeric(allSamples_1kb_bins$X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.5.1A_shifted_lessthan130_1MNorm) ~ as.numeric(allSamples_1kb_bins$X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.6.1p_shifted_lessthan130_1MNorm))


png('020918_1avs1p_Scatterplot_Norm.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_1kb_bins,
       aes(x = as.numeric(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.6.1p_shifted_lessthan130_1MNorm),y= (as.numeric(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.5.1A_shifted_lessthan130_1MNorm)/ as.numeric(oneAvsoneP_lm$coefficients[2])))) +
  scale_color_manual(name = 'Legand', values=c("grey70")) +
  geom_point(alpha = 0.3, colour = "grey44") +
  scale_size_continuous(range = c(1, 2)) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020918_1avs1p_Scatterplot_Norm') +
  #scale_x_log10() +
  #scale_y_log10() +
  xlab('1A') +
  ylab('1P') +
  #geom_smooth(method='lm', formula = y~x) +
  #ylim(0,300) + 
  #xlim(0,300) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  #annotate("text", x=100, y=300, label= st5_st14_Cor, colour = 'black') +
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

png('020918_1avs1p_Scatterplot.png', width = 4000, height = 2000, units = "px",  res=300) 
ggplot(allSamples_1kb_bins,
       aes(x = as.numeric(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.5.1A_shifted_lessthan130_1MNorm),y= as.numeric(X011617_081817_RemDUP_041217_Bowtie2_ME_JH_112315_1105.ATACSlice2.6.1p_shifted_lessthan130_1MNorm))) +
  scale_color_manual(name = 'Legand', values=c("grey70")) +
  geom_point(alpha = 0.3, colour = "grey44") +
  scale_size_continuous(range = c(1, 2)) +
  geom_density2d(colour="lightblue1") + 
  ggtitle('020918_1avs1p_Scatterplot') +
  #scale_x_log10() +
  #scale_y_log10() +
  xlab('1A') +
  ylab('1P') +
  #geom_smooth(method='lm', formula = y~x) +
  #ylim(0,300) + 
  #xlim(0,300) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  #annotate("text", x=100, y=300, label= st5_st14_Cor, colour = 'black') +
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


JustAPromoters_one = allSamples_yes_peaks_ATACSkew_1Avs1P %>%
  select(Chr, Start, End, Location, Type, Name, OneA, OneP, ATACSkewScore) %>%
  filter(Location %in%  c('Anterior', 'Mostly Ant')) %>%
  filter(Type == 'Promoter') %>%
  mutate(PositionalScore = (as.numeric(OneA) - as.numeric(OneP)) / (as.numeric(OneA) + as.numeric(OneP))) %>%
  distinct()

JustPPromoters_one = allSamples_yes_peaks_ATACSkew_1Avs1P %>%
  select(Chr, Start, End, Location, Type, Name, OneA, OneP, ATACSkewScore) %>%
  filter(Location %in%  c('Posterior', 'Mostly Post')) %>%
  filter(Type == 'Promoter') %>%
  mutate(PositionalScore = (as.numeric(OneA) - as.numeric(OneP)) / (as.numeric(OneA) + as.numeric(OneP))) %>%
  distinct()

JustAPPromoters_one = bind_rows(JustAPromoters_one, JustPPromoters_one)

#write.csv(JustAPEnhancers_one, "020318_2Reps_NewGenelist_PE_Analysis_APEnhancers_PositionalScore.csv")


JustAPPromoters_one$Location <- factor(JustAPPromoters_one$Location,
                                       levels = c('Anterior', 'Mostly Ant', 'Posterior', 'Mostly Post'),ordered = TRUE)

png('020918_PE_2reps_nodups_SingleHalves_AP_Promoter_PositionalScore_Bargraph_Position_0118newgeneList.png', width = 5000, height = 2000, units = "px",  res=300)
ggplot(JustAPPromoters_one,
       aes(x = reorder(as.factor(Name),- as.numeric(PositionalScore)),
           y = as.numeric(PositionalScore),
           fill = Location)) +
  #guides(fill= FALSE) +
  geom_col(position = "stack", colour = "black") +
  #ylim(-0.6, 0.6) +
  xlab("A-P Promoters") +
  ylab("PositionalScore") +
  ggtitle("020918_PE_2reps_nodups_SingleHalves_AP_Promoter_PositionalScore_Bargraph_Position_0118newgeneList") +
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
JustDorsalEnhancers_One = allSamples_yes_peaks_ATACSkew_1Avs1P %>%
  select(Location, Type, Name, OneA, OneP, ATACSkewScore) %>%
  filter(Location == "Dorsal") %>%
  filter(Type == 'Enhancer') %>%
  mutate(PositionalScore = (as.numeric(OneA) - as.numeric(OneP)) / (as.numeric(OneA) + as.numeric(OneP))) %>%
  distinct()

JustVentralEnhancers_One = allSamples_yes_peaks_ATACSkew_1Avs1P %>%
  select(Location, Type, Name, OneA, OneP, ATACSkewScore) %>%
  filter(Location == "Ventral") %>%
  filter(Type == 'Enhancer') %>%
  mutate(PositionalScore = (as.numeric(OneA) - as.numeric(OneP)) / (as.numeric(OneA) + as.numeric(OneP))) %>%
  distinct()


JustDVEnhancers_One = bind_rows(JustDorsalEnhancers_One, JustVentralEnhancers_One)

JustDVEnhancers_One$Location <- factor(JustDVEnhancers_One$Location,
                                   levels = c("Dorsal", "Ventral"),ordered = TRUE)

ggplot(JustDVEnhancers_One,
       aes(x = reorder(as.factor(Name),- as.numeric(PositionalScore)),
           y = as.numeric(PositionalScore),
           fill = Location)) +
  #guides(fill= FALSE) +
  geom_col(position = "stack", colour = "black") +
  #ylim(-0.6, 0.6) +
  xlab("D-V Enhancers") +
  ylab("PositionalScore") +
  ggtitle("020918_PE_2reps_nodups_SingleHalves_DV_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList") +
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
  png('020918_PE_2reps_nodups_SingleHalves_DV_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList.png', width = 5000, height = 2000, units = "px",  res=300)
dev.off()
######### ######### ######### DV Promoter Bargraphs 
JustDorsalPromoter_One = allSamples_yes_peaks_ATACSkew_1Avs1P %>%
  select(Location, Type, Name, OneA, OneP, ATACSkewScore) %>%
  filter(Location == "Dorsal") %>%
  filter(Type == 'Promoter') %>%
  mutate(PositionalScore = (as.numeric(OneA) - as.numeric(OneP)) / (as.numeric(OneA) + as.numeric(OneP))) %>%
  distinct()

JustVentralPromoter_One = allSamples_yes_peaks_ATACSkew_1Avs1P %>%
  select(Location, Type, Name, OneA, OneP, ATACSkewScore) %>%
  filter(Location == "Ventral") %>%
  filter(Type == 'Promoter') %>%
  mutate(PositionalScore = (as.numeric(OneA) - as.numeric(OneP)) / (as.numeric(OneA) + as.numeric(OneP))) %>%
  distinct()

JustDVPromoterOne = bind_rows(JustDorsalPromoter_One, JustVentralPromoter_One)

JustDVPromoterOne$Location <- factor(JustDVPromoterOne$Location,
                                  levels = c("Dorsal", "Ventral"),ordered = TRUE)

ggplot(JustDVPromoterOne,
       aes(x = reorder(as.factor(Name),- as.numeric(PositionalScore)),
           y = as.numeric(PositionalScore),
           fill = Location)) +
  #guides(fill= FALSE) +
  geom_col(position = "stack", colour = "black") +
  ylim(-1, 1) +
  xlab("D-V Promoters") +
  ylab("PositionalScore") +
  ggtitle("020918_2reps_nodups_SingleHalves_DV_Promoters_PositionalScore_Bargraph_Position_0118newgeneList") +
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
  png('020918_2reps_nodups_SingleHalves_DV_Promoters_PositionalScore_Bargraph_Position_0118newgeneList.png', width = 5000, height = 2000, units = "px",  res=300)
dev.off()