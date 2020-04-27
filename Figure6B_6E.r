############################################################################################
######### Figure6B_6E.r                                                             ########
######### Import output from randomregion output, list of AP patterned enhancers    ########
######### output all TF signal traces for each region                               ########
############################################################################################

allSINGLECELLSamples <-read.delim2("040518_sciATAC_Randomregv2_output_LinRegWhole.txt", sep = "\t", stringsAsFactors = FALSE)

allPeaksFile <- read.delim2("021518_2reps_REppeaks_overlap_REVISEDGENELIST.txt", sep ="\t", stringsAsFactors = FALSE, header = FALSE)

allPeaks = allPeaksFile %>%
  select(V1, V2, V3, V6, V11, V12, V13, V17, V21)
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

allSINGLECELLSamples$WholePeaks = as.numeric(allPeaksWhole$fold_enrichment)[match(allSINGLECELLSamples$Name, allPeaksWhole$name)]
allSINGLECELLSamples$AntPeaks = as.numeric(allPeaksAnt$fold_enrichment)[match(allSINGLECELLSamples$Name, allPeaksAnt$name)]
allSINGLECELLSamples$PostPeaks = as.numeric(allPeaksPost$fold_enrichment)[match(allSINGLECELLSamples$Name, allPeaksPost$name)]
allSINGLECELLSamples$NoPeak = as.character(allPeaksNone$fold_enrichment)[match(allSINGLECELLSamples$Name, allPeaksNone$name)]

allSamples_SINGLECELL_enhancers_ATACSkew_Peaks = allSINGLECELLSamples %>%
  filter(X011617.use == "yes") %>%
  filter(New.Location.Assignment %in%  c('Anterior', 'Posterior', 'Mostly Post')) %>%
  filter(Type == 'Enhancer') %>%
  mutate(Dotsize = 1.5) %>%
  distinct()

allSamples_SINGLECELL_enhancers_ATACSkew_Peaks$New.Location.Assignment <- factor(allSamples_SINGLECELL_enhancers_ATACSkew_Peaks$New.Location.Assignment,
                                                                                 levels = c('Anterior', 'Posterior', 'Mostly Post'),ordered = TRUE)
ggplot(allSamples_SINGLECELL_enhancers_ATACSkew_Peaks,
       aes(x = reorder(as.factor(Name),- as.numeric(ATACSkewScore)),
           y = as.numeric(ATACSkewScore),
           fill = New.Location.Assignment)) +
  geom_col(position = "stack", colour = "black") +
  ylim(-1, 1) +
  xlab("A-P Enhancers") +
  ylab("PositionalScore") +
  ggtitle("040118_sciATAC_AP_Enhancers_PositionalScore_Bargraph_Position_0118newgeneList_LINREGNORM") +
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
  png('040118_sciATAC_APEnhancers_Bargraph_LinRegNorm.png', width = 5000, height = 2000, units = "px",  res=300)
dev.off()

SignificantPvalue = allSINGLECELLSamples %>%
  select(Name, Type, New.Location.Assignment, PValue, ATACSkewScore) %>%
  arrange(as.numeric(PValue))
######### ######### ######### AP Promoters Bargraphs 
#########  Filter out A vs P enhancers to calculate the positional score

allSamples_SINGLECELL_Promoters_ATACSkew_Peaks = allSINGLECELLSamples %>%
  filter(X011617.use == "yes") %>%
  filter(New.Location.Assignment %in%  c('Anterior', 'Mostly Ant', 'Posterior', 'Mostly Post')) %>%
  filter(Type == 'Promoter') %>%
  mutate(Dotsize = 1.5) %>%
  distinct()

allSamples_SINGLECELL_Promoters_ATACSkew_Peaks$New.Location.Assignment <- factor(allSamples_SINGLECELL_Promoters_ATACSkew_Peaks$New.Location.Assignment,
                                                                                 levels = c('Anterior', 'Mostly Ant', 'Posterior', 'Mostly Post'),ordered = TRUE)
# AP Promoters Bargraph
ggplot(allSamples_SINGLECELL_Promoters_ATACSkew_Peaks,
       aes(x = reorder(as.factor(Name),- as.numeric(ATACSkewScore)),
           y = as.numeric(ATACSkewScore),
           fill = New.Location.Assignment)) +
  geom_col(position = "stack", colour = "black") +
  ylim(-1, 1) +
  xlab("A-P Promoters") +
  ylab("PositionalScore") +
  ggtitle("040518_sciATAC_AP_Promoters_PositionalScore_Bargraph_Position_0118newgeneList_LINREGNORM") +
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
        strip.text.y = element_text(size = 12)) +
  png('040518_sciATAC_APPromoters_Bargraph_LinRegNorm.png', width = 5000, height = 2000, units = "px",  res=300)
dev.off()

############# scATAC vs Halves Positional Skew Score comparison

sci_vs_Halves_skewScore = allSINGLECELLSamples %>%
  mutate(MYDATA_Skewscore = (as.numeric(X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole) - as.numeric(X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole)) / (as.numeric(X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Ant_shifted_lessthan130_10MNorm_linregwhole) + as.numeric(X011718_MERGED_2reps_081817_RemDUP_041217_Bowtie2_ME_JH_112315_20Post_shifted_lessthan130_10MNorm_linregwhole))) %>%
  select(Name, MYDATA_Skewscore, ATACSkewScore, New.Location.Assignment)

cor(sci_vs_Halves_skewScore$MYDATA_Skewscore, as.numeric(sci_vs_Halves_skewScore$ATACSkewScore))

png('040918_sci_vs_Halves_skewScore.png', width = 2000, height = 2000, units = "px",  res=300) 
ggplot(sci_vs_Halves_skewScore,
       aes(x = as.numeric(sci_vs_Halves_skewScore$MYDATA_Skewscore),y= as.numeric(sci_vs_Halves_skewScore$ATACSkewScore), fill = as.factor(New.Location.Assignment))) +
  geom_point(colour="grey5", alpha = 0.5) +
  ggtitle('040918_sci_vs_Halves_skewScore') +
  xlab('Halves ATACSkewScore') +
  ylab('Sc ATACSkewScore') +
  ylim(-1,1) + 
  xlim(-1,1) +
  coord_equal() +
  geom_abline(slope = 1, intercept=0, linetype = 'dotted') +
  theme(panel.border = element_rect(fill = NA),
        panel.background = element_rect(fill= NA),
        axis.title.x = element_text(vjust = 0, size = 12),
        axis.title.y = element_text(vjust = 1, size = 12),
        axis.text.x = element_text(size=15),
        axis.text.y  = element_text(size=15),
        plot.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15))
dev.off()