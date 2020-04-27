####################################################################
######### Figure5_and_S5.r                                 #########
######### Import output from wigfilesnapshots              #########
######### output all TF signal traces for each region      #########
####################################################################


library("plyr")
library("dplyr")
library("ggplot2")
library("grid")

#Set working directory
setwd("/021318_TF_Wigfilesplitter/")
RegionFile = read.delim2("011618_ReviewsRevised_RevisedAPDV_enhPromRegions.txt", sep = "\t", stringsAsFactors = FALSE, header = TRUE)
colnames(RegionFile) = c("chrom", "start", "end", "Type", "Location", "Name", "use", "length", "fbID", "source")
RegionFile = RegionFile %>% 
  filter(use == "yes")
  

#Make a list of the files
FileList = list.files("021318_TF_Wigfilesplitter/")

TestList = FileList[0:4]
LengthofList = as.numeric(n_distinct(FileList))
counter = 1

BicoidList = list()
HbList = list()
KRList = list()
GTList = list()
ZLD_st3List = list()
ZLD_st4List = list()
ZLD_st5List = list()
CADrList = list()
KNIList = list()
AntList = list()
PostList = list()
DnaseIList = list()
EnhancerNameList = list()

while(counter <= LengthofList){

  ExtendedRegionFile = read.delim2(paste(as.character(FileList[counter])), sep = "\t", stringsAsFactors = FALSE, header = FALSE)
  #Rename the columns
  colnames(ExtendedRegionFile) = c("chrom", "start", "end", 'whole', 'ant', 'post', 'dnaseI', 'Bcd',	'Hb',	'KR',	'GT',	'ZLD_st3', 'ZLD_st4', 'ZLD_st5', 'CADr', 'KNI', "chr1", "start1", "end1",'DnaseI2')

  ExtendedRegionFile[ExtendedRegionFile == 'nan'] <- 0 #remove the nan
  ExtendedRegionFile[ExtendedRegionFile == 'NaN'] <- 0 #remove the nan

  #Make the title the name of the file
  titleLabel = paste(as.character(FileList[counter]))

  #Extract out the region name from the file name
  NameofRegion = strsplit(titleLabel, '_WigFileSnapshot.txt')

  #Look up the name in the RegionFile and obtain the chr, start, and end coordinates
  RegionCoord = RegionFile %>%
    filter(Name == paste(as.character(NameofRegion[1])))

  RegionChr = as.character(RegionCoord$chrom)
  RegionStart = as.numeric(RegionCoord$start)
  RegionEnd = as.numeric(RegionCoord$end)

  #make the output file name
  fileLabel = paste(NameofRegion,"WigFileSnapshot.png", sep = "_")
  filelabelpath = paste("~/042418_TFWig/", fileLabel, sep = "")
  print("fileLabel")

  
   #strategy: I am going make a new graph for each column and then arange them in a 1x9 strip using grid
  #         Then I will merge all the grids into one panel for each region

  png(filelabelpath, width = 1500, height = 5000, units = "px",  res=300)
  grid.newpage()
  AntpostATAC = ggplot(ExtendedRegionFile,
                       aes(x = as.numeric(ExtendedRegionFile$start))) +
    annotate("rect", xmin = RegionStart[1], xmax = RegionEnd[1],   ymin = 0, ymax = Inf,   alpha = 0.2, fill = "slategray4") +
    geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$ant)), colour = "darkorange1", lwd = 1, fill = "darkorange1", alpha = 0.3)  +
    geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$post)), colour = "dodgerblue3", lwd = 1, fill = "dodgerblue3", alpha = 0.2)  +
    xlab(paste(RegionChr[1])) +
    #ylim(0,3000) +
    scale_y_continuous(limits = c(0,5000), breaks=c(0,5000)) +
    ggtitle(as.character(NameofRegion)) +
    #ylab("Normalized Wig Score") +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill= NA),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 1, size = 20),
          axis.text.x = element_blank(),
          axis.text.y  = element_text(size=20),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20),
          plot.margin = grid::unit(c(0.5, 0.5,
                                     0.5, 0.5), "cm"))
  DnaseI_plot = ggplot(ExtendedRegionFile,
                       aes(x = as.numeric(ExtendedRegionFile$start)))+
    annotate("rect", xmin = RegionStart[1], xmax = RegionEnd[1],   ymin = 0, ymax = Inf,   alpha = 0.2, fill = "slategray4") +
    geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$DnaseI2)), colour = "green", lwd = 1, fill = "green", alpha = 0.3, show.legend = TRUE)  +
    #xlab(paste(RegionChr[1])) +
    #ylim(0,3000) +
    scale_y_continuous(limits = c(0,3000), breaks=c(3000)) +
    #ylab("Normalized Wig Score") +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill= NA),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 1, size = 20),
          axis.text.x = element_blank(),
          axis.text.y  = element_text(size=20),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20))

  BCD_plot = ggplot(ExtendedRegionFile,
                    aes(x = as.numeric(ExtendedRegionFile$start)))+
    annotate("rect", xmin = RegionStart[1], xmax = RegionEnd[1],   ymin = 0, ymax = Inf,   alpha = 0.2, fill = "slategray4") +
    geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$Bcd)), colour = "indianred4", lwd = 1, fill = "indianred4", alpha = 0.3, show.legend = TRUE)  +
    #xlab(paste(RegionChr[1])) +
    #ylim(0,50) +
    scale_y_continuous(limits = c(0,50), breaks=c(0,50)) +
    #ylab("Normalized Wig Score") +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill= NA),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 1, size = 20),
          axis.text.x = element_blank(),
          axis.text.y  = element_text(size=20),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20))

  Hb_plot = ggplot(ExtendedRegionFile,
                   aes(x = as.numeric(ExtendedRegionFile$start)))+
    annotate("rect", xmin = RegionStart[1], xmax = RegionEnd[1],   ymin = 0, ymax = Inf,   alpha = 0.2, fill = "slategray4") +
    geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$Hb)), colour = "sienna2", lwd = 1, fill = "sienna2", alpha = 0.3, show.legend = TRUE)  +
    #ylab("Normalized Wig Score") +
    #ylim(0,150) +
    scale_y_continuous(limits = c(0,180), breaks=c(0,180)) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill= NA),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 1, size = 20),
          axis.text.x = element_blank(),
          axis.text.y  = element_text(size=20),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          axis.ticks.x=element_blank(),
          strip.text.y = element_text(size = 20))

  KR_plot = ggplot(ExtendedRegionFile,
                   aes(x = as.numeric(ExtendedRegionFile$start)))+
    annotate("rect", xmin = RegionStart[1], xmax = RegionEnd[1],   ymin = 0, ymax = Inf,   alpha = 0.2, fill = "slategray4") +
    geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$KR)), colour = "goldenrod2", lwd = 1, fill = "goldenrod2", alpha = 0.3, show.legend = TRUE)  +
    #ylab("Normalized Wig Score") +
    #ylim(0,150) +
    scale_y_continuous(limits = c(0,150), breaks=c(0,150)) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill= NA),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 1, size = 20),
          axis.text.x = element_blank(),
          axis.text.y  = element_text(size=20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          axis.ticks.x=element_blank(),
          strip.text.y = element_text(size = 20))

  GT_plot = ggplot(ExtendedRegionFile,
                   aes(x = as.numeric(ExtendedRegionFile$start)))+
    annotate("rect", xmin = RegionStart[1], xmax = RegionEnd[1],   ymin = 0, ymax = Inf,   alpha = 0.2, fill = "slategray4") +
    geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$GT)), colour = "darkseagreen4", lwd = 1, fill = "darkseagreen4", alpha = 0.3, show.legend = TRUE)  +
    #ylab("Normalized Wig Score") +
    #ylim(0,150) +
    scale_y_continuous(limits = c(0,150), breaks=c(0,150)) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill= NA),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 1, size = 20),
          axis.text.x = element_blank(),
          axis.text.y  = element_text(size=20),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          axis.ticks.x=element_blank(),
          strip.text.y = element_text(size = 20))

  ZLD3_plot = ggplot(ExtendedRegionFile,
                     aes(x = as.numeric(ExtendedRegionFile$start)))+
    annotate("rect", xmin = RegionStart[1], xmax = RegionEnd[1],   ymin = 0, ymax = Inf,   alpha = 0.2, fill = "slategray4") +
    geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$ZLD_st3)), colour = "steelblue", lwd = 1, fill = "steelblue", alpha = 0.3, show.legend = TRUE)  +
    #geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$ZLD_st4)), colour = "steelblue3", lwd = 1, fill = NA, alpha = 0.3, show.legend = TRUE)  +
    #geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$ZLD_st5)), colour = "steelblue1", lwd = 1, fill = NA, alpha = 0.3, show.legend = TRUE)  +
    #ylab("Normalized Wig Score") +
    #ylim(0,500) +
    scale_y_continuous(limits = c(0,500), breaks=c(0,500)) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill= NA),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 1, size = 20),
          axis.text.x = element_blank(),
          axis.text.y  = element_text(size=20),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          axis.ticks.x=element_blank(),
          strip.text.y = element_text(size = 20))
  ZLD4_plot = ggplot(ExtendedRegionFile,
                     aes(x = as.numeric(ExtendedRegionFile$start)))+
    annotate("rect", xmin = RegionStart[1], xmax = RegionEnd[1],   ymin = 0, ymax = Inf,   alpha = 0.2, fill = "slategray4") +
    #geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$ZLD_st3)), colour = "steelblue", lwd = 1, fill = NA, alpha = 0.3, show.legend = TRUE)  +
    geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$ZLD_st4)), colour = "steelblue3", lwd = 1, fill = "steelblue3", alpha = 0.3, show.legend = TRUE)  +
    #geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$ZLD_st5)), colour = "steelblue1", lwd = 1, fill = NA, alpha = 0.3, show.legend = TRUE)  +
    #ylab("Normalized Wig Score") +
    #ylim(0,500) +
    scale_y_continuous(limits = c(0,500), breaks=c(0,500)) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill= NA),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 1, size = 20),
          axis.text.x = element_blank(),
          axis.text.y  = element_text(size=20),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          axis.ticks.x=element_blank(),
          strip.text.y = element_text(size = 20))

  ZLD5_plot = ggplot(ExtendedRegionFile,
                     aes(x = as.numeric(ExtendedRegionFile$start)))+
    annotate("rect", xmin = RegionStart[1], xmax = RegionEnd[1],   ymin = 0, ymax = Inf,   alpha = 0.2, fill = "slategray4") +
    #geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$ZLD_st3)), colour = "steelblue", lwd = 1, fill = NA, alpha = 0.3, show.legend = TRUE)  +
    #geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$ZLD_st4)), colour = "steelblue3", lwd = 1, fill = NA, alpha = 0.3, show.legend = TRUE)  +
    geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$ZLD_st5)), colour = "steelblue1", lwd = 1, fill = "steelblue1", alpha = 0.3, show.legend = TRUE)  +
    #ylab("Normalized Wig Score") +
    #ylim(0,500) +
    scale_y_continuous(limits = c(0,500), breaks=c(0,500)) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill= NA),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 1, size = 20),
          axis.text.x = element_blank(),
          axis.text.y  = element_text(size=20),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          axis.ticks.x=element_blank(),
          strip.text.y = element_text(size = 20))

  Cad_plot = ggplot(ExtendedRegionFile,
                    aes(x = as.numeric(ExtendedRegionFile$start)))+
    annotate("rect", xmin = RegionStart[1], xmax = RegionEnd[1],   ymin = 0, ymax = Inf,   alpha = 0.2, fill = "slategray4") +
    geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$CADr)), colour = "mediumpurple4", lwd = 1, fill = "mediumpurple3", alpha = 0.3, show.legend = TRUE)  +
    #ylab("Normalized Wig Score") +
    #ylim(0,50) +
    scale_y_continuous(limits = c(0,50), breaks=c(0,50)) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill= NA),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 1, size = 20),
          axis.text.x = element_blank(),
          axis.text.y  = element_text(size=20),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          axis.ticks.x=element_blank(),
          strip.text.y = element_text(size = 20))

  Kni_plot = ggplot(ExtendedRegionFile,
                    aes(x = as.numeric(ExtendedRegionFile$start)))+
    annotate("rect", xmin = RegionStart[1], xmax = RegionEnd[1],   ymin = 0, ymax = Inf,   alpha = 0.2, fill = "slategray4") +
    geom_ribbon(aes(ymin = 0, ymax = as.numeric(ExtendedRegionFile$KNI)), colour = "navajowhite4", lwd = 1, fill = "navajowhite4", alpha = 0.3, show.legend = TRUE)  +
    #ylab("Normalized Wig Score") +
    #ylim(0,50) +
    scale_y_continuous(limits = c(0,50), breaks=c(0,50)) +
    xlab(paste(RegionChr[1])) +
    theme(panel.border = element_blank(),
          panel.background = element_rect(fill= NA),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(vjust = 1, size = 20),
          axis.text.x = element_text(size=16),
          axis.text.y  = element_text(size=20),
          plot.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20),
          plot.margin = grid::unit(c(1, 0.5,
                                     0.5, 0.5), "cm"))



  #grid.draw(rbind(ggplotGrob(AntpostATAC), ggplotGrob(DnaseI_plot), ggplotGrob(BCD_plot), ggplotGrob(Hb_plot), ggplotGrob(KR_plot), ggplotGrob(GT_plot), ggplotGrob(ZLD3_plot), ggplotGrob(ZLD4_plot), ggplotGrob(ZLD5_plot), ggplotGrob(Cad_plot), ggplotGrob(Kni_plot), size = "last"))
  g1 = ggplotGrob(AntpostATAC)
  g2 = ggplotGrob(DnaseI_plot)
  g3 = ggplotGrob(BCD_plot)
  g4 = ggplotGrob(Hb_plot)
  g5 = ggplotGrob(KR_plot)
  g6 = ggplotGrob(GT_plot)
  g7 = ggplotGrob(ZLD3_plot)
  g8 = ggplotGrob(ZLD4_plot)
  g9 = ggplotGrob(ZLD5_plot)
  g10 = ggplotGrob(Cad_plot)
  g11 = ggplotGrob(Kni_plot)
  g = rbind(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, size = "last")
  #g = rbind(ggplotGrob(AntpostATAC), ggplotGrob(DnaseI_plot), ggplotGrob(BCD_plot), ggplotGrob(Hb_plot), ggplotGrob(KR_plot), ggplotGrob(GT_plot), ggplotGrob(ZLD3_plot), ggplotGrob(ZLD4_plot), ggplotGrob(ZLD5_plot), ggplotGrob(Cad_plot), ggplotGrob(Kni_plot), size = "last")
  g$widths = unit.pmax(g1$widths, g2$widths, g3$widths, g4$widths, g5$widths, g6$widths, g7$widths, g8$widths, g9$widths, g10$widths, g11$widths)
  grid.draw(g)
  #print(p)
  dev.off()
  counter = counter + 1
  }

####################################################################
######### 020318 TF Heatmap                 ############ 
####################################################################

library("plyr")
library("dplyr")
library("ggplot2")
library("grid")
library("ComplexHeatmap")
library("circlize")
#read in file
RegionFile = read.csv2("020318_2Reps_NewGenelist_PE_Analysis_APEnhancers_PositionalScore.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE)

#open region file and order by positional score and only do enhancers
RegionFileEnhOnly = RegionFile %>%
  select(Chr, Start, End, Name, Type, Location, PositionalScore) %>%
  filter(Type == "Enhancer") %>%
  filter(Location %in% c('Anterior', 'Posterior', 'Mostly Post')) %>%
  arrange(desc(as.numeric(PositionalScore)))
colnames(RegionFileEnhOnly) = c("chrom", "start", "end", "Name", "Type", "Location","Positional_score")

#find the length of the list
Length = as.numeric(length(RegionFileEnhOnly$chrom))

#make the names a list
Namelist = as.list(as.character(RegionFileEnhOnly$Name))


# Loop that sums the values in the window
counter = 1

#make a bunch of empty lists -- one for each row of the heatmap
BicoidList_win = list()
HbList_win = list()
KRList_win = list()
GTList_win = list()
ZLD_st3List_win = list()
ZLD_st4List_win = list()
ZLD_st5List_win = list()
CADrList_win = list()
KNIList_win = list()
# AntList_win = list()
# PostList_win = list()
# DnaseIList_win = list()
EnhancerNameList = list()
#WholeList_win = list()

while(counter <= Length){
  
  # Make the extended region file name from the name list
  NameofRegion = paste((Namelist[counter]) , 'WigFileSnapshot.txt', sep = "_")
  
  #load the extended region file
  ExtendedRegionFile = read.delim2(paste(NameofRegion), sep = "\t", stringsAsFactors = FALSE, header = FALSE)
  #Rename the columns
  colnames(ExtendedRegionFile) = c("chrom", "start", "end", 'whole', 'ant', 'post', 'dnaseI', 'Bcd',	'Hb',	'KR',	'GT',	'ZLD_st3', 'ZLD_st4', 'ZLD_st5', 'CADr', 'KNI')
  
  #remove the nan's 
  ExtendedRegionFile[ExtendedRegionFile == 'nan'] <- 0 #remove the nan
  ExtendedRegionFile[ExtendedRegionFile == 'NaN'] <- 0 #remove the nan
  
  #  # Need to find the max peak height for each TF for that region
  # # 'Bcd',	'Hb',	'KR',	'GT',	'ZLD_st3', 'ZLD_st4', 'ZLD_st5', 'CADr', 'KNI'
  EnhancerNameList = append(EnhancerNameList, as.character((Namelist[counter])))
  
  bcdwin = sum(as.numeric(ExtendedRegionFile$Bcd))
  BicoidList_win = append(BicoidList_win, bcdwin)
  
  Hbwin = sum(as.numeric(ExtendedRegionFile$Hb))
  HbList_win = c(HbList_win, Hbwin )
  
  KRwin = sum(as.numeric(ExtendedRegionFile$KR))
  KRList_win = c(KRList_win, KRwin )
  
  GTwin = sum(as.numeric(ExtendedRegionFile$GT))
  GTList_win = c(GTList_win, GTwin)
  
  ZLD_st3win = sum(as.numeric(ExtendedRegionFile$ZLD_st3))
  ZLD_st3List_win = c(ZLD_st3List_win, ZLD_st3win )
  
  ZLD_st4win = sum(as.numeric(ExtendedRegionFile$ZLD_st4))
  ZLD_st4List_win = c(ZLD_st4List_win, ZLD_st4win )
  
  ZLD_st5win = sum(as.numeric(ExtendedRegionFile$ZLD_st5))
  ZLD_st5List_win = c(ZLD_st5List_win, ZLD_st5win )
  
  CADrwin = sum(as.numeric(ExtendedRegionFile$CADr))
  CADrList_win = c(CADrList_win, CADrwin )
  
  KNIwin = sum(as.numeric(ExtendedRegionFile$KNI))
  KNIList_win = c(KNIList_win, KNIwin)
  
  # antwin = sum(as.numeric(ExtendedRegionFile$ant))
  # AntList_win = c(AntList_win, antwin)
  # 
  # postmax = max(as.numeric(ExtendedRegionFile$post))
  # PostList = c(PostList, postmax)
  # 
  # DnaseImax = max(as.numeric(ExtendedRegionFile$dnaseI))
  # DnaseIList = c(DnaseIList, DnaseImax)
  # 
  # WholeMax = max(as.numeric(ExtendedRegionFile$whole))
  # WholeList = c(WholeList, WholeMax)
  
  #Extend the counter
  counter = counter + 1
}

# Then I need to make a matrix with regions as the rows and the columns are TFs 
maxMatrix = matrix(c(EnhancerNameList, BicoidList_win, CADrList_win, KNIList_win, GTList_win, HbList_win, KRList_win, ZLD_st3List_win,  ZLD_st4List_win, ZLD_st5List_win), nrow=length(EnhancerNameList))
colnames(maxMatrix) = c("EnhancerNameList", "bicoid", "caudal", "knirps", "giant", "hunchback", "kruppel", "zelda stage 3", "zelda stage 4", "zelda stage 5")
maxMatrix = as.data.frame(maxMatrix)

#Then I need to transpose the matrix so that the rows are TFs and the columns are regions but keeping the first row as the headers
transposed_maxMatrix = setNames(data.frame(t(maxMatrix[,-1])), maxMatrix[,1])
Max_datamatrix = data.matrix(transposed_maxMatrix)

#Scale by TF 
maxs = apply(Max_datamatrix,1,max)
mins = apply(Max_datamatrix,1,min)
Max_datamatrixScaled = scale(t(Max_datamatrix), center = mins, scale = maxs - mins)


png('020318_XYscript_NewList_TFbinding_heatmap_1kbwin_Colclust.png', width = 4000, height = 2000, units = "px",  res=300) 
ha = HeatmapAnnotation(df = data.frame(PositionalScore = as.numeric(RegionFileEnhOnly$Positional_score)), col = list(PositionalScore = colorRamp2(c(-0.3,0, 0.46),c("dodgerblue2", "white", "darkorange"))))
Heatmap(t(Max_datamatrixScaled), top_annotation = ha, cluster_columns = TRUE, cluster_rows = FALSE, col = colorRampPalette(c("white", "black"))(500), column_title = "1kb windows") 
dev.off()