#Author: Mehdin Masinovic

####Info####
"
 In this file, we present how the simulation of a DNA sequence works. 
 It refers to the thesis chapter 4.2 (\"Validation of STRAH ’s output\")
 
 It includes:
 
 1. Loading necessary libraries and functions.
 2. Generating the simulated sequence.
 3. Applying STRAH on the sequence.
 4. Transforming the result.
 5. Plotting it.
 
 In order to analyze the sequence with STRAH, it is not necessary to download STRAH with all its dependencies. 
 Instead, I provided the two most important functions that just have to be sourced into R.
"

####Libraries & Working Directory####
library(seqinr)
library(reticulate)
library(Biostrings)
library(dplyr)
library(plyr)
library(ggplot2)
library(cowplot)
library(grid)

setwd("~/Simulate-Sequence")

####Sourcing necessary functions#####
"
 The sequence generating functions are a set of self-written R functions. 
 They are saved in seperate files, therefore, we source them here:
"

source("simulation_plots.R")
source("create_seq.R")
source("create_field.R")
source('generate_sequence.R')
source('create_fake_dsb_map.R')
reticulate::source_python("fuzzyMatchNM.py")

"
 We use two of the most crucial STRAH functions for the analysis of the simulated sequence:
"
source('STR_analysis.R')
source('STR_detection.R')

#Function to transform result into dataframe
obtain_table <- function(data, pathtodsbmap){
  genome <- data
  
  #What is the total HS-zone length for each chromosome?
  dsb_length <- readRDS(pathtodsbmap)
  hs_length <- c()
  for(i in unique(dsb_length$chrom)){hs_length <- c(hs_length, sum(((dsb_length[,3]-dsb_length[,2])+1000)[which(dsb_length$chrom == i)]))}

  ####We prepare our dataframe for transforming####
  #We subset the zones and repeat length, and combine them into a df.
  Zone <- genome$Zone
  Length_As <- genome$`Length of STR stretch in bp`
  testdf <- data.frame(Length_As, Zone)
  colnames(testdf) <- c("Length_As","Zone")
  #head(testdf)
  
  #We group the dataframe by zone and count the repeats found for each zone.
  test = NULL
  test = dplyr::group_by(testdf, Zone, Length_As) %>%
    dplyr::summarise(n=length(Length_As))
  #head(test)
  
  
  #We add a new column into our dataframe for the zone length.
  test = cbind(test, zone_length = rep(0, nrow(test)))
  #We start with the HS-zone.
  test[which(test$Zone == 1),]$zone_length <- rep(sum(hs_length), length(test[which(test$Zone == 1),]$zone_length))
  #print(test[which(test$Zone == 1),])
  
  #We continue with the surrounding zones.
  for(i in 2:6){
    test[which(test$Zone == i),]$zone_length <- rep(nrow(dsb_length)*2000, length(test[which(test$Zone == i),]$zone_length))
  }
  #print(test[which(test$Zone == 2),])
  
  #We end with the outside zone
  oz <- 200000 - nrow(dsb_length)*2000 * 5 - sum(hs_length)
  #only use this line if not comparing to old version.
  test[which(test$Zone == 7),]$zone_length <- rep(oz, length(test[which(test$Zone == 7), ]$zone_length))
  #print(test[which(test$Zone == 7),])
  
  ####We transform our data for plotting####
  #Data for Panel B1 (Poly-A Density):
  test = cbind(test, panelB1dens = test$n * 0)
  #Replace 7 by 6 if you want to compare to paper
  for ( i in 1 : 7 ) {
    test[test$Zone == i,]$panelB1dens <- (test[test$Zone == i,]$n)/test[test$Zone == i,]$zone_length
  }
  
  #Data for Panel B2 (Density of As):
  test = cbind(test, panelB2dens = test$n * 0)
  #Replace 7 by 6 if you want to compare to paper
  for ( i in 1 : 7 ) {
    test[test$Zone == i,]$panelB2dens <- (test[test$Zone == i,]$n * test[test$Zone == i,]$Length_As)/test[test$Zone == i,]$zone_length
  }
  
  #Data for Panel C (Poly-A density per base pair). 
  "For panel C, we will reuse the densities of figure A."
  
  ##We collapse our data:
  #We add a new column for our densities.
  test = cbind(test, densities = (test$n)/(test$zone_length))
  head(test)
  
  #We collapse the densities for all repeats >= 10 and save it into the row for polyA-Length == 10.
  collapsed <- plyr::ddply(test[test$Length_As>=10,], .(Zone), summarize, n = sum(n), panelB1dens = sum(panelB1dens), panelB2dens = sum(panelB2dens), densities=sum(densities))
  test[test$Length_As==10,][,c(3,5:7)] <- collapsed[,-1]
  test = test[test$Length_As <= 10,]
  #print(test, n= Inf)
  
  #Data for Panel A:
  #We add a new column for our renormalized densities.
  test = cbind(test, renormdens = test$n * 0)
  #We renormalize each base pair group by its sum of densities.
  for(i in 6:10){
    test[test$Length_As == i,]$renormdens <- test[test$Length_As == i,]$densities / sum(test[test$Length_As == i,]$densities)
  }
  #Data for Panel D (Poly-A densities per poly-A tract length (i) displayed as fractions):
  test = cbind(test, panelDdens = test$n * 0)
  for ( i in 6:10 ) {
    
    test[test$Length_As==i,]$panelDdens <- test[test$Length_As == i,]$densities / test[test$Length_As == i,]$densities[1]
  }
  colnames(test)[2] <- "Length_Polys"
  return(test)
  
}



####Simulating a sequence####

#We create a sequence full with mononucleotide repeats (poly-As)
sequence <- generate_sequence(lengthofseq = 200000, nr_HS = 10, le_HS = 2000,
                              nr_zones = 5, le_zones = 2000, dens_As = c(Z7 = 0.25, rep(0.25,5) , Z1 = 0.25), 
                              STR= "A", setofAs= c(seq(6,10,1)), set_distr = c(0.2,0.2,0.2,0.2,0.2))

#We extract starting & ending positions + the sequence string.
positions <- unlist(sequence[2])
sequence <- unlist(sequence[1])

#We write it out as a fasta file to load it into create_fake_dsb_map.
write.fasta(sequence, "Simulation", file.out = "SimSeq.txt")
dsb_map_fake <- create_fake_dsb_map(sequence = "SimSeq.txt")
saveRDS(dsb_map_fake, "dsb_map_fake.Rds")

#We replace Z with T.
sequence <- stringr::str_replace_all(sequence, "Z", "T")
fuzzyMatch("Z", sequence) #No Z contained
write.fasta(sequence, "Simulation", file.out = "SimSeq.txt")


####Applying STRAH, analyzing the result, saving plot####
#We apply STRAH.
result <- STR_analysis(seqName = "SimSeq.txt", nr.STRs = 6, nr.mismatch = 0, chrs = "chrsim", STR = "A", lens.grey = 0:5*1000,
                       start.position = 1, end.position = 200000, reverse.comp = FALSE, addToHs = 500,
                       dsb_map = dsb_map_fake, species = "Simulation")

result <- lapply(result, function(x) unlist(x, use.names = F))
#saveRDS(result, "Simulation-Result-Correct.Rds")

#We transform the result obtained from STRAH and plot it.
check <- obtain_table(result, "dsb_map_fake.Rds")
simpA <- plotfigA(dataframe = check, letter = "A")
simpC <- plotfigC(dataframe = check, letter = "A")


#We combine the plots same as in the thesis.
legendA <- cowplot::get_legend(simpA + theme(legend.box.margin = margin(0, 0, 0, 0)))
prow <- cowplot::plot_grid(simpA + theme(legend.position = "none"), 
                           simpC + theme(legend.position = "none"),
                           align = 'vh',
                           labels = c("A", "B"),
                           hjust = -1,
                           nrow = 1)
onebyone <- cowplot::plot_grid(prow, legendA, rel_widths = c(1.5,.2))
#grid.newpage(); onebyone


#We save the plot.
ggsave(
  filename="Simulation-Plot-Correct.png",
  plot = onebyone,
  device = "png",
  path = NULL,
  scale = 1,
  width = 38,
  height = 16,
  units = c("cm"),
  limitsize = TRUE,
)




