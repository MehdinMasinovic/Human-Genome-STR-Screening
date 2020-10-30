#Author: Mehdin Masinovic

####Info####
"
 While updating STRAH to be able to detect repeats with a repeat unit greater than two, we detected a coding error within the code of STRAH. 
 The consequence of the error is explained in the thesis chapter 8.4 (\"Correction\").
 In the following code, I want to list how I obtained the plots and describe the workflow. 
 
 I am including the incorrect and correct version of STR_analysis() for future reference and reproducibility. 
"

####Libraries & working directory####

rm(list=ls())
setwd("~/Correction")
library(ggplot2)
library(cowplot)
library(grid)
source("create_plots.R")
source("STR_detection.R")

####Loading functions####
applying_STRAH_old <- function(letter){
        library(BSgenome.Hsapiens.UCSC.hg19)
        source("STR_analysis_old.R")
        
        
        system(paste0("mkdir old_hg19_", letter))
        setwd(paste0("old_hg19_", letter))
        
        
        pos_matrix <- readRDS("position_matrix_hg19.Rds")
        for(i in 1:nrow(pos_matrix)) {
                starts <- Sys.time()
                res <- STR_analysis(nr.STRs = 6, nr.mismatch = 0, STR = letter, lens.grey = 0:5*1000,  reverse.comp = FALSE, pos_matrix = pos_matrix[i,])
                ends <- Sys.time()
                assign(as.character(pos_matrix[i,1]), res)
                assign(paste0(as.character(pos_matrix[i,1]), "_time"), ends - starts)
                
                
                saveRDS(get(as.character(pos_matrix[i,1])), paste0(pos_matrix[i,1],"_res.Rds"))
                saveRDS(get(paste0(as.character(pos_matrix[i,1]), "_time")), paste0(as.character(pos_matrix[i,1]), "_time.Rds"))
                rm(res, ends, starts)
        }
}

applying_STRAH_new <- function(letter){
        library(BSgenome.Hsapiens.UCSC.hg19)
        source("STR_analysis_new.R")
        
        system(paste0("mkdir new_hg19_", letter))
        setwd(paste0("new_hg19_", letter))
        
        pos_matrix <- readRDS("position_matrix_hg19.Rds")
        for(i in 1:nrow(pos_matrix)) {
                starts <- Sys.time()
                res <- STR_analysis(nr.STRs = 6, max.nr.STRs = 200, nr.mismatch = 0, STR = letter, addToHs = 500, lens.grey = 0:5*1000,  reverse.comp = FALSE, pos_matrix = pos_matrix[i,], species=BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
                ends <- Sys.time()
                assign(as.character(pos_matrix[i,1]), res)
                assign(paste0(as.character(pos_matrix[i,1]), "_time"), ends - starts)
                
                
                saveRDS(get(as.character(pos_matrix[i,1])), paste0(pos_matrix[i,1],"_res.Rds"))
                saveRDS(get(paste0(as.character(pos_matrix[i,1]), "_time")), paste0(as.character(pos_matrix[i,1]), "_time.Rds"))
                rm(res, ends, starts)
        }
}

combine_and_save <- function(old, new, panel, STR){

        prow <- cowplot::plot_grid(old + theme(legend.position = "none"), 
                                   new + theme(legend.position = "none", axis.title.y = element_blank()),
                                   align = 'vh',
                                   labels = c("A", "B"),
                                   hjust = -1,
                                   nrow = 1)
        legend <- cowplot::get_legend(old + theme(legend.box.margin = margin(0, 0, 0, 0)))
        onebyone <- cowplot::plot_grid(prow, legend, rel_widths = c(3,.2))

        ggsave(
                filename=paste0("Poly-", STR, "-density-panel-", panel, "-old-new.png"),
                plot = onebyone,
                device = "png",
                path = NULL,
                scale = 1,
                width = 34,
                height = 16,
                units = c("cm"),
                limitsize = TRUE,
        )
}

####Applying STRAH###
"
 For anyone curious on how both version (before and after correction) were run, I have included the code below.
 In case someone does not want to run the STRAH application, it is possible to just skip to the plotting section below.
"

#Old STRAH
applying_STRAH_old(letter = "A")

#New STRAH
applying_STRAH_new(letter = "A")

####Analysing the result: Poly-As####
"The result of the above run can be transformed using the data transformation steps from 
 the \"Applying-STRAH\" folder.
 In this script, we already provide the transformed data.

"

hg19new2kA = readRDS("/home/roots/Desktop/projects/Project_STRAH/STRAH/my_codes/strah_writing_thesis/github_STRAH/Correction/hg19_new_2k_A.Rds")
hg19old2kA = readRDS("/home/roots/Desktop/projects/Project_STRAH/STRAH/my_codes/strah_writing_thesis/github_STRAH/Correction/hg19_old_2k_A.Rds")

colnames(hg19new2kA)[2] <- "Length_Polys"
colnames(hg19old2kA)[2] <- "Length_Polys"

#Panel A
oldpA <- plotfigA(dataframe = hg19old2kA, letter = "A")
newpA <- plotfigA(dataframe = hg19new2kA, letter = "A")
combine_and_save(old = oldpA, new = newpA, panel = "A", STR = "A")

#Panel B1
oldBp1 <- plotfigB1panel(dataframe = hg19old2kA, letter = "A")
newBp1 <- plotfigB1panel(dataframe = hg19new2kA, letter = "A")
combine_and_save(old = oldBp1, new = newBp1, panel = "B1", STR = "A")

#Panel B2
oldBp2 <- plotfigB2panel(dataframe = hg19old2kA, letter = "A")
newBp2 <- plotfigB2panel(dataframe = hg19new2kA, letter = "A")
combine_and_save(old = oldBp2, new = newBp2, panel = "B2", STR = "A")

#Panel C
oldpC <- plotfigC(dataframe = hg19old2kA, letter = "A")
newpC <- plotfigC(dataframe = hg19new2kA, letter = "A")
combine_and_save(old = oldpC, new = newpC, panel = "C", STR = "A")

#Panel D
oldpD <- plotfigD(dataframe = hg19old2kA, letter = "A")
newpD <- plotfigD(dataframe = hg19new2kA, letter = "A")
combine_and_save(old = oldpD, new = newpD, panel = "D", STR = "A")


####Analysing the result: Poly-Ts####


hg19new2kT = readRDS("hg19_new_2k_T.Rds")
hg19old2kT = readRDS("hg19_old_2k_T.Rds")

colnames(hg19new2kT)[2] <- "Length_Polys"
colnames(hg19old2kT)[2] <- "Length_Polys"

#Panel A
oldpA <- plotfigA(dataframe = hg19old2kT, letter = "T")
newpA <- plotfigA(dataframe = hg19new2kT, letter = "T")
combine_and_save(old = oldpA, new = newpA, panel = "A", STR = "T")

#Panel B1
oldBp1 <- plotfigB1panel(dataframe = hg19old2kT, letter = "T")
newBp1 <- plotfigB1panel(dataframe = hg19new2kT, letter = "T")
combine_and_save(old = oldBp1, new = newBp1, panel = "B1", STR = "T")

#Panel B2
oldBp2 <- plotfigB2panel(dataframe = hg19old2kT, letter = "T")
newBp2 <- plotfigB2panel(dataframe = hg19new2kT, letter = "T")
combine_and_save(old = oldBp2, new = newBp2, panel = "B2", STR = "T")

#Panel C
oldpC <- plotfigC(dataframe = hg19old2kT, letter = "T")
newpC <- plotfigC(dataframe = hg19new2kT, letter = "T")
combine_and_save(old = oldpC, new = newpC, panel = "C", STR = "T")

#Panel D
oldpD <- plotfigD(dataframe = hg19old2kT, letter = "T")
newpD <- plotfigD(dataframe = hg19new2kT, letter = "T")
combine_and_save(old = oldpD, new = newpD, panel = "D", STR = "T")









