#Author: Mehdin Masinovic

####Info####
"
 The data analysis of all short tandem repeat types in the thesis results into ~15GB of data. 
 For this reason, I am providing not the datasets that result from applying STRAH, but rather how STRAH was applied (with what arguments, what repeat types, etc.)
 Here, I show how I applied STRAH on the full set of distinct dinucleotide repeats. 
 The number of distinct combinations is 6, therefore, the process of running STRAH is automated via a function (applying_STRAH()).
 
 I provide the functions that are specifically used in my analysis for reproducibility (in case STRAH will be updated in the future), see STR_analysis() and STR_detection().
"

####Libraries & Working Directory####
library(STRAH)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
source("STR_analysis.R")
source("STR_detection.R")

setwd("~/Applying-STRAH")

####Loading function#####

applying_STRAH <- function(letter){
  
  #Create a folder for every repeat type that is analyzed.
  system(paste0("mkdir exact_hg38_", letter))
  setwd(paste0("exact_hg38_", letter))
  
  #Load position matrix.
  pos_matrix <- readRDS("position_matrix_hg38.Rds")
  
  #For every chromosome in the position matrix, apply STRAH and detect all occurences of the repeat.
  for(i in 1:nrow(pos_matrix)) {
    #Record time
    starts <- Sys.time()
    #Apply Strah
    res <- STR_analysis(nr.STRs = 3, max.nr.STRs = 200, nr.mismatch = 0, STR = letter, addToHs = 500, lens.grey = 0:5*1000,  reverse.comp = FALSE, pos_matrix = pos_matrix[i,], species=BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
    ends <- Sys.time()
    #Assign result and time to variable
    assign(as.character(pos_matrix[i,1]), res)
    assign(paste0(as.character(pos_matrix[i,1]), "_time"), ends - starts)
    #Save result and time into file
    saveRDS(get(as.character(pos_matrix[i,1])), paste0(pos_matrix[i,1],"_res.Rds"))
    saveRDS(get(paste0(as.character(pos_matrix[i,1]), "_time")), paste0(as.character(pos_matrix[i,1]), "_time.Rds"))
    #Remove result and time from environment for the chromosome analyze (to save space)
    rm(res, ends, starts)
  }
}

####Apply STRAH####
#We sequentially run all repeat types. STRAH is not parallelized, but can be run simultaneously on several threads, which will speed up the run.
mySTR = "AC"
applying_STRAH(letter = mySTR)
mySTR = "AG"
applying_STRAH(letter = mySTR)
mySTR = "AT"
applying_STRAH(letter = mySTR)
mySTR = "CG"
applying_STRAH(letter = mySTR)
mySTR = "CT"
applying_STRAH(letter = mySTR)
mySTR = "GT"
applying_STRAH(letter = mySTR)







