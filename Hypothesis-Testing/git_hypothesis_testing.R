#Author: Mehdin Masinovic

####Info####
"
 In my thesis, the density of all analysed short tandem repeats is tested on whether there are differences among zones. 
 In this script, I want to describe:
 
 1. How I conducted the hypothesis testing
 2. Why I chose Kruskal Wallis and Dunn's test (Test for normality: Kolmogorov-Smirnov Test, Test for homogeneity of variances: Fligner-Killeen Test)
 3. How I used Kruskal Wallis and Dunn's test
 
 
 Similarly to other folders on this GitHub repository, I am only going to present the hypothesis testing using some representatives, not the whole data set. 
 
 First, I read into the data frames and apply a Kruskal Wallis test + Dunns test for post-hoc comparisons. 
 
 I chose these tests, because the data is non-normally distributed and the homogeneity of variance is valid. 
 I test these two via a Kolmogornov-Smirnov and a Fligner-Killeen test, respectively. 
 
 The result is obtained using obtain_pvals() and obtain_pvals_motifs(), which is a list of the p-value and Dunn's test.
 
 To transform them into tables as they are seen in the thesis, I wrote the functions get_pval_df() and get_Dunns_test_df() to transform the data. 
 
 
 "

####Libraries & working directory####
rm(list=ls())
setwd("~/Hypothesis-Testing/")
library(plyr)
library(dplyr)
library(FSA)
library(ggplot2)
library(tidyverse)

####Loading functions#####
obtain_pvals <- function(paths, strs, cutdata = T){
        count = 1
        stats <- list()
        for(i in paths){
                orig_data = readRDS(i)
                #orig_data = orig_data[-which(orig_data$Zone == 7),]
                colnames(orig_data)[2] <- "Length_As"
                
                if(cutdata){
                        sums <- ddply(orig_data, .(Length_As), summarise, Sum = sum(n))
                        print(paste0("STR: ", strs[count]))
                        
                        print(paste0("Limiting value: ", sums[1,2]*0.01))
                        remains <- diff(which(sums$Sum > (sums[1,2]*0.01)))
                        
                        if(length(which(remains != 1))==0){
                                limitnr <- sums$Length_As[length(remains) + 1]
                        } else {
                                limitnr <- sums$Length_As[which(diff(which(sums$Sum > (sums[1,2]*0.01))) != 1)]
                        }
                        print(paste0("Limiting number: ", limitnr))
                        my_data = orig_data[which(orig_data$Length_As <= limitnr),]
                } else {my_data = orig_data}
                
                print(ks.test(my_data$renormdens, "pnorm"))
                print(fligner.test(renormdens ~ as.factor(Zone), data = my_data))
                
                my_data$Zone <- my_data$Zone - 1
                my_data$Zone[which(my_data$Zone == 0)] <- "HS"
                my_data$Zone[which(my_data$Zone == 6)] <- "OS"
                stats[count] <- list(c(list(strs[count]), list(kruskal.test(renormdens~as.factor(Zone), data = my_data)),
                                       list(dunnTest(renormdens~as.factor(Zone), data = my_data))))
                count = count + 1
        }
        
        return(stats)
}

obtain_pvals_motifs <- function(paths, strs){
        count = 1
        stats <- list()
        motifdata <- data.frame()
        
        for(i in paths){
                orig_data = readRDS(i)
                #orig_data = orig_data[-which(orig_data$Zone == 7),]
                motifdata <- rbind(motifdata, cbind(STR = strs[count], orig_data))
                count = count + 1
        }  
        print(ks.test(motifdata$densities, "pnorm"))
        print(fligner.test(densities ~ as.factor(Zone), data = motifdata))
        
        motifdata$Zone <- motifdata$Zone - 1
        motifdata$Zone[which(motifdata$Zone == 0)] <- "HS"
        motifdata$Zone[which(motifdata$Zone == 6)] <- "OS"
        stats <- list(c(list(strs[count]), list(kruskal.test(densities~as.factor(Zone), data = motifdata)),
                        list(dunnTest(densities~as.factor(Zone), data = motifdata))))
        
        return(stats)
}

get_pval_df <- function(mylist){
        pvals <- unlist(lapply(mylist, function(l) c(l[[1]], c(l[[2]][[1]], l[[2]][[3]]))    )     ) 
        names(pvals) <- NULL
        
        df_pvals <- data.frame(STR = pvals[c(T, F, F)],
                               ChiSquared = pvals[c(F, T, F)],
                               PValue = pvals[c(F, F, T)]
        )
        
        df_pvals[,2] <- as.numeric(as.character(df_pvals[,2]))
        df_pvals[,3] <- as.numeric(as.character(df_pvals[,3]))
        
        
        return(df_pvals)
}

get_dunnsTest_df <- function(mylist){
        
        duns <- lapply(mylist, function(l){if(l[[2]][[3]] < 0.05){
                return(data.frame(STR = l[[1]],
                                  overall.pval = l[[2]][[3]],
                                  l[[3]][[2]][,c(1,2,3,4)]))}     
        })
        
        duns2 <- lapply(compact(duns), function(l){
                l = l[which(l[,5]<=0.05),]
                if(nrow(l) == 0) {
                        return(l)
                } else {
                        l[,5] <- format.pval(pv = l[,5], digits = 1, eps = 0.001, nsmall = 3)
                        l[,4] <- round(l[,4], 2)
                        return(l[,-2])
                }
        })
        
        
        return(compact(duns2))
}

####Hypothesis testing####

#Monos
monos <- c("A", "C")
monos_path <- paste0("df_hg38_2k_", monos, ".Rds")
monos_uncut <- obtain_pvals(paths = monos_path, strs = monos, cutdata = F)
saveRDS(object = get_pval_df(monos_uncut), file = "kw_monos.Rds")
saveRDS(object = get_dunnsTest_df(monos_uncut), file = "dunns_monos.Rds")

#Dis
dis = c("AT", "CG")
dis_path <- paste0("df_hg38_2k_", dis, ".Rds")
dis_uncut <- obtain_pvals(paths = dis_path, strs = dis, cutdata = F)
saveRDS(object = get_pval_df(dis_uncut), file = "kw_dis.Rds")
saveRDS(object = get_dunnsTest_df(dis_uncut), file = "kw_dis.Rds")

#Tris
tris <- c("CCG", "ATT")
tris_path <- paste0("df_hg38_2k_", tris, ".Rds")
tris_uncut <- obtain_pvals(paths = tris_path, strs = tris, cutdata = F)
saveRDS(object = get_pval_df(tris_uncut), file = "kw_tris.Rds")
saveRDS(object = get_dunnsTest_df(tris_uncut), file = "kw_tris.Rds")

#Tetras
tetras <- c("AACT", "AAGA", "AATC", "ACAA")
tetras_path <- paste0("df_hg38_2k_", tetras, ".Rds")
tetras_kw <- obtain_pvals_motifs(paths = tetras_path, strs = tetras)
saveRDS(object = get_pval_df(tetras_kw), file = "kw_tetras.Rds")
saveRDS(object = get_dunnsTest_df(tetras_kw), file = "kw_tetras.Rds")

#Motifs
motifs <- c("AAAAC", "AAAAG")
motifs_path <- paste0("df_hg38_2k_", motifs, ".Rds")
motifs_kw <- obtain_pvals_motifs(paths = motifs_path, strs = motifs)
saveRDS(object = get_pval_df(motifs_kw), file = "kw_motifs.Rds")
saveRDS(object = get_dunnsTest_df(motifs_kw), file = "kw_motifs.Rds")











