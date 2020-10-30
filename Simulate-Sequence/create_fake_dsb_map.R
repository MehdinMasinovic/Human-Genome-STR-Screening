#Author: Mehdin Masinovic

####Info####
"
 When using STRAH, a DNA sequence and the coordinates of recombination hotspots have to be provided. 
 As part of the process of simulating a DNA sequence, it is also necessary to create this position matrix. 
 In this program, we will create a fake double-stranded-break map (DSB map) that will provide the coordinates of all hotspots inside the simulated DNA sequence.
"
####Function#####
create_fake_dsb_map <- function(sequence = "", chrs = "chrsim"){

  sequence <- read.fasta(sequence, as.string = T, seqonly = T)
  
  #We check for our "Z" at the beginning and the end of a hotspot zone. 
  
  #We load our python file:
  
  source_python('/home/roots/Desktop/projects/Project_STRAH/somecodeandwsp/fuzzyMatchNMtest.py')
  str <- "Z"
  
  positionsnm <- (fuzzyMatch(str, paste(unlist(sequence)))) + 1 #We add one to each position because pyhton's indexing starts from 0. 
  
  startpos <- positionsnm[seq(1,length(positionsnm), 2)]
  endpos <- positionsnm[seq(2, length(positionsnm), 2)]
  
  
  dsb_map_fake <- as.data.frame(cbind(chrom=chrs,startpos,endpos))
  
  dsb_map_fake[,1] <- as.character( dsb_map_fake[, 1] )
  dsb_map_fake[,2] <- as.numeric(as.character( dsb_map_fake[, 2] ))
  dsb_map_fake[,3] <- as.numeric(as.character( dsb_map_fake[, 3] ))

  return(dsb_map_fake)
}


