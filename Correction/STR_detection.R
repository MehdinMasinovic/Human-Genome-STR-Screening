#Author: Mehdin Masinovic

####Info####
"
 STRAH is a set of functions that work together in order to detect short tandem repeats in a DNA sequence and relate them to g.e. recombination hotspots via a position matrix. 
 The main two functions of STRAH are STR_analysis() and STR_detection(). 
 Downloading STRAH together with its dependencies requires a relatively high amount of storage, which is why I decided to provide the main functions of STRAH seperately. 
"
####Function#####

STR_detection = function(seqName, chrs, start.position = NA, end.position = NA, bed_file,
                         pos_matrix, nr.STRs, max.nr.STRs,  nr.mismatch = 0, reverse.comp = F, STR = "A",
                         species=BSgenome.Hsapiens.UCSC.hg19::Hsapiens, translated_regions=F, output_file) {

  min.nr = nr.STRs
  start.position = max(start.position, 1, na.rm = T)
  df <- ""
  if(!missing(seqName)) {
    if(is(seqName, "DNAStringSet")) {
      sequence = seqName
    } else {
      sequence = Biostrings::readDNAStringSet(seqName, format = "fasta")
    }
  }
  else if(!missing(bed_file)){
    bed <- read.table(bed_file,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    chrs <- bed[,1]
    start.position <- bed[,2] + 1
    end.position <- bed[,3]
    if (translated_regions == TRUE){
      original_region <- bed[,4]
    }
    sequence <- Biostrings::DNAStringSet(getflank2(species = species, chrs = chrs, start.position = start.position, end.position = end.position))
  }
  else if(!missing(pos_matrix)){
    chrs <- pos_matrix[,1]
    start.position <- pos_matrix[,2]
    end.position <- pos_matrix[,3]
    sequence <- Biostrings::DNAStringSet(getflank2(species = species, chrs = chrs, start.position = start.position, end.position = end.position))
  }
  else {
    sequence <- Biostrings::DNAStringSet(getflank2(species = species, chrs = chrs, start.position = start.position, end.position = end.position))
    }
  if(reverse.comp) {
    sequence <- Biostrings::reverseComplement(s)
    STR <- toString(Biostrings::reverseComplement(Biostrings::DNAString(STR)))
  }


  seq_name_list <- list()
  nr.matches <- list()
  start.position_list <- list()
  nr.STRs.c_list <- list()
  matches_list <- list()
  original_region_list <- list()

  for (s in 1:(length(sequence))){
    Sys.sleep(0.2)
    if(missing(seqName)) {
      if(missing(bed_file) && missing(pos_matrix)){
        name <- paste0(chrs[s], ":", formatC(start.position[s], format = "fg"), "-", formatC(end.position[s], format = "fg"))
      }
      else{
        name <- paste0(chrs[s], ":", formatC(start.position[s], format = "fg"), "-", formatC(end.position[s], format = "fg"))
      }
    } else {
      name = names(seqName)
    }
    if(missing(seqName)) {
      message(paste0(name, " of ", toString(species), " is under study!"),"\r",appendLF=TRUE)
    } else {
        message(paste0(name, " is under study!"),"\r",appendLF=TRUE)
    }
if(nchar(STR) == 0){
    tp = unlist(strsplit(Biostrings::toString(sequence[s,]), split  ="")) == STR
    #motif <- paste(rep(STR, nr.STRs), collapse = "")
    #motif <- paste(STR, "{", nr.STRs, ",}")
    # tp_2 <- gregexpr(paste(motif, "+", sep=""), Biostrings::toString(s))

    #tp_2 <- Biostrings::matchPattern(motif, Biostrings::toString(s), max.mismatch = nr.mismatch)
    #print(tp_2)

    # tp_2 <- agrep( paste0(motif, '+'), Biostrings::toString(s),
    #                max = list(sub = nr.mismatch),  ignore.case = T, fixed=FALSE, useBytes = TRUE)
    # print(tp_2)

    if (length(which(tp)) == 0 | (length(which(tp))-(nr.STRs-1)) <= 0){
      message(paste0("No STRs are contained in ", name),"\r",appendLF=TRUE)
      message("","\r",appendLF=TRUE)
      next
    }
    tp2 = which(
      sapply(1:(length(tp)-(nr.STRs-1)), function(i) {

        sum(tp[i:(i+(nr.STRs-1))])
      }
      )>= (nr.STRs-nr.mismatch))
    start.pos = tp2[tp[tp2]]
    if(length(start.pos)==0) {
      message(paste0("No STRs are contained in ", name),"\r",appendLF=TRUE)
      next
    }

    add = c(ave(start.pos, cumsum(c(F, diff(start.pos) > 1)), FUN=seq_along) - 1,0)
    nr.STRs.c = add[which(add==0)-1]+nr.STRs


    #### remove duplicates based on the difference in the starting position ####
    ind.double = which(diff(start.pos) == 1)+1
    if (length(ind.double) != 0){ # only if there are more than one start positions of STRs
      start.pos <- start.pos[-ind.double]
    }
    
} else {
  tsequence <- Biostrings::toString(sequence[s,])
  bigdf <- data.frame(start = 0, end= 0, width = 0)
  boundary <- (max.nr.STRs+1)*nchar(STR)
  if(max.nr.STRs > nr.STRs){max.nr.STRs = max.nr.STRs + 1; boundary = max.nr.STRs*nchar(STR)}
  if(max.nr.STRs < nr.STRs){message(paste0("Maximum number of STRs to search for is smaller than the minimum! ", name),"\r",appendLF=TRUE); next}
  
  for(i in max.nr.STRs:nr.STRs){
    
    pattern <- paste0(rep(STR, i), collapse="")
    matches <- Biostrings::matchPattern(pattern,tsequence, max.mismatch = nr.mismatch)
    smalldf <- data.frame(IRanges::ranges(matches))
    
    if(nrow(smalldf) == 0){next} 
    else{
      duplicates = sapply(1:nrow(smalldf), function(j) {
        return(any((smalldf$start[j]>=bigdf$start) & (smalldf$end[j] <= bigdf$end)))
      })
      
      bigdf <- rbind(bigdf, smalldf[!duplicates,])
    }
  }

  bigdf <- bigdf[order(bigdf$start),]
  bigdf <- bigdf[!(bigdf$width >= boundary),]
  
  
  start.pos <- bigdf$start[-1]
  end_pos_STR.c <- bigdf$end[-1]
  nr.STRs.c <- bigdf$width[-1]/nchar(STR)
  
  if(length(start.pos)==0) {
    message(paste0("No STRs are contained in ", name),"\r",appendLF=TRUE)
    next
  }
  
  }

    
    
    
    matches = sapply(1:length(start.pos), function(k) { 
      
      
      if(nchar(STR) == 1){
    # if STR ends at last position of interval
    if ((start.pos[k]+nr.STRs.c[k]-1) >= Biostrings::width(sequence[s,])){
      end_pos_STR <- width(sequence[s,]) # end of interval is end of STR
    }
    else{
      end_pos_STR <- start.pos[k]+nr.STRs.c[k]-1
    }

      temp = unlist(Biostrings::subseq(sequence[s,], start.pos[k], end_pos_STR))
      ### we do not want to end a stretch with anything unequal to STR so as long as the last character is unequal to it, this character will be removed ####
      while(Biostrings::toString(temp[nchar(Biostrings::toString(temp))]) != STR) {
        nr.STRs.c[k] = nr.STRs.c[k]-1
        temp = unlist(Biostrings::subseq(temp, 1, nchar(Biostrings::toString(temp))-1))
      }
      # if STR starts at first position of interval
      } else {end_pos_STR <- end_pos_STR.c[k]}
      
      if (start.pos[k] <= 5){
        start <- start.pos[k]
      }
      else{
        start <- start.pos[k]-5
      }
      # if STR ends with last or less than 5nt before end position of interval
      if ((end_pos_STR+5)>= Biostrings::width(sequence[s,])){
        end <- Biostrings::width(sequence[s,])
      }
      else{
        end <- end_pos_STR+5
      }
      temp.2 = unlist(Biostrings::subseq(sequence[s,], start, end))

      #print(temp.2)
      return(Biostrings::toString(temp.2))
    })
    print(matches)
    message("","\r",appendLF=TRUE)

    seq_name_list <- append(seq_name_list, name)
    nr.matches<- append(nr.matches, length(start.pos))
    nr.STRs.c_list <- append(nr.STRs.c_list,list(nr.STRs.c))
    start.position_list <- append(start.position_list, list(start.pos))
    matches_list <- append(matches_list, list(matches))

    if (translated_regions){
      original_region_list <- append(original_region_list, original_region[s])
      header <- c("chr", "start_STR", "end_STR", "len_STR", "sequence_STR", "translated chr_start_stop", "untranslated chr_start_stop")
    }
    else{
      header <- c("chr", "start_STR", "end_STR", "len_STR", "sequence_STR", "chr_start_stop")
    }

    flush.console()
    if(!missing(output_file)){
      if (is.data.frame(df) == FALSE){
        write.table(rbind(header),paste0(output_file, ".bed"), sep = "\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
      }
      if (translated_regions){
        df <- data.frame(chr_str = rep(chrs[s], length(start.pos)),start_str = start.pos-1, end_str = start.pos + nr.STRs.c-1, len_str = nr.STRs.c, sequence_str = matches, chr_name = rep(name, length(start.pos)), chr_name_original = rep(original_region[s], length(start.pos)))
      }
      else{
        df <- data.frame(chr_str = rep(chrs[s], length(start.pos)),start_str = start.pos-1, end_str = start.pos + nr.STRs.c-1, len_str = nr.STRs.c, sequence_str = matches, chr_name = rep(name, length(start.pos)))
      }
      write.table(df, paste0(output_file, ".bed"), col.names=FALSE, sep="\t", row.names = FALSE, append=TRUE, quote=FALSE)

    }
  }
  if (translated_regions){
    output <- list("Sequence name (translated)" = seq_name_list, "Sequence name (untranslated)" = original_region_list, "Reverse complement" = reverse.comp, "Number of allowed mismatches" = nr.mismatch, "Minimum length" = min.nr, "Number of matches" = nr.matches,
                   "Length of STR stretch in bp" = nr.STRs.c_list, "Start positions" = start.position_list, "Matched segments" = matches_list)
  }
  else{
    output <- list("Sequence Name" = seq_name_list, "Reverse complement" = reverse.comp, "Number of allowed mismatches" = nr.mismatch, "Minimum length" = min.nr, "Number of matches" = nr.matches,
                   "Length of STR stretch in bp" = nr.STRs.c_list, "Start positions" = start.position_list, "Matched segments" = matches_list)
  }

    output <- lapply(output,FUN=function(x) {
    if(length(x) == length(unlist(seq_name_list)) & length(x) != 1){
      names(x) <- unlist(seq_name_list)
    }
    return(x)
  })
  closeAllConnections()
  return(output)
  }

