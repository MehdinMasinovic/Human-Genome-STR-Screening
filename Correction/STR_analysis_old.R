#Author: Mehdin Masinovic

####Info####
"
 STRAH is a set of functions that work together in order to detect short tandem repeats in a DNA sequence and relate them to g.e. recombination hotspots via a position matrix. 
 The main two functions of STRAH are STR_analysis() and STR_detection(). 
 For reproducibility, I am providing the STR_analysis() code before the error correction was made. 
"
####Function#####


STR_analysis = function(nr.STRs = 10, nr.mismatch = 0, chrs = "", STR = "A", lens.grey = lens.grey, start.position = NA,
                        end.position = NA, reverse.comp = FALSE, bed_file = "", pos_matrix = "", output_file="") {
 ind = length.As = pos.As = seq_name_list = nr_matches = original_region_list = list()
 df <- ""
 if(bed_file != ""){
   bed <- read.table(bed_file,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
   chrs <- bed[,1]
   start.position <- bed[,2] +1
   end.position <- bed[,3]
 }
 else if(is.data.frame(pos_matrix) | is.matrix(pos_matrix)){
   chrs <- pos_matrix[,1]
   start.position <- pos_matrix[,2]
   end.position <- pos_matrix[,3]
 }

 index_chr_no_str <- vector(mode="logical", length=0)
 index_str = 1
 for(index in 1:length(chrs)) {
  if(bed_file == "" && is.data.frame(pos_matrix) == FALSE && is.matrix(pos_matrix) == FALSE){
     assign(paste("results.", chrs[index], sep = ""),
            STR_detection(chrs = chrs[index], nr.mismatch = nr.mismatch, nr.STRs = nr.STRs, STR = STR,
                                 start.position = start.position, end.position = end.position,
                                 translated_regions=F, output_file = "")) # start.position = 1, end.position = chr.lengths[index],
   }
   else{
     assign(paste("results.", chrs[index], sep = ""),
            STR_detection(chrs = chrs[index], nr.mismatch = nr.mismatch, nr.STRs = nr.STRs, STR = STR,
                                 start.position = start.position[index], end.position = end.position[index],
                                 translated_regions=F, output_file="")) # start.position = 1, end.position = chr.lengths[index],
   }

  temp2 <- which(unlist(get(paste("results.", chrs[index], sep = ""))[[6]])>= nr.STRs)
  if(length(temp2)==0) {
    index_chr_no_str <- append(index_chr_no_str,FALSE) # add index of chr where no STR is present
    next
  }
  else{
    index_chr_no_str <- append(index_chr_no_str,TRUE)
  }
  ind[[index_str]] = temp2
  length.As[[index_str]] = unname(unlist(get(paste("results.", chrs[index], sep = ""))[[6]][temp2]))
  pos.As[[index_str]] = unname(unlist(get(paste("results.", chrs[index], sep = ""))[[7]][temp2]))
  nr_matches[[index_str]] = length(length.As[[index_str]])
  index_str = index_str + 1
 }

if(length(which(index_chr_no_str == FALSE)) == length(chrs)){
  stop("Analysis aborted since no STRs were found in the regions!")
}
 else{ # only some chrs contain contain STRs or all
   length_chr <- which(index_chr_no_str == T)
 }
### We manually filtered the DSB Map to contain only those where AA-type is present ####
 # dsb_map = read.csv("DSB-map_Pratto_man.csv", header = T, sep = ";")
 pos.chr = sapply(length_chr, function(i) { # choose only those chromosomes where STRs are present
 #pos.chr = sapply(1:length(chrs), function(i) {
   which(STRAH::dsb_map$chrom == as.character(chrs[i]))
   })

 if(length(length_chr) == 1) {
   pos.chr <- list(c(pos.chr))
 }
 if (is.data.frame(pos.chr) == TRUE | is.matrix(pos.chr) == TRUE){
   pos.chr <- lapply(apply(pos.chr, 2, list), unlist)
 }
 start.dsb.int    = lapply(1:length(pos.chr), function(j) {STRAH::dsb_map[pos.chr[[j]],2]})
 end.dsb.int      = lapply(1:length(pos.chr), function(j) {STRAH::dsb_map[pos.chr[[j]],3]})
#### For being within it has to be larger than the start minus 500bp and smaller than the end position plus 500bp. ####
#### The greyzone is restructured in 5 segments of length 2 kb increasing in the distance to the hotspot ####
#### Moreover it is not allowed to contain a hotspot ####
 within = lapply(1:length(pos.chr), function(j) {
   sapply(1:length(ind[[j]]), function(i) {
     return(any(pos.As[[j]][i] >= (start.dsb.int[[j]]-500) & pos.As[[j]][i] <= (end.dsb.int[[j]]+500)))
     })
   })
 Sys.sleep(0.2)
  for(len.grey in 2:length(lens.grey)) {
    assign(paste("greyzone",lens.grey[len.grey-1]/1000,lens.grey[len.grey]/1000,sep="_"), lapply(1:length(pos.chr), function(j) {
      !within[[j]] & sapply(1:length(ind[[j]]), function(i) {
        return(any((pos.As[[j]][i] < (start.dsb.int[[j]]-500-lens.grey[len.grey-1]) & (pos.As[[j]][i] >= (start.dsb.int[[j]]-500-lens.grey[len.grey]))) | (pos.As[[j]][i] <= (end.dsb.int[[j]]+500+lens.grey[len.grey-1]) & pos.As[[j]][i] > (end.dsb.int[[j]]+500+lens.grey[len.grey]))))
        })
      }))
    message(paste("Number of greyzones finished: ", len.grey-1, " (of ", length(lens.grey)-1, ")", sep = ""))


 }## end lens.grey
 message("")
 flush.console()

 name.greyzones = unlist(strsplit(paste("greyzone", lens.grey[-length(lens.grey)]/1000, lens.grey[-1]/1000, collapse = " ", sep = "_"), split = " "))
 code.zone = lapply(1:length(get(name.greyzones[1])), function(x){
   return(sapply(1:length(get(name.greyzones[1])[[x]]), function(y) {
    if(within[[x]][y]) {
      return(1)
      } # "within"
    else for(k in 1:length(name.greyzones)) {
      if(get(name.greyzones[k])[[x]][[y]]) {return(k+1)}
    }
    return(length(name.greyzones)+2)}))
   })

 counter = 1
 for (s in length_chr) {
   if(bed_file == "" && is.data.frame(pos_matrix) == FALSE && is.matrix(pos_matrix) == FALSE){
     name_file <- paste0(chrs[s], ":", start.position, "-", end.position)
  }
  else{
   name_file <- paste0(chrs[s], ":", start.position[s], "-", end.position[s])
  }
  seq_name_list <- append(seq_name_list, name_file)
  header <- c("chr", "start_STR", "end_STR", "len_STR", "zones","chr_start_stop")


  if(output_file != ""){
    if (is.data.frame(df) == FALSE){
      write.table(rbind(header),output_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
    }
    df <- data.frame(chr_str = rep(chrs[s], length(pos.As[[counter]])),start_str = pos.As[[counter]], end_str = pos.As[[counter]] + length.As[[counter]] -1, len_str = length.As[[counter]], zones = code.zone[[counter]], chr_name = rep(name_file, length(pos.As[[counter]])))

    write.table(df, output_file, col.names=FALSE, sep="\t", row.names = FALSE, append=TRUE, quote=FALSE)
   }
  counter = counter + 1
 }

  output <- list("Sequence Name" = seq_name_list,
                  "Reverse Complement" = reverse.comp, "Number of allowed Mismatches" = nr.mismatch, "Minimum Length" = nr.STRs,
                  "Number of Matches" = nr_matches, "Length of STR stretch in bp" = length.As, "Start positions" = pos.As, "Zone" = code.zone)

  output <- lapply(output,FUN=function(x) {
   if(length(x) == length(unlist(seq_name_list)) & length(x) != 1){
     names(x) <- unlist(seq_name_list)
   }
 return(x)
 })

 return(output)
 closeAllConnections()

  }

