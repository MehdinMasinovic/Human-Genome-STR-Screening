#Author: Mehdin Masinovic

####Info####
"
 In the set of functions to simulate a DNA sequence, three functions are working together: generate_sequence(), create_seq(), create_field().
 generate_sequence() is the main function that calls the other two. In a nutshell, generate_sequence() prepares the sequence arguments, 
 it passes those to create_seq(), which creates the actual DNA sequence, and create_field() puts the sequence fragments into order. 
"
####Function#####

create_field <- function(le_zones, le_HS, oz, dens_As, setofAs, set_distr, nr_zones, STR, rest_n, startpos, endpos, vec){
  
  final_seg <- ""
  
  for ( i in c(1:(nr_zones+1), nr_zones:1)){
    
    if(( i >=1 ) && ( i <= 5 )){
      segment <- create_seq(lengthseq = le_zones/2, dens_A = dens_As[i+1], setofAs = setofAs, set_distr = set_distr, STR = STR, rest_n = rest_n)
      final_seg <- c(final_seg, segment)
      startpos <- endpos + 1; endpos <- endpos + nchar(segment); vec[which(vec == 0)[c(1,2)]] <- c(startpos,endpos)
    }
    if(i == (nr_zones + 1)){
      segment1 <- create_seq(lengthseq = le_HS/2, dens_A = dens_As[i+1], setofAs = setofAs, set_distr = set_distr, STR = STR, rest_n = rest_n) # not generic
      segment2 <- create_seq(lengthseq = le_HS/2, dens_A = dens_As[i+1], setofAs = setofAs, set_distr = set_distr, STR = STR, rest_n = rest_n)
      segment <- paste(segment1, segment2, collapse = "", sep = "")
      substr(segment, 500, 500) <- "Z"
      substr(segment, nchar(segment)-500, nchar(segment)-500) <- "Z"
      final_seg <- c(final_seg, segment)
      startpos <- endpos + 1; endpos <- endpos + nchar(segment); vec[which(vec == 0)[c(1,2)]] <- c(startpos,endpos)
    }
  }
  
  segment <- create_seq(lengthseq = oz, dens_A = dens_As[1], setofAs = setofAs, set_distr = set_distr, STR = STR, rest_n = rest_n)
  final_seg <- c(final_seg, segment)
  final_seg <- paste(final_seg, collapse = "")
  startpos <- endpos + 1; endpos <- endpos + nchar(segment); vec[which(vec == 0)[c(1,2)]] <- c(startpos,endpos)
  
  return(list(final_seg, startpos, endpos, vec))
}
