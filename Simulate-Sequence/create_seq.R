#Author: Mehdin Masinovic

####Info####
"
 In the set of functions to simulate a DNA sequence, three functions are working together: generate_sequence(), create_seq(), create_field().
 generate_sequence() is the main function that calls the other two. In a nutshell, generate_sequence() prepares the sequence arguments, 
 it passes those to create_seq(), which creates the actual DNA sequence, and create_field() puts the sequence fragments into order. 
"
####Function#####

create_seq <- function(lengthseq = NULL, dens_A = NULL, setofAs = NULL, set_distr = NULL, STR = "A", rest_n = c("C", "G", "T")){
  

  amount_STR <- floor(lengthseq * dens_A) #amount of our STR, g.e. "A"
  amount_rest <- lengthseq - amount_STR #amount of remaining nucleotides, g.e. "G", "C", and "T"
  
  sets <- amount_STR * set_distr
  
  amount_STR <- floor(sets/setofAs)
  
  if(sum(amount_STR) == 0){return(paste(sample(rest_n, lengthseq, replace = T), collapse =""))}
  if(any(amount_STR == 0)){setofAs <- setofAs[-which(amount_STR == 0)]; amount_STR <- amount_STR[-which(amount_STR == 0)]}
  pool <- c()
  for (i in 1:length(amount_STR)){
    pool <- c(pool, replicate(expr = paste(rep(STR, setofAs[i]), collapse = ""), n = amount_STR[i]))
  }
  
  chosen_pool <- sample(pool, length(pool), replace = F)
  
  
  
  if(length(chosen_pool) > amount_rest){
    warn <- TRUE
  } else { warn <- FALSE }
  
  while( length(chosen_pool) > amount_rest ) {
    
    for (i in 1: length(setofAs)){
      
      
      chosen_pool <- chosen_pool[-which(nchar(chosen_pool) == setofAs[i])[1]]
      
      amount_STR <- amount_STR - setofAs[i]
      amount_rest <- amount_rest + setofAs[i]
      
      sets <- amount_STR * set_distr
      
      amount_STR <- floor(sets/setofAs)
      
      
      if(sum(amount_STR) <= amount_rest){
        break
      }
      
    }
    
  }
  
  if(warn){
    warning(paste("Your chosen percentage for A's will not allow for at least one letter inbetween the desired repeats. The function will continue on the value which suits your requirements the best. The maximum density of A for a length of ", lengthseq, "is smaller equal", amount_STR/lengthseq))
  }
  
  
  
  inbetween <- (floor(amount_rest / length(chosen_pool)+1)) - 1
  
  
  final_seg <- c()
  for(i in 1:length(chosen_pool)){
    final_seg <- c(final_seg, paste(sample(rest_n, inbetween, replace = T), collapse = ""), chosen_pool[i])
  }
  final_seg <- paste(final_seg, collapse = "")
  end <- paste(sample(rest_n, floor(lengthseq - nchar(final_seg)), replace = T), collapse ="")
  final_seg <- paste(final_seg, end, collapse = "", sep = "")
  
  return(final_seg)
}
