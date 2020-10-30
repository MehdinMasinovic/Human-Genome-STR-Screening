#Author: Mehdin Masinovic

####Info####
"
 In the set of functions to simulate a DNA sequence, three functions are working together: generate_sequence(), create_seq(), create_field().
 generate_sequence() is the main function that calls the other two. In a nutshell, generate_sequence() prepares the sequence arguments, 
 it passes those to create_seq(), which creates the actual DNA sequence, and create_field() puts the sequence fragments into order. 
"
####Function#####

generate_sequence <- function(lengthofseq = 150000, nr_HS = 10, le_HS = 2000,
                              nr_zones = 5, le_zones = 2000, dens_As = c(OZ = 0.5, rep(0.2,5) , HS = 0.8), 
                              STR= "A", setofAs= c(6,8,10), set_distr = c(0.5, 0.25, 0.25))
{
 

  if(sum(set_distr) != 1){
    stop("The distribution of your set \"set_distr\" of repeats is not summing up to 1.")
  }
  nucl <- c("A", "C", "T", "G")
  rest_n <- nucl[-which(nucl == STR)]
  
  names(dens_As)[c(-1, -length(dens_As))] <- paste("Z", (length(dens_As[c(-1, -length(dens_As))])+1):2, sep="") #we name our zones
  
  ###########################################
  
  Z1 = nr_HS * le_HS
  for(i in 2:(nr_zones+1)){
    assign(paste("Z", i, sep = ""), le_zones*nr_HS)
    
  }
  Z7 = lengthofseq - Z1 - nr_zones*le_zones*nr_HS
  
  count = 1
  sequencelist <- list()
  for (i in names(dens_As)){
    sequencelist[[count]] <- create_seq(get(i), dens_As[count], setofAs, set_distr, STR, rest_n)
    count = count + 1
  }
  
  ###########################################
  
  field_len <- (le_HS + nr_zones * le_zones) #length of one field. 
  
  oz <- floor((lengthofseq - (field_len * nr_HS))/ (nr_HS + 1))   #length of outside zone (oz) around fields. 
  
  all_len <- floor(field_len + oz) #length of field + oz
  
  
  if( lengthofseq < (all_len * nr_HS + oz) ){
    while(lengthofseq < (all_len * nr_HS)){
      nr_HS <- nr_HS - 1
    }
    warning(paste("It is not possible to create the number of Hotspot-Zones (HSZs) with the sequence and zone length provided.", "\n", "We reduce the number of HSZs to", nr_HS, "."))
  }
  
  pos <- function(vec, startpos, endpos){
    vec[which(vec == 0)[c(1,2)]] <- c(startpos,endpos)
    return(vec)
  }
  
  final_seq <- ""
  vec <- integer((nr_HS*(nr_zones+1)+ nr_HS + 6)*4)
  
  startpos = start_hs = start_zo = start_oz = 1
  end_oz = oz
  end_zo = le_zones/2
  end_hs = le_HS
  
  
  for(i in 1:nr_HS){
    final_seq = paste(final_seq, substr(sequencelist[[1]], start_oz, end_oz), collapse = "",sep = "")
    vec = pos(vec, startpos,endpos <- nchar(final_seq))
    startpos = endpos + 1
    for(j in c(2:(nr_zones + 2), (nr_zones + 1):2)){
      if(j == 7){
        segment = substr(sequencelist[[j]], start_hs, end_hs)
        substr(segment, 500, 500) <- "Z"
        substr(segment, nchar(segment)-500, nchar(segment)-500) <- "Z"
        
        
        final_seq = paste(final_seq, segment, collapse = "",sep = "")
        start_zo = start_zo + le_zones/2
        end_zo = end_zo + le_zones/2
        vec = pos(vec, startpos,endpos <- nchar(final_seq))
        startpos = endpos + 1
        }else{
        final_seq = paste(final_seq, substr(sequencelist[[j]], start_zo, end_zo), collapse = "",sep = "")
        vec = pos(vec, startpos,endpos <- nchar(final_seq))
        startpos = endpos + 1
      }
    }
      
    start_oz <- start_oz + oz
    end_oz <- end_oz + oz
    
    start_hs <- start_hs + le_HS
    end_hs <- end_hs + le_HS
    
    start_zo <- start_zo + le_zones/2
    end_zo <- end_zo + le_zones/2
    
    
    if(i == nr_HS){
      final_seq = paste(final_seq, substr(sequencelist[[1]], start_oz, end_oz), collapse = "",sep = "")
      vec = pos(vec, startpos,endpos <- nchar(final_seq))
      startpos = endpos + 1
      
    }
  }
  vec = vec[which(vec != 0)]
  if(nchar(final_seq)< lengthofseq){
    final_seq = paste(final_seq, paste(sample(c(rest_n), (lengthofseq - nchar(final_seq)), replace = T), collapse = "", sep = ""), collapse = "", sep = "")
    vec[length(vec)]<-nchar(final_seq)
    }
  
  
  return(list(final_seq, vec))
  
  ########################################################
}














