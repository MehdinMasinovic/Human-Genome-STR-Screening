#Author: Mehdin Masinovic

####Info####
"
 This file contains the main plots that were used in my thesis. 
"

####plotfigb1b2####
plotfigB1panel <- function(dataframe, letter){
  my_colors <- c("#DC143C", "#4E84C4", "#8FBC8F", "#9932CC", "#FF8C00", "#FFFF00", "#008080")
  my_zones <- c("HS", as.character(1:5), "OS")


  plot_panelB1 <- ggplot(data = dataframe) + labs(x = "Zone", y = paste0("Poly-", letter, " density"), fill = "Zone") +
    theme(title = element_text(color = "black"),
          axis.title = element_text(color = "black", size = 14),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          legend.key.size = unit(0.7, "cm")) +
    scale_fill_manual(values = my_colors[1:length(unique(dataframe$Zone))], labels=my_zones[1:length(unique(dataframe$Zone))]) 
  
  p1 = plot_panelB1 + geom_boxplot(aes(x=as.factor(dataframe$Zone), y=dataframe$panelB1dens, fill = as.factor(dataframe$Zone)),
                                   outlier.colour = "black", outlier.shape = 16,
                                   outlier.size = 2, notch = F) + scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + scale_x_discrete(labels= my_zones) 
  
  return(p1)
}

plotfigB2panel <- function(dataframe, letter){
  my_colors <- c("#DC143C", "#4E84C4", "#8FBC8F", "#9932CC", "#FF8C00", "#FFFF00", "#008080")
  my_zones <- c("HS", as.character(1:5), "OS")
  
  plot_panelB2 <- ggplot(data = dataframe) + labs(x = "Zone", y = paste0("Density of ",letter,"'s"), fill = "Zone") +
    theme(title = element_text(color = "black"),
          axis.title = element_text(color = "black", size = 14),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          legend.key.size = unit(0.7, "cm")) +
    scale_fill_manual(values = my_colors[1:length(unique(dataframe$Zone))], labels=my_zones[1:length(unique(dataframe$Zone))])
  
  
  p2 = plot_panelB2 + geom_boxplot(aes(x=as.factor(dataframe$Zone), y=dataframe$panelB2dens, fill = as.factor(dataframe$Zone)),
                                   outlier.colour = "black", outlier.shape = 16,
                                   outlier.size = 2, notch = F) + scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + scale_x_discrete(labels= my_zones) 
  
  return(p2)
  
}

####plotfigA####
plotfigA <- function(dataframe, letter) {
  my_colors <- c("#DC143C", "#4E84C4", "#8FBC8F", "#9932CC", "#FF8C00", "#FFFF00", "#008080")
  my_zones <- c("HS", as.character(1:5), "OS")
  
  plot_panelA <- ggplot(data = dataframe) + labs(y = paste0("Normalized poly-",letter," density per zone"), fill="Zone") +
    theme(title = element_text(color = "black"),
          axis.title = element_text(color = "black", size = 14),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),
          legend.key.size = unit(0.7, "cm")) +
    scale_fill_manual(values = my_colors[1:length(unique(dataframe$Zone))], labels=my_zones[1:length(unique(dataframe$Zone))]) + 
    scale_x_continuous(name = paste0("Poly-",letter," length"), breaks = min(unique(dataframe$Length_Polys/nchar(letter))):max(unique(dataframe$Length_Polys/nchar(letter))), 
                       labels = c(paste0(unique(dataframe$Length_Polys/nchar(letter)))[1:(length(unique(dataframe$Length_Polys/nchar(letter)))-1)], paste0("\u2265",unique(dataframe$Length_Polys/nchar(letter))[length(unique(dataframe$Length_Polys/nchar(letter)))]))
    )
  plot_labels <- paste0("n=",aggregate(n ~ Length_Polys, data = dataframe, FUN = sum)$n)
  pA <- plot_panelA + geom_bar(aes(x=dataframe$Length_Polys/nchar(letter), y=dataframe$renormdens, fill = as.factor(dataframe$Zone)), 
                               stat = "identity", position = "stack")  + 
    annotate(geom = "text", x = unique(dataframe$Length_Polys/nchar(letter))[1]:unique(dataframe$Length_Polys/nchar(letter))[length(unique(dataframe$Length_Polys/nchar(letter)))], y = 1-nchar(plot_labels)/100, label = plot_labels, 
             color = "black", angle = 90, family = "Arial")
  
  return(pA)
}
####plotfigC####
plotfigC <- function(dataframe, letter){
  my_colors <- c("#DC143C", "#4E84C4", "#8FBC8F", "#9932CC", "#FF8C00", "#FFFF00", "#008080")
  my_zones <- c("HS", as.character(1:5), "OS")
  
  plot_panelC <- ggplot(data = dataframe) + labs(y = paste0("Poly-",letter," density per base pair"), fill = "Zone") +
    theme(title = element_text(color = "black"),
          axis.text = element_text(size = 12),
          axis.title = element_text(color = "black", size = 14),
          legend.key.size = unit(0.7, "cm",),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),) +
    scale_fill_manual(values = my_colors[1:length(unique(dataframe$Zone))], labels=my_zones[1:length(unique(dataframe$Zone))]) + 
    scale_x_continuous(name = paste0("Poly-",letter," length"), breaks = min(unique(dataframe$Length_Polys/nchar(letter))):max(unique(dataframe$Length_Polys/nchar(letter))), 
                       labels = c(paste0(unique(dataframe$Length_Polys/nchar(letter)))[1:(length(unique(dataframe$Length_Polys/nchar(letter)))-1)], paste0("\u2265",unique(dataframe$Length_Polys/nchar(letter))[length(unique(dataframe$Length_Polys/nchar(letter)))]))
    )
  
  
  pC <- plot_panelC + geom_bar(aes(x=dataframe$Length_Polys/nchar(letter), y=dataframe$densities, fill = as.factor(dataframe$Zone)), 
                               stat = "identity", position = "dodge") + scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) 
  return(pC)
}

####plotfigD####
plotfigD <- function(dataframe, letter){
  my_colors <- c("#DC143C", "#4E84C4", "#8FBC8F", "#9932CC", "#FF8C00", "#FFFF00", "#008080")
  my_zones <- c("HS", as.character(1:5), "OS")
  
  plot_panelD <- ggplot(data = dataframe) + labs(y = paste0("Proportions of poly-", letter," density"), fill = "Zone") +
    theme(title = element_text(color = "black"),
          axis.text = element_text(size = 12),
          axis.title = element_text(color = "black", size = 14),
          legend.key.size = unit(0.7, "cm",),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 17),) +
    scale_fill_manual(values = my_colors[1:length(unique(dataframe$Zone))], labels=my_zones[1:length(unique(dataframe$Zone))]) + 
    scale_x_continuous(name = paste0("Poly-",letter," length"), breaks = min(unique(dataframe$Length_Polys/nchar(letter))):max(unique(dataframe$Length_Polys/nchar(letter))), 
                       labels = c(paste0(unique(dataframe$Length_Polys/nchar(letter)))[1:(length(unique(dataframe$Length_Polys/nchar(letter)))-1)], paste0("\u2265",unique(dataframe$Length_Polys/nchar(letter))[length(unique(dataframe$Length_Polys/nchar(letter)))]))
    ) + scale_y_continuous(breaks = seq(0, 1.25, by=0.25)) +
    coord_cartesian(ylim=c(0, 1.25))
  
  
  pD <- plot_panelD + geom_bar(aes(x=dataframe$Length_Polys/nchar(letter), y=dataframe$panelDdens, fill = as.factor(dataframe$Zone)), 
                               stat = "identity", position = "dodge")
  return(pD)
}