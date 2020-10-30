#Author: Mehdin Masinovic

####Info####
"
 In the folder \"Hypothesis-Testing\" I explain how the Kruskal-Wallis test was used
 to check for differences among zones, including a post-hoc comparison using a Dunn's test. 
 
 In order to compare densities among zones for short tandem repeats, 
 it is important to distinguish between untransformed and transformed densities. 
 
 In this chapter, I want to explain why I choose to test transformed densities. 
"

####Libraries & working directory####

rm(list=ls())
setwd("~/Why-Transformed-Densities")
library(plyr)
library(dplyr)
library(FSA)
library(ggplot2)

source("create_plots.R")

####Loading functions####


####Explaining test####
"In this example, I want to illustrate why it is important to transform the data before it is tested. 
 I create a simple example with dummy values that will emulate repeat densities stratified by repeat length. "

Zone <- as.factor(rep(1:7,7))
densities <- c(24.555, 22.888, 22.887, 22.886, 22.885, 22.884, 22.883,
               20.666, 18.339, 18.338, 18.337, 18.336, 18.335, 18.334,
               16.333, 14.005, 14.004, 14.003, 14.002, 14.001, 14.000,
               11.999, 9.669, 9.668, 9.667, 9.666, 9.665, 9.664,
               7.666, 5.339, 5.338, 5.337, 5.336, 5.335, 5.334,
               3.333, 3.005, 3.004, 3.003, 3.002, 3.001, 3.000,
               0.999, 0.666, 0.665, 0.664, 0.663, 0.662, 0.661)

Length_Polys <- rep(6:12, each = 7)
simdf <- data.frame(Zone, densities, Length_Polys)

simp <- plotfigC(simdf, "A")
simp

ggsave(
  filename=paste0("Simulated-poly-A-densities.png"),
  plot = simp,
  device = "png",
  path = NULL,
  scale = 1,
  width = 20,
  height = 16,
  units = c("cm"),
  limitsize = TRUE,
)

"In the plot seen when running plotfigC(), we can imagine that the repeat densities of repeat length 6 are obviously higher
 than the repeat densities of length 7. The repeat densities of length 7 are higher than 8, 8 higher than 9, and so on. 
 We expect this behaviour because we assume that it is less probably for a repeat to be consecutively repeated 15 times than 14, etc.
 Expecting a behaviour is not proof. Therefore, I added a folder in the GitHub repository that lists the repeat densities stratified
 by repeat length for all 70 short tandem repeats with a repeat unit smaller 4 (mono-, di-, and trinucleotide repeats).
 In all cases, the following holds: The longer the repeat, the less frequent it is. 
 

 In my dummy example, I created a set of fake densities which should emulate a normal scenario. 
 Here, the hotspot zone is different from all other zones in really every repeat length, and
 we don't see any difference in the surrounding five zones + the outside zone. 
 
 A kruskal-wallis test, therefore, is expected to tell you that there are significant differences, not knowing where they actually are. 
 Now, we run a kruskal-wallis test using the fake densities:
"


x <- simdf$densities
g <- as.factor(simdf$Zone)
k <- nlevels(g)
n <- length(x)
r <- rank(x)

mydf <- data.frame(DENS = x, ZONE = g, RANK = r)


mydf[order(-mydf$RANK),]

"
 How does a Kruskal-Wallis test work?
 The Kruskal-Wallis test is a non-parametric test that tests whether there is a difference in distribution in two samples.
 The way it does it is to \"sort\" your data and give it ranks, i.e. it ranks your data from lowest to highest.
 Then, it computes the statistic of the test by summing the ranks of one group, taking the squared value of it, and then divide it
 by the number of observations of the group. 
 This statistic is then used in the formula suggested in the paper 
 ((12 * STATISTIC / (n * (n + 1)) - 3 * (n + 1)) / (1 - sum(TIES^3 - TIES) / (n^3 - n)))
 
 and ultimately determines how strong the difference will be in the groups.  
 
 Above this text, I printed the ranked values in a dataframe. Intruigingly, the values are ranked but the zones always remain in the same order:
 1,2,3,4,5,6,7,1,2,3,4,5,6,7, and so on. This is the reason why the kruskal-wallis test is bound to be insignificant:
 The repeat density of length 6 is at a different scale than length 7, length 7 on a different scale than length 8, and so on. 
 As seen below:
"

TIES <- table(x)
STATISTIC <- sum(tapply(r, g, "sum")^2 / tapply(r, g, "length"))
## keep as n+1 to avoid (implausible) integer overflows
STATISTIC <- ((12 * STATISTIC / (n * (n + 1)) - 3 * (n + 1)) /
                (1 - sum(TIES^3 - TIES) / (n^3 - n)))
PARAMETER <- k - 1L
PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
names(STATISTIC) <- "Kruskal-Wallis chi-squared"
names(PARAMETER) <- "df"
PVAL #insignificant


"
 Instead, I decided to use the normalized repeat densities, where every density of a repeat length is divided by the 
 sum of all densities in a repeat length. Here, I would actually see how different the zones are among one repeat length proportionally speaking. 
 The code is shown below.
"


simdf$renormdens <- rep(0, 49)
for(i in 6:12){
  simdf[simdf$Length_Polys == i,]$renormdens <- simdf[simdf$Length_Polys == i,]$densities / sum(simdf[simdf$Length_Polys == i,]$densities)
}
simdf$densities <- simdf$renormdens

simp <- plotfigC(simdf, "A")
simp

ggsave(
  filename=paste0("Simulated-poly-A-densities-normalized.png"),
  plot = simp,
  device = "png",
  path = NULL,
  scale = 1,
  width = 20,
  height = 16,
  units = c("cm"),
  limitsize = TRUE,
)

"
 Now, we see that the densities across all repeat lengths are brought on to the same scale, and what we are testing are whether the
 hotspot zone is proportionally larger than the remaining zones. 
 If we consider the ranking of the Kruskal Wallis, we will now see that the hotspot zone will receive the highest rank among all zones
 (as we would expect). See below.
"

x <- simdf$renormdens
g <- as.factor(simdf$Zone)
k <- nlevels(g)
n <- length(x)
r <- rank(x)

mydf <- data.frame(DENS = x, ZONE = g, RANK = r)
mydf[order(-mydf$RANK),] #Hotspot gets highest ranks, therefore there is a strong difference in the statistic value. 

TIES <- table(x)
STATISTIC <- sum(tapply(r, g, "sum")^2 / tapply(r, g, "length"))
## keep as n+1 to avoid (implausible) integer overflows
STATISTIC <- ((12 * STATISTIC / (n * (n + 1)) - 3 * (n + 1)) /
                (1 - sum(TIES^3 - TIES) / (n^3 - n)))
PARAMETER <- k - 1L
PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
names(STATISTIC) <- "Kruskal-Wallis chi-squared"
names(PARAMETER) <- "df"
PVAL # significant difference

"The test now tells us that there are significant differences, which we expect as long as we made an extreme example of differences among zones.
 Now, we can see where those differences exist using a Dunns test:"

FSA::dunnTest(renormdens ~ as.factor(Zone), data = simdf) 

"Here, the hotspot zone is the only zone that results in significant differences among zones. 
"

####Cut off values or not?####

"When considering differences among zones in normalized densities, there is one more topic that has to be considered: 
 Whether repeats with high repeat lengths should be considered in the analysis. 
 As we are now considering proportional differences, it is important to question whether repeats of high repeat length 
 contribute equally to the repeat type than repeats of low repeat length. 
 
 We can take one more example to illustrate this topic using actual repeat densities of the poly-G repeat in the human genome:"

polyG <- readRDS("/home/roots/Desktop/projects/Project_STRAH/STRAH/my_codes/strah_writing_thesis/github_STRAH/Why-Transformed-Densities/hg38_2k_G.Rds")
colnames(polyG)[2] <- "Length_Polys"

polyGp <- plotfigC(polyG, "G")

ggsave(
  filename=paste0("Poly-G-densities.png"),
  plot = polyGp,
  device = "png",
  path = NULL,
  scale = 1,
  width = 20,
  height = 16,
  units = c("cm"),
  limitsize = TRUE,
)

"Here, we see a similar behaviour like in our previous example: 
 1. The higher the repeat length, the less frequent it becomes
 2. Hotspot zone is always higher than the remaining zones.

 In this example, a Kruskal-Wallis test tells you there are no differences with a p-value extremely close to 1.

 Now, we use the transformed values:"


polyG$densities <- polyG$renormdens
polyGp <- plotfigC(polyG, "G")


ggsave(
  filename=paste0("Poly-G-densities-normalized.png"),
  plot = polyGp,
  device = "png",
  path = NULL,
  scale = 1,
  width = 20,
  height = 16,
  units = c("cm"),
  limitsize = TRUE,
)

"Now, we see how zones differ proportionally speaking in every repeat-length. 
 The problem now is that the Kruskal-Wallis test will not use the information of the repeat length in the testing,
 it will not remember that the poly-G repeats of length 6 are six times more frequent than length 7, 
 29 times more frequent than length 8, 81 times more frequent than length 9, and so on. 
 
 Therefore, one could argue that the densities with lower repeat length should have a higher weight in the testing,
 or maybe that repeats with a very high repeat length should not be considered at all or tested seperately. 
 
 When cutting off data or giving it less weight, nevertheless, you are assuming that the differences in higher repeat lengths
 are \" less important \" or that they don't represent the poly-G repeats as much as repeats with lower repeat length. 
 
 For my thesis, I considered both cases. In the \"Hypothesis-Testing\" chapter, i provide an option \"cutdata=T/F\", 
 where the user can choose to cut of the lowest 1% of the dataset. In my testing, not cutting or cutting off data resulted 
 in the same differences, meaning that cutting off data has not changed the interpretation of the result. 
 
 Therefore, I chose to present the transformed, uncut testing of the data, because I have the \"possibility\" to interpret
 the data the same way without making assumptions that might not be true in the case of short tandem repeats. 
"














