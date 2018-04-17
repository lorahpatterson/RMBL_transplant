# April 2018
# Principal Components Analysis - 2017 RMBL turf transplant data
# Lorah Patterson

library(dplyr)
library(tidyr)
library(plyr)

setwd("~/Documents/R_data/PCA_turfs_April2018/")
turf_traits <- read.csv("Leaf_Traits_2017_040818.csv", header = TRUE)
turf_abund <- read.csv("Abundance_2017_JS_LP.csv", header = TRUE)

### GOAL 1: Calculate weighted trait means per block.

  # 1) Calculate the mean of each trait per species/block/site. Note: In 2017, we collected trait data by block, not by plot.

mean_la <- turf_traits %>% group_by(Site,Block,Taxon) %>% summarise(mean(Leaf_Area_cm2))  # Mean of Leaf Area per species/block/site. 
mean_sla <- turf_traits %>% group_by(Site,Block,Taxon) %>% summarise(mean(SLA))  # Mean of SLA per species/block/site.
mean_ldmc <-turf_traits %>% group_by(Site,Block,Taxon) %>% summarise(mean(LDMC)) # Mean of LDMC per species/block/site.
mean_thick <-turf_traits %>% group_by(Site,Block,Taxon) %>% summarise(mean(Thickness_avg_or_single_mm)) # Mean of thickness per species/block/site.

mean_traits <- left_join(mean_la, mean_sla, by=c("Site", "Block", "Taxon")) 
mean_traits2 <- left_join(mean_traits, mean_ldmc, by=c("Site", "Block", "Taxon"))
mean_traits3 <- left_join(mean_traits2, mean_thick, by=c("Site", "Block", "Taxon")) #Join all of the mean trait data.

names(mean_traits3)[names(mean_traits3)=='mean(Leaf_Area_cm2)'] <- 'mean_leafarea'
names(mean_traits3)[names(mean_traits3)=='mean(SLA)'] <- 'mean_SLA'
names(mean_traits3)[names(mean_traits3)=='mean(LDMC)'] <- 'mean_LDMC'
names(mean_traits3)[names(mean_traits3)=='mean(Thickness_avg_or_single_mm)'] <- 'mean_thick'  # Change the column names


  # 2) Merge the abundance data to the mean trait data.  Note: We do not have abundance data for every species in the trait data. NAs will be inserted when this is the case. We measured abundance for one plot per block, so this one plot will have to substitute for a whole block.

turf_abund$Block <- as.factor(turf_abund$Block) # converts Block from integer to factor so it will join
mean_traits3$Block <- as.factor(mean_traits3$Block) # converts Block from integer to factor so it will join

traits_abundance <- left_join(mean_traits3, turf_abund, by = c("Site","Block","Taxon")) #join mean trait data to abundance data


  # 3) Create weights for each site/block/species combo: Divide count per species (i.e. abundance) by total number of plants counted in that block. For example, if we counted 10 of species A out of 100 total plants in Block 1, the weight would be 0.10 for species A in Block 1. Weights per block should add up to 1.

totalblock <- turf_abund %>% group_by(Site, Block) %>% summarise(sum(Count))# Calculate the total number of plants per "block" (In reality, we only counted abundance for one plot in each block of five plots.)
names(totalblock)[names(totalblock)=='sum(Count)'] <- 'sum_count' # Rename column

tot_traits_abundance <- left_join(traits_abundance, totalblock, by= c("Site", "Block"))# Join total count per block to abundance data

tot_traits_abundance <- mutate(tot_traits_abundance, weight = (Count / sum_count)) # Calculates weight (count of species per plot divided by count of total species per block)


  # 4) Calulate the weighted trait mean per species per block. Multiply the mean of each trait per species by the weight. Then sum these numbers per block and divide by the sum of the weights per block (which is 1, so I ignored this step). 

tot_traits_abundance <- mutate(tot_traits_abundance, weighted_leafarea = (mean_leafarea * weight)) # Multiply the mean of each trait by the weight
tot_traits_abundance <- mutate(tot_traits_abundance, weighted_SLA = (mean_SLA * weight))
tot_traits_abundance <- mutate(tot_traits_abundance, weighted_LDMC = (mean_LDMC * weight))
tot_traits_abundance <- mutate(tot_traits_abundance, weighted_thick = (mean_thick * weight))

# Couldn't figure out how to do this next part in R, so did it in excel. :(
#weighted_trait_means <- tot_traits_abundance %>% group_by(Site, Block) %>% summarize(tot_traits_abundance,leafarea= sum(weighted_leafarea)) #Sum the weighted traits for each block.
# Add a column concatenating Site and Block number and make this column 1. Delete Site and Block individual columns.
# write.csv(tot_traits_abundance, file = "tot_traits_abundance.csv")

tot_weight_traits4 <- read.csv("tot_weight_traits4.csv", header = TRUE, row.names=1)  # Import data from excel. Make sure row names are column 1.
traits_matrix <- data.matrix(tot_weight_traits4, rownames.force = TRUE) # Convert data frame into matrix. Samples should be rows and traits should be columns for PCA.



### GOAL 2: Principal Components Analysis (PCA) of traits across sites (Source: Joshua Starmer StatQuest https://www.youtube.com/watch?v=0Jp4gsfOLMs)

  # 1) Run the prcomp function on the traits_matrix
pca_traits <- prcomp((traits_matrix), scale=TRUE) 
pca_traits # PC1 accounts for most variation in the original data, PC2 accounts for second most, etc.

  # 2) simple plots using base graphics
plot(pca_traits$x[,1], pca_traits$x[,2]) # 2D PCA graph. PC1 is on the X axis and is the 1st column in pca_traits; PC2 is on the Y axis and is the 2nd column in pca_traits.
pca_traits_var <- pca_traits$sdev^2 # Calculates how much variation in the original data each principal component accounts for.
pca_traits_var_per <- round(pca_traits_var/sum(pca_traits_var)*100,1) # Calculate percentages
barplot(pca_traits_var_per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

  # 3) fancier plots using ggplot2 
pca_data <- data.frame(Sample=rownames(pca_traits$x), X=pca_traits$x[,1], Y=pca_traits$x[,2])
pca_data

ggplot(data=pca_data, aes(x=X, y=Y, label=Sample)) + geom_text() + xlab(paste("PC1 -", pca_traits_var_per[1], "%", sep="")) + ylab(paste("PC2 - ", pca_traits_var_per[2], "%", sep="")) + theme_bw() + ggtitle("Principal Components Analysis")

  # 4) "Loading scores" show us which traits on PC1 push samples to the left side of the graph (negative score) or right side (positive score)

loading_scores_traits <-pca_traits$rotation[,2] # Calc laoding scores for PC1 (ie column 1)
trait_scores <- abs(loading_scores_traits) # Sort based on the magnitude of the number, not from high to low
trait_scores_ranked <- sort(trait_scores, decreasing=TRUE) # Sort the magnitude from high to low
top_traits <- names(trait_scores_ranked[1:4]) # list of traits from high to low magnitude
pca_traits$rotation[top_traits,1] # This shows the scores and the positive or negative values. LDMC, thickness, and SLA push the samples right (positive score) and leaf area pushes the samples left (neg score) along PC1


### GOAL 3: create data matrix of cover of species per plot (maybe just do dominant species)