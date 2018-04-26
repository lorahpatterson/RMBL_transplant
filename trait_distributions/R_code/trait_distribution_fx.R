#Fx to create trait distribution
#Brian Maitner 4/17/2018

#Inputs:
  #Number of replicated outputs
  # Species abundance dataframe (2 columns, first column: species name, second column: abundance)
  # Trait data frame (2 or more columns: first column: species name, columns 2+ : traits)

#Output
  #Matrix with nrows = number of replicates, ncols = total abundance

trait_distributions<-function(number_replicates, abundance_data, trait_data){
  
  output<-list()
  
  for(t in 2:ncol(trait_data)){
    trait_t<-colnames(trait_data)[t]  
    
    out_t<-NULL
    for(n in 1:number_replicates){
      rep_n<-NULL  
      
      for( i in 1:nrow( abundance_data)){
        
        species_i<-as.character(abundance_data[i,1])    
        abund_i<-abundance_data[i,2]
        traits_i<-na.omit(trait_data[which(trait_data$Taxon==species_i),trait_t])
        
        #Dont do anything if there's no trait data
        if(length(traits_i)!=0){rep_n<-c(rep_n,sample(x = traits_i,size = abund_i,replace = T))}
        
      }# i abundance loop
      
      out_t<-rbind(out_t,rep_n)
      
    }#n replicates loop
    
    
    output[[t-1]]<-out_t
    
  }#t traits loop  
  names(output)<-colnames(trait_data)[2:ncol(trait_data)]
  
  return(output)  
  
}# trait_distribution function

#################################################
# Trait distributions by site for RMBL Transplant
# Lorah Patterson
# 4-25-18

library(tidyverse)
library(dplyr)

# load data
trait <- read.csv("Leaf_Traits_2017_040818.csv", header = TRUE)
abundance <- read.csv("Abundance_2017_JS_LP.csv", header = TRUE)

# Get species abundance dataframe in shape (2 columns, first column: species name, second column: abundance)
abundance_CBT <- abundance %>% filter(Site == "CBT") %>% group_by(Taxon) %>% summarise(sum(Count))
abundance_UM <- abundance %>% filter(Site == "Upper Montane") %>% group_by(Taxon) %>% summarise(sum(Count))
abundance_Pf <- abundance %>% filter(Site == "Pfeiler") %>% group_by(Taxon) %>% summarise(sum(Count))
abundance_Mon <- abundance %>% filter(Site == "Monument") %>% group_by(Taxon) %>% summarise(sum(Count))

# Get trait dataframe in shape (2 or more columns: first column: species name, columns 2+ : traits)
trait_CBT <- trait %>% filter(Site == "CBT") %>% select(Taxon, Leaf_Area_cm2, SLA, LDMC, Thickness_avg_or_single_mm)
trait_UM <- trait %>% filter(Site == "Upper Montane") %>% select(Taxon, Leaf_Area_cm2, SLA, LDMC, Thickness_avg_or_single_mm)
trait_Pf <- trait %>% filter(Site == "Pfeiler") %>% select(Taxon, Leaf_Area_cm2, SLA, LDMC, Thickness_avg_or_single_mm)
trait_Mon <- trait %>% filter(Site == "Monument") %>% select(Taxon, Leaf_Area_cm2, SLA, LDMC, Thickness_avg_or_single_mm)

# Run the trait distributions function
dist_CBT <- trait_distributions(10,abundance_CBT, trait_CBT)




