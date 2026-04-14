# Author: Mayra Beatriz Mendoza Velazquez
# Title: Exploration of reconstruction methods

library(readr)
library(mlBioNets)
library(dplyr)
library(tibble)


# Load csv archives 
metadata <- as.data.frame(read_tsv(file = "01_RawData/metadata_clean.tsv"))
metadata[is.na(metadata)] <- 0
f_clean <- as.data.frame(read_tsv("01_RawData/f_clean.tsv"))
rzcompositiondata <- read.csv("01_RawData/rzcomposition.csv")

#### Function ####
# Function for selecting specific columns based on the community sample 
community_isolation <- function(rcommunities, sample, sample_col, rcommunities_col, df, arrangev, interest_column, dfwvals, composition){

   library(dplyr)
# create an empty list   
community_list <- list()
k <- 1

# subset the commuinity values based on the day and community name 

for (i in 1:length(rcommunities)){
  for (x in 1:length(sample)){
    
    community_list[[k]] <- subset(df, sample_col %in% c(0, sample[x])  & rcommunities_col == rcommunities[i]) %>% 
      arrange(arrangev) %>% 
      pull(interest_column) %>% 
      as.character()
    
    k <- k + 1
  }
}

# create an empty list for the abundances 
abundances_tables <- list()

# Sbased on the column names, select the values from f_clean 
#
for (id in seq_along(community_list)) {
  abundances_tables[[id]] <- dfwvals %>%
    dplyr::select(1 | all_of(community_list[[id]])) %>%
    column_to_rownames(dfwvas[,1])
  
}

abnds_filtered <- list()

for (g in seq_len(ncol(composition))) {
  
  # select the community name / for arranging the list 
  comm_name <- colnames(composition)[g]
  
  # select the specific bacterial id for the selected community 
  spp <- composition[, g]
  spln <- 1:length(spp)
  # indexes [to select in the list of data.frames]
  idx <- (2*g - 1):(2*g)
  
  for (j in seq_along(idx)) {
    
    k <- idx[j]
    
    abnds_filtered[[paste0(comm_name, "_", sample[j])]] <- abundances_tables[[k]] %>% 
      filter(row.names(abundances_tables[[k]]) %in% spp)
    
    
  }
}

return(abnds_filtered)

}

## Function testing 
comms_rhiz <- unique(metadata$community) # to use all the community values for subsetting the df 
temps_rhiz <- c(28, 32) # without the 0 temperature, because it is already being used in the function 

# Using function 
rz_communities <- community_isolation(rcommunities = comms_rhiz , sample = temps_rhiz , sample_col = "temp", 
                                      rcommunities_col = "community", df = metadata, arrangev = "day", interest_column = "label", dfwvals = f_clean)


## FUNCTION DOESN'T WORK ## 

# isolated "for loops" to subset the data.frames #

rz_list_not_filtered <- list()
k <- 1
# subset the commuinity values based on the day and community name 

for (i in 1:length(comms_rhiz)){
  for (x in 1:length(temps_rhiz)){
    
    rz_list_not_filtered [[k]] <- subset(metadata, temp %in% c(0, temps_rhiz[x])  & community == comms_rhiz[i]) %>% 
      arrange(day) %>% 
      pull(label) %>% 
      as.character()
    
    k <- k + 1
  }
}
rz_list_not_filtered
# create an empty list for the abundances 
abundances_rz_not_filtered <- list()

# Based on the column names, select the values from f_clean 
for (id in seq_along(rz_list_not_filtered)) {
  abundances_rz_not_filtered[[id]] <- f_clean %>%
    dplyr::select(1 | all_of(rz_list_not_filtered[[id]]) )  %>%
    column_to_rownames(var = "row.names")
}

# filter the list of data.frames based on the community composition 
rz_filtered <- list()

for (g in seq_len(ncol(rzcompositiondata))) {
  
  # select the community name / for arranging the list 
  comm_name <- colnames(rzcompositiondata)[g]
  
  # select the specific bacterial id for the selected community 
  spp <- rzcompositiondata[, g]
  spln <- 1:length(spp)
  # indexes [to select in the list of data.frames]
  idx <- (2*g - 1):(2*g)
  
  for (j in seq_along(idx)) {
    
    k <- idx[j]
    
    rz_filtered[[paste0(comm_name, "_", temps_rhiz[j])]] <- abundances_rz_not_filtered[[k]] %>% 
      filter(row.names(abundances_rz_not_filtered[[k]]) %in% spp)

    
  }
}

rz_filtered 

##### Testing network inference algortihms ####

#### SPARCc ####

# R1
R1_28_inference <- net_inference(taxa_abs = t(rz_filtered[[1]]), method = "sparcc", p = 0.5)
plot(R1_28_inference)

R1_32_inference <- net_inference(taxa_abs = t(rz_filtered[[2]]), method = "sparcc", p = 0.5)
plot(R1_32_inference)

# R2
R2_28_inference <- net_inference(taxa_abs = t(rz_filtered[[3]]), method = "sparcc", p = 0.5)
plot(R2_28_inference)

R2_32_inference <- net_inference(taxa_abs = t(rz_filtered[[4]]), method = "sparcc", p = 0.5)
plot(R2_32_inference )

# R3
R3_28_inference <- net_inference(taxa_abs = t(rz_filtered[[5]]), method = "sparcc", p = 0.5)
plot(R3_28_inference)

R3_32_inference <- net_inference(taxa_abs = t(rz_filtered[[6]]), method = "sparcc", p = 0.5)
plot(R3_32_inference)

# R4
R4_28_inference <- net_inference(taxa_abs = t(rz_filtered[[7]]), method = "sparcc", p = 0.5)
plot(R4_28_inference)

R4_32_inference <- net_inference(taxa_abs = t(rz_filtered[[8]]), method = "sparcc", p = 0.5)
plot(R4_32_inference)

# R5
R5_28_inference <- net_inference(taxa_abs = t(rz_filtered[[9]]), method = "sparcc", p = 0.5)
plot(R5_28_inference)

R5_32_inference <- net_inference(taxa_abs = t(rz_filtered[[10]]), method = "sparcc", p = 0.5)
plot(R5_32_inference)

# R6
R6_28_inference <- net_inference(taxa_abs = t(rz_filtered[[11]]), method = "sparcc", p = 0.5)
plot(R6_28_inference)

R6_32_inference <- net_inference(taxa_abs = t(rz_filtered[[12]]), method = "sparcc", p = 0.5)
plot(R6_32_inference)

# R7
R7_28_inference <- net_inference(taxa_abs = t(rz_filtered[[13]]), method = "sparcc", p = 0.5)
plot(R7_28_inference )

R7_32_inference <- net_inference(taxa_abs = t(rz_filtered[[14]]), method = "sparcc", p = 0.5)
plot(R7_32_inference)

# R8
R8_28_inference <- net_inference(taxa_abs = t(rz_filtered[[15]]), method = "sparcc", p = 0.5)
plot(R8_28_inference)

R8_32_inference <- net_inference(taxa_abs = t(rz_filtered[[16]]), method = "sparcc", p = 0.5)
plot(R8_32_inference)

# R9
R9_28_inference <- net_inference(taxa_abs = t(rz_filtered[[17]]), method = "sparcc", p = 0.5)
plot(R9_28_inference)

R9_32_inference <- net_inference(taxa_abs = t(rz_filtered[[18]]), method = "sparcc", p = 0.5)
plot(R9_32_inference)

# R10
R10_28_inference <- net_inference(taxa_abs = t(rz_filtered[[19]]), method = "sparcc", p = 0.5)
plot(R10_28_inference)

R10_32_inference <- net_inference(taxa_abs = t(rz_filtered[[20]]), method = "sparcc", p = 0.5)
plot(R10_32_inference)

# R11
R11_28_inference <- net_inference(taxa_abs = t(rz_filtered[[21]]), method = "sparcc", p = 0.5)
plot(R11_28_inference)

R11_32_inference <- net_inference(taxa_abs = t(rz_filtered[[22]]), method = "sparcc", p = 0.5)
plot(R11_32_inference)

# R12
R12_28_inference <- net_inference(taxa_abs = t(rz_filtered[[23]]), method = "sparcc", p = 0.5)
plot(R12_28_inference)

R12_32_inference <- net_inference(taxa_abs = t(rz_filtered[[24]]), method = "sparcc", p = 0.5)
plot(R12_32_inference)
