#DATA WRANGLE SRA MALINSKY 2017 DATA: BIOPROJECT PRJEB1254 and PRJEB15289
#MOORE LAB
#OCTOBER 2026

#Load required packages for data wrangling

library(dplyr)
library(tidyverse)

#Data downloaded from Bioproject numbers on NCBI from Malinsky et al. (2017) paper
#Load data into Rscript

bioproject1 <- read.csv("C:/Users/crazy/Desktop/Moore Lab/Code/cichild_compartive_genomics-main/SraRunTable_PRJEB1254.csv", header = TRUE)

bioproject2 <- read.csv("C:/Users/crazy/Desktop/Moore Lab/Code/cichild_compartive_genomics-main/SraRunTable_PRJEB15289.csv", header = TRUE)

#Join the two datasets together

bioproject <- bind_rows(bioproject1, bioproject2)

#Look at the number of unique species and accession number

unique_ids <- bioproject %>% 
  group_by(Organism, BioSample) %>% 
  tally()
#519 unique species and BioSamples

#Look at the total number of species

bioproject %>% 
  summarise(total_species = n_distinct(Organism))
#127 species

# Organize data by alphabetical order in Organism column

bioproject <- bioproject %>% 
  arrange(Organism)

#Keep only one representative from each species

bioproject_species_unique <- bioproject %>%
  group_by(Organism) %>%
  slice_head(n = 1) %>%
  ungroup()

#Generate file with just accession numbers

accession_number <- bioproject_species_unique$Run
writeLines(as.character(accession_number), "C:/Users/crazy/Desktop/Moore Lab/Code/cichild_compartive_genomics-main/sra_accession_number.txt")
