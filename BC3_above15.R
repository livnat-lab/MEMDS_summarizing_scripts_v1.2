# This script reports frequency of read families with more than 15 different BC3_groups in them
# Based on the data in the "mutationFrequncyPerBarcode" files

# Define libraries and I/O
rm(list=ls(all=TRUE))
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

cat("Starting", date(), "\n")
x <- Sys.time()

args <- commandArgs(TRUE)
cat(c("Run args", args), sep = "\n")

infile_dir <- args[1]
outfile <- args[2]

#infile_dir <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/tables.BC3cutoff1.BC_dist"
#outfile <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/BC3_types_above15.txt"

############################################################
infile_vec <- grep("BC_dist.txt$", dir(infile_dir), value=T)
cat("Detected", length(infile_vec), "input files in", infile_dir, "\n")

# Initialize output data frame
high_BC3 <- data.frame(Sample = character(), "Families_15+_reads" = integer(), BC3_types_above15 = double(), check.names = FALSE, stringsAsFactors = FALSE)

# Iterate over input files
for(infile in infile_vec){
 # Read the input
 input <- paste(infile_dir, infile, sep="/")
 
 tab <- read_tsv(input, quote = "", col_types = "iidddid") # Remember to check that col_types match the input
 tab <- tab %>% rowwise() %>% mutate(BC_above15 = Family_count*BC3_above15_freq)
 
 sample_name <- sub(".bwa.sorted.bam.BC_dist.txt", "", infile)
 
 # Count number of read families with >15 reads and frequency of families with >15 different BC3_groups among them
 tab2 <- tab[tab$Family_size > 15, ]
 if(nrow(tab2) > 0) {
   large_families <- sum(unlist(tab2$Family_count))
   high_BC3_types <- sum(unlist(tab2$BC_above15))/large_families
 }
 else {
   large_families <- high_BC3_types <- 0
 }
 
 # Add analysis results to the output table
 high_BC3[nrow(high_BC3)+1,] <- c(sample_name, large_families, high_BC3_types)
} 

# Write output
write_delim(high_BC3, outfile, delim = "\t")

############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")