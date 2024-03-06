# This script reports proportion of primary barcodes within group of same-size read families that are found within Hamming distance of one from primary barcodes of larger size families
# # Based on the data in the "tables_consensus" files 

# Define libraries and I/O
rm(list=ls(all=TRUE))
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(foreach)
library(doParallel)


cat("Starting", date(), "\n")
x <- Sys.time()

args <- commandArgs(TRUE)
cat(c("Run args", args), sep = "\n")

infile <- args[1]
outfile <- args[2]

# infile <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Apl1_try/try.txt"
# outfile <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Apl1_try/Hamming_try2.txt"

############################################################
## Use this part to run parallel computing of Hamming distance
#cores <- as.integer(args[3])
#cores <- 8
#cl <- makeCluster(cores) #not to overload your computer
#registerDoParallel(cl)

############################################################
# Helper functions
# Find proportion of primary barcodes within Hamming distance of one from primary barcodes of larger families
## V1 -analyzed families' barcodes; V2 - larger families' barcodes
hamming_1_count <- function(v1, v2){
  bases <- c("A", "C", "G", "T")
  
  if(length(v1) == 0 | length(v2) == 0) {
    return(NA)
  }
  
  k1 <- length(v1)
  N <- nchar(v1)[1]
  
  # Create all possible permutations of the analyzed barcodes within Hamming distance of one from original ones 
  v1_expanded <- rep(v1, 4*N)
  
  for (i in 1:N) {
    substr(v1_expanded[((i-1)*4*k1 + 1):((i-1)*4*k1 + k1)], i, i) <- "A"
    substr(v1_expanded[((i-1)*4*k1 + k1 + 1):((i-1)*4*k1 + 2*k1)], i, i) <- "C"
    substr(v1_expanded[((i-1)*4*k1 + 2*k1 + 1):((i-1)*4*k1 + 3*k1)], i, i) <- "G"
    substr(v1_expanded[((i-1)*4*k1 + 3*k1 + 1):((i-1)*4*k1 + 4*k1)], i, i) <- "T"
  }
  
  # Return number of intersections between the permutations and actual barcodes of larger families
  # Intersect already returns unique values
  return(length(intersect(v1_expanded, v2))) 
}

# Estimate expected proportion of primary barcodes within Hamming distance of one from larger families' barcodes 
hamming_exp_count <- function(nA, nB, hamm_length){
  a <- (-3*hamm_length*nB/(4^hamm_length - nA)) * (1 - (nA/(4^hamm_length)))
  exp_count <- nA * (1 - exp(a))
  return(exp_count) # intersect already returns unique values
}
############################################################
# Read the input file; extract primary barcode sequences and family size data
tab <- read_tsv(infile)
bc5_tab <- tab %>% group_by(Barcode) %>% summarise(size = total_count, .groups = "drop")
sizes <- sort(unique(bc5_tab$size))

n_families <- n_remaining <- n_families_hamming_1 <- n_families_expected <- c()

# Iterate over size groups of read families and compare each groups to the larger size families to find Hamming distance
for (i in sizes) {
  v1 <- bc5_tab %>% filter(size == i) %>% pull(Barcode)
  v2 <- bc5_tab %>% filter(size > i) %>% pull(Barcode)
  
  n_families <- c(n_families, length(v1))
  n_remaining <- c(n_remaining, length(v2))
  n_families_hamming_1 <- c(n_families_hamming_1, hamming_1_count(v2, v1))
  n_families_expected <- c(n_families_expected, hamming_exp_count(length(v1), length(v2), (nchar(v1)[1] - 4)))
}

############################################################
# Create output data frame from the Hamming analysis results
out_tab <- tibble(Family_size = sizes, Family_num = n_families, Families_gt_size = n_remaining,
                  Hamming_1 = n_families_hamming_1, Hamming_gt1 = (n_families - n_families_hamming_1),
                  Hamming_1_fract = (n_families_hamming_1/n_families), Hamming_gt1_fract = (1 - (n_families_hamming_1/n_families)),
                  Hamming_1_exp_fract = (n_families_expected/n_families), Obs_to_exp = (n_families_hamming_1/n_families_expected))

# Write the output
write_delim(out_tab, outfile, delim="\t")

############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")