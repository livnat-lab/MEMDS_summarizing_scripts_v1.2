# This script checks relative frequency of individual nucleotides (A,C,T,G) within primary barcodes (BC5) of read families
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

# Input parameters, see example below
infile <- args[1]
outfile <- args[2]
BC5_len <- as.numeric(args[3])
BC_id <- args[4]
read_pos <- args[5]
id_length <- args[6]

# Example of input parameters
# infile <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/mut_example.txt"
# outfile <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/BC5_nucl_dist2.txt"
# BC5_len <- as.numeric("18") # primary barcode length
# BC_id <- "ATCGCT,TGTC" # sample identifiers, as listed in parameter table (BC5, BC3)
# read_pos <- "15,-6"# first position of identifier seq, as listed in parameter table (BC5, BC3)
# id_length <- "4" # identifier length (BC5)
############################################################
# Helper functions
# Return primary barcode identifier without the constant sample ID sequence
check_BC5_id <- function(seq5, BC5_ids, BC5_pos, id_lengths) {
  flag <- 0
  
  for (i in seq(length(BC5_ids))) {
    id <- substr(BC5_ids[i], 1, 1+id_lengths[i]-1)
    pos <- BC5_pos[i]
    detected_id <- substr(seq5, pos, pos+nchar(id)-1)
    
    if(detected_id == id) {
      return(str_sub(seq5, 1, pos-1))
      flag <- flag+1
      break
    }
  } 
  
  if (flag == 0) {return(as.character(0))}
}
############################################################
# Read the input file
tab <- read_tsv(infile, quote = "", col_types = "ccdddddccddccccdd")
tab <- tab %>% rowwise () %>% mutate(total_BC3_types = length(unique(unlist(strsplit(barcode3_all, ";")))), .after = barcode3_all)

# Check that mutation summary table doesn't produce non-matching values for family sizes and BC3_type counts, when looking at mutations in the same family
tab_check <- tab %>% select (ID, Total_count, total_BC3_types)
tab_check <- tab_check %>% group_by(ID)
check_results <- tab_check %>% summarise(Totals_diff = any(diff(Total_count) != 0),
                                         BC3_diffs = any(diff(total_BC3_types) != 0))

if (sum(check_results$Totals_diff) > 0) {
  stop("Mismatch in total read counts was found for family entries: ", 
       paste(check_results$ID[which(check_results$Totals_diff == TRUE)], collapse = ", "))
}

if (sum(check_results$BC3_diffs) > 0) {
  stop("Mismatch in total BC3 groups was found for family entries: ", 
       paste(check_results$ID[which(check_results$BC3_diffs == TRUE)], collapse = ", "))
}
############################################################
# Parse sample identifier data
BC_ids <- unlist(strsplit(BC_id, ","))
BC5_ids <- unlist(strsplit(BC_ids[1], "|", fixed = "TRUE"))

# Parse first position of identifier sequence in the barcodes, relative to the first position of the read
id_pos <- unlist(strsplit(read_pos, ","))
BC5_pos <- unlist(strsplit(id_pos[1], "|", fixed = "TRUE"))
BC5_pos <- as.numeric(BC5_pos)

# Parse length of sample identifier sequence
id_lengths <- unlist(strsplit(id_length, "|", fixed = "TRUE"))
id_lengths <- as.numeric(id_lengths)

############################################################
# Select BC5 sequence information from "mutationFrequncyPerBarcode" input table and retain only single entry per read family 
tab_work <- tab %>% select (ID, Total_count)
tab_work <- distinct(tab_work, ID, .keep_all = TRUE)

# Remove the sample ID part of the barcode that remains constant in all barcodes 
tab_work <- tab_work %>% rowwise() %>% mutate (ID2 = check_BC5_id(ID, BC5_ids, BC5_pos, id_lengths))
tab_work <- tab_work[which(tab_work$ID2 != 0), ]

# Add information on the individual nucleotide frequencies within primary barcode
tab_work <- tab_work %>% rowwise() %>% mutate(A_freq = str_count(ID2, regex("A", ignore_case = TRUE))/nchar(ID2),
                                              C_freq = str_count(ID2, regex("C", ignore_case = TRUE))/nchar(ID2),
                                              T_freq = str_count(ID2, regex("T", ignore_case = TRUE))/nchar(ID2),
                                              G_freq = str_count(ID2, regex("G", ignore_case = TRUE))/nchar(ID2))

# Summarize information on nucleotide frequencies across BC5 of same-size families (mean,median)
dist_table <- tab_work %>% group_by(Total_count) %>% summarise(Family_count = n(), 
                                                               Mean_A = mean(A_freq), Median_A = median(A_freq), SD_A = sd(A_freq),
                                                               Mean_C = mean(C_freq), Median_C = median(C_freq), SD_C = sd(C_freq),
                                                               Mean_T = mean(T_freq), Median_T = median(T_freq), SD_T = sd(T_freq),
                                                               Mean_G = mean(G_freq), Median_G = median(G_freq), SD_G = sd(G_freq)
                                                               )
dist_table <- rename(dist_table, Family_size = Total_count)

# Output the data
write_delim(dist_table, outfile, delim = "\t")

############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")