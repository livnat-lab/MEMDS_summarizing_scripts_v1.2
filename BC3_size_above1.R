# This script reports how prevalent are BC3 groups with more than one read in them in the read families
# Reported are mean and median values across same-size read families
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

infile <- args[1]
outfile <- args[2]

#infile <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/mut_example.txt"
#outfile <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/BC3_size_above1.txt"

############################################################
# Helper functions
# Report frequency of BC3 groups with more than read in them out of total BC3 groups in the read family
BC3_above1 <- function(barcode3_all_counts){
  count_list <- as.numeric(unlist(strsplit(barcode3_all_counts, ";")))
  types_above1 <- length(count_list[which(count_list > 1)])/length(count_list)
  return(types_above1)
}

# Report proportion of reads contained in BC3 groups with more than read in them out of total reads in the read family
BC3_reads_above1 <- function(barcode3_all_counts){
  count_list <- as.numeric(unlist(strsplit(barcode3_all_counts, ";")))
  reads_above1 <- sum(count_list[which(count_list > 1)])/sum(count_list)
  return(reads_above1)
}

############################################################
# Read the input file
tab <- read_tsv(infile, quote = "", col_types = "ccdddddccddccccdd")
tab <- tab %>% rowwise () %>% mutate(total_BC3_types = length(unique(unlist(strsplit(barcode3_all, ";")))), .after = barcode3_all)

# Check that mutation summary table doesn't produce non-matching values for family sizes and BC#_type counts, when looking at mutations in the same family
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
# Select BC3 group information from "mutationFrequncyPerBarcode" input table and retain only single entry per read family 
tab_work <- tab %>% select (ID, Total_count, barcode3_all, barcode3_all_counts, total_BC3_types)
tab_work <- distinct(tab_work, ID, .keep_all = TRUE)

# Add information on BC3 groups with more than one read in them per read family
tab_work <- tab_work %>% rowwise() %>% mutate(types_1 = BC3_above1(barcode3_all_counts), reads_1 = BC3_reads_above1(barcode3_all_counts))

# Summarize information on BC3 groups with more than one read in them across same-size read families (mean, median)
dist_table <- tab_work %>% group_by(Total_count) %>% summarise(Family_count = n(), Mean_types_above1 = mean(types_1),
                                                               Median_types_above1 = median(types_1), SD_types_above1 = sd(types_1), Sep = "\t",
                                                               Mean_reads_above1 = mean(reads_1), Median_reads_above1 = median(reads_1),
                                                               SD_reads_above1 = sd(reads_1))
dist_table <- rename(dist_table, Family_size = Total_count)

# Output the results
write_delim(dist_table[order(dist_table$Family_size, decreasing = TRUE), ], outfile, delim = "\t")

############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")