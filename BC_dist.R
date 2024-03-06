# This script summarizes information on distribution of primary and secondary barcodes (family sizes and distribution and BC3 group counts)
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
#outfile <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/BC_example_dist.txt"
############################################################
# Read the input file and extract primary barcode IDs (BC5)
tab <- read_tsv(infile, quote = "", col_types = "ccdddddccddccccdd")
tab <- tab %>% rowwise () %>% mutate(total_BC3_types = length(unique(unlist(strsplit(barcode3_all, ";")))), .after = barcode3_all)
BC5s <- unique(tab$ID)

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
# Generate distribution table for BC data
# Select BC3 group information from "mutationFrequncyPerBarcode" input table and retain only single entry per read family 
tab_work <- tab %>% select (ID, Total_count, barcode3_all, total_BC3_types)
tab_work <- distinct(tab_work, ID, .keep_all = TRUE) %>% group_by(Total_count)

# Summarize family size distribution and BC3 group counts across same-size families (mean, median, max)
dist_table <- tab_work %>% summarise(Family_count = n(), Mean_BC3_types = mean(total_BC3_types), 
                                     Median_BC3_types = median(total_BC3_types), SD_BC3_types = sd(total_BC3_types), 
                                     Max_BC3_types = max(total_BC3_types), BC3_above15_freq = sum(total_BC3_types > 15)/n())
dist_table <- rename(dist_table, Family_size = Total_count)

# Output the data
write_delim(dist_table, outfile, delim = "\t")

############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")