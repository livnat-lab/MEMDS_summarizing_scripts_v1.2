# This script checks size evenness of BC3 groups within read families.
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
#outfile <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/Eveness_out.txt"

############################################################
# Helper functions
# Calculate Shannon eveness of BC3 groups sizes within read family
Shannon_eveness <- function(total_BC3_types, barcode3_all_counts){
  
  counts_vector <- as.numeric(unlist(strsplit(barcode3_all_counts, ";")))
  freqs_vector <- counts_vector/sum(counts_vector)
  shannon_diversity <- -1*sum(freqs_vector * log(freqs_vector))
  eveness <- shannon_diversity/log(as.numeric(total_BC3_types))
  
  return(eveness)
}

############################################################
# Read the input file and count number of different BC3 types per read family
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

# Calculate Shannon evenness per read family
tab_work <- tab_work %>% rowwise() %>% mutate(Eveness = Shannon_eveness(total_BC3_types, barcode3_all_counts))
tab_work$Eveness[tab_work$total_BC3_types == 1] <- 1

# Check that there were no computational errors
if (any(is.nan(tab_work$Eveness))) { 
  stop("Evenness results couldn't be computed for families: ", 
       paste(tab_work$ID[which(is.nan(tab_work$Eveness))], collapse = ", "))
  }

# Summarize BC3 group size evenness across same-size families (mean, median, min)
dist_table <- tab_work %>% group_by(Total_count) %>% summarise(Family_count = n(), Mean_eveness = mean(Eveness),
                                                               Median_eveness = median(Eveness), SD_eveness = sd(Eveness),
                                                               Min_eveness = min(Eveness))
dist_table <- rename(dist_table, Family_size = Total_count)

# Output the data
write_delim(dist_table[order(dist_table$Family_size, decreasing = TRUE), ], outfile, delim = "\t")

############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")