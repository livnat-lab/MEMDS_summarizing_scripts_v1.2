# This script gathers information on BC3 group sizes within read families and outputs them from largest to smallest
# The script reports information on N largest BC3_groups, averaged across families of same size. The number of reported groups is user-defined
# The information is reported on all BC3 groups and on BC3 groups having more than one read in them only.
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
outprefix <- args[2]
col_report <- as.numeric(args[3])

#infile <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/mut_example.txt"
#outprefix <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/"
#col_report <- as.numeric("5")

############################################################
# Helper functions
# Calculate frequencies of BC3 groups within read family, from largest to smallest
top_freq <- function(Total_count, barcode3_all_counts, col_report){
  # Gather BC3 groups size info
  count_list <- as.numeric(unlist(strsplit(barcode3_all_counts, ";")))
  
  # If user-defined number of BC3 groups to report is higher than their number in the family, add zeroes to pad the array
  if (length(count_list) < col_report) {
    count_list <- append(count_list, numeric(col_report - length(count_list)))
  }
  
  # Report BC3_group frequencies within read family, sorted from largest to smallest
  count_freq <- sort((count_list/Total_count), decreasing = TRUE)
  return(paste(count_freq[1:col_report], collapse = ";"))
}

# Calculate frequencies of BC3 groups within read family, from largest to smallest, excluding groups with only one read
top_freq_above1 <- function(barcode3_all_counts, col_report){
  # Gather BC3 groups size info, filtering groups with size of one read
  count_list <- as.numeric(unlist(strsplit(barcode3_all_counts, ";")))
  count_list <- count_list[count_list > 1]
  
  # If user-defined number of BC3 groups to report is higher than their number in the family, add zeroes to pad the array
  if (length(count_list) < col_report) {
    count_list <- append(count_list, numeric(col_report - length(count_list)))
  }
  
  # Report relative frequency of each BC3_group within read family, out of all BC3_groups with more than one read in them
  totals <- sum(count_list)
  if (totals > 0) {count_freq <- as.numeric(count_list/totals)}
  else {count_freq <- numeric(col_report)}
  count_freq <- sort(count_freq, decreasing = TRUE)
  
  zz <- paste(append(count_freq[1:col_report], totals), collapse = ";")  
  return(zz)
}

# Report number of BC3 groups with more than one read in them and total number of reads, excluding BC3 groups with one read per read family
BC3_above1 <- function(barcode3_all_counts, cutoff) {
  # Gather BC3 groups size info, filtering groups with size of one read
  ordered_counts <- sort(as.numeric(unlist(strsplit(barcode3_all_counts, ";"))), decreasing = TRUE)
  effective_counts <- ordered_counts[ordered_counts > cutoff]
  
  # Calculate family size and number of BC3 groups w/o BC3 groups with one read 
  effective_size <- sum(effective_counts)
  effective_types <- length(effective_counts)
  
  # Report the results
  effective_counts_str <- paste(effective_counts, collapse = ";")
  effective_counts_str <- ifelse(nchar(effective_counts_str) == 0, 0, effective_counts_str)
  
  zz <- paste(effective_size, effective_types, effective_counts_str, sep = "-")
  return(zz)
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
tab_work1 <- tab_work2 <- distinct(tab_work, ID, .keep_all = TRUE)

# Add BC3 group size frequencies per read family ("top_freq" function)
tab_work1 <- tab_work1 %>% rowwise() %>% mutate(top_BC3 = top_freq(Total_count, barcode3_all_counts, col_report))
tab_work1 <- tab_work1 %>% separate(top_BC3, paste0("BC3_", 1:col_report), sep = ";", convert = TRUE)

# Add BC3 group size frequencies per read family, excluding BC3 groups with only one read in them ("top_freq_above1" function)
tab_work2 <- tab_work2 %>% rowwise() %>% mutate(top_BC3 = top_freq_above1(barcode3_all_counts, col_report))
tab_work2 <- tab_work2 %>% separate(top_BC3, c(paste0("BC3_", 1:col_report), "Effective_size"), sep = ";", convert = TRUE)
tab_work2 <- tab_work2[which(tab_work2$Effective_size > 0), ]

# Summarize frequencies of top-sized BC3 groups across same-size read families
for (i in 1:col_report) {
  col_choice <- paste0("BC3_", i)
  
  # Create temporary table with average frequency of BC3 group i
  dist_tmp1 <- tab_work1 %>% group_by(Total_count) %>% summarise(Family_count = n(), 
                                                                !! paste0(col_choice, "_mean") := mean(!! sym(col_choice)),
                                                                !! paste0(col_choice, "_SD") := sd(!! sym(col_choice)))
  
  dist_tmp2 <- tab_work2 %>% group_by(Total_count) %>% summarise(Family_count = n(), Effective_size_mean = mean(Effective_size),
                                                                   Effective_size_SD = sd(Effective_size),
                                                                   !! paste0(col_choice, "_mean") := mean(!! sym(col_choice)),
                                                                   !! paste0(col_choice, "_SD") := sd(!! sym(col_choice)))
  
  # Add data from temporary table to the large output table
  if (i == 1) {
    dist_table1 <- dist_tmp1
    dist_table2 <- dist_tmp2
  }
  else{
    cols_take1 <- which(!(colnames(dist_tmp1) %in% c("Family_count")))
    cols_take2 <- which(!(colnames(dist_tmp2) %in% c("Family_count", "Effective_size_mean", "Effective_size_SD")))
    
    dist_table1 <- full_join(dist_table1, dist_tmp1[, cols_take1], by = "Total_count")
    dist_table2 <- full_join(dist_table2, dist_tmp2[, cols_take2], by = "Total_count")
  }
}

# Rename "Total count" column to "Family size", to provide more accurate description of the data
dist_table1 <- rename(dist_table1, Family_size = Total_count)
dist_table2 <- rename(dist_table2, Family_size = Total_count)

# Output the results
outfile1 <- paste0(outprefix,"BC3_freq.txt")
outfile2 <- paste0(outprefix,"BC3_above1_freq.txt")

write_delim(dist_table1[order(dist_table1$Family_size, decreasing = TRUE), ], outfile1, delim = "\t")
write_delim(dist_table2[order(dist_table2$Family_size, decreasing = TRUE), ], outfile2, delim = "\t")

############################################################
# The second part of the script creates a table with ordered counts of BC3 groups (from largest to smallest) in the top sized 5000 families.
tab_top <- distinct(tab_work[order(tab_work$Total_count, decreasing = TRUE), ], ID, .keep_all = TRUE)

# Extract data on the 5000 top sized families from the input table, or report the whole table if there are less than 5000 families
row_take <- ifelse(nrow(tab_top) > 5000, 5000, nrow(tab_top))
tab_top <- tab_top[1:row_take, ]
tab_top <- tab_top %>% rowwise() %>% mutate(BC3_counts_ordered = paste(sort(as.numeric(unlist(strsplit(barcode3_all_counts, ";"))), decreasing = TRUE), collapse = ";"))

# Add information on BC3 groups with more than one read in them ("BC3_above1" function)
tab_top <- tab_top %>% rowwise() %>% mutate(BC3_effective = BC3_above1(barcode3_all_counts, 1))
tab_top <- tab_top %>% separate(BC3_effective, c("Effective_size", "Effective_types", "BC3_counts_above1"), sep = "-", convert = TRUE)
tab_top <- rename(tab_top, Family_size = Total_count, Total_BC3_types = total_BC3_types)

# Output the results
outfile3 <- paste0(outprefix,"top_5000.txt")
tab_top_out <- tab_top[, !(names(tab_top) %in% c("barcode3_all", "barcode3_all_counts"))]
write_delim(tab_top_out, outfile3, delim = "\t")
############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")