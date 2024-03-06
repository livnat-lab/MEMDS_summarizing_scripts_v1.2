# This script checks incidence of relabeling events within the read families 
# First it uses "wrongId.barcodes" file to identify reads containing relabeling oligo
# Next, it groups relabeled reads in "wrongId.barcodes" file by primary barcode and to form read clusters (sharing same BC5). These clusters are sorted by size and the script checks how many read groups are found in each size category.
# After that it takes "mutationFrequncyPerBarcode" file and checks what proportion of "good families" that are used to analyze mutations in the samples share primary barcode with the relabeled reads. The families are grouped by size. 

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
oligo_seq <- args[3]
BC5_len <- as.numeric(args[4])
BC_id <- args[5]
read_pos <- args[6]
id_length <- args[7]
mut_input <- args[8]

# Example of input parameters
# cat(c(infile, outfile, oligo_seq, BC5_len, BC_id, read_pos, id_length, mut_input), sep = "\n")
# 
# infile <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/Apl1_exp_idx1.assembled.filtered.fastq.trimmed.wrongId.barcodes"
# outfile <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/relabeling_dist.txt"
# oligo_seq <- "TTTTACACT" # relabeling oligo
# BC5_len <- as.numeric("18") # primary barcode length
# BC_id <- "GATCCT,CTAA" # sample identifiers, as listed in parameter table (BC5, BC3)
# read_pos <- "15,-6"# first position of identifier seq, as listed in parameter table (BC5, BC3)
# id_length <- "4" # identifier length (BC5)
# mut_input <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/Apl1_exp_idx1.APL1.bwa.sorted.bam.mutationFrequncyPerBarcode.txt"
############################################################
# Helper functions
# Check if primary barcode identifier fits expected sequence (user-defined, changes per sample)
check_BC5_id <- function(seq5, BC5_ids, BC5_pos, id_lengths) {
  flag <- 0
  
  for (i in seq(length(BC5_ids))) {
    id <- substr(BC5_ids[i], 1, 1+nchar(id_lengths[i])-1)
    pos <- BC5_pos[i]
    detected_id <- substr(seq5, pos, pos+nchar(id)-1)
    
    if(detected_id == id) {flag <- flag+1}
  } 
  
  return(ifelse(flag > 0, 1, 0))
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
# Read "wrongId.barcodes" file and extract reads with control oligo
tab <- read_tsv(infile, quote = "", col_types = "cc")
names(tab) <- c("seq5", "seq3")

tab <- tab[which(tab$seq3 == oligo_seq), 1:2]
tab <- tab %>% group_by(seq5) %>% mutate(Family_size = n())

# If no relabeling occurred - write out empty table and exit
if (nrow(tab) == 0) {
  cat("No relabeling was found", file = outfile)
  stop ("No relabeling was found\n")
}

# Evaluate relabeled read group sizes (read group - relabeled reads sharing same BC5)
tab_wrongID <- distinct(tab, seq5, .keep_all = TRUE)
dist_wrongID <- tab_wrongID %>% group_by(Family_size) %>% summarize(Family_count = n())

# Make a separate data frame of reads with good primary barcode (containing correct identifier)
tab_cleaned <- tab_wrongID[which(nchar(tab_wrongID$seq5) == BC5_len), ]
tab_cleaned <- tab_cleaned %>% rowwise() %>% mutate (good_BC5 = check_BC5_id(seq5, BC5_ids, BC5_pos, id_lengths))
tab_cleaned <- tab_cleaned[which(tab_cleaned$good_BC5 > 0), ]

# Evaluate relabeled read group sizes among reads with good primary barcode only
# Add the data to the relabeled read group size distribution table
if (nrow(tab_cleaned) > 0) {
  dist_wrongID_clean <- tab_cleaned %>% group_by(Family_size) %>% summarize(Family_count = n())
  dist_final <- full_join(dist_wrongID, dist_wrongID_clean, by = "Family_size", suffix = c("_all", "_cleaned"))
  dist_final[is.na(dist_final)] <- 0

# If no reads with good primary barcodes are found, set "Family_count_cleaned" column values to zero  
} else {
  dist_final <- cbind(dist_wrongID, numeric(nrow(dist_wrongID)))
  names(dist_final) <- c("Family_size", "Family_count_all", "Family_count_cleaned")
}

# Add percent distribution of  relabeled reads size group sizes
dist_final <- dist_final %>% rowwise() %>% mutate(Family_perc_all = (Family_count_all * 100)/sum(dist_final$Family_count_all),
                                                  Family_perc_cleaned = (Family_count_cleaned * 100)/sum(dist_final$Family_count_cleaned))

# Write out relabeled read size group distribution table
write_delim(dist_final, outfile, delim = "\t")
############################################################
# Read "mutationFrequncyPerBarcode" to identify read families related to the relabeled reads
mut_tab <- read_tsv(mut_input, quote = "", col_types = "ccdddddccddccccdd")
mut_tab <- distinct(mut_tab, ID, .keep_all = TRUE)

# Identify read families related to the relabeled reads
common_tab <- mut_tab[which(mut_tab$ID %in% tab_cleaned$seq5), c("ID", "Total_count")]

# Find how many families (total and relabeled read-related) are found in each size group
mut_dist <- common_tab %>% group_by(Total_count) %>% summarize(Family_count = n())
all_dist <- mut_tab %>% group_by(Total_count) %>% summarize(Family_count = n())

# Merge size group count results - total families and relabeled read-related
mut_dist <- merge(mut_dist, all_dist, by = "Total_count", all = FALSE, suffixes = c("_labeled", "_all"))
mut_dist <- rename(mut_dist, Family_size_mut = Total_count)

# Add percent distribution of size group sizes
mut_dist <- mut_dist %>% rowwise() %>% mutate(Family_perc_labeled = (Family_count_labeled * 100)/sum(mut_dist$Family_count_labeled),
                                              Group_frequency = Family_count_labeled/Family_count_all)

# Write the output 
cat("\n", file = outfile, append = TRUE)
write_delim(mut_dist, outfile, delim = "\t", col_names = TRUE, append = TRUE)

############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")