# This script calculates ratio of plasmid to WT families in the analyzed data at different cut-off criteria thresholds

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

results_dir <- args[1]
outfile_path <- args[2]
mut_profile0 <- args[3]

#results_dir <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline2_R_processing/cons-count_summary"
#outfile_path <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline2_R_processing/cons-count_summary"
#mut_profile0 <- "72A-,73A-,74G-,79-AGA,80TG;74-TTTC,77G-,78C-,79T-,80T-"
############################################################
# Define empty output table
ratio_table <- data.frame(stringsAsFactors = FALSE)

# Parse substitution count files
infile_dirs <- grep("BC3above[0-9]$", dir(results_dir), value=T)
cat("Detected", length(infile_dirs), "input files in", results_dir, "\n")

infile_counter <- 0

# Iterate over input data
for(dir in infile_dirs){
  
  infile_counter <- infile_counter + 1
  cat(infile_counter, "/", length(infile_dirs), ":", dir, "\n")
  
  # Check if "others" folder is analyzed, adjust folder name accordingly
  others_flag <- FALSE
  dir2 <- dir
  
  if (grepl("others", dir)) {
    others_flag <- TRUE
    dir2 <- sub("others.", "", dir2)
  }

  # Extract information on the used cutoff criteria
  cutoff_names <- c("Name",	"MutFreq",	"ReadCount",	"BC3_groups",	"BC3_above") 
  cutoff_data <- numeric(length(cutoff_names))
  names(cutoff_data) <- cutoff_names
  
  x1 <- unlist(str_split(dir2, "\\.", 3))
  x2 <- unlist(str_split(x1[3], "\\_"))
  
  cutoff_data["Name"] <- paste(x1[1], x1[2], sep="_")
  if(others_flag) {cutoff_data["Name"] <- paste(cutoff_data["Name"], "others", sep="_")}
  
  cutoff_data["MutFreq"] <- as.numeric(sub("mutFreq", "", unlist(x2[grepl("mutFreq", x2)])))
  cutoff_data["ReadCount"] <- as.numeric(sub("readCount", "", unlist(x2[grepl("readCount", x2)])))
  cutoff_data["BC3_groups"] <- as.numeric(sub("BC3WithMut", "", unlist(x2[grepl("BC3WithMut", x2)])))
  cutoff_data["BC3_above"] <- as.numeric(sub("BC3above", "", unlist(x2[grepl("BC3above", x2)])))
  
  #####  
  # Read in consensus count summary file to extract counts of analyzed families
  tab_summary <- read_lines(paste(results_dir, dir, "Mut_totals_summary.txt", sep="/"))
  
  mut_profiles <- unlist(str_split(mut_profile0, ";"))
  mut_profiles <- gsub(",", ";", mut_profiles)
  
  fields <- c("WT\t", mut_profiles)
  count <- c()
  
  for (i in seq(length(fields))) {
    data <- unlist (strsplit(tab_summary[grepl(fields[i], tab_summary)], "\t"))
    if (is.null(data) & i == 1) {data <- c(0,0)}
    if (is.null(data) & i > 1) {data <- c(0,0,0)}
    count <- append(count,as.numeric(data[-1]))
  }
  
  #####
  # Calculate the WT-plasmid ratios
  if (count[1] > 0) {ratios <- count[-1]/count[1]}
  else {ratios <- rep(NA, length(count[-1]))}
  
  # Create an output table with the WT-plasmid ratios
  ratio_table <- rbind(ratio_table, c(t(cutoff_data), as.numeric(t(count)), as.numeric(t(ratios))), stringsAsFactors=FALSE)
  
  count_header <- c()
  for (i in seq(length(fields[-1]))) {
    headers <- paste(c("P", "exact_P"), i, sep ="")
    count_header <- c(count_header, headers)
  }
  
  ratio_header <- paste(count_header, "WT_ratio", sep = "-")
  
  if (nrow(ratio_table) == 1) {
    names(ratio_table) <- c(cutoff_names, "WT", count_header, ratio_header) 
  }
}

############################################################
# Output plasmid-WT ratio table
outfile <- paste(outfile_path, paste0(cutoff_data["Name"], "_plasmid-WT_ratio.txt"), sep = "/")
write_delim(ratio_table, outfile, delim = "\t")

############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")