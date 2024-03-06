# This script re-calculates error and mutation rates within studied reads omitting user-defined positions and mutation types
# The goal of the script is to allow rate re-calculations without the need to run the whole mutation count summarizing script anew

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

results_dir <- args[2]
param_file <- args[3]

#results_dir <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline2_R_processing/cons-count_summary"
#param_file <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline2_R_processing/cons-count_summary/summarize_params.txt"

############################################################
# Open parameter table and extract parameters for mutation profile analysis
gene <- args[1]
#gene <- "APL1"

param_tbl <- read_tsv(param_file, col_types = cols(.default = "c", filtering_pos = "d"))
param_tbl <- subset(param_tbl, ref_name == gene)

############################################################
# Parse sequence of interest range
ranges <- unlist(str_split(param_tbl$seq_range, ";"))
seq_range <- vector()

for (range0 in ranges) {
  seq_pos <- unlist(strsplit(range0, "-"))
  seq_range <- append(seq_range, seq_pos[1]:seq_pos[2])
}

seq_range <- sort(unique(seq_range))

# Parse restriction site range
RS_ranges <- unlist(str_split(param_tbl$RS_range, ";"))
BSU_range <- vector ()

for (RS_range0 in RS_ranges) {
  RS_pos <- unlist(strsplit(RS_range0, "-"))
  BSU_range <- append(BSU_range, RS_pos[1]:RS_pos[2])
}

BSU_range <- sort(unique(BSU_range))

# Parse read positions to skip while counting substitutions
skip_pos <- args[5]
#skip_pos <- "60"

excl_pos <- unlist(strsplit(skip_pos, ","))
excl_pos <- sort(excl_pos)
if (excl_pos != "None" & excl_pos != "NONE") {excl_pos <- as.numeric(excl_pos)}

# Parse substitutions to exclude from counts (e.g CA,GA...)
skip_sub <- args[6]
#skip_sub <- "CA,GA,GT"

excl_sub <- unlist(strsplit(skip_sub, ","))
excl_sub <- sort (excl_sub)

# Parse specific variants to exclude from counts (e.g 39TC,...)
skip_var <- args[7]
#skip_var <- "NONE"

excl_var <- unlist(strsplit(skip_var, ","))
excl_var <- sort (excl_var)

# Define data frames to collect substitution counts
err_rate_all <- data.frame(stringsAsFactors = FALSE)
mut_rate_all <- data.frame(stringsAsFactors = FALSE)
############################################################
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

  # Extract information on cutoff criteria used to obtain given mutation counts
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
  # Read in cons-count summary file to extract counts of analyzed families
  tab_summary <- read_lines(paste(results_dir, dir, "Mut_totals_summary.txt", sep="/"))
  fields <- c("WT\t", "Total no. of families excluding rejected, ambiguous, plasmids, and having indels",
              "Total no. of families with indels \\(excluding rejected, ambiguous and plasmids\\)")
  
  family_count <- c(0,0,0)
  count <- c()
  
  for (i in seq(length(fields))) {
    data <- unlist (strsplit(tab_summary[grepl(fields[i], tab_summary)], "\t"))
    if (is.null(data)) {data <- c(0,0)}
    count[i] <- as.numeric(data[2])
  }
  
  family_count[1] <- count[1]
  family_count[2] <- count[2]
  family_count[3] <- count[2]+count[3]
  
  ######
  # Read in summary file of substitution counts per position
  infile <- paste(results_dir, dir, "Mut_per_position_all.txt", sep="/")
  tab <- read_tsv(infile, col_types = c("ddcdddd"))
  
  # Compute sequence length for error and mutation rate calculations
  length_NRS <- as.numeric(length(seq_range) - length(BSU_range)) # Flank length - for error-rate
  length_RS <- as.numeric(length(BSU_range)) # ROI length - for mutation rate
  
  # Remove excluded positions from total sequence length for rate calculations
  for (position in excl_pos) {
    if (position %in% BSU_range) {length_RS <- length_RS - 1}
    if (position %in% seq_range && !(position %in% BSU_range)) {length_NRS <- length_NRS - 1}
  }
  
  # Set counts at skipped positions to zero, to remove them from rate calculations
  tab[tab$read_position %in% excl_pos, 4:ncol(tab)] <- 0

  # Set counts of skipped variants to zero, to remove them from rate calculations
  if (excl_var != "NONE" && excl_var != "None") {
    for (var0 in excl_var) {
      read_pos <- str_extract(var0, "^[0-9]+")
      if (is.na(read_pos)) {stop("Excluded variant ", var0, " lacks position")}
      
      nucl1 <- str_sub(var0, nchar(read_pos)+1, nchar(read_pos)+1)
      nucl2 <- str_sub(var0, nchar(read_pos)+2)
      if (nchar(nucl2) != 1)  {stop("Excluded variant ", var0, " should be a point substitution")}
       
      tab[(tab$sequence == nucl1) & (tab$read_position == read_pos), nucl2] <- 0
    }
  }
  
  # Set counts of skipped substitutions to zero, to remove them from rate calculations
  if (excl_sub != "NONE" && excl_sub != "None") {
    for (sub0 in excl_sub) {
      nucl1 <- str_sub(sub0, 1, 1)
      nucl2 <- str_sub(sub0, 2)
      if (nchar(nucl2) != 1)  {stop("Excluded variant ", var0, " should be a point substitution")}
      
      tab[(tab$sequence == nucl1), nucl2] <- 0
    }
  }  
  
  # Sum surviving substitution counts, for flanks and ROI
  count_sub_NRS <- sum(tab[!(tab$read_position %in% BSU_range), 4:ncol(tab)])
  count_sub_RS <- sum(tab[tab$read_position %in% BSU_range, 4:ncol(tab)])
  
  # Calculate error and mutation rates (mutation rate is corrected by error rate)
  err_rate <- count_sub_NRS/(length_NRS*family_count)
  mut_rate <- count_sub_RS/(length_RS*family_count)
  mut_rate <- mut_rate - err_rate
  
  # Add rate calculation results to the output data frame
  err_rate_all <- rbind(err_rate_all, c(t(cutoff_data), family_count[3], as.numeric(t(err_rate))), stringsAsFactors=FALSE)  
  mut_rate_all <- rbind(mut_rate_all, c(t(cutoff_data), family_count[3], as.numeric(t(mut_rate))), stringsAsFactors=FALSE) 
  
  # Define output column names 
  if (nrow(err_rate_all) == 1) {
    names(err_rate_all) <- names(mut_rate_all) <- c(cutoff_names, "Family_counts", "WT", "WT+subs", "WT+subs+indel") 
  }
}

############################################################
# Output rate tables
## Create the output directory
outfile_path <- args[4]
#outfile_path <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline2_R_processing/cons-count_summary/re_err"

out_excep <- paste("excl_pos", paste(excl_pos, collapse = "_"),"subs", paste(excl_sub, collapse = "_"), 
                   "muts", paste(excl_var, collapse = "_"), sep="_")
outdir <- paste (outfile_path, out_excep, sep = "/")

if(!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }

## Write the output
outfile_err <- paste(outdir, paste0(cutoff_data["Name"], "_re-err.txt"), sep = "/")
write_delim(err_rate_all, outfile_err, delim = "\t")

outfile_mut <- paste(outdir, paste0(cutoff_data["Name"], "_re-mut.txt"), sep = "/")
write_delim(mut_rate_all, outfile_mut, delim = "\t")

############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")