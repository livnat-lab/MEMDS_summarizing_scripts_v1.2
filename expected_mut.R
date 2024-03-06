# This script outputs expected substitution frequency within ROI, based on substitution frequencies outside it. 

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
# Read-in nucleotides to check; use nucleotide alphabet to create substitutions' names
bases <- unlist(strsplit(param_tbl$bases, ","))
sub_names <- c()
for (b1 in bases) {
  for (b2 in bases) {
    mut <- paste0(b1,b2,collapse="")
    if (b1 != b2) {
      sub_names <- c(sub_names, mut)
    }
  }
}

# Define matrix to count different substitution types
base_counts <- matrix(numeric(length(sub_names)*2), nrow = 2, ncol = length(sub_names), dimnames = list(rows = c("NRS", "RS"), sub_names))

# Parse analyzed sequence range
ranges <- unlist(str_split(param_tbl$seq_range, ";"))
seq_range <- vector()

for (range0 in ranges) {
  seq_pos <- unlist(strsplit(range0, "-"))
  seq_range <- append(seq_range, seq_pos[1]:seq_pos[2])
}

seq_range <- sort(unique(seq_range))

# Parse ROI (range of interest) range
RS_ranges <- unlist(str_split(param_tbl$RS_range, ";"))
BSU_range <- vector ()

for (RS_range0 in RS_ranges) {
  RS_pos <- unlist(strsplit(RS_range0, "-"))
  BSU_range <- append(BSU_range, RS_pos[1]:RS_pos[2])
}

BSU_range <- sort(unique(BSU_range))

# Parse read positions to skip while counting substitutions
skip_pos <- args[5]
#skip_pos <- "72,73,74,81,82,83"

excl_pos <- unlist(strsplit(skip_pos, ","))
excl_pos <- sort(excl_pos)
if (excl_pos != "None" & excl_pos != "NONE") {excl_pos <- as.numeric(excl_pos)}

# Parse substitutions to exclude from counts (e.g CA,GA...)
skip_sub <- args[6]
#skip_sub <- "NONE"

excl_sub <- unlist(strsplit(skip_sub, ","))
excl_sub <- sort (excl_sub)

sub_names <- sub_names[!(sub_names %in% excl_sub)] # Remove from substitution names excluded substitutions

# Parse specific variants to exclude from counts (e.g 39TC,...)
skip_var <- args[7]
#skip_var <- "None"

excl_var <- unlist(strsplit(skip_var, ","))
excl_var <- sort (excl_var)

# Define data frames to collect substitution counts
counts_all <- data.frame(stringsAsFactors = FALSE)
############################################################
# Parse substitution count files
infile_dirs <- grep("BC3above[0-9]$", dir(results_dir), value=T)
cat("Detected", length(infile_dirs), "input files in", results_dir, "\n")

infile_counter <- 0

# Iterate over the input data
for(dir in infile_dirs){
  
  infile_counter <- infile_counter + 1
  cat(infile_counter, "/", length(infile_dirs), ":", dir, "\n")
  
  # Check if "others" folder is analyzed, adjust folder name accordingly for parsing 
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
  # Read in consensus count summary file to extract counts of analyzed families
  tab_summary <- read_lines(paste(results_dir, dir, "Mut_totals_summary.txt", sep="/"))
  fields <- c("Total no. of families excluding rejected, ambiguous, plasmids, and having indels",
              "Total no. of families with indels \\(excluding rejected, ambiguous and plasmids\\)")

  family_count <- 0

  for (field in fields) {
    data <- unlist (strsplit(tab_summary[grepl(field, tab_summary)], "\t"))
    if (is.null(data)) {data <- c(0,0)}
    family_count <- (family_count + as.numeric(data[2]))
  }
  
  ######
  # Read in summary file of substitution counts per position
  infile <- paste(results_dir, dir, "Mut_per_position_all.txt", sep="/")
  tab <- read_tsv(infile, col_types = c("ddcdddd"))
  
  # Initialize substitution count vector
  sub_counts_RS <- sub_counts_NRS <- rep("NA",length(sub_names))
  names(sub_counts_RS) <- names(sub_counts_NRS) <- sub_names
  
  # Count number of times each nucleotide appears in the reference sequence to normalize variant counts
  for (b1 in bases) {
    for (b2 in bases[bases != b1]) {
      base_counts["NRS", paste0(b1,b2)] <- nrow(tab[(tab$sequence == b1) & !(tab$read_position %in% BSU_range), "sequence"])
      base_counts["RS", paste0(b1,b2)] <- nrow(tab[(tab$sequence == b1) & (tab$read_position %in% BSU_range), "sequence"])
    }
  }
  
  # Set mutation counts at excluded positions to zero, to remove them from calculation of expected mutations
  tab[tab$read_position %in% excl_pos, 4:ncol(tab)] <- 0
  
  # Remove nucleotides at excluded positions from counts of nucleotide types in the reference sequence
  for (b1 in bases) {
    for (b2 in bases[bases != b1]) {
      base_counts["NRS", paste0(b1,b2)] <- base_counts["NRS", paste0(b1,b2)] - nrow(tab[(tab$read_position %in% excl_pos) & (tab$sequence == b1) & !(tab$read_position %in% BSU_range), "sequence"])
      base_counts["RS", paste0(b1,b2)] <- base_counts["RS", paste0(b1,b2)] - nrow(tab[(tab$read_position %in% excl_pos) & (tab$sequence == b1) & (tab$read_position %in% BSU_range), "sequence"])
    }
  }
  
  # Set substitution counts of skipped variants to zero
  if (excl_var != "NONE" && excl_var != "None") {
    for (var0 in excl_var) {
      # Extract position of excluded variant
      read_pos <- str_extract(var0, "^[0-9]+")
      if (is.na(read_pos)) {stop("Excluded variant ", var0, " lacks position")}
      
      # Extract reference and variant nucleotide from excluded variant name
      nucl1 <- str_sub(var0, nchar(read_pos)+1, nchar(read_pos)+1)
      nucl2 <- str_sub(var0, nchar(read_pos)+2)
      if (nchar(nucl2) != 1)  {stop("Excluded variant ", var0, " should be a point substitution")}
       
      # Set excluded variant counts to zero
      tab[(tab$sequence == nucl1) & (tab$read_position == read_pos), nucl2] <- 0
      
      # Remove excluded variant reference nucleotide(s) from counts of nucleotide types in the reference sequence for given variant
      # Overlap check with excluded positions is used to ensure that same nucleotide is not removed twice from the counts
      if (!(read_pos %in% excl_pos)) {
        base_counts["NRS", paste0(nucl1,nucl2)] <- base_counts["NRS", paste0(nucl1,nucl2)] - nrow(tab[(tab$read_position == read_pos) & (tab$sequence == nucl1) & !(tab$read_position %in% BSU_range), "sequence"])
        base_counts["RS", paste0(nucl1,nucl2)] <- base_counts["RS", paste0(nucl1,nucl2)] - nrow(tab[(tab$read_position == read_pos) & (tab$sequence == nucl1) & (tab$read_position %in% BSU_range), "sequence"])
      } 
      
    }
  }
  
  ######
  # Count substitutions observed in the read families, for each possible nucleotide combination
  for (b1 in bases) {
    for (b2 in bases) {

      if (b1 != b2) {
        mut <- paste0(b1,b2,collapse="")
        if (mut %in% excl_sub) {next} # Skip excluded substitutions
        
        if(nrow(tab[(tab$sequence == b1) & (tab$read_position %in% BSU_range), b2]) > 0) {
          sub_counts_RS[mut] <- sum(tab[(tab$sequence == b1) & (tab$read_position %in% BSU_range), b2])
        }
        
        if (nrow(tab[(tab$sequence == b1) & !(tab$read_position %in% BSU_range), b2]) > 0) {
          sub_counts_NRS[mut] <- sum(tab[(tab$sequence == b1) & !(tab$read_position %in% BSU_range), b2])*(base_counts["RS", mut]/base_counts["NRS", mut])
        }
      }
    }
  }
  
  # Organize column names for the output table
  sub_counts_all <- unlist(strsplit(paste(sub_counts_NRS, sub_counts_RS, sep = "\t"), "\t"))
  sub_names_tbl <- rep(sub_names, each = 2)
  sub_names_tbl[seq(1, length(sub_names_tbl), 2)] <- paste(sub_names_tbl[seq(1, length(sub_names_tbl), 2)], "_e", sep ="")
  sub_names_tbl[seq(2, length(sub_names_tbl), 2)] <- paste(sub_names_tbl[seq(2, length(sub_names_tbl), 2)], "_c", sep ="")
  
  # Create an output table with substitution counts and frequencies
  counts_all <- rbind(counts_all, c(t(cutoff_data), family_count, as.numeric(t(sub_counts_all))), stringsAsFactors=FALSE)  
  
  if (nrow(counts_all) == 1) {
    names(counts_all) <- c(cutoff_names, "Family_counts", sub_names_tbl) 
  }
}

############################################################
# Output expected mutation count tables
## Create the output directory
outfile_path <- args[4]
#outfile_path <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline2_R_processing/cons-count_summary/Subs2"

out_excep <- paste("excl_pos", paste(excl_pos, collapse = "_"),"subs", paste(excl_sub, collapse = "_"), 
                   "muts", paste(excl_var, collapse = "_"), sep="_")
outdir <- paste (outfile_path, out_excep, sep = "/")

#if(dir.exists(outdir)) {unlink(outdir, recursive = TRUE)}
if(!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }

## Write the output
outfile_all <- paste(outdir, paste0(cutoff_data["Name"], "_expected_subs.txt"), sep = "/")
write_delim(counts_all, outfile_all, delim = "\t")

############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")