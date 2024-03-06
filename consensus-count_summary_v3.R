# This script organizes information from "consensus count" files produced by the pipeline. These files count how many families have a particular mutation profile when different cutoff criteria are applied to them
# Takes as input all files whose names end with "cons-count.txt" in the specified input directory 
# Each input file must include:
#   - two tab-separated columns with "consensus" and "barcodes_count" headers
#   - several (possibly zero) rows starting with "rejected", counting families that failed to pass cutoff criteria
#   - WT family count line (one per file)
#   - all other count lines
#"x**********************************x"#

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

infile_dir <- args[2]
param_file <- args[3]

#infile_dir <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline2_R_processing/cons-count_summary"
#param_file <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline2_R_processing/cons-count_summary/summarize_params.txt"

############################################################
# Helper functions
# Function counting number of common mutations between defined plasmid mutation profile and given list of mutations
n_matches_with_plasmid <- function(v, k){
  return(length(intersect(v, plasmid_list[[k]])))
}

# A function adjusting reported mutation position, to skip zero when moving from before TSS (translation start site) to after it
coding_adjustment <- function(v,k) {
  if (v < 0) {
    k[k >= 0] <- k[k >= 0] + 1
  }
  if (v > 0) {
    k[k <= 0] <- k[k <= 0] - 1
  }
  if(v == 0) {
    stop ("Error: coding position can't be zero")
  }
  
  return(k)
}
############################################################
# Open parameter table and extract parameters for mutation profile analysis
#gene <- "HBD"
gene <- args[1]

param_tbl <- read_tsv(param_file, col_types = cols(.default = "c", filtering_pos = "d"))
param_tbl <- subset(param_tbl, ref_name == gene)

# Define parameters to check, based on user-supplied parameter table
# Nucleotide alphabet
bases <- unlist(strsplit(param_tbl$bases, ","))

# Reference sequence
WT_string <- param_tbl$ref_seq

# Mutational profile of plasmids (used to calculate enrichment factor)
# Additional mutational profiles to be excluded from general counts. Designated as "plasmids", since plasmids are excluded from general counts
plasmid_vec <- unlist(strsplit(param_tbl$plasmid_mut_profile, ";"))

if (plasmid_vec[1] != toupper("None")) {
  plasmid_list <- list()
  
  for (i in seq(length(plasmid_vec))) {
    variants <- unlist(str_split(plasmid_vec[i], ","))
    plasmid_list[[i]] <- variants
  }
}

# Create vectors of positions and substitutions to exclude from mutation summary counts
excluded_pos <- unlist(strsplit(param_tbl$pos_excluded, ","))
excluded_pos <- as.numeric(excluded_pos)
excluded_var <- unlist(strsplit(param_tbl$subs_excluded, ","))

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

# Adjust reported mutation positions relative to TSS
coding_diff <- unlist(strsplit(param_tbl$coding_start, ","))
coding_diff <- as.numeric(coding_diff)
coding_range0 <- coding_diff[1] + (coding_diff[2] * (seq_range - 1))
coding_range <- coding_adjustment(coding_diff[1], coding_range0)

# Organize specific variants (e.g. 39TC) we want to exclude from further counts
excl_mut_names <- unlist(strsplit(param_tbl$var_excluded, ","))
excl_mut_names_RS <- excl_mut_names_NRS <- c()

for (name in excl_mut_names) {
  excl_pos <- as.numeric(str_extract(name, "^[0-9]+"))
  
  if (excl_pos %in% BSU_range) {
    excl_mut_names_RS <- c(excl_mut_names_RS, name)  
  }else if (!(excl_pos %in% BSU_range) && !(is.na(excl_pos))){
    excl_mut_names_NRS <- c(excl_mut_names_NRS, name)
  }
}
############################################################
# Gather "consensus count" files from the input directory
infile_vec <- grep("cons-count.txt$", dir(infile_dir), value=T)
cat("Detected", length(infile_vec), "input files in", infile_dir, "\n")

infile_counter <- 0

# Iterate over input files
for(infile in infile_vec){
  
  infile_counter <- infile_counter + 1
  cat(infile_counter, "/", length(infile_vec), ":", infile, "\n")
  
  WT_vec <- unlist(strsplit(WT_string, ""))[seq_range]
  
  # Read input file, add column identifying each mutation profile as plasmid (T/F) and filter rejected sequence counts
  tab <- read_tsv(paste(infile_dir, infile, sep="/"), col_types = list(col_character(), col_double()))
  tab <- cbind(tab, Plasmid = "F", stringsAsFactors = FALSE)
  tab <- tab %>% filter(!grepl("^rejected.", consensus))   # removing all lines starting with "rejected."
  
  # Gather information on WT family count
  WT_count <- tab %>% 
    filter(consensus == "WT") %>% 
    pull(barcodes_count)
  
  if(length(WT_count) > 1){
    stop("problem with format of ", infile, ": has more than one WT line") 
  }
  
  if(length(WT_count) == 0){
    WT_count <- 0
  }

  tab <- tab %>% filter(consensus != "WT")  # Filter WT counts, to analyze families with mutations
  
  # Initialize variables to count plasmids and families with ambiguous positions
  plasmid_counts <- numeric(length(plasmid_list))
  plasmid_counts2 <- numeric(length(plasmid_list))
  N_count <- 0   

  # Organize special position counts
  pos_39_44_49_names <- unlist(strsplit(param_tbl$pos_names, ","))
  pos_39_44_49_counts <- numeric(length(pos_39_44_49_names))
  names(pos_39_44_49_counts) <- pos_39_44_49_names
  
  excl_mut_counts_RS <- numeric(length(excl_mut_names_RS))
  excl_mut_counts_NRS <- numeric(length(excl_mut_names_NRS))
  
  names(excl_mut_counts_RS) <- excl_mut_names_RS
  names(excl_mut_counts_NRS) <- excl_mut_names_NRS
  ####
  # Initialize tables for counting point mutations and indels observed in the read families
  substitution_count_mat <- matrix(0, length(bases), length(bases), dimnames = list(bases, bases))
  
  substitution_by_position_table <- tibble(read_position = seq_range,
                                           coding_sequence_position = coding_range,
                                           sequence = WT_vec)
  substitution_by_position_table[,bases] <- 0
  
  substitution_count_mat_excluded <- substitution_count_mat_BSU_excl_49 <- substitution_count_mat_outside_BSU <- substitution_count_mat
  substitution_count_mat_BSU <- substitution_count_mat_outside_BSU_all <- substitution_count_mat
  indel_table <- tibble(indel = character(0), count = numeric(0))
  
  # Parse barcode count table row by row to extract mutation profile data
  for(i in seq_len(nrow(tab))){
    variant <- unlist(strsplit(tab$consensus[i], ";"))
    barcode_count <- tab %>% slice(i) %>% pull(barcodes_count)
    
    # Count ambiguous positions (cases whereby it is unclear if certain position contains mutation)
    if(any(grepl("N", variant))){  # ambiguous
      N_count <- N_count + barcode_count
    }
    # Count families containing only filtering position mutation as 'WT'. Filtering position is used for diagnostic purpose only and is outside our region of interest, thus it is not ocunted towards analyzed mutations  
    else if(length(variant) == 1 && substr(variant, 1, nchar(param_tbl$filtering_pos)) == param_tbl$filtering_pos){
      WT_count <- WT_count + barcode_count
      tab$consensus[i] <- "WT"
    }
    
    # Analyze all remaining case
    else{
      if (plasmid_vec[1] != toupper("None")) {
        # Identify plasmid reads among read families
        is_plasmid <- F
        for(k in 1:length(plasmid_list)){
          match_size = ifelse(length(plasmid_list[[k]]) > 4, (length(plasmid_list[[k]]) - 1), length(plasmid_list[[k]]))
          match_size2 = length(plasmid_list[[k]])
          
          # Match mutation profile of a read against plasmid signature, allowing up to one mismatch if the signature contains more than 4 different mutations
          if(n_matches_with_plasmid(variant, k) >= match_size){
            plasmid_counts[k] <- plasmid_counts[k] + barcode_count
            is_plasmid <- T
          }
          
          # Match mutation profile of a read against plasmid signature, with no mismatches allowed
          if(n_matches_with_plasmid(variant, k) >= match_size2){
            plasmid_counts2[k] <- plasmid_counts2[k] + barcode_count
          }
        }
        
        # Mark reads originating from plasmids
        if(is_plasmid) {
          tab$Plasmid[i] <- "T"
        }
      }
      # Count and summarize mutations at "good reads" (not plasmids, rejected, etc.)
      if(!is_plasmid){
        filtered <- paste0("^", param_tbl$filtering_pos)
        variant <- variant[!grepl(filtered, variant)]  # Removing filtering position mutations
        
        if(any(grepl(filtered, variant))) {
          tab$consensus[i] <- paste(variant, collapse = ";")
        }
        
        # Parse "indel" variants, e.g.: variant <- unlist(strsplit("15AG;39C-;46CA", ";"))
        if(any(grepl("-", variant))){  
          
          # Count deletions, collapsing consecutive deletions into a single variant (e.g. 15G-, 16G- -> 15GG-)
          deletions <- sort(variant[grepl("-$", variant)])
          if(length(deletions) > 0){
            positions <- as.numeric(str_extract(deletions, "^[0-9]+"))
            from <- str_sub(deletions, nchar(positions) + 1, nchar(positions) + 1)
            contiguous_positions_list <- split(positions, cumsum(c(1, diff(positions) != 1)))
            lengths <- sapply(contiguous_positions_list, length)
            
            for(j in seq_along(lengths)){  
              indices <- (sum(lengths[1:j]) - lengths[j] + 1):sum(lengths[1:j])
              combined_deletion <- paste(positions[indices[1]], paste(from[indices], collapse=""), "-", sep="")
              indel_table <- indel_table %>%
                bind_rows(tibble(indel = combined_deletion, count = barcode_count))
            }
          }
          
          # Count insertions
          insertions <- variant[grepl("[0-9]-", variant)]
          indel_table <- indel_table %>%
            bind_rows(tibble(indel = insertions, count = barcode_count))
        }
        
        # Parse point mutations
        for(mutation in variant){ 
          # Count mutations at special positions, as defined by user
          if(any(mutation == pos_39_44_49_names)){
            pos_39_44_49_counts[mutation] <- pos_39_44_49_counts[mutation] + barcode_count
          }
          
          # Count excluded variants (e.g. 39TC) inside ROI
          if(any(mutation == excl_mut_names_RS)){
            excl_mut_counts_RS[mutation] <- excl_mut_counts_RS[mutation] + barcode_count
          }
          
          # Count excluded variants (e.g. 39TC) in flanking sequences
          if(any(mutation == excl_mut_names_NRS)){
            excl_mut_counts_NRS[mutation] <- excl_mut_counts_NRS[mutation] + barcode_count
          }
          
          # Extract point mutation info
          position <- as.numeric(str_extract(mutation, "^[0-9]+"))
          from <- str_sub(mutation, nchar(position) + 1, nchar(position) + 1)
          to <- str_sub(mutation, nchar(position) + 2)
          
          # Check that "from" part of the mutation matches expected reference sequence
          if(from != "-" && from != substr(WT_string, position, position)){
            stop("'from' of mutation ", mutation, ", line ", i+1, " of ", infile, ", is different from wildtype (should be ", substr(WT_string, position, position), ")")
          }
          
          # Internal check to validate that point mutation info was parsed correctly
          if(paste(position, from, to, sep="") != mutation){
            stop("a problem")
          }
          
          # Parse only point mutations that contain nucleotides in the user-defined alphabet
          if(from %in% bases && to %in% bases){
            
            # Add counts of all identified mutations to the relevant table(s)
            substitution_count_mat[from, to] <- substitution_count_mat[from, to] + barcode_count
            substitution_by_position_table[substitution_by_position_table$read_position == position, to] <- substitution_by_position_table[substitution_by_position_table$read_position == position, to] + barcode_count
            
            # Add counts of mutations in the ROI to the relevant tables: all mutations and no excluded positions
            if(position %in% BSU_range && !(position %in% excluded_pos)){
              substitution_count_mat_BSU_excl_49[from, to] <- substitution_count_mat_BSU_excl_49[from, to] + barcode_count
            }
            if(position %in% BSU_range){
              substitution_count_mat_BSU[from, to] <- substitution_count_mat_BSU[from, to] + barcode_count
            }
            
            # Add counts of mutations in the flanking sequences to the relevant tables: all mutations and no excluded positions
            if(!(position %in% BSU_range) && !(position %in% excluded_pos)){
              substitution_count_mat_outside_BSU[from, to] <- substitution_count_mat_outside_BSU[from, to] + barcode_count
            }
            if(!(position %in% BSU_range)){
              substitution_count_mat_outside_BSU_all[from, to] <- substitution_count_mat_outside_BSU_all[from, to] + barcode_count
            }
            
            # Add counts of mutations at excluded positions (both ROI and flanks) to theexcluded position counts table
            if(position %in% excluded_pos){
              substitution_count_mat_excluded[from, to] <- substitution_count_mat_excluded[from, to] + barcode_count
            }
            
          }
        }
      }
    }
  }
  
  # Remove the WT lines again, after converting mutations at filtering positions to 'WT'
  tab <- tab %>% filter(consensus != "WT")  
  ############################################################ 
  # Verify consistency of mutation counts in different count tables
  for(b1 in bases){
    for(b2 in setdiff(bases, b1)){
      if(substitution_count_mat[b1, b2] != substitution_by_position_table %>% filter(sequence == b1) %>% pull(b2) %>% sum)
        stop("total substitution matrix is inconsistent with substitution by position table")
      if(substitution_count_mat[b1, b2] != substitution_count_mat_outside_BSU[b1, b2] + 
                                           substitution_count_mat_BSU_excl_49[b1, b2] +
                                           substitution_count_mat_excluded[b1, b2])
        stop("substitution matrices are inconsistent")
    }
  }
  
  ############################################################
  # Write output files
  # Define output directory
  outfile_path <- args[4]
  #outfile_path <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline2_R_processing/cons-count_summary"
  
  summary_dir <- gsub(".bwa.sorted.bam|.no_overlap|.cons-count.txt", "",infile)
  outfile_dir <- paste (outfile_path, summary_dir, sep = "/")
  
  if(dir.exists(outfile_dir)) {unlink(outfile_dir, recursive = TRUE)}
  dir.create(outfile_dir, recursive = TRUE);
  
  #outfile_prefix <- str_sub(infile, end = -5)
  
  # Write file 1 - general summary of mutation counts, grouped by various criteria
  outfile_1_full_path <- paste(outfile_dir, "Mut_totals_summary.txt", sep="/")
  
  ## Output total family counts
  cat ("#### Total_family_counts ####\n", file=outfile_1_full_path) #First line of the output, no need to append
  
  cat("Total no. of families (excluding rejected)\t", WT_count + tab %>% pull(barcodes_count) %>% sum, "\n", sep="", file=outfile_1_full_path, append=T) 
  
  cat("WT\t", WT_count, "\n", sep="", file=outfile_1_full_path, append=T)
  
  if (plasmid_vec[1] != toupper("None")) {
    for(k in seq_along(plasmid_list)){
      cat("Plasmid ", k, " (", paste(plasmid_list[[k]], collapse=";"),  ")\t", plasmid_counts[k], "\t", plasmid_counts2[k], "\n", sep="", file=outfile_1_full_path, append=T)
    }
  }
  
  cat("Ambiguous variants\t", N_count, "\n", sep="", file=outfile_1_full_path, append=T)
  
  # Counts w/o ambiguous and plasmids
  n_families_WT_or_subs_only <- WT_count + 
    tab %>% 
    filter(!grepl("N|-", consensus) & !Plasmid == "T") %>% 
    pull(barcodes_count) %>% 
    sum
  
  cat("Total no. of families excluding rejected, ambiguous, plasmids, and having indels\t", n_families_WT_or_subs_only, "\n", sep="", file=outfile_1_full_path, append=T)
  
  cat("Total no. of families excluding rejected, WT, ambiguous, plasmids, and having indels\t", n_families_WT_or_subs_only - WT_count, "\n", sep="", file=outfile_1_full_path, append=T)
  
  # Indel counts, w/ and w/o plasmids
  n_families_with_indels <- tab %>% 
    filter(grepl("-", consensus)) %>% 
    pull(barcodes_count) %>% 
    sum
  
  n_families_good_indel <- tab %>% 
    filter(grepl("-", consensus) & !grepl("N", consensus) & !Plasmid == "T") %>% 
    pull(barcodes_count) %>% 
    sum
  
  cat("Total no. of families with indels (excluding rejected, including ambiguous and plasmids)\t", n_families_with_indels, "\n", sep="", file=outfile_1_full_path, append=T)
  
  cat("Total no. of families with indels (excluding rejected, ambiguous and plasmids)\t", n_families_good_indel, "\n", sep="", file=outfile_1_full_path, append=T)
  
  ## Output counts at special positions and counts of excluded variants
  cat ("\n#### Special_position_counts ####\n", file=outfile_1_full_path, append=T)
  cat(paste(names(pos_39_44_49_counts), " substitutions\t", pos_39_44_49_counts, "\n", sep=""), sep="", file=outfile_1_full_path, append=T)
  
  cat ("\n#### Excluded_variant_counts ####\n", file=outfile_1_full_path, append=T)
  if (length(names(excl_mut_counts_NRS)) > 0) {
    cat(paste(names(excl_mut_counts_NRS), " substitutions\t", excl_mut_counts_NRS, "\n", sep=""), sep="", file=outfile_1_full_path, append=T)
  }  
  
  if (length(names(excl_mut_counts_RS)) > 0) {
    cat(paste(names(excl_mut_counts_RS), " substitutions\t", excl_mut_counts_RS, "\n", sep=""), sep="", file=outfile_1_full_path, append=T)
  }
  
  if (sum(length(names(excl_mut_counts_NRS)), length(names(excl_mut_counts_RS))) == 0) {
    cat(paste ("NONE substitutions\t", 0, "\n", sep=""), sep="", file=outfile_1_full_path, append=T)
  }
  
  ## Output substitution counts, for each possible combination of nucleotides
  cat ("\n#### Counts_per_substitution ####\n", file=outfile_1_full_path, append=T)
  for(b1 in bases){
    for(b2 in setdiff(bases, b1)){
      cat(paste(b1, b2, " substitutions", sep=""), "\t", substitution_count_mat[b1, b2], "\n",  sep="", file=outfile_1_full_path, append=T)
    }
  }
  
  cat ("\n#### Counts_per_substitution_ROI ####\n", file=outfile_1_full_path, append=T)
  RS_excl <- ifelse(length(excluded_pos[excluded_pos %in% BSU_range]), excluded_pos[excluded_pos %in% BSU_range], "None")
  RS_excl <- paste(RS_excl, collapse = ", ")
  
  for(b1 in bases){
    for(b2 in setdiff(bases, b1)){
      cat(paste(b1, b2, " substitutions inside ROI", sep=""), "\t", substitution_count_mat_BSU[b1, b2], "\t", paste ("excluding pos. ", RS_excl, sep=""),"\t", substitution_count_mat_BSU_excl_49[b1, b2], "\n", sep="", file=outfile_1_full_path, append=T)
    }
  }
  
  cat ("\n#### Counts_per_substitution_outside_ROI ####\n", file=outfile_1_full_path, append=T)
  NRS_excl <- ifelse(length(excluded_pos[!(excluded_pos %in% BSU_range)]) && excluded_pos[1] != 0, excluded_pos[!(excluded_pos %in% BSU_range)], "None")
  NRS_excl <- paste(NRS_excl, collapse = ", ")
  
  for(b1 in bases){
    for(b2 in setdiff(bases, b1)){
      cat(paste(b1, b2, " substitutions outside ROI ", sep=""), "\t", substitution_count_mat_outside_BSU_all[b1, b2], "\t", paste("excluding pos. ", NRS_excl), "\t", substitution_count_mat_outside_BSU[b1,b2], "\n", sep="", file=outfile_1_full_path, append=T)
    }
  }
  
  ## Output total number of mutations found in the read families
  tot_subs_excl <- sum(substitution_count_mat)
  RS_subs_excl <- sum(substitution_count_mat_BSU_excl_49)
  NRS_subs_excl <- sum(substitution_count_mat_outside_BSU)
  
  ### Calculate total number of mutations, w/o excluded variants
  if (excluded_var[1] != toupper("None")) {
    for (mut in excluded_var) {
      nucl <- unlist(strsplit(mut, ""))
      if (length(nucl) != 2) {
        stop( "Err in excluded mutations: ", mut, " is not a point mutation\n")
      }
      
      tot_subs_excl <- tot_subs_excl - substitution_count_mat[nucl[1], nucl[2]]
      RS_subs_excl <- RS_subs_excl - substitution_count_mat_BSU_excl_49[nucl[1], nucl[2]]
      NRS_subs_excl <- NRS_subs_excl - substitution_count_mat_outside_BSU[nucl[1], nucl[2]]
    }
  }
  
  ### Output total counts
  cat ("\n#### Total_substitution_counts ####\n", file=outfile_1_full_path, append=T)
  cat("Total no. of substitutions\t", sum(substitution_count_mat), "\n", sep="", file=outfile_1_full_path, append=T)
  cat("Total no. of substitutions excluding muts: ", param_tbl$subs_excluded, "\t", tot_subs_excl, "\n", sep="", file=outfile_1_full_path, append=T)
  
  cat("Total no. of substitutions inside ROI excluding pos. ", RS_excl, "\t", sum(substitution_count_mat_BSU_excl_49), "\n", sep="", file=outfile_1_full_path, append=T)
  cat("Total no. of substitutions inside ROI excluding pos. ", RS_excl, ", muts: ", param_tbl$subs_excluded, "\t", RS_subs_excl, "\n", sep="", file=outfile_1_full_path, append=T)
  
  cat("Total no. of substitutions outside ROI excluding pos. ", NRS_excl, "\t", sum(substitution_count_mat_outside_BSU), "\n", sep="", file=outfile_1_full_path, append=T)
  cat("Total no. of substitutions outside ROI excluding pos. ", NRS_excl, ", muts: ", param_tbl$subs_excluded, "\t", NRS_subs_excl, "\n", sep="", file=outfile_1_full_path, append=T)
  
  ## Calculate error rate and expected mutation rate inside ROI (substitution rate in ROI minus error rate)
  ## Err_rate = Subs_count_flanks/(Family_number * flank_length); Mut_rate = Subs_count_ROI/(Family_number * ROI_length)
  cat ("\n#### Error+mut rate ####\n", file=outfile_1_full_path, append=T)
  # Calculate good families
  n_families_WT_subs_indel <- WT_count + 
    tab %>% 
    filter(!grepl("N", consensus) & !Plasmid == "T") %>% 
    pull(barcodes_count) %>% 
    sum
  
  # Calculate number of positions in ROI and flanking sequences 
  length_NRS <- as.numeric(length(seq_range) - length(BSU_range))
  length_RS <- as.numeric(length(BSU_range))
  
  # Remove excluded positions from total position counts
  for (pos0 in excluded_pos) {
    if (pos0 %in% BSU_range) {length_RS <- length_RS - 1}
    if (pos0 %in% seq_range && !(pos0 %in% BSU_range)) {length_NRS <- length_NRS - 1}
  }
  
  # Calculate error rate, using different family counts for normalization: WT, WT+subs, WT+subs+indel
  err_WT <- (NRS_subs_excl - sum(excl_mut_counts_NRS))/(WT_count * length_NRS)
  err_WT_subs_indel <- (NRS_subs_excl - sum(excl_mut_counts_NRS))/(n_families_WT_subs_indel * length_NRS)
  err_WT_subs <- (NRS_subs_excl - sum(excl_mut_counts_NRS))/(n_families_WT_or_subs_only * length_NRS)
  
  # Calculate mutation rate, with error rate correction, using different family counts for normalization: WT, WT+subs, WT+subs+indel
  mut_WT <- ((RS_subs_excl - sum(excl_mut_counts_RS))/(WT_count * length_RS)) -  err_WT
  mut_WT_subs_indel <- ((RS_subs_excl - sum(excl_mut_counts_RS))/(n_families_WT_subs_indel * length_RS)) - err_WT_subs_indel
  mut_WT_subs <- ((RS_subs_excl - sum(excl_mut_counts_RS))/(n_families_WT_or_subs_only * length_RS)) - err_WT_subs
  
  # Output: err_rate\tmut_rate
  cat("Error+mut rate WT:\t", err_WT, "\t", mut_WT, "\n", sep="", file=outfile_1_full_path, append=T)
  cat("Error+mut rate WT+subs:\t", err_WT_subs, "\t", mut_WT_subs, "\n", sep="", file=outfile_1_full_path, append=T)
  cat("Error+mut rate WT+subs+indel:\t", err_WT_subs_indel, "\t", mut_WT_subs_indel, "\n", sep="", file=outfile_1_full_path, append=T)
  
  ## Output indel counts, for all types of indels in the input file
  cat ("\n#### Indel_counts ####\n", file=outfile_1_full_path, append=T)
  indel_table <- indel_table %>%  # combining identical indels
    group_by(indel) %>% 
    summarize(count = sum(count)) %>% 
    arrange(desc(count))
  
  if(nrow(indel_table) > 0){
    for(j in 1:nrow(indel_table))
      cat(indel_table$indel[j], "\t", indel_table$count[j], "\n", sep="", file=outfile_1_full_path, append=T)
  }
  
  
  # Write file 2 - count point mutations by type and position in the analyzed read
  
  #outfile_2 <- paste(outfile_prefix, "summary_2.txt", sep="_")
  outfile_2_full_path <- paste(outfile_dir, "Mut_per_position_all.txt", sep="/")
  write_delim(substitution_by_position_table, outfile_2_full_path, delim="\t")
  
  
  # Write file 3 - count excluded point mutations (e.g. CA,GA; as specified in the parameter table) by type and position in the analyzed read
  ## Initialize exception table
  tab_excep <- substitution_by_position_table
  tab_excep[ ,4:ncol(tab_excep)] <- 0
  
  ## Add counts of the excluded variants
  if (excluded_var[1] != toupper("None")) {
    for (mut in excluded_var) {
      nucl <- unlist(strsplit(mut, ""))
      if (length(nucl) != 2) {
        stop( "Err in excluded mutations: ", mut, " is not a point mutation\n")
      }
      
      tab_excep[tab_excep$sequence %in% nucl[1], nucl[2]] <- substitution_by_position_table[substitution_by_position_table$sequence %in% nucl[1], nucl[2]]
    }
  }
  
  ## Write the output
  #outfile_3 <- paste(outfile_prefix, "summary_3.txt", sep="_")
  outfile_3_full_path <- paste(outfile_dir, "Mut_per_position_exceptions.txt", sep="/")
  write_delim(tab_excep, outfile_3_full_path, delim="\t")
  
  # File 4 - count point mutations by type and position in the analyzed read, without excluded mutations
  ## Initialize output table
  tab_no_excep <- substitution_by_position_table
  
  ## Remove counts of the excluded variants
  if (excluded_var[1] != toupper("None")) {
    for (mut in excluded_var) {
      nucl <- unlist(strsplit(mut, ""))
      if (length(nucl) != 2) {
        stop( "Err in excluded mutations: ", mut, " is not a point mutation\n")
      }
      
      tab_no_excep[tab_no_excep$sequence %in% nucl[1], nucl[2]] <- 0
    }
  }
  
  ## Write the output
  #outfile_4 <- paste(outfile_prefix, "summary_4.txt", sep="_")
  outfile_4_full_path <- paste(outfile_dir, "Mut_per_position_no_exceptions.txt", sep="/")
  write_delim(tab_no_excep, outfile_4_full_path, delim="\t")
}

############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")