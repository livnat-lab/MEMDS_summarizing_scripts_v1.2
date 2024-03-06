# This script output plots depicting mutations per position in the analyzed reads, as a percent out of total observed mutations
# Plots are produced for each set of cutoffs analyzed

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

#results_dir <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3"
#param_file <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline2_R_processing/cons-count_summary/summarize_params.txt"

############################################################
# Open parameter table and extract parameters for mutation profile analysis
gene <- args[1]
#gene <- "APL1"

param_tbl <- read_tsv(param_file, col_types = cols(.default = "c", filtering_pos = "d"))
param_tbl <- subset(param_tbl, ref_name == gene)

############################################################
# Parse read positions to skip while counting substitutions
skip_pos <- args[5]
#skip_pos <- "None"

excl_pos <- unlist(strsplit(skip_pos, ","))
excl_pos <- sort(excl_pos)
if (excl_pos != "None" & excl_pos != "NONE") {excl_pos <- as.numeric(excl_pos)}

# Parse substitutions to exclude from counts (e.g CA,GA...)
skip_sub <- args[6]
#skip_sub <- "CA,GA"

excl_sub <- unlist(strsplit(skip_sub, ","))
excl_sub <- sort (excl_sub)

# Parse specific variants to exclude from counts (e.g 39TC,...)
skip_var <- args[7]
#skip_var <- "NONE"

excl_var <- unlist(strsplit(skip_var, ","))
excl_var <- sort (excl_var)

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

  # Extract sample name from the input data path
  x1 <- unlist(str_split(dir2, "\\.", 3))
  x2 <- unlist(str_split(x1[3], "\\_"))
  
  data_name <- paste(x1[1], x1[2], sep="_")
  if(others_flag) {data_name <- paste(cutoff_data["Name"], "others", sep="_")}

  #####
  # Read in consensus count summary file to extract counts of analyzed families
  # Uncomment the lines below if we want to add normalization of mutation counts by number of families at different cutoffs
  # tab_summary <- read_lines(paste(results_dir, dir, "Mut_totals_summary.txt", sep="/"))
  # fields <- c("Total no. of families excluding rejected, ambiguous, plasmids, and having indels",
  #             "Total no. of families with indels \\(excluding rejected, ambiguous and plasmids\\)")
  # 
  # count <- c()
  # 
  # for (i in seq(length(fields))) {
  #   data <- unlist (strsplit(tab_summary[grepl(fields[i], tab_summary)], "\t"))
  #   if (is.null(data)) {data <- c(0,0)}
  #   count[i] <- as.numeric(data[2])
  # }
  # 
  # family_count <- count[1]+count[2]
  
  #####
  # Read in counts file and extract substitution count
  infile <- paste(results_dir, dir, "Mut_per_position_all.txt", sep="/")
  tab <- read_tsv(infile, col_types = c("ddcdddd"))
  
  # Set counts at skipped positions to zero
  tab[tab$read_position %in% excl_pos, 4:ncol(tab)] <- 0

  # Set counts of skipped variants to zero
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
  
  # Set counts of skipped substitutions to zero
  if (excl_sub != "NONE" && excl_sub != "None") {
    for (sub0 in excl_sub) {
      nucl1 <- str_sub(sub0, 1, 1)
      nucl2 <- str_sub(sub0, 2)
      if (nchar(nucl2) != 1)  {stop("Excluded variant ", var0, " should be a point substitution")}
      
      tab[(tab$sequence == nucl1), nucl2] <- 0
    }
  }  
  
  #####
  # Calculate mutation counts per position as a percentage of total mutation counts
  z <- ifelse(sum(tab[, 4:ncol(tab)]) > 0, 100/sum(tab[, 4:ncol(tab)]), 0)
  tab[,4:ncol(tab)] <- tab[, 4:ncol(tab)] * z
  tab <- as.data.frame(tab)
  
  ############################################################
  # Prepare data table for plotting
  rownames(tab) <- paste(tab$read_position, tab$sequence, sep ="")
  plot_table <- t(tab[,4:ncol(tab)])
  plot_table <- plot_table[nrow(plot_table):1, ]
  
  # Create output directory for the plots
  outfile_path <- args[4]
  #outfile_path <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/Plots"
  
  out_excep <- paste("excl_pos", paste(excl_pos, collapse = "_"),"subs", paste(excl_sub, collapse = "_"), 
                     "muts", paste(excl_var, collapse = "_"), sep="_")
  
  outdir <- paste (outfile_path, out_excep, paste(x2, collapse = "_"), sep = "/")
  if(!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }
  
  outfile <- paste(outdir, paste0(data_name, ".variant_plot.pdf"), sep = "/")
  
  # Define plot parameters and create barplot of normalized mutation counts per read position
  out_width <- ncol(plot_table)*12/52 # Customize plot width
  pdf(file = outfile, width = out_width, height = 5) # Define output format and plot height
  
  op <- par(mar = c(5,5,4,2) + 0.1)
  y_max = max(colSums(plot_table))
  
  barplot(plot_table, 
          main = "Substitution per position",
          xlab = "Reference position", ylim = range(pretty(c(0, y_max))),
          col = c("Red", "Yellow", "Green", "Blue"), las = 2, cex.names = .85,
          xpd = TRUE,
          legend.text = rownames(plot_table),
          args.legend = list(x = "right", title = "Variant", inset = c(-0.03, 0), xpd = TRUE),
          beside = FALSE)  
  title(ylab = "Percent total variants", cex.lab = 1, line = 3)
    
  dev.off()
}

############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")