# This script summarizes calculated error- and mutation- rates in the read families at different cut-off criteria thresholds

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

#results_dir <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3/Plots"
#outfile_path <- "/media/evgeni/Data_drive/Evgeni/Non-random_mutation/Pipeline/Pipeline-master/try_BC3"
############################################################
# Define empty output table
err_table <- mut_table <- data.frame(stringsAsFactors = FALSE)

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
  # Read in consensus count summary file to family sizes
  tab_summary <- read_lines(paste(results_dir, dir, "Mut_totals_summary.txt", sep="/"))
  fields <- c("Total no. of families excluding rejected, ambiguous, plasmids, and having indels",
              "Total no. of families with indels \\(excluding rejected, ambiguous and plasmids\\)")
  
  family_count <- 0
  
  for (i in seq(length(fields))) {
    data <- unlist (strsplit(tab_summary[grepl(fields[i], tab_summary)], "\t"))
    if (is.null(data)) {data <- c(0,0)}
    family_count <- family_count + as.numeric(data[2])
  }
  
  # Read in consensus count summary file to extract error and mutation rates
  tab_summary <- read_lines(paste(results_dir, dir, "Mut_totals_summary.txt", sep="/"))
  fields <- c('Error\\+mut rate WT\\:', 'Error\\+mut rate WT\\+subs\\:', 'Error\\+mut rate WT\\+subs\\+indel\\:')
  count_err <- count_mut <- c()
  
  for (i in seq(length(fields))) {
    data <- unlist (strsplit(tab_summary[grepl(fields[i], tab_summary)], "\t"))
    if (is.null(data)) {data <- c(0,0,0)}
    count_err <- append(count_err,as.numeric(data[2]))
    count_mut <- append(count_mut,as.numeric(data[3]))
  }
  
  ##### 
  # Create summary output tables
  err_table <- rbind(err_table, c(t(cutoff_data), family_count, as.numeric(t(count_err))), stringsAsFactors=FALSE)
  mut_table <- rbind(mut_table, c(t(cutoff_data), family_count, as.numeric(t(count_mut))), stringsAsFactors=FALSE)
  
  if (nrow(err_table) == 1) {
    names(err_table) <- names(mut_table) <- c(cutoff_names, "Family_number", "WT", "WT+subs", "WT+subs+indel") 
  }
}

############################################################
# Output summary tables
out_err <- paste(outfile_path, paste0(cutoff_data["Name"], "_error.txt"), sep = "/")
out_mut <- paste(outfile_path, paste0(cutoff_data["Name"], "_muts.txt"), sep = "/")

write_delim(err_table, out_err, delim = "\t")
write_delim(mut_table, out_mut, delim = "\t")

############################################################
y <- Sys.time()

cat (difftime(y,x, units = "sec"), " secs\n")
cat("Done", date(), "\n")