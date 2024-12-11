## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
# This snippet is run at the beginning of a snakemake run to setup the env
# Helps to load the workspace if the script is run independently or debugging
if(exists("snakemake")){
  INPUT <- snakemake@input
  OUTPUT <- snakemake@output
  WILDCARDS <- snakemake@wildcards
  THREADS <- snakemake@threads

  # setup logger if log file is provided
  if(length(snakemake@log)>0) 
    sink(
      file = snakemake@log[[1]], 
      append = FALSE, 
      type = c("output", "message"), 
      split = TRUE
  )

  # Assuming that this script is named after the rule
  # Saves the workspace to "resources/"process_combotherapy"
  file.path("resources", paste0(snakemake@rule, ".RData")) |> 
    save.image()
}else{
  # If the snakemake object does not exist, load the workspace
  file.path("resources", "process_combotherapy.RData") |>
    load()
}
library(log4r)
library(data.table)
library(reshape2)
library(readxl)
logger <- create.logger(logfile = stdout(), level = "DEBUG")
info(logger, "Starting process...")

###############################################################################
# Load INPUT
###############################################################################
info(logger, "Loading sample metadata from input file...")
sample <- readxl::read_excel(INPUT$samplemetadata, sheet = 1) |>
  data.table::as.data.table()
info(logger, "Sample metadata loaded successfully.")


info(logger, "Loading input files...")
raw_zip <- INPUT$raw

# temporary directory
temp_dir <- file.path("temp", "combotherapy")
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
info(logger, paste("Created temporary directory:", temp_dir))

# unzip files
unzip(raw_zip, exdir = temp_dir)
info(logger, paste("Unzipped files to:", temp_dir))

# read in files
files <- fs::dir_ls(temp_dir, glob = "*.xlsx", recurse = TRUE)
info(logger, paste("Found files:", length(files)))

###############################################################################
# Main Script
###############################################################################
prefix_1 <- "TableToDecrease_"
prefix_2 <- "Synergy_input_"
prefixes <- paste(prefix_1, prefix_2, sep = "|")

info(logger, "Processing files...")
raw_combo_dt <- lapply(files, function(file){
  base_name <- fs::path_file(file)

  sampleid <- gsub(prefixes, "", base_name) |>
    gsub("\\..*", "",x= _) |>
    gsub("_.*", "",x= _)
  
  df <- readxl::read_excel(file, sheet = 1) |>
    setDT()
  df$sampleid <- sampleid
  data.table::setnames(
    df, 
    c("Drug1", "Conc1", "Drug2", "Conc2", "Response"),
    c("treatment1id", "treatment1conc", "treatment2id", "treatment2conc", "viability")
  )
  info(logger, paste("Processed file:", file, "with sampleid:", sampleid))
  df
}) |> data.table::rbindlist(use.names = TRUE, fill = TRUE)
info(logger, "All files processed and combined into a data table.")

# map sample names
sample$toMap <- gsub("-", "", gsub(" ", "", sample$sampleid)) 
sample$toMap[sample$toMap == "HUT78"] = "Hut78"
sample$toMap[sample$toMap == "MOLT4"] = "Molt4"
sample$toMap[sample$toMap == "OCILy13.2"] = "OciLy132"
sample$toMap[sample$toMap == "Peer"] = "PEER"

raw_combo_dt$sampleid <- sample$sampleid[match(raw_combo_dt$sampleid, sample$toMap)]

###############################################################################
# Save OUTPUT 
###############################################################################
info(logger, "Saving output...")
data.table::fwrite(raw_combo_dt, OUTPUT$combo_processed)

# delete temporary directory
info(logger, "Cleaning up temporary directory...")
unlink(temp_dir, recursive = TRUE)
info(logger, "Temporary directory cleaned up.")