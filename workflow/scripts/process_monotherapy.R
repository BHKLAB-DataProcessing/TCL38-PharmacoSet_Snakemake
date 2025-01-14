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
      type = c("output"), 
      split = TRUE
  )

  # Assuming that this script is named after the rule
  # Saves the workspace to "resources/"process_monotherapy"
  file.path("resources", paste0(snakemake@rule, ".RData")) |> 
    save.image()
}else{
  # If the snakemake object does not exist, load the workspace
  file.path("resources", "process_monotherapy.RData") |>
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
# load in cell line metadata
info(logger, "Loading sample metadata from input file...")
sample_metadata <- read.csv(INPUT$samplemetadata) |> as.data.frame()
rownames(sample_metadata) <- sample_metadata$sampleid <- sample_metadata$cellosaurus.accession

# load in drug sensitivity data and format drug and sample names
info(logger, "Loading drug sensitivity data...")

dss <- readxl::read_excel(INPUT$raw, sheet = 1) |> setDT()
data.table::setnames(dss,
                     c("HuT78", "Karpas299", "Karpas384", "PEER"),
                     c("Hut78", "Karpas299", "Karpas384", "Peer"))

colnames(dss)[2:39] <- sample_metadata$sampleid[match(colnames(dss)[2:39], sample_metadata$Cell.Line)]

###############################################################################
# Main Script
###############################################################################
info(logger, "Processing drug sensitivity data...")

# create data.table object
info(logger, "Creating data.table object...")

raw_single_dt <- reshape2::melt(dss) |> setDT()
colnames(raw_single_dt) <- c("treatmentid", "sampleid", "dss")

###############################################################################
# Save OUTPUT 
###############################################################################


data.table::fwrite(raw_single_dt, OUTPUT$mono_processed)
info(logger, "Saving processed data...")
data.table::fwrite(raw_single_dt, OUTPUT$mono_processed)
info(logger, "Process completed successfully.")