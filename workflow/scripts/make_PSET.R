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
  # Saves the workspace to "resources/"make_PSET"
  file.path("resources", paste0(snakemake@rule, ".RData")) |> 
    save.image()
}else{
  # If the snakemake object does not exist, load the workspace
  file.path("resources", "make_PSET.RData") |>
    load()
}
library(log4r)
library(data.table)
suppressMessages(library(PharmacoGx))
logger <- create.logger(logfile = stdout(), level = "DEBUG")
info(logger, "Starting process...")

###############################################################################
# Load INPUT
###############################################################################
info(logger, "Loading input files...")
tre <- readRDS(INPUT$tre)
info(logger, "Treatment response experiment loaded successfully.")
mae <- readRDS(INPUT$mae)
info(logger, "Monotherapy activity experiment loaded successfully.")


sampleMetadata <- data.table::fread(INPUT$processed_sample_metadata)
info(logger, "Sample metadata loaded successfully.")

treatmentMetadata <- data.table::fread(INPUT$processed_treatment_metadata)
info(logger, "Treatment metadata loaded successfully.")

###############################################################################
# Main Script
###############################################################################


## ------------------------------------------------------------------------- ##
# SANITY CHECKS

# Need to make sure that the sampleMetadata and the colnames of the summarized
# experiments are the same

allSamples <- sapply(names(mae), function(se) colnames(mae[[se]])) |>
    unlist() |>
    unique()

# Check if all the samples in the MultiAssayExperiment are in the sampleMetadata
stopifnot(all(allSamples %in% sampleMetadata$sampleid))

treSamples <- tre@colData$sampleid |> 
    unique() 

if (!all(treSamples %in% sampleMetadata$sampleid)) {
  missing <- treSamples[!treSamples %in% sampleMetadata$sampleid]
  error(logger, paste("The following samples are missing from the sampleMetadata:", paste(missing, collapse = ", ")))
  stop("Not all samples in the TreatmentResponseExperiment are in the sampleMetadata")
}

allTreatments <- unique(tre@rowData$treatmentid)

# Check if all the treatments in the TreatmentResponseExperiment are in the treatmentMetadata
# this should be true since treatmentMetadata is created from the TreatmentResponseExperiment
stopifnot(all(allTreatments %in% treatmentMetadata$treatmentid))

## ------------------------------------------------------------------------- ##

treatmentMetadata <- treatmentMetadata[!duplicated(treatmentid),]
treatment_df <-  data.frame(
    treatmentMetadata, 
    row.names = treatmentMetadata$treatmentid)

sampleMetadata <- sampleMetadata[!duplicated(sampleid),]
sample_df <- data.frame(
    sampleMetadata, 
    row.names = sampleMetadata$sampleid
)


## ------------------------------------------------------------------------- ##

name <- "TCL38"

tryCatch({
  capture.output({
    pset <- PharmacoGx::PharmacoSet2(
      name = name,
      treatment = treatment_df,
      sample = sample_df,
      molecularProfiles = mae,
      treatmentResponse = tre,
      curation = list(
        sample = sample_df,
        treatment = treatment_df,
        tissue = data.frame()
      ),
      datasetType = "sensitivity"
    )
  }, type = "message", file = textConnection("stdout_msgs"))
}, error = function(e) {
  error(logger, sprintf("Error creating PharmacoSet: %s", e$message))
  stop(e)
}, warning = function(w) {
  warn(logger, sprintf("Warning during PharmacoSet creation: %s", w$message))
})

###############################################################################
# Save OUTPUT 
###############################################################################

info(logger, "Saving output files...")
saveRDS(pset, OUTPUT$pset)
info(logger, "Output files saved successfully.")