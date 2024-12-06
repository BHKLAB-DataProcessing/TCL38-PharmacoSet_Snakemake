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
  # Saves the workspace to "resources/"make_TRE"
  file.path("resources", paste0(snakemake@rule, ".RData")) |> 
    save.image()
}else{
  # If the snakemake object does not exist, load the workspace
  file.path("resources", "make_TRE.RData") |>
    load()
}
library(log4r)
library(data.table)
logger <- create.logger(logfile = stdout(), level = "DEBUG")
info(logger, "Starting process...")
###############################################################################
# Load INPUT
###############################################################################
info(logger, "Loading input files...")
mono_processed <- data.table::fread(INPUT$mono_processed)
combo_processed <- data.table::fread(INPUT$combo_processed)
info(logger, "Input files loaded successfully.")

###############################################################################
# Main Script
###############################################################################

# create drug index for single-agent response
info(logger, "Creating drug index for single-agent response...")
Index <- data.frame(treatmentid = unique(mono_processed$treatmentid), 
                    index = 1:length(unique(mono_processed$treatmentid)))
info(logger, "Drug index created.")

# format individual response data
info(logger, "Formatting individual response data...")
mono_processed$TreatmentIndex <- Index$index[match(mono_processed$treatmentid, Index$treatmentid)]
mono_processed$treatment2id <- mono_processed$treatment1conc <- mono_processed$treatment2conc <- mono_processed$ConcUnit <- mono_processed$viability <- combo_processed$dss <- NA
combo_processed$PairIndex <- combo_processed$PairIndex + 173

mono_processed$Testing <- "Single-Agent"
combo_processed$Testing <- "Combo"

data.table::setnames(mono_processed, c("treatmentid"), c("treatment1id"))
data.table::setnames(combo_processed, c("PairIndex"), c("TreatmentIndex"))

combined_raw_treatment <- rbind(mono_processed, combo_processed)
info(logger, "Individual response data formatted.")

order_col <- c("TreatmentIndex", "Testing", "sampleid", "treatment1id", "treatment2id", "treatment1conc", "treatment2conc", "dss", "viability", "ConcUnit")
combined_raw_treatment <- combined_raw_treatment[, ..order_col]
info(logger, "Combined raw treatment data created.")

# capitalize sampleid
combined_raw_treatment$sampleid_capitalized <- toupper(combined_raw_treatment$sampleid)

# create treatment response experiment
info(logger, "Creating treatment response experiment...")
tremapper <- CoreGx::TREDataMapper(rawdata = combined_raw_treatment)

guess <- CoreGx::guessMapping(
  tremapper,
  list(
    rowDataMap = c("TreatmentIndex", "Testing", "treatment1id", "treatment2id", "treatment1conc", "treatment2conc"),
    colDataMap = c("sampleid"),
    assayMap = c("TreatmentIndex", "Testing", "treatment1id","treatment2id", "treatment1conc", "treatment2conc", "sampleid", "dss", "viability")
  ),
  subset = TRUE
)

CoreGx::rowDataMap(tremapper) <- guess$rowDataMap
CoreGx::colDataMap(tremapper) <- guess$colDataMap
CoreGx::assayMap(tremapper) <- list(sensitivity = guess$assayMap)

(tre <- CoreGx::metaConstruct(tremapper))
info(logger, "Treatment response experiment created.")

###############################################################################
# Save OUTPUT 
###############################################################################
info(logger, "Saving output...")
saveRDS(tre, OUTPUT$tre)
info(logger, "Output saved successfully.")
