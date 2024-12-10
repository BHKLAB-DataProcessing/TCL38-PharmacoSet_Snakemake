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
  # Saves the workspace to "resources/"make_sampleMetadata"
  file.path("resources", paste0(snakemake@rule, ".RData")) |> 
    save.image()
}else{
  # If the snakemake object does not exist, load the workspace
  file.path("resources", "make_sampleMetadata.RData") |>
    load()
}
library(log4r)
library(data.table)
library(readxl)
logger <- create.logger(logfile = stdout(), level = "DEBUG")
info(logger, "Logger initialized.")

# Set this for mapping 
options("mc.cores" = THREADS)
info(logger, paste("Number of threads set to:", THREADS))

options("log_level" = "INFO")   # AnnotationGx logging level

###############################################################################
# Load INPUT
###############################################################################
info(logger, "Loading sample metadata from input file...")
sample <- readxl::read_excel(INPUT$sample_metadata, sheet = 1) |>
  data.table::as.data.table()
info(logger, "Sample metadata loaded successfully.")

###############################################################################
# Main Script
###############################################################################
info(logger, "Mapping sample IDs to accessions...")
mapped <- AnnotationGx::mapCell2Accession(sample$sampleid)
info(logger, "Mapping completed.")

missing <- mapped[is.na(mapped$accession), query]
warn(logger, paste("The following samples are missing from the annotation database:", paste(missing, collapse = ", ")))

### MANUAL FIXES ###
warn(logger, "Applying manual fixes...")
mapped[query =="MyLa"]$accession <- "CVCL_5341"
info(logger, "Manual fixes applied.")

# old columns in `sample`
# sampleid  Cell.Line     Subtype Maturity     TP53 Mutation
# Rename
# sampleid  Cell.Line     Subtype Maturity     TP53.Mutation
colnames(sample) <- gsub(" ", "_", colnames(sample))
info(logger, "Column names in sample metadata standardized.")

data.table::setnames(mapped, c("query", "accession"), c("sampleid", "cellosaurus.accession"))
info(logger, "Mapped data column names standardized.")

# Merge the mapped data with the sample metadata
info(logger, "Merging mapped data with sample metadata...")
merged <- merge(
  sample, 
  mapped[, .(sampleid,cellosaurus.accession)],
  by = "sampleid", 
  all.x = TRUE
)
info(logger, "Merge completed.")

info(logger, "Annotating cell accessions...")
annotated_accessions <- AnnotationGx::annotateCellAccession(
    accessions = merged$cellosaurus.accession,
)
info(logger, "Annotation completed.")

info(logger, paste("Number of annotated accessions:", nrow(annotated_accessions)))

info(logger, "Number of unique categories: ")
annotated_accessions[, .N, by = "category"]

info(logger, "Number of unique sexOfCell: ")
annotated_accessions[, .N, by = "sexOfCell"]

annotated_accessions[, synonyms := sapply(synonyms, function(x) paste(x, collapse = "; "))]
annotated_accessions[, diseases := sapply(diseases, function(x) paste(x, collapse = "; "))]

annotated_accessions[, c("crossReferences", "hierarchy", "comments") := NULL]

names(annotated_accessions) <- paste0("cellosaurus.", names(annotated_accessions))
annotated_accessions <- unique(annotated_accessions)
info(logger, "Annotated accessions data cleaned and standardized.")

info(logger, "Merging annotated accessions with merged data...")
final_annotated <- merge(
    merged, 
    annotated_accessions,
    by= "cellosaurus.accession",
    all.x = TRUE
) |> unique()
info(logger, "Final merge completed.")

###############################################################################
# Save OUTPUT 
###############################################################################

info(logger, "Saving output...")
data.table::fwrite(final_annotated, OUTPUT$processed_sample_metadata)
info(logger, "Output saved successfully.")