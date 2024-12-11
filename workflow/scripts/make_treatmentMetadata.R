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
  # Saves the workspace to "resources/"make_treatmentMetadata"
  file.path("resources", paste0(snakemake@rule, ".RData")) |> 
    save.image()
}else{
  # If the snakemake object does not exist, load the workspace
  file.path("resources", "make_treatmentMetadata.RData") |>
    load()
}
library(log4r)
library(data.table)
suppressMessages(library(CoreGx))
logger <- create.logger(logfile = stdout(), level = "DEBUG")
info(logger, "Starting process...")

# Set this for mapping 
options("mc.cores" = THREADS)
info(logger, paste("Number of threads set to:", THREADS))

options("log_level" = "INFO")   # AnnotationGx logging level

###############################################################################
# Load INPUT
###############################################################################


syn_data <- readRDS(INPUT$treatment_syn_rds) |> data.table::as.data.table()

extracted_meta <- syn_data[
	, {
		# Split the Treatment column into Drug1 and Drug2
		treatment_split <- tstrsplit(Treatment, " - ", fixed = TRUE)
    list(
      treatmentid = c(treatment_split[[1]], treatment_split[[2]]),
      targets = c(TargetsD1, TargetsD2)
    )
	}
] |> unique()
extracted_meta

mono_processed <- data.table::fread(INPUT$mono_processed)
combo_processed <- data.table::fread(INPUT$combo_processed)

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
combined_raw_treatment$sampleid <- combined_raw_treatment$sampleid

t1 <- combined_raw_treatment$treatment1id |> unique()
t2 <- combined_raw_treatment$treatment2id |> unique()

# combine into single unique list
treatments <- c(t1, t2) |> unique() 

# remove na values
treatments <- treatments[!is.na(treatments)]

# Create a data.table with the treatment ids
treatment_df <- data.table(treatmentid = treatments)

treatment_metadata <- merge(
  treatment_df, extracted_meta, by = "treatmentid", all.x = TRUE
)

###############################################################################

(compound_nameToCIDS <- AnnotationGx::mapCompound2CID(
    treatments,
    first = TRUE
))

na_rows <- compound_nameToCIDS[is.na(cids),]

warn(logger, paste("The following treatments are missing from the PubChem database:", paste(na_rows$name, collapse = ", ")))

compound_nameToCIDS <- compound_nameToCIDS[
  !is.na(compound_nameToCIDS$name) & !duplicated(compound_nameToCIDS$name) & !is.na(cids),
]
names(compound_nameToCIDS) <- c("treatmentid", "pubchem.CID")


# 3.0 Get properties from PubChem
# -------------------------------
properties=c('Title', 'MolecularFormula', 'InChIKey', 'CanonicalSMILES')
info(
  logger,
  "Getting the following properties from PubChem for ", nrow(compound_nameToCIDS), " compounds\n\t",
  paste(properties, collapse= "\n\t")
)

(pubchemProperties <- AnnotationGx::getPubchemCompound(
    ids = compound_nameToCIDS$pubchem.CID, 
    from = 'cid', 
    to = 'property', 
    properties= properties
))
names(pubchemProperties) <- paste0("pubchem.", names(pubchemProperties))


pubchemProperties[, pubchem.CID := as.character(pubchem.CID)]
compound_nameToCIDS[, pubchem.CID := as.character(pubchem.CID)]
# set each pubchem.CID column to character
pubchem_annotated <- merge(compound_nameToCIDS, pubchemProperties, by= "pubchem.CID", all.x = TRUE)
show(pubchem_annotated)


treatment_metadata_pubchem <- merge(
  treatment_metadata, 
  pubchem_annotated,
  by = "treatmentid",
  all.x = TRUE
)

# save checkpoint to restart from
file.path("resources", paste0(snakemake@rule, ".RData")) |> 
    save.image()

###############################################################################
# Annotate the Rest
###############################################################################
annotatedCIDs <- treatment_metadata_pubchem

unichem_sources <- AnnotationGx::getUnichemSources(T)
data.table::setkey(unichem_sources, Name)

sources_of_interest <- c("chembl", "drugbank", "chebi", "phamgkb", "lincs", "clinicaltrials", "nih_ncc", "fdasrs", "pharmgkb", "rxnorm")

sourceID <- unichem_sources[Name == "pubchem", SourceID]

info(logger, "\n\nAnnotating with unichem...")
unichem_mappings <- parallel::mclapply(annotatedCIDs$pubchem.CID, function(x){
  tryCatch({
    result <- AnnotationGx::queryUnichemCompound(type = "sourceID", compound = x, sourceID = sourceID)

    subset <- result$External_Mappings[Name %in% sources_of_interest, .(compoundID, Name)]
    # make Name the column names and the values the compoundID 
    subset$cid <- x
    data.table::dcast(subset, cid ~ Name, value.var = "compoundID", fun.aggregate = list)
  }, error = function(e) {
    error(logger, "Error: ", e$message)
    })
  } 
  ) |> data.table::rbindlist(fill = T)
show(unichem_mappings)


## ------------------------------------------------------------------------- ##

info(logger, "Unlisting unichem_mappings")
# for each column, if its a list then make it a string with a comma separator
for(col in names(unichem_mappings)){
  if(is.list(unichem_mappings[[col]])){
    unichem_mappings[[col]] <- sapply(unichem_mappings[[col]], function(x) paste(x, collapse = ","))
  }
}

# Rename columns like drugbank to unichem.DrugBank etc, using the unichem_sources "NameLabel" column 
names(unichem_mappings) <- paste(
    "unichem", 
    unichem_sources[names(unichem_mappings), gsub(" ", "_", NameLabel)], 
    sep = "."
    )

all_annotated_treatmentMetadata <- merge(
    annotatedCIDs, 
    unichem_mappings, 
    by.x = "pubchem.CID", 
    by.y = "unichem.NA",        # dw about name being NA itll get removed  by this merge
    all.x = T
)




###############################################################################
# Save OUTPUT 
###############################################################################
info(logger, "Saving output...")
data.table::fwrite(all_annotated_treatmentMetadata, OUTPUT$treatment_metadata)
info(logger, "Output saved successfully.")