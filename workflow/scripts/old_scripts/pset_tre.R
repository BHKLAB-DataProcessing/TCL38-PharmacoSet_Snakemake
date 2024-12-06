set.seed(123)

library(CBWWorkshop2024)
data("NCI_ALMANAC_sample_metadata")
data("NCI60_molecular_data")

data("NCI_ALMANAC_raw")
data("NCI_ALMANAC_treatment_metadata")
# CREATING DUMMY TREATMENT RESPONSE EXPERIMENT
#

# Need: "dose" "viability" "sampleid" "treatmentid"

subset_treatmentmeta <- NCI_ALMANAC_treatment_metadata[, c("treatmentid", "cid", "inchikey")]

# grab only rows where "treatmentid" is less than 12 characters

subset_treatmentmeta <- subset_treatmentmeta[nchar(subset_treatmentmeta$treatmentid) < 12, ]

treatmentids <- subset_treatmentmeta$treatmentid

sampleids <- NCI_ALMANAC_raw$sampleid[1:50]

# Create a long table of treatment response data for each pair of treatmentid and sampleid
# Each row will represent a single treatment response measurement
# There will be the following doses
subset_raw_tr <- NCI_ALMANAC_raw[
  treatment2id == "",
  viability,
  by = .(treatment1id, treatment1dose, sampleid, PANEL, CONC1, tech_rep, bio_rep)
][
  treatment1id %in% treatmentids & sampleid %in% sampleids
]
data.table::setnames(
  subset_raw_tr, 
  c("treatment1id", "treatment1dose"), c("treatmentid", "treatmentdose")
)


tremapper <- CoreGx::TREDataMapper(rawdata = subset_raw_tr)

guess <- CoreGx::guessMapping(
  tremapper,
  list(
    rowDataMap = c("treatmentdose", "treatmentid", "tech_rep"),
    colDataMap = c("sampleid", "bio_rep"),
    assayMap = c("treatmentdose", "treatmentid", "sampleid", "tech_rep", "bio_rep")
  ),
  subset = TRUE
)

CoreGx::rowDataMap(tremapper) <- guess$rowDataMap
CoreGx::colDataMap(tremapper) <- guess$colDataMap
CoreGx::assayMap(tremapper) <- list(sensitivity = guess$assayMap)

(dummy_tre <- CoreGx::metaConstruct(tremapper))