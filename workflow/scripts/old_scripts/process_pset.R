setwd("~/TCL38")
#setwd("/home/bioinf/bhklab/julia/projects/TCL38")

suppressMessages(library(data.table))
suppressMessages(library(S4Vectors))
suppressMessages(library(CoreGx))
suppressMessages(library(PharmacoGx))
suppressMessages(library(readxl))
suppressMessages(library(reshape2))

# read in existing PSet object
load("old-pset-2024-10-08/PSet_TCL38_23.rds")
#tcl38 <- readRDS("PSet_TCL38_23.rds")

# read in RDS objects from Aleksandr
var <- readRDS("RDSobjects/VariantCallsTCL.RDS")
rna <- readRDS("RDSobjects/TCL_RNA.RDS")
atac <- readRDS("RDSobjects/atac2_transformed_rnagene_peaks_all_TSS.RDS")


### ===== Process Monotherapy Response ===== ###

# load in drug sensitivity data and format drug and sample names
dss <- read_excel("2024-10-09-drugscreening/DSS/DSS_all_corrected_TC38.xlsx", sheet = 1)
dss$DRUG <- gsub(" ", "_", dss$DRUG)
colnames(dss) <- toupper(gsub("-", "", colnames(dss)))
colnames(dss)[colnames(dss) == "OCILY13.2"] <- "OCILY132"

# create treatment response experiment
raw_single_tre <- melt(dss) |> setDT()
colnames(raw_single_tre) <- c("treatmentid", "sampleid", "dss")


### ===== Process Combination Therapy Response ===== ###

# read in files
files <- list.files("2024-10-26-combo", recursive = T, full.name = T, pattern = "xlsx")

# set file prefixes
p1 <- "2024-10-26-combo/ViabilityCombinations/Decrease_input_CTG_26_cell_lines/TableToDecrease_"
p2 <- "2024-10-26-combo/ViabilityCombinations/Synergy_input_"

# create dataframe to store combined response
raw_combo_tre <- data.frame(matrix(nrow = 0, ncol = 8))

# process individual files
for (file in files) {

  # get sampleid
  sampleid <- gsub("\\..*", "", (gsub("_.*", "", gsub(p1, "", gsub(p2, "", file)))))

  # format response
  df <- read_excel(file, sheet = 1) |> setDT()
  df$sampleid <- sampleid
  data.table::setnames(
    df, 
    c("Drug1", "Conc1", "Drug2", "Conc2", "Response"), c("treatment1id", "treatment1conc", "treatment2id", "treatment2conc", "viability")
  )

  # merge sample response
  raw_combo_tre <- rbind(raw_combo_tre, df)
}


### ===== Create Treatment Response Experiment ===== ###

# create drug index for single-agent response
Index <- data.frame(treatmentid = unique(raw_single_tre$treatmentid), 
                    index = 1:length(unique(raw_single_tre$treatmentid)))

# format individual response data
raw_single_tre$TreatmentIndex <- Index$index[match(raw_single_tre$treatmentid, Index$treatmentid)]
raw_single_tre$treatment2id <- raw_single_tre$treatment1conc <- raw_single_tre$treatment2conc <- raw_single_tre$ConcUnit <- raw_single_tre$viability <- raw_combo_tre$dss <- NA
raw_combo_tre$PairIndex <- raw_combo_tre$PairIndex + 173

raw_single_tre$Testing <- "Single-Agent"
raw_combo_tre$Testing <- "Combo"

data.table::setnames(raw_single_tre, c("treatmentid"), c("treatment1id"))
data.table::setnames(raw_combo_tre, c("PairIndex"), c("TreatmentIndex"))

# combined response data
combined_tre <- rbind(raw_single_tre, raw_combo_tre)

order_col <- c("TreatmentIndex", "Testing", "sampleid", "treatment1id", "treatment2id", "treatment1conc", "treatment2conc", "dss", "viability", "ConcUnit")
combined_tre <- combined_tre[, ..order_col]

# create treatment response experiment
tremapper <- CoreGx::TREDataMapper(rawdata = combined_tre)

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


### ===== Treatment Metadata ===== ###

# read in targets
treatment_meta <- unique(combined_tre[, c("TreatmentIndex", "Testing", "treatment1id","treatment2id"), with = FALSE])

### ===== Sample Metadata ===== ###

sample_meta <- read_excel("metadata/sample-metadata.xlsx", sheet = 1) |> as.data.frame()
rownames(sample_meta) <- sample_meta$sampleid


### ====== Create SE object for rnaseq ====== ###

rna <- rna[,order(colnames(rna))]

# create elementMetadata
elementMetadata <- data.frame(gene_name = rownames(rna))
rownames(elementMetadata) <- rownames(rna)

# create colData
colData <- data.frame(sampleid = colnames(rna))
rownames(colData) <- colnames(rna)
colData$batchid <- 1

# create se object
rnaSE <- SummarizedExperiment::SummarizedExperiment(
  assays = list(assay = rna), # counts data
  rowData = elementMetadata, # feature data
  colData = colData # sample info
)
rnaSE@metadata$annotation <- "rnaseq"


### ====== Create SE object for variant calls ====== ###

var <- as.data.frame(var[,order(colnames(var))])

# remove duplicate RPL21 gene
var <- var[-which(var$Gene.refGene == "RPL21" & var$Chr == "10"),]

# create elementMetadata
elementMetadata <- data.frame(gene_name = var$Gene.refGene, chr = var$Chr)
rownames(elementMetadata) <- rownames(var) <- var$Gene.refGene
var$Gene.refGene <- var$Chr <- NULL

# create colData
colnames(var) <- toupper(gsub("\\.", "", colnames(var))) # reformat sample names to match other dataframes
colData <- data.frame(sampleid = colnames(var))
rownames(colData) <- colnames(var)
colData$batchid <- 1

# create se object
varSE <- SummarizedExperiment::SummarizedExperiment(
  assays = list(assay = var), # variant calls
  rowData = elementMetadata, # feature data
  colData = colData # sample info
)
rnaSE@metadata$annotation <- "variantcalls"


### ====== Create SE object for atacseq ====== ###

atac <- atac[,order(colnames(atac))]

# create elementMetadata
elementMetadata <- data.frame(gene_name = rownames(atac))
rownames(elementMetadata) <- rownames(atac)

# create colData
colData <- data.frame(sampleid = colnames(atac))
rownames(colData) <- colnames(atac)
colData$batchid <- 1

# create se object
atacSE <- SummarizedExperiment::SummarizedExperiment(
  assays = list(assay = atac), # atac-seq counts data
  rowData = elementMetadata, # feature data
  colData = colData # sample info
)
atacSE@metadata$annotation <- "atacseq"


### ====== Create SE object for acgh ====== ###

acgh <- tcl38@molecularProfiles$acgh
assay <- acgh@assays@data$assay

# create elementMetadata
elementMetadata <- acgh@elementMetadata

# create colData
colData <- acgh@colData

# reformat sample names to match other dataframes
colnames(assay) <- toupper(gsub("_", "", colnames(assay)))
assay <- assay[,order(colnames(assay))]
colData$sampleid <- toupper(gsub("_", "", colData$cell))
rownames(colData) <- colData$sampleid
colData <- colData[order(rownames(colData)),]
colData$batchid <- 1

# create se object
acghSE <- SummarizedExperiment::SummarizedExperiment(
  assays = list(assay = assay), # atac-seq counts data
  rowData = elementMetadata, # feature data
  colData = colData # sample info
)
acghSE@metadata$annotation <- "acgh"


### ====== Create MAE object for Molecular Profiles ====== ###

tcl38_mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = list(
    "rnaseq" = rnaSE, 
    'acgh'= acghSE, 
    "variantcall" = varSE, 
    "atacseq" = atacSE
  ),
  colData = sample_meta
)

save(tcl38_mae, tre, treatment_meta, sample_meta, file = "temp-toJermiah.RData")

# create PSet
(TCL38_pset <- PharmacoGx::PharmacoSet2(
  name = "TCL38",
  treatment = treatment_meta,
  sample = sample_meta,
  molecularProfiles = tcl38_mae,
  treatmentResponse = tre,
  curation = list(sample = data.frame(), treatment = data.frame(), tissue = data.frame())
))

## save tcl38 pset
saveRDS(tcl38, file = "PSet_TCL38_23-dss.rds")