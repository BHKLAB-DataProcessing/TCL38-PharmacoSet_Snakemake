## Script to compile latest version of TCL38 PSet

setwd("~/TCL38")

# load libraries
suppressMessages(library(data.table))
suppressMessages(library(S4Vectors))
suppressMessages(library(CoreGx))
suppressMessages(library(PharmacoGx))
suppressMessages(library(readxl))
suppressMessages(library(reshape2))

# read in RDS objects from Aleksandr
var <- readRDS("RDSobjects/VariantCallsTCL.RDS")
rna <- readRDS("RDSobjects/TCL_RNA.RDS")
atac <- readRDS("RDSobjects/atac2_transformed_rnagene_peaks_all_TSS.RDS")

# read in RDS objects from Minoru
acgh <- readRDS("RDSobjects/aCGH.RDS")
acghEM <- readRDS("RDSobjects/aCGH-elementmeta.RDS")

# read in cell line metadata
metadata <- read_excel("metadata/sample-metadata.xlsx", sheet = 1) |> as.data.frame()
metadata$sampleid <- toupper(gsub("-", "", gsub("\\.", "", metadata$sampleID)))
metadata <- metadata[order(metadata$sampleid),]
rownames(metadata) <- metadata$sampleid



## === Process monotherapy drug response data === ###

# load in drug sensitivity data and format drug and sample names
dss <- read_excel("2024-10-09-drugscreening/DSS/DSS_all_corrected_TC38.xlsx", sheet = 1)
dss$DRUG <- gsub(" ", "_", dss$DRUG)
colnames(dss) <- toupper(gsub("-", "", colnames(dss)))
colnames(dss)[colnames(dss) == "OCILY13.2"] <- "OCILY132"

# create info sheet
#x <- tcl38@treatmentResponse$info
samples <- colnames(dss)[-c(1)]
drugs <- dss$DRUG
info <- data.frame(sampleid = rep(samples, each = length(drugs)), 
                   treatmentid = rep(drugs, length(samples)))
rownames(info) <- paste0(info$treatmentid, "_", info$sampleid)

# create profiles matrix
#x <- tcl38@treatmentResponse$profiles
profiles <- melt(dss)
rownames(profiles) <- paste0(profiles$DRUG, "_", profiles$variable)
profiles$DRUG <- profiles$variable <- NULL
colnames(profiles) <- "dss"



## ====== Creating SE object for RNA-Seq ====== ##

rna <- rna[,order(colnames(rna))]

# create elementMetadata
elementMetadata <- data.frame(gene_name = rownames(rna))
rownames(elementMetadata) <- rownames(rna)

# create se object
rnaSE <- SummarizedExperiment(
  assays = list(assay = rna), # counts data
  rowData = elementMetadata, # feature data
  colData = metadata, # sample info
  metadata = list(
    annotation = "rnaseq"
  )
)


### ====== Create SE object for variant calls ====== ###

var <- as.data.frame(var[,order(colnames(var))])

# remove duplicate RPL21 gene
var <- var[-which(var$Gene.refGene == "RPL21" & var$Chr == "10"),]

# create elementMetadata
elementMetadata <- data.frame(gene_name = var$Gene.refGene, chr = var$Chr)
rownames(elementMetadata) <- rownames(var) <- var$Gene.refGene
var$Gene.refGene <- var$Chr <- NULL

# create se object
colnames(var) <- toupper(gsub("\\.", "", colnames(var))) # reformat sample names to match other dataframes

varSE <- SummarizedExperiment(
  assays = list(assay = var), # variant calls
  rowData = elementMetadata, # feature data
  colData = metadata, # sample info
  metadata = list(
    annotation = "variantcalls"
  )
)


### ====== Create SE object for atacseq ====== ###

atac <- atac[,order(colnames(atac))]

# create elementMetadata
elementMetadata <- data.frame(gene_name = rownames(atac))
rownames(elementMetadata) <- rownames(atac)

# create se object
atacSE <- SummarizedExperiment(
  assays = list(assay = atac), # atac-seq counts data
  rowData = elementMetadata, # feature data
  colData = metadata, # sample info
  metadata = list(
    annotation = "atacseq"
  )
)


### ====== Create SE object for acgh ====== ###

colnames(acgh) <- toupper(gsub("_", "", colnames(acgh)))
acgh <- acgh[,order(colnames(acgh))]

# create se object
acghSE <- SummarizedExperiment(
  assays = list(assay = acgh), # atac-seq counts data
  rowData = acghEM, # feature data
  colData = metadata, # sample info
  metadata = list(
    annotation = "acgh"
  )
)


### ====== Create MAE object of all molecular profiles ====== ###

tcl38_mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = list(
    "rnaseq" = rnaSE,
    "variantcalls" = varSE,
    "atacseq" = atacSE,
    "acgh" = acghSE
  ),
  colData = metadata
)


### ====== Create TCL38 PharmacoSet ====== ###

tcl38_pset <- PharmacoGx::PharmacoSet2(
  name = "TCL38",
  treatment = tr_meta,
  sample = metadata,
  molecularProfiles = tcl38_mae,
  treatmentResponse = tcl38_tre,
  curation = list(sample = data.frame(), treatment = data.frame(), tissue = data.frame())
)

saveRDS(tcl38_pset, file = "PSet_TCL38.rds")