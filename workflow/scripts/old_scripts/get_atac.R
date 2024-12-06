# Run this script after get_variant.R

setwd("~/TCL38")

suppressMessages(library(data.table))
suppressMessages(library(S4Vectors))
suppressMessages(library(CoreGx))
suppressMessages(library(PharmacoGx))

# read in existing PSet object
load("PSet_TCL38_23.rds")

# load in TMM counts
counts <- readRDS("ATACseq/TMM_counts_atac_TCL38.RDS")

# remove controls
counts <- counts[-which(rownames(counts) %in% rownames(counts)[grep("chrUn", rownames(counts))]),]
counts <- counts[-which(rownames(counts) %in% rownames(counts)[grep("random", rownames(counts))]),]

# read in metadata
#anno <- fread("ATACseq/TCL38_annotation.tsv")
#anno <- anno[anno$cell_line %in% colnames(counts),]

# map anno sample names to sample names in counts matrix, then change sample names to subrun
mapping <- c("NYKS" = "NK-YS", "CCRFHSB2" = "CCRF-H-SB2", "HUT78" = "Hut78", "Mac2A" = "Mac2a", 
             "SUDHL1" = "Su-DHL-1")
for (i in 1:length(colnames(counts))) {
    cell = colnames(counts)[i]
    if (cell %in% names(mapping)) {colnames(counts)[i] <- unname(mapping[cell])}
}

# create colData object
colData <- data.frame(sampleid = colnames(counts))
rownames(colData) <- colData$sampleid

# create elementMetadata object
elementMetadata <- as.data.frame(do.call(rbind, strsplit(rownames(counts), "_")))
colnames(elementMetadata) <- c("seqnames", "start", "end")
rownames(elementMetadata) <- rownames(counts)
elementMetadata$peak <- rownames(counts)

# create summarized experiment object for variant calls
atSE <- SummarizedExperiment(
  assays = list(assay = counts), # TMM in data frame
  rowData = elementMetadata, # variant Data
  colData = colData # sample info
)
atSE@metadata$annotation <- "atacseq"

# add in new molecular profile
tcl38@molecularProfiles <- list("rnaseq" = tcl38@molecularProfiles$rnaseq, 'acgh'= tcl38@molecularProfiles$acgh, "variantcall" = tcl38@molecularProfiles$variantcall, "atacseq" = atSE)


save(tcl38, file = "PSet_TCL38_23.rds")


##### compare samples
setwd("~/TCL38")
suppressMessages(library(data.table))
suppressMessages(library(S4Vectors))
suppressMessages(library(CoreGx))
suppressMessages(library(PharmacoGx))

load("PSet_TCL38_23.rds")

atac <- toupper(gsub("-", "", tcl38@molecularProfiles$atacseq@colData$sampleid))


rna <- toupper(gsub("-", "", unique(tcl38@molecularProfiles$rnaseq@colData$cell)))
var <- toupper(gsub("-", "", unique(tcl38@molecularProfiles$variantcall@colData$cell)))
acgh <- toupper(gsub("_", "", unique(tcl38@molecularProfiles$acgh@colData$cell)))
# acgh = "OCI-Ly13-2", rna = "OCI-Ly13.2"
acgh[-which(acgh %in% rna)] <- "OCILY13.2"


atac <- atac[order(atac)]
rna <- rna[order(rna)]
var <- var[order(var)]
acgh <- acgh[order(acgh)]

### FIX MAPPING