## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
  INPUT <- snakemake@input
  OUTPUT <- snakemake@output
  WILDCARDS <- snakemake@wildcards
  THREADS <- snakemake@threads

  if(length(snakemake@log)>0) 
    sink(
      file = snakemake@log[[1]], 
      append = FALSE, 
      type = c("output", "message"), 
      split = TRUE
  )

  file.path("resources", paste0(snakemake@rule, ".RData")) |> 
    save.image()
}else{
  file.path("resources", "get_MAE.RData") |>
    load()
}
library(log4r)
logger <- create.logger(logfile = stdout(), level = "DEBUG")
info(logger, "Starting process...")
suppressMessages(library(data.table))
suppressMessages(library(S4Vectors))
suppressMessages(library(readxl))
suppressMessages(library(reshape2))

###############################################################################
# Load INPUT
###############################################################################
info(logger, "Loading input files...")
variantcall <- readRDS(INPUT$variantcall)
rnaseq <- readRDS(INPUT$rnaseq)
atac <- readRDS(INPUT$atac)
acgh <- readRDS(INPUT$acgh)

sample_metadata <- read_excel(INPUT$samplemetadata, sheet=1) |> as.data.frame()
rownames(sample_metadata) <- sample_metadata$sampleid
info(logger, "Input files loaded successfully")

###############################################################################
# Main Script
###############################################################################

info(logger, "Creating SummarizedExperiment object for RNA-seq data...")
rna <- rnaseq[,order(colnames(rnaseq))]
elementMetadata <- data.frame(gene_name = rownames(rna))
rownames(elementMetadata) <- rownames(rna)
colData <- data.frame(sampleid = colnames(rna))
rownames(colData) <- colnames(rna)
colData$batchid <- 1

rnaSE <- SummarizedExperiment::SummarizedExperiment(
  assays = list(assay = rna),
  rowData = elementMetadata,
  colData = colData,
  metadata = list(annotation = "rnaseq")
)
info(logger, "RNA-seq SummarizedExperiment object created")

info(logger, "Creating SummarizedExperiment object for variant calls...")
var <- as.data.frame(variantcall[,order(colnames(variantcall))])
var <- var[-which(var$Gene.refGene == "RPL21" & var$Chr == "10"),]
elementMetadata <- data.frame(gene_name = var$Gene.refGene, chr = var$Chr)
rownames(elementMetadata) <- rownames(var) <- var$Gene.refGene
var$Gene.refGene <- var$Chr <- NULL
colData <- data.frame(sampleid = colnames(var))
rownames(colData) <- colnames(var)
colData$batchid <- 1

varSE <- SummarizedExperiment::SummarizedExperiment(
  assays = list(assay = var),
  rowData = elementMetadata,
  colData = colData,
  metadata = list(annotation = "variantcalls")
)
info(logger, "Variant calls SummarizedExperiment object created")

info(logger, "Creating SummarizedExperiment object for ATAC-seq data...")
atac <- atac[,order(colnames(atac))]
elementMetadata <- data.frame(gene_name = rownames(atac))
rownames(elementMetadata) <- rownames(atac)
colData <- data.frame(sampleid = colnames(atac))
rownames(colData) <- colnames(atac)
colData$batchid <- 1

atacSE <- SummarizedExperiment::SummarizedExperiment(
  assays = list(assay = atac),
  rowData = elementMetadata,
  colData = colData,
  metadata = list(annotation = "atacseq")
)
info(logger, "ATAC-seq SummarizedExperiment object created")

info(logger, "Creating SummarizedExperiment object for aCGH data...")
assay <- acgh@assays@data$assay
elementMetadata <- acgh@elementMetadata
colData <- acgh@colData
assay <- assay[,order(colnames(assay))]
# rownames(colData) <- colData$sampleid

# colData <- colData[order(rownames(colData)),]
colData$batchid <- 1

acghSE <- SummarizedExperiment::SummarizedExperiment(
  assays = list(assay = assay),
  rowData = elementMetadata,
  colData = colData,
  metadata = list(annotation = "acgh")
)
info(logger, "aCGH SummarizedExperiment object created")

info(logger, "Creating MultiAssayExperiment object...")
tcl38_mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = list(
    "rnaseq" = rnaSE, 
    'acgh'= acghSE, 
    "variantcall" = varSE, 
    "atacseq" = atacSE
  ),
  colData = sample_metadata
)
info(logger, "MultiAssayExperiment object created successfully")

###############################################################################
# Save OUTPUT 
###############################################################################
info(logger, "Saving MultiAssayExperiment object...")
saveRDS(tcl38_mae, OUTPUT$mae)
info(logger, "Process completed successfully")