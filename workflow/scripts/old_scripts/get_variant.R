setwd("~/TCL38")

library(stringr)
library(dplyr)
suppressMessages(library(data.table))
suppressMessages(library(S4Vectors))
suppressMessages(library(CoreGx))
suppressMessages(library(PharmacoGx))
library(readxl)
df <- as.data.frame(read_excel("202211_VariantCall/Annovar_merge/dt_merge.xlsx", sheet = 1))

# filter for variant calls to make @molecularProfiles$variantcall@assays@data@assay
assay <- df[,which(names(df) %in% names(df)[grep("^A0062001", names(df))])]
rownames(assay) <- df$Coord

# match sample order with rnaseq
assay <- assay[,order(gsub(".*_", "", colnames(assay)))]


# filter for variant call info to make @molecularProfiles$variantcall@elementMetadata
qc <- c(grep("^QUAL", names(df), value = TRUE), grep("^FILTER", names(df), value = TRUE))
elementMetadata <- df[,which(names(df) %in% c("Coord", "Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "ExonicFunc.refGene", "AAChange.refGene", qc))]


# create colData object
tcl38 <- readRDS("20230406/PSet_TCL38_4.2.rds")
mData <- tcl38@molecularProfiles$rnaseq@colData

colData <- data.frame(vc_sample_name = colnames(assay), sample_name = mData$sample_name, cell = mData$cell, num_id = gsub(".*_", "", colnames(assay)), sampleid = paste0(rownames(mData)))
rownames(colData) <- colData$sampaleid
colnames(assay) <- colData$sampleid

# create summarized experiment object for variant calls
vcSE <- SummarizedExperiment(
  assays = list(assay = assay), # variant calls in data frame
  rowData = elementMetadata, # variant Data
  colData = colData # sample info
)
vcSE@metadata$annotation <- "variantcall"

# add in new molecular profile
tcl38@molecularProfiles <- list("rnaseq" = tcl38@molecularProfiles$rnaseq, 'acgh'= tcl38@molecularProfiles$acgh, "variantcall" = vcSE)


save(tcl38, file = "PSet_TCL38_23.rds")








# what is everything else?

misc <- df[,-which(names(df) %in% names(df)[grep("^A0062001", names(df))])]
misc <- misc[,-which(names(misc) %in% c("Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene"))]
#misc <- misc[,-which(names(misc) %in% c("Chr", "Start", "End", "Ref", "Alt", "Func.refGene", "Gene.refGene", "ExonicFunc.refGene", "AAChange.refGene"))]
misc[1:5,1:15]

# check CHROM
chrom <- misc[,which(names(misc) %in% names(misc)[grep("^CHROM", names(misc))])]
chrom[1:5,1:5]
table(misc$Chr == chrom)
table(misc$Chr == chrom[,1])
table(misc$Chr == chrom[,5])
misc$Chr == chrom[,5]
head(chrom[,5])

# check POS !!!!!!! not perfect match
pos <- misc[,which(names(misc) %in% names(misc)[grep("^POS", names(misc))])]
pos[1:5,1:5]
table(misc$Start == pos)
table(misc$End == pos)

# check REF and ALT
ref <- misc[,which(names(misc) %in% names(misc)[grep("^REF", names(misc))])]
alt <- misc[,which(names(misc) %in% names(misc)[grep("^ALT", names(misc))])]
table(misc$Ref == ref)
table(misc$Alt == alt)