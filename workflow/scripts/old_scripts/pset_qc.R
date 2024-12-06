#qc
setwd("~/TCL38")
suppressMessages(library(data.table))
suppressMessages(library(S4Vectors))
suppressMessages(library(CoreGx))
suppressMessages(library(PharmacoGx))
load("PSet_TCL38_23.rds")

# colData matches sample dataframe
v_colData <- tcl38@molecularProfiles$variantcall@colData
r_colData <- tcl38@molecularProfiles$rnaseq@colData
a_colData <- tcl38@molecularProfiles$acgh@colData

table(rownames(v_colData) %in% rownames(tcl38@sample))
table(rownames(r_colData) %in% rownames(tcl38@sample))
table(rownames(a_colData) %in% tcl38@sample$sampleid)

# colData matches assays
vc <- tcl38@molecularProfiles$variantcall@assays@data$assay
rs <- tcl38@molecularProfiles$rnaseq@assays@data$exprs
ac <- tcl38@molecularProfiles$acgh@assays@data$assay

table(colnames(vc) %in% rownames(v_colData))
table(colnames(rs) %in% rownames(r_colData))
table(colnames(ac) %in% rownames(a_colData))


# elementMetadata matches assays
v_emData <- tcl38@molecularProfiles$variantcall@elementMetadata
r_emData <- tcl38@molecularProfiles$rnaseq@elementMetadata
a_emData <- tcl38@molecularProfiles$acgh@elementMetadata

table(rownames(vc) %in% v_emData$Coord)
table(rownames(rs) %in% r_emData$gene_id)
table(rownames(ac) %in% paste0(a_emData$cellid, "_", a_emData$abberation_num))


