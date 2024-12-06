setwd("~/TCL38")

library(stringr)
library(readxl)
library(PharmacoGx)

options(stringsAsFactors = FALSE)

x <- read_excel("2024-10-09-drugscreening/raw/SNK6_all_drugs/220927_170_01 Timonen custom Annotations corrected.xlsx", sheet = 3)

input_dir <- "2024-10-09-drugscreening/raw/Re-screen"
# output_dir <- '/Users/minoru/Code/bhklab/JAKSTAT/TCL38_PSet-snakemake/data/processed'

sens_dir <- file.path(input_dir, "sensitivity")
dir.create(sens_dir)
unzip(zipfile = file.path(input_dir, "sensitivity.zip"), exdir = sens_dir)

files <- list.files(input_dir, full.name = T)
dfs <- list()
for (file in files) {
  if (str_detect(file, "CTG")) {
    df <- read.csv(file, skip = 12, col.names = F)
    df <- df[!df$ProductId %in% c("empty", "BzCl", "DMSO"), ]
    name <- str_replace(file, "Raw_Data_", "")
    name <- str_replace(name, ".xlsx", "")
    dfs[[name]] <- data.frame(df)
  }
}

all_drugs <- c()
all_cells <- c()
merged_df <- data.frame(matrix(ncol = length(colnames(dfs[[1]]))))
colnames(merged_df) <- colnames(dfs[[1]])
for (name in names(dfs)) {
  merged_df <- rbind(merged_df, dfs[[name]])
}
merged_df <- merged_df[!is.na(merged_df$DWell), ]
all_drugs <- sort(unique(merged_df$ProductId))
all_cells <- sort(unique(merged_df$screen_id))

cols <- c("cell", "drug", "Dose1", "Dose2", "Dose3", "Dose4", "Dose5")
dose_df <- data.frame(matrix(ncol = length(cols), nrow = 0))
viability_df <- data.frame(matrix(ncol = length(cols), nrow = 0))
colnames(dose_df) <- cols
colnames(viability_df) <- cols

for (cell in all_cells) {
  for (drug in all_drugs) {
    subset <- merged_df[merged_df$screen_id == cell & merged_df$ProductId == drug, ]
    subset <- subset[order(subset$Concentration), ]

    added_dose <- setNames(as.list(c(cell, drug, subset$Concentration)), cols)
    dose_df <- rbind(dose_df, added_dose)

    added_viability <- setNames(as.list(c(cell, drug, subset$inhibition_percent)), cols)
    viability_df <- rbind(viability_df, added_viability)
  }
}

# Convert screen_id into standardized cell names
samples <- read.csv(file.path(output_dir, "rnaseq_samples.csv"))
cells_standardized <- sort(unique(samples$cell))

modified <- str_to_lower(all_cells)
modified <- str_replace_all(modified, "[\\W_]", "")
cell_map <- as.data.frame(list(screenid = all_cells, modified = modified))
cell_map$standard_name <- NA
for (cell in cells_standardized) {
  comp <- str_to_lower(cell)
  comp <- str_replace_all(comp, "[\\W_]", "")
  found <- modified[str_detect(modified, comp)]
  print(found)
  if (length(found) > 0) {
    cell_map[cell_map$modified %in% found, ]$standard_name <- cell
  }
}
# manually add missing standard cell name and make modifications:
cell_map[cell_map$screenid == "SUH-DHL1", ]$standard_name <- "Su_DHL_1"
cell_map[cell_map$screenid == "CTG_Phera_NK_YS_w", ]$standard_name <- "NK_YS_w"
cell_map[cell_map$screenid == "CTG_Phera_NK_YS_wo", ]$standard_name <- "NK_YS_wo"
cell_map[cell_map$screenid == "CTG_Phera_TLBR2_w", ]$standard_name <- "TLBR2_w"
cell_map[cell_map$screenid == "CTG_Phera_TLBR2_wo", ]$standard_name <- "TLBR2_wo"

dose_df$cell <- unlist(lapply(dose_df$cell, function(cell) {
  cell_map[cell_map$screenid == cell, c("standard_name")]
}))
viability_df$cell <- unlist(lapply(viability_df$cell, function(cell) {
  cell_map[cell_map$screenid == cell, c("standard_name")]
}))

dose_df$drug <- str_replace_all(dose_df$drug, "\\s", "_")
viability_df$drug <- str_replace_all(viability_df$drug, "\\s", "_")

doses <- c("Dose1", "Dose2", "Dose3", "Dose4", "Dose5")
# invert the inhibition rate
for (col in doses) {
  viability_df[[col]] <- 100 - as.numeric(viability_df[[col]])
  dose_df[[col]] <- as.numeric(dose_df[[col]]) / 1000 # get micro molar concentration
}

rownames(dose_df) <- paste0(dose_df$drug, "_", dose_df$cell)
rownames(viability_df) <- paste0(viability_df$drug, "_", viability_df$cell)
raw.sensitivity <- cbind(dose_df, viability_df)
raw.sensitivity[, c("cell", "drug")] <- NULL
raw.sensitivity[, c("cell", "drug")] <- NULL

# Make the merged dataframe into an array
raw.sensitivity <- array(
  c(
    as.matrix(raw.sensitivity[, 1:length(doses)]),
    as.matrix(raw.sensitivity[, (length(doses) + 1):(2 * length(doses))])
  ),
  c(nrow(raw.sensitivity), length(doses), 2),
  dimnames = list(
    rownames(raw.sensitivity),
    colnames(raw.sensitivity[, 1:length(doses)]),
    c("Dose", "Viability")
  )
)

# sensitivity_info
sensitivity_info <- data.frame(matrix(data = NA, ncol = 4, nrow = length(rownames(dose_df))))
colnames(sensitivity_info) <- c("cellid", "drugid", "chosen.min.range", "chosen.max.range")
sensitivity_info$cellid <- dose_df$cell
sensitivity_info$drugid <- dose_df$drug
sensitivity_info$chosen.min.range <- dose_df$Dose1
sensitivity_info$chosen.max.range <- dose_df$Dose5
rownames(sensitivity_info) <- paste0(dose_df$drug, "_", dose_df$cell)

# sensitivity_profile
sensitivity_profile <- data.frame(matrix(data = NA, ncol = 5, nrow = 0))
colnames(sensitivity_profile) <- c(
  "aac_recomputed",
  "ic50_recomputed",
  "HS",
  "E_inf",
  "EC50"
)
calculated_profiles <- PharmacoGx:::.calculateFromRaw(raw.sensitivity)

for (row in rownames(raw.sensitivity[, , "Dose"])) {
  sensitivity_profile[row, ] <- c(
    calculated_profiles$AUC[row],
    calculated_profiles$IC50[row],
    calculated_profiles$pars[[row]]$HS,
    calculated_profiles$pars[[row]]$E_inf,
    calculated_profiles$pars[[row]]$EC50
  )
}

saveRDS(raw.sensitivity, file.path(output_dir, "raw_sensitivity.rds"))
saveRDS(sensitivity_info, file.path(output_dir, "sensitivity_info.rds"))
saveRDS(sensitivity_profile, file.path(output_dir, "sensitivity_profile.rds"))