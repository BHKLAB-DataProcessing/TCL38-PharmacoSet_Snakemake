from pathlib import Path

RAWDATA = Path("rawdata")
PROCDATA = Path("procdata")
METADATA =  Path("metadata")
RESULTS = Path("results")

rule make_PSET:
  input:
    tre = PROCDATA / "TCL38_TRE.RDS",
    mae = PROCDATA / "TCL38_MAE.RDS",
    processed_sample_metadata = PROCDATA / "sample-metadata.csv",
    processed_treatment_metadata = PROCDATA / "treatment-metadata.csv"
  output:
    pset = RESULTS / "TCL38_PSet.RDS"
  log:
    "logs/make_PSET.log"
  script:
    "workflow/scripts/make_PSET.R"

rule make_treatmentMetadata:
  input:
    mono_processed = PROCDATA / "mono_treatment_response.csv",
    combo_processed = PROCDATA / "combo_treatment_response.csv",
    treatment_syn_rds = METADATA / "ZIPsynergyTCL38.RDS"
  output:
    treatment_metadata = PROCDATA / "treatment-metadata.csv"
  log:
    "logs/make_treatmentMetadata.log"
  conda:
    "workflow/envs/annotationgx.yaml"
  threads:
    4
  script:
    "workflow/scripts/make_treatmentMetadata.R"

rule make_sampleMetadata:
  input:
    sample_metadata = METADATA / "sample-metadata.xlsx"
  output:
    processed_sample_metadata = PROCDATA / "sample-metadata.csv"
  log:
    "logs/make_sampleMetadata.log"
  conda:
    "workflow/envs/annotationgx.yaml"
  threads:
    4
  script:
    "workflow/scripts/make_sampleMetadata.R"

rule make_TRE:
  input:
    mono_processed = PROCDATA / "mono_treatment_response.csv",
    combo_processed = PROCDATA / "combo_treatment_response.csv",
  output:
    tre = PROCDATA / "TCL38_TRE.RDS"
  log:
    "logs/make_TRE.log"
  script:
    "workflow/scripts/make_TRE.R"

rule process_combotherapy:
  input:
    raw = RAWDATA / "ViabilityCombinationsTCL38.zip"
  output:
    combo_processed = PROCDATA / "combo_treatment_response.csv"
  log:
    "logs/process_combotherapy.log"
  script:
    "workflow/scripts/process_combotherapy.R"

rule process_monotherapy:
  input:
    raw = RAWDATA / "2024-10-09-drugscreening" / "DSS" / "DSS_all_corrected_TC38.xlsx",
  output:
    mono_processed = PROCDATA / "mono_treatment_response.csv"
  log:
    "logs/process_monotherapy.log"
  script:
    "workflow/scripts/process_monotherapy.R"

rule get_MAE:
  input:
    samplemetadata= METADATA / "sample-metadata.xlsx",
    variantcall= RAWDATA / "RDSobjects" / "VariantCallsTCL.RDS",
    rnaseq= RAWDATA / "RDSobjects" / "TCL_RNA.RDS",
    atac= RAWDATA / "RDSobjects" / "atac2_transformed_rnagene_peaks_all_TSS.RDS",
    acgh= RAWDATA / "RDSobjects" / "aCGH.RDS",
  output:
    mae= PROCDATA / "TCL38_MAE.RDS"
  log:
    "logs/get_MAE.log"
  script:
    "workflow/scripts/get_MAE.R"