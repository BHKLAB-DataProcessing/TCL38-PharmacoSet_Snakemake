from pathlib import Path

RAWDATA = Path("rawdata")
PROCDATA = Path("procdata")
METADATA =  Path("metadata")

rule make_PSET:
  input:
    tre = PROCDATA / "TCL38_TRE.RDS",
    mae = PROCDATA / "TCL38_MAE.RDS",
  output:
    pset = PROCDATA / "TCL38_pset.RDS"
  log:
    "logs/make_PSET.log"
  script:
    "workflow/scripts/make_PSET.R"

rule make_TRE:
  input:
    mono_processed = PROCDATA / "mono_treatment_response.csv",
    combo_processed = PROCDATA / "combo_treatment_response.csv",
    # sample_metadata = METADATA / "sample-metadata.xlsx",
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