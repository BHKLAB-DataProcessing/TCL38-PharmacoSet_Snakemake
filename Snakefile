from pathlib import Path

RAWDATA = Path("rawdata")
PROCDATA = Path("procdata")
METADATA =  Path("metadata")
RESULTS = Path("results")

configfile: "config/tcl38.yaml"
annotationGx_conda = "workflow/envs/annotationgx.yaml"
pharmacoGx_conda = "workflow/envs/pharmacoset.yaml"
treatmentResponse_conda = "workflow/envs/treatmentResponse.yaml"
r_essential_conda = "workflow/envs/r-essential.yaml"

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
  conda:
    pharmacoGx_conda
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
    annotationGx_conda
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
    annotationGx_conda
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
  conda:
    treatmentResponse_conda
  script:
    "workflow/scripts/make_TRE.R"

rule process_combotherapy:
  input:
    samplemetadata= METADATA / "sample-metadata.xlsx",
    raw = RAWDATA / "ViabilityCombinationsTCL38.zip"
  output:
    combo_processed = PROCDATA / "combo_treatment_response.csv"
  log:
    "logs/process_combotherapy.log"
  script:
    "workflow/scripts/process_combotherapy.R"

rule process_monotherapy:
  input:
    raw = RAWDATA / "DSS_all_corrected_TC38.xlsx",
  output:
    mono_processed = PROCDATA / "mono_treatment_response.csv"
  log:
    "logs/process_monotherapy.log"
  conda:
    r_essential_conda
  script:
    "workflow/scripts/process_monotherapy.R"

rule get_MAE:
  input:
    samplemetadata= METADATA / "sample-metadata.xlsx",
    variantcall= RAWDATA / "VariantCallsTCL.RDS",
    rnaseq= RAWDATA / "TCL_RNA.RDS",
    atac= RAWDATA / "atac2_transformed_rnagene_peaks_all_TSS.RDS",
    acgh= RAWDATA / "aCGH.RDS",
  output:
    mae= PROCDATA / "TCL38_MAE.RDS"
  log:
    "logs/get_MAE.log"
  conda:
    treatmentResponse_conda
  script:
    "workflow/scripts/get_MAE.R"

rule download_treatment_data:
  output:
    mono = RAWDATA / "DSS_all_corrected_TC38.xlsx",
    combo = RAWDATA / "ViabilityCombinationsTCL38.zip"
  params:
    mono = config["molecularData"]["treatmentResponse"][1]["data"][1]["url"],
    combo = config["molecularData"]["treatmentResponse"][1]["data"][2]["url"]
  run:
    import asyncio
    import aiohttp

    async def download_file(session, url, dest):
        async with session.get(url) as response:
            with open(dest, 'wb') as f:
                while True:
                    chunk = await response.content.read(1024*64)
                    if not chunk:
                        break
                    f.write(chunk)

    async def download_all():
        async with aiohttp.ClientSession() as session:
            tasks = [
                download_file(session, params.mono, output.mono),
                download_file(session, params.combo, output.combo)
            ]
            await asyncio.gather(*tasks)

    asyncio.run(download_all())

rule download_molecular_profiles:
  output:
    variantcall = RAWDATA / "VariantCallsTCL.RDS",
    rnaseq = RAWDATA / "TCL_RNA.RDS",
    atacseq = RAWDATA / "atac2_transformed_rnagene_peaks_all_TSS.RDS",
    acgh = RAWDATA / "aCGH.RDS"
  params:
    variantcall = config["molecularData"]["rna"]["variantcall"][1]["url"],
    rnaseq = config["molecularData"]["rna"]["rnaseq"][1]["url"],
    atacseq = config["molecularData"]["dna"]["atacseq"][1]["url"],
    acgh = config["molecularData"]["dna"]["acgh"][1]["url"]
  run:
    import asyncio
    import aiohttp

    async def download_file(session, url, dest):
        async with session.get(url) as response:
            with open(dest, 'wb') as f:
                while True:
                    chunk = await response.content.read(1024*64)
                    if not chunk:
                        break
                    f.write(chunk)

    async def download_all():
        async with aiohttp.ClientSession() as session:
            tasks = [
                download_file(session, params.variantcall, output.variantcall),
                download_file(session, params.rnaseq, output.rnaseq),
                download_file(session, params.atacseq, output.atacseq),
                download_file(session, params.acgh, output.acgh)
            ]
            await asyncio.gather(*tasks)

    asyncio.run(download_all())

rule download_metadata:
  output:
    samplemetadata = METADATA / "sample-metadata.xlsx",
    treatment_syn_rds = METADATA / "ZIPsynergyTCL38.RDS"
  params:
    samplemetadata = config["metadata"]["sampleMetadata"][1]["url"],
    treatment_syn_rds = config["metadata"]["treatmentMetadata"][1]["url"]
  run:
    import asyncio
    import aiohttp

    async def download_metadata():
      async with aiohttp.ClientSession() as session:
        async with session.get(params.samplemetadata) as response:
          with open(output.samplemetadata, "wb") as f:
            f.write(await response.read())
        async with session.get(params.treatment_syn_rds) as response:
          with open(output.treatment_syn_rds, "wb") as f:
            f.write(await response.read())

    asyncio.run(download_metadata())



rule download_all:
  input:
    mono = RAWDATA / "DSS_all_corrected_TC38.xlsx",
    combo = RAWDATA / "ViabilityCombinationsTCL38.zip",
    variantcall = RAWDATA / "VariantCallsTCL.RDS",
    rnaseq = RAWDATA / "TCL_RNA.RDS",
    atacseq = RAWDATA / "atac2_transformed_rnagene_peaks_all_TSS.RDS",
    acgh = RAWDATA / "TCL_ACGH.RDS",
    samplemetadata = METADATA / "sample-metadata.xlsx",
    treatment_syn_rds = METADATA / "ZIPsynergyTCL38.RDS"
  shell:
    "echo 'Downloaded all data'"


rule clean:
  shell:
    "rm -rf procdata results rawdata metadata"