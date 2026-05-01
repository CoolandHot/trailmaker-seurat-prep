# Trailmaker Seurat Data Preparation Pipeline

This repository contains a pipeline to prepare single-cell RNA-seq datasets (specifically Smart-seq2 GBM IDHwt data) into a Seurat object `.rds` file that is fully compatible with the Parse Biosciences Trailmaker platform.

## Scripts

- `run_prepare_trailmaker.sh`: A shell wrapper that automatically decompresses the input dataset and invokes the R script.
- `prepare_trailmaker.R`: The core R script that constructs the Seurat object, ensures proper Assay class formatting, performs dimensionality reduction (PCA & UMAP), and spoofs the object version to bypass platform validation constraints.

## Trailmaker Seurat Object Requirements

This pipeline explicitly adheres to the Trailmaker upload specifications. For a Seurat `.rds` object to be successfully uploaded, it must meet the following criteria:

- **Format:** Seurat objects must be in the `.rds` format and **must be strictly in the Seurat v4 format**. (v5 formats will fail validation. This script natively handles v5-to-v4 downgrading/spoofing).
- **Size Limit:** There is a maximum size limit of **15GB**. (If over 15GB, remove non-essential assays).
- **Default Reduction:** The default dimensionality reduction must be named exactly `umap` or `tsne`. (If it includes umap/tsne like `ref.umap`, it is automatically renamed. If it doesn't contain these names, the upload will fail).
- **Slots and Metadata required:**
  - `scdata$samples`: Sample assignment. If absent, the data is treated as a single-sample experiment.
  - `scdata[['RNA']]@counts`: Must contain raw feature counts.
  - `scdata@reductions`: Must contain the embeddings for PCA, as well as either UMAP or tSNE.
- **Auto-detection:**
  - Cluster metadata located in `scdata@meta.data` is auto-detected.
  - Sample-level metadata in `scdata@meta.data` that groups samples in `scdata$samples` is auto-detected for downstream analysis.

## Usage

### Running via R script directly
You can use the standalone R script. It requires three arguments: the input TSV file, the output RDS file path, and the sample name.

```bash
Rscript prepare_trailmaker.R <input_tsv> <output_rds> <sample_name>
```

### Running via Wrapper Script
Alternatively, use the shell script to automatically unzip the TSV and invoke the R script.

```bash
chmod +x run_prepare_trailmaker.sh
./run_prepare_trailmaker.sh
```

### Official Documentation
For more detailed information, please refer to the official [Trailmaker User Guide](https://support.parsebiosciences.com/hc/en-us/articles/27076682137236-Trailmaker-User-Guide#h_01HZ4VDNQ8HQNH0BC12VX54PJR).
