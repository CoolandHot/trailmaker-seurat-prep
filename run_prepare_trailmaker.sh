#!/bin/bash
# Exit immediately if a command exits with a non-zero status
set -e

# File paths
WORKSPACE_DIR="/scratch/hh01116/trailmaker_scr"
GZ_FILE="$WORKSPACE_DIR/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv.gz"
UNZIPPED_FILE="$WORKSPACE_DIR/GSM3828672_Smartseq2_GBM_IDHwt_processed_TPM.tsv"
RDS_FILE="$WORKSPACE_DIR/GSM3828672_Smartseq2_GBM_IDHwt_Trailmaker.rds"
R_SCRIPT="$WORKSPACE_DIR/prepare_trailmaker.R"
SAMPLE_NAME="GSM3828672_GBM_IDHwt"

# Path to Rscript in the micromamba gene_r environment
RSCRIPT_BIN="/vol/research/brainTumorST/mamba_baseEnv/envs/gene_r/bin/Rscript"

echo "Checking if unzipped TSV exists..."
if [ ! -f "$UNZIPPED_FILE" ]; then
    echo "Unzipping data file..."
    # -d for decompress, -k to keep the original .gz file
    gzip -dk "$GZ_FILE"
else
    echo "Unzipped file already exists, skipping unzip."
fi

echo "Running R script to generate Seurat object..."
"$RSCRIPT_BIN" "$R_SCRIPT" "$UNZIPPED_FILE" "$RDS_FILE" "$SAMPLE_NAME"

echo "Pipeline finished successfully!"
echo "Your Trailmaker-ready .rds file is located at:"
echo "$RDS_FILE"
