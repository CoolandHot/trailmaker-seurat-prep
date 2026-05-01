# Note: Seurat v3 and v4 share the exact same 'Assay' class structure.
# Setting assay.version = 'v3' instructs Seurat v5 to use the v3/v4 'Assay' format instead of the newer 'Assay5' format.
options(Seurat.object.assay.version = 'v3')

library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript prepare_trailmaker.R <input_tsv> <output_rds> <sample_name>")
}

input_file <- args[1]
output_file <- args[2]
sample_name <- args[3]

cat("Reading input file:", input_file, "\n")
# Read TSV file using data.table for speed
raw_data <- data.table::fread(input_file, data.table = FALSE)

# Assume the first column contains gene names
rownames(raw_data) <- raw_data[, 1]
raw_data <- raw_data[, -1]

# Convert to dgCMatrix so Seurat creates a v3 Assay instead of v5 Assay
raw_data <- as(as.matrix(raw_data), "dgCMatrix")

cat("Creating Seurat object...\n")
# Create Seurat object. Trailmaker requires scdata[['RNA']]@counts.
# We place the processed TPM data in the counts slot as well to fulfill requirements.
scdata <- CreateSeuratObject(counts = raw_data, assay = "RNA")

cat("=== Dataset Information ===\n")
print(scdata)
cat("Preview of Feature Names (Genes):\n")
print(head(rownames(scdata)))
cat("Preview of Cell Names (Barcodes):\n")
print(head(colnames(scdata)))
cat("===========================\n")

# Trailmaker requirement: scdata$samples assignment
scdata$samples <- sample_name

cat("Running normalisation, feature selection, and scaling...\n")
# Although it's TPM, we apply NormalizeData (log1p) to ensure the data slot is populated
# and distributions are standard for downstream Seurat processing
scdata <- NormalizeData(scdata)
scdata <- FindVariableFeatures(scdata, selection.method = "vst", nfeatures = 2000)
scdata <- ScaleData(scdata)

cat("Running PCA...\n")
# Trailmaker requirement: scdata@reductions contains the embeddings for pca
scdata <- RunPCA(scdata, features = VariableFeatures(object = scdata), verbose = FALSE)

cat("Running UMAP...\n")
# Trailmaker requirement: scdata@reductions contains either umap or tsne
# default dimensionality reduction must be named exactly umap or tsne
scdata <- RunUMAP(scdata, dims = 1:20, reduction.name = "umap", n.components = 3L, verbose = FALSE)

cat("Spoofing Seurat object version to v4.4.0 for Trailmaker compatibility...\n")
scdata@version <- package_version("4.4.0")

cat("Saving Seurat object to:", output_file, "\n")
saveRDS(scdata, file = output_file)
cat("Done!\n")
