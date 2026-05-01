library(Seurat)

# Ensure Seurat uses v3 assays by default for backward compatibility with Trailmaker
options(Seurat.object.assay.version = "v3")

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

# Convert to v4 format immediately after creation to ensure Trailmaker compatibility
scdata@assays$RNA <- as(object = scdata@assays$RNA, Class = "Assay")

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
# Ensure v4 compatibility for PCA reduction
scdata[["pca"]] <- as(object = scdata[["pca"]], Class = "DimReduc")

cat("Running UMAP...\n")
# Trailmaker requirement: scdata@reductions contains either umap or tsne
# default dimensionality reduction must be named exactly umap or tsne
scdata <- RunUMAP(scdata, dims = 1:20, reduction.name = "umap", n.components = 3L, verbose = FALSE)
# Ensure v4 compatibility for UMAP reduction
scdata[["umap"]] <- as(object = scdata[["umap"]], Class = "DimReduc")

# Remove all other assays and reductions (keep only RNA, pca, umap/tsne)
scdata@assays <- scdata@assays["RNA"]
scdata@reductions <- scdata@reductions[intersect(names(scdata@reductions), c("pca", "umap", "tsne"))]

cat("Saving Seurat object to:", output_file, "\n")
saveRDS(scdata, file = output_file)
cat("Done!\n")
