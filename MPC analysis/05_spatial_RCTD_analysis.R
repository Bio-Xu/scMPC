# ============================================================
# Spatial transcriptomic deconvolution using RCTD
# ============================================================

library(Seurat)
library(spacexr)
library(Matrix)
library(dplyr)

# ----------------------------
# User settings
# ----------------------------
spatial_data_dir <- "path/to/10x_spatial_sample"
sc_reference_rds <- "path/to/scRNA_reference.rds"
output_dir <- "results/spatial_RCTD"

# Metadata column containing reference cell-type labels.
reference_label_col <- "cellchat"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1. Load and preprocess the spatial dataset
# ----------------------------
stRNA <- Load10X_Spatial( data.dir = spatial_data_dir, slice = "image")
stRNA <- SCTransform( stRNA, assay = "Spatial", verbose = FALSE)
stRNA <- RunPCA( stRNA, assay = "SCT",  verbose = FALSE)
stRNA <- FindNeighbors( stRNA, reduction = "pca", dims = 1:15)
stRNA <- FindClusters( stRNA,verbose = FALSE)
stRNA <- RunUMAP( stRNA,  reduction = "pca", dims = 1:15)

# ----------------------------
# 2. Build the SpatialRNA object
# ----------------------------
st_count <- GetAssayData( stRNA, assay = "Spatial", layer = "counts")
coords <- GetTissueCoordinates(stRNA)
st_loc <- coords[, 1:2, drop = FALSE]

# Ensure the spatial counts and coordinates contain the same spots.
common_spots <- intersect(colnames(st_count), rownames(st_loc))
st_count <- st_count[, common_spots, drop = FALSE]
st_loc <- st_loc[common_spots, , drop = FALSE]
puck <- SpatialRNA( st_loc, st_count)

# ----------------------------
# 3. Build the single-cell reference
# ----------------------------
scRNA <- readRDS(sc_reference_rds)
sc_count <- GetAssayData( scRNA, assay = "RNA", layer = "counts")
cell_types <- as.factor(scRNA@meta.data[[reference_label_col]])
names(cell_types) <- colnames(scRNA)
reference <- Reference(sc_count,cell_types)
saveRDS( reference,file.path(output_dir, "RCTD_reference.rds"))

# ----------------------------
# 4. RCTD deconvolution
# ----------------------------
myRCTD <- create.RCTD( puck, reference, max_cores = 4)
myRCTD <- run.RCTD( myRCTD,  doublet_mode = "doublet")
saveRDS( myRCTD,  file.path(output_dir, "RCTD_result.rds"))

# ----------------------------
# 5. Export deconvolution results
# ----------------------------
stRNA <- AddMetaData( stRNA,  metadata = myRCTD@results$results_df)
weights <- as.data.frame(myRCTD@results$weights)

# normalize_weights() follows the original RCTD analysis workflow.
norm_weights <- normalize_weights(weights)

write.csv( norm_weights, file.path(output_dir, "RCTD_deconvolution_weights.csv"), row.names = TRUE)

# Add normalized cell-type weights as a Seurat assay for downstream
# spatial visualization or spatial co-occurrence analyses.
norm_weights_matrix <- as.matrix(t(norm_weights))
stRNA[["norm_weights"]] <- CreateAssayObject(data = norm_weights_matrix)

saveRDS( stRNA, file.path(output_dir, "spatial_seurat_with_RCTD.rds"))

