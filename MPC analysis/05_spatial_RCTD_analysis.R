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

# The scRNA-seq reference originates from GSE161529. The "cellchat" metadata
# column contains the final cell-type annotations obtained by integrating
# cell-type annotation with MPC subpopulation identification in this study.
# MPC and NPC are retained as separate labels in this column and are supplied
# to RCTD as distinct reference cell types. This enables mapping the study-
# defined cell populations to spatial transcriptomic spots.
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

# ----------------------------
# 6. Spatial co-localization analysis using mistyR
# ----------------------------
library(mistyR)

mistyr_output_dir <- file.path(output_dir, "mistyR_colocalization")
mistyr_juxta_threshold <- 5
mistyr_para_scale <- 10
dir.create(mistyr_output_dir, recursive = TRUE, showWarnings = FALSE)

# Use normalized RCTD cell-type weights as spatial features.
features <- as.data.frame(norm_weights)
common_spots <- Reduce(intersect, list(
  rownames(features), rownames(st_loc), colnames(stRNA)
))
if (length(common_spots) < 10) {
  stop("Insufficient matching spots for mistyR analysis.")
}
features <- features[common_spots, , drop = FALSE]
coords <- as.data.frame(st_loc[common_spots, , drop = FALSE])
coords[] <- lapply(coords, as.numeric)
colnames(coords) <- c("x", "y")
if (anyNA(coords)) stop("Spatial coordinates contain missing values.")

# Remove invariant or nonnumeric features before model fitting.
features[] <- lapply(features, as.numeric)
keep <- vapply(features, function(x) {
  all(is.finite(x)) && stats::sd(x) > 0
}, logical(1))
features <- features[, keep, drop = FALSE]
if (ncol(features) < 2) {
  stop("At least two variable cell-type features are required.")
}

set.seed(1234)

# Intra: same-spot features; juxta: immediate spatial neighborhood;
# para: broader spatial neighborhood.
misty_views <- mistyR::create_initial_view(features, coords)
misty_views <- mistyR::add_juxtaview(
  misty_views, neighbor.thr = mistyr_juxta_threshold
)
misty_views <- mistyR::add_paraview(
  misty_views, l = mistyr_para_scale
)
mistyR::run_misty(misty_views, results.folder = mistyr_output_dir)

write.csv(features, file.path(mistyr_output_dir, "input_cell_type_weights.csv"))
write.csv(coords, file.path(mistyr_output_dir, "input_spatial_coordinates.csv"))
writeLines(capture.output(sessionInfo()),
           file.path(mistyr_output_dir, "R_sessionInfo.txt"))
