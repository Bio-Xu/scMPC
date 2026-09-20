# ============================================================
# Monocle2 trajectory analysis after MPC identification
# The trajectory root is selected using inferCNV-derived CNV burden.
# ============================================================

library(Seurat)
library(monocle)
library(dplyr)

# ----------------------------
# User settings
# ----------------------------
seurat_rds <- "path/to/seurat_object.rds"
infercnv_rds <- "path/to/infercnv_final_object.rds"
gene_order_file <- "path/to/hg38_gencode_v27.txt"
output_dir <- "results/Monocle2"
patient_id <- "GSE195_2"
patient_col <- "PatientID"
celltype_col <- "Tissue_Cell_Type"

# Epithelial populations used in the original trajectory analysis.
epithelial_types <- c("Met-MC_Epi", "Nor_Epi", "Pri-MC_Epi")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helper: inferCNV-based CNV burden
# ----------------------------
CNVburden <- function(cnvscore, gene_info) {
  step1 <- abs(cnvscore - 1)
  gene_info <- gene_info[rownames(gene_info) %in% rownames(step1), , drop = FALSE]
  step1 <- step1[rownames(gene_info), , drop = FALSE]
  gene_length <- gene_info[, 4] - gene_info[, 3]
  burden <- vapply(colnames(step1),function(cell_id) {sum(step1[, cell_id] * gene_length) / sum(gene_length)},numeric(1) )
  data.frame(cnv_burden = burden,row.names = names(burden))
}

# ----------------------------
# 1. Prepare patient epithelial cells
# ----------------------------
seurat_obj <- readRDS(seurat_rds)
meta <- seurat_obj@meta.data
keep_cells <- rownames(meta)[meta[[patient_col]] == patient_id & meta[[celltype_col]] %in% epithelial_types]
seurat_sub <- subset(seurat_obj, cells = keep_cells)
seurat_sub <- FindVariableFeatures(seurat_sub, selection.method = "vst", nfeatures = 2000)

# ----------------------------
# 2. Build the Monocle2 CellDataSet
# ----------------------------
count_mat <- GetAssayData(seurat_sub, assay = "RNA", layer = "counts")
gene_anno <- data.frame( gene_short_name = rownames(count_mat), row.names = rownames(count_mat))
pd <- new("AnnotatedDataFrame", data = seurat_sub@meta.data)
fd <- new("AnnotatedDataFrame", data = gene_anno)
cds <- newCellDataSet( count_mat, phenoData = pd, featureData = fd)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
ordering_genes <- VariableFeatures(seurat_sub, nfeatures = 2000)
ordering_genes <- intersect(ordering_genes, rownames(cds))
cds <- setOrderingFilter(cds, ordering_genes)

# ----------------------------
# 3. DDRTree trajectory inference
# ----------------------------
cds <- reduceDimension( cds, max_components = 2, reduction_method = "DDRTree")
cds <- orderCells(cds)

# ----------------------------
# 4. Select the root state using CNV burden
# ----------------------------
infercnv_obj <- readRDS(infercnv_rds)
common_cells <- intersect( colnames(infercnv_obj@expr.data), colnames(exprs(cds)))
cnvscore <- infercnv_obj@expr.data[, common_cells, drop = FALSE]
gene_info <- read.table( gene_order_file, header = FALSE, row.names = 1, sep = "\t",  check.names = FALSE)
cnv_burden <- CNVburden(cnvscore, gene_info)
cnv_burden$CellID <- rownames(cnv_burden)
state_info <- pData(cds)
state_info$CellID <- rownames(state_info)
state_info <- merge(state_info,cnv_burden, by = "CellID")
state_mean <- aggregate( cnv_burden ~ State, data = state_info, FUN = mean)
root_state <- state_mean$State[which.min(state_mean$cnv_burden)]
cds <- orderCells(cds, root_state = root_state)

# ----------------------------
# 5. Save trajectory results
# ----------------------------
saveRDS(cds, file.path(output_dir, paste0(patient_id, "_monocle2_cds.rds")))

trajectory_meta <- pData(cds)
trajectory_meta$CellID <- rownames(trajectory_meta)

write.table( trajectory_meta, file = file.path(output_dir, paste0(patient_id, "_trajectory_metadata.tsv")),sep = "\t", row.names = FALSE, quote = FALSE)
