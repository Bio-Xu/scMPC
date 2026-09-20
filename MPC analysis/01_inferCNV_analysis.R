# ============================================================
# inferCNV analysis after MPC identification
# ============================================================

library(Seurat)
library(infercnv)
library(dplyr)

# ----------------------------
# User settings
# ----------------------------
seurat_rds <- "path/to/seurat_object.rds"
gene_order_file <- "path/to/hg38_gencode_v27.txt"
output_dir <- "results/inferCNV"
patient_id <- "GSE195_5"

patient_col <- "PatientID"
celltype_col <- "Tissue_Cell_Type"

# Set the malignant epithelial label used in the corresponding dataset.
malignant_label <- "Pri-MC_Epi"

# Immune-cell populations used as normal references in the original analysis.
reference_cell_types <- c(
  "class switched memory B cell",
  "effector memory CD8-positive, alpha-beta T cell",
  "mature NK T cell",
  "IgG plasma cell",
  "activated CD4-positive, alpha-beta T cell",
  "regulatory T cell",
  "IgA plasma cell",
  "natural killer cell",
  "macrophage"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1. Prepare patient-level data
# ----------------------------
seurat_obj <- readRDS(seurat_rds)
meta <- seurat_obj@meta.data
patient_cells <- rownames(meta)[meta[[patient_col]] == patient_id]
patient <- subset(seurat_obj, cells = patient_cells)
cell_counts <- table(patient@meta.data[[celltype_col]])
reference_groups <- intersect(names(cell_counts[cell_counts > 5]),reference_cell_types)
keep_groups <- c(reference_groups, malignant_label)
keep_cells <- rownames(patient@meta.data)[patient@meta.data[[celltype_col]] %in% keep_groups]
patient <- subset(patient, cells = keep_cells)

# ----------------------------
# 2. Create inferCNV object
# ----------------------------
raw_counts <- GetAssayData(patient, assay = "RNA", layer = "counts")
annotations <- data.frame(
  cell_type = patient@meta.data[[celltype_col]],
  row.names = rownames(patient@meta.data)
)

gene_order <- read.table(gene_order_file,header = FALSE, row.names = 1,sep = "\t",check.names = FALSE)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = raw_counts,
  annotations_file = annotations,
  gene_order_file = gene_order,
  min_max_counts_per_cell = c(100, Inf),
  ref_group_names = reference_groups
)

saveRDS(infercnv_obj, file.path(output_dir, "infercnv_input_object.rds"))

# ----------------------------
# 3. Run inferCNV
# ----------------------------
options(scipen = 100)

infercnv_result <- infercnv::run(
  infercnv_obj = infercnv_obj,
  out_dir = output_dir,
  cutoff = 0.1,
  analysis_mode = "subclusters",
  tumor_subcluster_partition_method = "random_trees",
  HMM = TRUE,
  HMM_type = "i6",
  HMM_report_by = "subcluster",
  denoise = TRUE,
  num_threads = 4,
  plot_probabilities = FALSE
)

saveRDS(infercnv_result, file.path(output_dir, "infercnv_final_object.rds"))

# ----------------------------
# 4. Export tumor-subcluster membership
# This table can be used for downstream clonal-tree annotation.
# ----------------------------
if (malignant_label %in% names(infercnv_result@tumor_subclusters$subclusters)) {
  malignant_groups <- infercnv_result@tumor_subclusters$subclusters[[malignant_label]]

  group_table <- do.call(
    rbind,
    lapply(names(malignant_groups), function(group_name) {
      idx <- malignant_groups[[group_name]]
      data.frame(
        cell_group_name = group_name,
        cell_index = idx,
        cell_id = colnames(infercnv_result@expr.data)[idx],
        stringsAsFactors = FALSE
      )
    })
  )

  write.table(group_table,file = file.path(output_dir, "tumor_subcluster_membership.tsv"),sep = "\t",row.names = FALSE,quote = FALSE )

  # Uphyloplot2 uses a two-column *.cell_groupings file.
  write.table( group_table[, c("cell_group_name", "cell_index")],file = file.path(output_dir, "tree.cell_groupings"),sep = "\t", row.names = FALSE,quote = FALSE)
}
