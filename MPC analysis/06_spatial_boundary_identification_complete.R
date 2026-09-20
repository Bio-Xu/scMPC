# ============================================================
# 06. Spatial tumor-boundary identification
# ============================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(Cottrazm)
  library(infercnv)
  library(reticulate)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(ggpubr)
  library(readr)
  library(dendextend)
  library(ape)
})
set.seed(1234)
options(scipen = 100)

# -------------------- User settings --------------------
visium_dir <- "path/to/Visium_sample"                       # contains filtered_feature_bc_matrix(.h5) and spatial/
gene_order_file <- "path/to/gencode_v38_gene_pos.txt"      # inferCNV gene-position file
tissue_positions_file <- file.path(visium_dir, "spatial", "tissue_positions_list.csv")
output_dir <- "results/spatial_boundary"
sample_name <- "BRCA"
conda_env <- "TumorBoundary"                                # Cottrazm morphology Python environment
cluster_resolution <- 1.5
infercnv_assay <- "Spatial"
infercnv_threads <- 24
cnv_k <- 8
malignant_cnv_labels <- NULL                               # e.g. c("1","5","6","7","8"); NULL = top 2 median CNV-score groups
max_iterations <- 6
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
infercnv_dir <- file.path(output_dir, "InferCNV")
dir.create(infercnv_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------- Helpers --------------------
get_layer <- function(x, assay, layer = "counts") {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") SeuratObject::LayerData(x, assay = assay, layer = layer)
  else Seurat::GetAssayData(x, assay = assay, slot = layer)
}
prepare_spatial_neighbors <- function(position_file) {
  position <- read.csv(position_file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(position) != 6) stop("tissue_positions_list.csv should contain 6 columns.")
  position <- position[, -2, drop = FALSE]
  colnames(position) <- c("spot_id","row","col","imagerow","imagecol")
  position$spot.ids <- seq_len(nrow(position))
  dists <- compute_interspot_distances(position = position, scale.factor = 1.05)
  df_j <- find_neighbors(position = position, radius = dists$radius, method = "manhattan")
  names(df_j) <- position$spot_id
  lookup_vec <- setNames(position$spot_id, position$spot.ids)
  df_j <- lapply(df_j, function(x) {
    y <- lookup_vec[as.character(x)]
    y[is.na(y)] <- as.character(x[is.na(y)])
    unname(y)
  })
  list(position = position, neighbors = df_j)
}
get_normal_cluster <- function(TumorST) {
  x <- tapply(TumorST$NormalScore, TumorST$seurat_clusters, mean, na.rm = TRUE)
  names(sort(x, decreasing = TRUE))[1]
}

# -------------------- 1. Read Visium data --------------------
STPreProcess <- function(InDir, Sample, OutDir) {
  spatial_path <- file.path(InDir, "spatial")
  image <- Read10X_Image(image.dir = spatial_path)
  DefaultAssay(image) <- "Spatial"
  h5 <- file.path(InDir, "filtered_feature_bc_matrix.h5")
  if (file.exists(h5)) {
    TumorST <- Load10X_Spatial(data.dir = InDir, filename = basename(h5), assay = "Spatial", slice = "image", image = image)
  } else {
    matrix_dir <- file.path(InDir, "filtered_feature_bc_matrix")
    if (!dir.exists(matrix_dir))
      stop("Cannot find filtered_feature_bc_matrix.h5 or filtered_feature_bc_matrix/.")
    counts <- Read10X(data.dir = matrix_dir)
    TumorST <- CreateSeuratObject(counts = counts, assay = "Spatial", project = Sample)
    image <- image[colnames(TumorST)]
    TumorST[["image"]] <- image
  }
  if ("image" %in% names(TumorST@images)) {
    sf <- TumorST@images$image@scale.factors
    if (!is.null(sf$hires)) TumorST@images$image@scale.factors$lowres <- sf$hires
  }
  TumorST[["Mito.percent"]] <- PercentageFeatureSet(TumorST, pattern = "^MT-")
  qc_dir <- file.path(OutDir, "QC")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  pdf(file.path(qc_dir, "Vlnplot.pdf"), width = 7, height = 4)
  print(VlnPlot(TumorST, features = c("nFeature_Spatial","nCount_Spatial","Mito.percent"), pt.size = 0, ncol = 3))
  dev.off()

  pdf(file.path(qc_dir, "SpatialFeaturePlot.pdf"), width = 12, height = 4)
  print(SpatialFeaturePlot(TumorST, features = c("nFeature_Spatial","nCount_Spatial","Mito.percent"), ncol = 3))
  dev.off()
  write.csv(TumorST@meta.data[, c("nCount_Spatial","nFeature_Spatial","Mito.percent")],
            file.path(qc_dir, "QCData.csv"), row.names = TRUE)
  TumorST
}

# -------------------- 2. Morphology-adjusted clustering + NormalScore --------------------
STModiCluster <- function(InDir, Sample, OutDir, TumorST, res = 1.5, conda_env = "TumorBoundary") {
  if (!is.null(conda_env) && nzchar(conda_env)) reticulate::use_condaenv(conda_env, required = TRUE)
  py_script <- system.file("python/Rusedtile.py", package = "Cottrazm")
  if (!nzchar(py_script)) stop("Cottrazm python/Rusedtile.py was not found.")
  reticulate::source_python(py_script)
  Adjusted_expr_mtx <- ME_normalize(inDir = InDir, outDir = paste0(normalizePath(OutDir, mustWork = FALSE), "/"), sample = Sample)
  spatial_counts <- get_layer(TumorST, "Spatial", "counts")
  expected_spots <- colnames(spatial_counts)
  expected_genes <- rownames(spatial_counts)
  if (is.null(dim(Adjusted_expr_mtx))) stop("ME_normalize did not return a matrix.")
  if (nrow(Adjusted_expr_mtx) == length(expected_spots) && ncol(Adjusted_expr_mtx) == length(expected_genes)) {
    rownames(Adjusted_expr_mtx) <- expected_spots
    colnames(Adjusted_expr_mtx) <- expected_genes
  } else {
    mtx_file <- file.path(OutDir, paste0(Sample, "_raw_SME_normalizeA.mtx"))
    if (!file.exists(mtx_file)) stop("Unexpected ME_normalize dimensions and fallback matrix was not found: ", mtx_file)
    Adjusted_expr_mtx <- Matrix::readMM(mtx_file)
    rownames(Adjusted_expr_mtx) <- expected_spots
    colnames(Adjusted_expr_mtx) <- expected_genes
  }
  morph_counts <- Matrix::Matrix(t(Adjusted_expr_mtx), sparse = TRUE)
  morph_obj <- CreateSeuratObject(counts = morph_counts, assay = "RNA")
  morph_obj <- subset(morph_obj, cells = rownames(TumorST@meta.data))
  TumorST[["Morph"]] <- morph_obj[["RNA"]]
  TumorST <- NormalizeData(TumorST, assay = "Morph", verbose = FALSE)
  TumorST <- FindVariableFeatures(TumorST, mean.function = ExpMean, dispersion.function = LogVMR, assay = "Morph", verbose = FALSE)
  TumorST <- ScaleData(TumorST, assay = "Morph", verbose = FALSE)
  TumorST <- RunPCA(TumorST, npcs = 50, assay = "Morph", verbose = FALSE)
  TumorST <- FindNeighbors(TumorST, reduction = "pca", dims = 1:50, assay = "Morph", verbose = FALSE)
  TumorST <- RunUMAP(TumorST, dims = 1:50, assay = "Morph", seed.use = 1234, verbose = FALSE)
  TumorST <- FindClusters(TumorST, resolution = res, algorithm = 1, graph.name = "Morph_snn", random.seed = 1234, verbose = FALSE)
  cluster_col <- paste0("Morph_snn_res.", res)
  if (!cluster_col %in% colnames(TumorST@meta.data)) cluster_col <- tail(grep("^Morph_snn_res\\.", colnames(TumorST@meta.data), value = TRUE), 1)
  TumorST$seurat_clusters <- factor(TumorST@meta.data[[cluster_col]])
  normal_features <- c("PTPRC","CD2","CD3D","CD3E","CD3G","CD5","CD7","CD79A","MS4A1","CD19")
  morph_data <- get_layer(TumorST, "Morph", "data")
  normal_features <- intersect(normal_features, rownames(morph_data))
  if (length(normal_features) == 0) stop("None of the NormalScore marker genes were found in the Morph assay.")
  TumorST$NormalScore <- Matrix::colMeans(morph_data[normal_features, , drop = FALSE])
  NormalCluster <- get_normal_cluster(TumorST)
  message("NormalCluster = ", NormalCluster)
  pdf(file.path(OutDir, paste0(Sample, "_Spatial_SeuratCluster.pdf")), width = 7, height = 7)
  print(SpatialDimPlot(TumorST, group.by = "seurat_clusters", pt.size.factor = 1, alpha = 0.8) +
        ggtitle(paste0("Resolution = ", res)))
  dev.off()
  pdf(file.path(OutDir, paste0(Sample, "_UMAP_SeuratCluster.pdf")), width = 7, height = 7)
  print(DimPlot(TumorST, group.by = "seurat_clusters") +
        ggtitle(paste0("Resolution = ", res)))
  dev.off()
  pdf(file.path(OutDir, paste0(Sample, "_NormalScore.pdf")), width = 7, height = 4)
  print(VlnPlot(TumorST, features = "NormalScore", pt.size = 0, group.by = "seurat_clusters") +
      geom_boxplot(width = 0.2, outlier.shape = NA) +
      NoLegend())
  dev.off()
  annotation <- data.frame(CellID = rownames(TumorST@meta.data), DefineTypes = as.character(TumorST$seurat_clusters))
  write.table(annotation, file.path(OutDir, "InferCNV", "CellAnnotation.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# -------------------- 3. Spatial inferCNV --------------------
STCNV <- function(TumorST, assay = "Spatial", OutDir, Sample, gene_order_file, num_threads = 24) {
  matrix <- as.matrix(get_layer(TumorST, assay, "counts"))
  annotation_file <- file.path(OutDir, "InferCNV", "CellAnnotation.txt")
  NormalCluster <- get_normal_cluster(TumorST)
  if (!file.exists(gene_order_file)) stop("gene_order_file not found: ", gene_order_file)
  infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix = matrix, annotations_file = annotation_file,
    delim = "\t", gene_order_file = gene_order_file, ref_group_names = NormalCluster)

  cnv_outdir <- file.path(OutDir, "InferCNV", paste0("output_", assay))
  infercnv_obj <- infercnv::run(infercnv_obj, cutoff = 0.1, out_dir = cnv_outdir, cluster_by_groups = FALSE,
    analysis_mode = "subclusters", denoise = TRUE, HMM = TRUE, tumor_subcluster_partition_method = "random_trees",
    HMM_type = "i6", BayesMaxPNormal = 0, num_threads = num_threads, write_expr_matrix = TRUE, write_phylo = TRUE)
}

# -------------------- 4. CNVLabel + cnv_score --------------------
STCNVScore <- function(TumorST, assay = "Spatial", OutDir, Sample, k = 8) {
  cnv_outdir <- file.path(OutDir, "InferCNV", paste0("output_", assay))
  dendro_file <- list.files(cnv_outdir, pattern = "HMM_predHMMi6.*observations_dendrogram\\.txt$", full.names = TRUE)
  obs_file <- list.files(cnv_outdir, pattern = "HMM_predHMMi6.*observations\\.txt$", full.names = TRUE)
  if (length(dendro_file) == 0) stop("inferCNV observation dendrogram file was not found in: ", cnv_outdir)
  if (length(obs_file) == 0) stop("inferCNV HMM observation file was not found in: ", cnv_outdir)
  cell_groupings <- ape::read.tree(dendro_file[1])
  infercnv_label <- dendextend::cutree(cell_groupings, k = k)
  label_df <- data.frame(CNVLabel = as.character(infercnv_label), row.names = names(infercnv_label), stringsAsFactors = FALSE)
  missing_spots <- setdiff(rownames(TumorST@meta.data), rownames(label_df))
  if (length(missing_spots) > 0) label_df <- rbind(label_df, data.frame(CNVLabel = "Normal", row.names = missing_spots))
  TumorST$CNVLabel <- label_df[rownames(TumorST@meta.data), "CNVLabel"]
  cnv_table <- read.table(obs_file[1], header = TRUE, check.names = TRUE)
  cnv_score_table <- as.matrix(cnv_table)
  storage.mode(cnv_score_table) <- "numeric"
  scores <- colSums(abs(cnv_score_table - 3), na.rm = TRUE)
  names(scores) <- gsub("\\.", "-", names(scores))
  TumorST$cnv_score <- unname(scores[rownames(TumorST@meta.data)])
  TumorST$cnv_score[TumorST$CNVLabel == "Normal"] <- 0
  pdf(file.path(OutDir, paste0(Sample, "_cnv_label.pdf")), width = 7, height = 7)
  print(SpatialDimPlot(TumorST, group.by = "CNVLabel", pt.size.factor = 1, alpha = 0.6))
  dev.off()
  pdf(file.path(OutDir, paste0(Sample, "_reduction_cnvlabel.pdf")), width = 7, height = 7)
  print(DimPlot(TumorST, group.by = "CNVLabel"))
  dev.off()
  pdf(file.path(OutDir, paste0(Sample, "_cnv_observation_vlnplot.pdf")), width = 7, height = 4)
  print(ggplot(TumorST@meta.data, aes(x = CNVLabel, y = cnv_score, fill = CNVLabel)) +
    geom_violin(alpha = 0.5) + geom_boxplot(width = 0.25, outlier.size = 0.4) + labs(y = "CNV score") +
    theme_classic() + NoLegend())
  dev.off()
  TumorST
}

# -------------------- 5. Boundary definition --------------------
BoundaryDefine <- function(TumorST, position_file, MalLabel = NULL, max_iterations = 6) {
  required <- c("seurat_clusters","NormalScore","CNVLabel","cnv_score")
  missing <- setdiff(required, colnames(TumorST@meta.data))
  if (length(missing) > 0) stop("Missing metadata: ", paste(missing, collapse = ", "))
  if (!"umap" %in% names(TumorST@reductions)) stop("UMAP reduction is required.")
  UMAPembeddings <- as.data.frame(Embeddings(TumorST, "umap")[, 1:2, drop = FALSE])
  colnames(UMAPembeddings) <- c("x","y")
  spatial_info <- prepare_spatial_neighbors(position_file)
  position <- spatial_info$position
  df_j <- spatial_info$neighbors
  if (is.null(MalLabel)) {
    med <- tapply(TumorST$cnv_score, TumorST$CNVLabel, median, na.rm = TRUE)
    med <- med[names(med) != "Normal" & is.finite(med)]
    MalLabel <- names(sort(med, decreasing = TRUE))[seq_len(min(2, length(med)))]
  }
  MalLabel <- as.character(MalLabel)
  message("Malignant CNV labels: ", paste(MalLabel, collapse = ", "))
  MalCellID <- rownames(TumorST@meta.data)[as.character(TumorST$CNVLabel) %in% MalLabel]
  NormalCluster <- get_normal_cluster(TumorST)
  NormalCellID <- rownames(TumorST@meta.data)[as.character(TumorST$seurat_clusters) == NormalCluster]
  cnv_tab <- table(as.character(TumorST$CNVLabel), as.character(TumorST$seurat_clusters))
  cluster_sizes <- table(as.character(TumorST$seurat_clusters))
  valid_labels <- intersect(MalLabel, rownames(cnv_tab))
  ClusterID <- colnames(cnv_tab)[vapply(colnames(cnv_tab), function(z) sum(cnv_tab[valid_labels, z, drop = FALSE]) > 0.5 * cluster_sizes[z], logical(1))]
  if (length(ClusterID) == 0) stop("No morphology cluster passed the >50% malignant-CNV criterion.")
  CiMal <- tibble(cluster = ClusterID) %>%
    mutate(sub_MalCellID = map(cluster, ~intersect(MalCellID, rownames(TumorST@meta.data)[as.character(TumorST$seurat_clusters) == .x]))) %>%
    filter(lengths(sub_MalCellID) > 0) %>%
    mutate(sub_CiMal = map(sub_MalCellID, ~colMeans(UMAPembeddings[.x, , drop = FALSE], na.rm = TRUE)))
  CiNormal <- colMeans(UMAPembeddings[NormalCellID, , drop = FALSE], na.rm = TRUE)
  MalCellIDsi <- map2(CiMal$sub_MalCellID, CiMal$sub_CiMal, function(ids, center) {
    ids[vapply(ids, function(id) {
      pos <- as.numeric(UMAPembeddings[id, ])
      rt <- sqrt(sum((pos - center)^2))
      rn <- sqrt(sum((pos - CiNormal)^2))
      rt < (1/3) * rn
    }, logical(1))]
  }) %>% unlist(use.names = FALSE) %>% unique()
  MalCellIDL <- vapply(MalCellIDsi, function(id) {
    nb <- df_j[[id]]
    if (is.null(nb) || sum(nb %in% MalCellIDsi) == 0) id else NA_character_
  }, character(1)) %>% na.omit() %>% as.character()
  BdyCellID <- character(0)
  nbrs_of_MalL <- nbrs(df_j = df_j, MalCellIDAdd = MalCellIDL, CellIDRaw = c(MalCellIDsi, NormalCellID, BdyCellID))
  ClusterL <- do.call(rbind, lapply(names(nbrs_of_MalL), function(cell_id) {
    cluster_id <- as.character(TumorST@meta.data[cell_id, "seurat_clusters"])
    center_index <- match(cluster_id, CiMal$cluster)
    if (is.na(center_index)) return(NULL)
    malignant_center <- CiMal$sub_CiMal[[center_index]]
    do.call(rbind, lapply(nbrs_of_MalL[[cell_id]], function(neighbor_id) {
      pos <- as.numeric(UMAPembeddings[neighbor_id, ])
      rt <- sqrt(sum((pos - malignant_center)^2))
      rn <- sqrt(sum((pos - CiNormal)^2))
      data.frame(CellID = neighbor_id, Location = ifelse(rt < (1/3) * rn, "Mal", "Bdy"), stringsAsFactors = FALSE)
    }))
  }))
  if (!is.null(ClusterL) && nrow(ClusterL) > 0) {
    vote <- as.data.frame.matrix(table(ClusterL$CellID, ClusterL$Location))
    if (!"Mal" %in% colnames(vote)) vote$Mal <- 0
    ClusterL <- data.frame(CellID = rownames(vote), Location = ifelse(vote$Mal > 0, "Mal", "Bdy"), stringsAsFactors = FALSE)
  } else ClusterL <- data.frame(CellID = character(0), Location = character(0))
  MalCellID <- unique(c(MalCellIDsi, ClusterL$CellID[ClusterL$Location == "Mal"]))
  BdyCellID <- unique(ClusterL$CellID[ClusterL$Location == "Bdy"])
  MalCellIDN <- MalCellID
  Clustern <- rbind(data.frame(CellID = NormalCellID, Location = "Normal"), data.frame(CellID = MalCellID, Location = "Mal"),
                    data.frame(CellID = BdyCellID, Location = "Bdy"))
  TumorSTn <- TumorST
  TumorSTn$LabelNew <- Clustern$Location[match(rownames(TumorSTn@meta.data), Clustern$CellID)]
  n <- 1
  repeat {
    if (length(MalCellIDN) < 3 || n > max_iterations) break
    nbrs_of_Mal <- nbrs(df_j = df_j, MalCellIDAdd = MalCellIDN, CellIDRaw = c(MalCellID, NormalCellID, BdyCellID))
    new_neighbors <- unique(unlist(nbrs_of_Mal, use.names = FALSE))
    if (length(new_neighbors) < 3) break
    cells_to_keep <- intersect(unique(c(new_neighbors, MalCellID, NormalCellID, BdyCellID)), colnames(TumorST))
    TumorSTn <- subset(TumorST, cells = cells_to_keep)
    TumorSTn$Label <- Clustern$Location[match(rownames(TumorSTn@meta.data), Clustern$CellID)]
    ClusterAdd <- ClusterUpdate(x = n, position = position, df_j = df_j, UMAPembeddings = UMAPembeddings,
      MalCellIDN = MalCellIDN, BdyCellID = BdyCellID, NormalCellID = NormalCellID, MalCellID = MalCellID)
    Clustern <- rbind(Clustern, ClusterAdd)
    TumorSTn$LabelNew <- Clustern$Location[match(rownames(TumorSTn@meta.data), as.character(Clustern$CellID))]
    malignant_levels <- c("Mal", paste0("Mal", seq_len(n)))
    MalCellID <- rownames(TumorSTn@meta.data)[TumorSTn$LabelNew %in% malignant_levels]
    MalCellIDN <- rownames(TumorSTn@meta.data)[TumorSTn$LabelNew == paste0("Mal", n)]
    BdyCellID <- rownames(TumorSTn@meta.data)[TumorSTn$LabelNew == "Bdy"]
    n <- n + 1
  }
  message("Boundary expansion stopped after ", n - 1, " iteration(s).")
  TumorSTn
}

# -------------------- 6. Final Mal/Bdy/nMal annotation --------------------
BoundaryPlot <- function(TumorSTn, TumorST, position_file, OutDir, Sample) {
  spatial_info <- prepare_spatial_neighbors(position_file)
  df_j <- spatial_info$neighbors
  if (!"LabelNew" %in% colnames(TumorSTn@meta.data)) stop("TumorSTn does not contain LabelNew.")
  Mal_barcode <- rownames(TumorSTn@meta.data)[grepl("Mal", TumorSTn$LabelNew)]
  Bdy_barcode <- rownames(TumorSTn@meta.data)[grepl("Bdy", TumorSTn$LabelNew)]
  Normal_Bdy_barcode <- unique(unlist(nbrs(df_j = df_j, MalCellIDAdd = Mal_barcode, CellIDRaw = c(Bdy_barcode, Mal_barcode)), use.names = FALSE))
  Normal_Bdy_barcode <- setdiff(Normal_Bdy_barcode, Mal_barcode)
  nMal_barcode <- setdiff(rownames(TumorST@meta.data), unique(c(Mal_barcode, Bdy_barcode, Normal_Bdy_barcode)))
  annotation <- rbind(data.frame(barcode = Mal_barcode, Location = "Mal"),
    data.frame(barcode = unique(c(Bdy_barcode, Normal_Bdy_barcode)), Location = "Bdy"),
    data.frame(barcode = nMal_barcode, Location = "nMal"))
  TumorST$Location <- annotation$Location[match(rownames(TumorST@meta.data), annotation$barcode)]
  TumorST$Location <- factor(TumorST$Location, levels = c("Mal","Bdy","nMal"))
  write.table(data.frame(spot_id = rownames(TumorST@meta.data), Location = TumorST$Location),
    file.path(OutDir, paste0(Sample, "_boundary_annotation.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  pdf(file.path(OutDir, paste0(Sample, "_BoundaryDefine.pdf")), width = 7, height = 7)
  print(SpatialDimPlot(TumorST, group.by = "Location", cols = c("#CB181D","#1f78b4","#fdb462")))
  dev.off()
  readr::write_rds(TumorST, file.path(OutDir, paste0(Sample, "_BoundaryDefine.rds.gz")), compress = "gz")
  TumorST
}

# -------------------- Run --------------------
TumorST <- STPreProcess(InDir = visium_dir, Sample = sample_name, OutDir = output_dir)
saveRDS(TumorST, file.path(output_dir, paste0(sample_name, "_01_preprocessed.rds")))
TumorST <- STModiCluster(InDir = visium_dir, Sample = sample_name, OutDir = output_dir, TumorST = TumorST,
                         res = cluster_resolution, conda_env = conda_env)
saveRDS(TumorST, file.path(output_dir, paste0(sample_name, "_02_morph_clustered.rds")))
STInferCNV <- STCNV(TumorST = TumorST, assay = infercnv_assay, OutDir = output_dir, Sample = sample_name,
                    gene_order_file = gene_order_file, num_threads = infercnv_threads)
saveRDS(STInferCNV, file.path(output_dir, paste0(sample_name, "_03_infercnv_object.rds")))
TumorST <- STCNVScore(TumorST = TumorST, assay = infercnv_assay, OutDir = output_dir, Sample = sample_name, k = cnv_k)
saveRDS(TumorST, file.path(output_dir, paste0(sample_name, "_04_with_CNV_scores.rds")))
TumorSTn <- BoundaryDefine(TumorST = TumorST, position_file = tissue_positions_file,
                           MalLabel = malignant_cnv_labels, max_iterations = max_iterations)
TumorST <- BoundaryPlot(TumorSTn = TumorSTn, TumorST = TumorST, position_file = tissue_positions_file,
                        OutDir = output_dir, Sample = sample_name)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
message("Done. Final boundary labels are stored in TumorST$Location (Mal/Bdy/nMal).")
