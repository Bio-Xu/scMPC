# ============================================================
# MPC transcriptional-regulatory analysis
#
# using all malignant cells from each individual patient following the standard pySCENIC workflow:
#   1) out_SCENIC.loom( Output according to the pySCENIC standard analysis procedure)
#   2) adj.sample.csv (TF-target importance scores)
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(SCopeLoomR)
library(SCENIC)

# ----------------------------
# User settings
# ----------------------------
mpc_label_rds <- "path/to/PriMCLable.rds"
scenic_dir <- "path/to/SCENIC_results"
output_dir <- "results/regulatory_analysis"

patient_col <- "PatientID"
label_col <- "MPCLable"
cell_id_col <- "CellID"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helper: read the pySCENIC loom output
# ----------------------------
ReadSCENICloom <- function(loom_file) {
  loom <- open_loom(loom_file)
  on.exit(close_loom(loom), add = TRUE)
  regulon_mat <- get_regulons(loom, column.attr.name = "Regulons")
  regulon_list <- regulonsToGeneLists(regulon_mat)
  regulon_auc <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
  list(regulon.list = regulon_list,regulon.AUC = regulon_auc)
}

# ----------------------------
# 1. Load MPC/NPC labels
# ----------------------------
labels <- readRDS(mpc_label_rds)
if (cell_id_col %in% colnames(labels)) {
  rownames(labels) <- labels[[cell_id_col]]
}

labels[[patient_col]] <- gsub("GSE", "", labels[[patient_col]])
patients <- unique(labels[[patient_col]])

# ----------------------------
# 2. Identify MPC-enriched regulons in each patient
# RSS is calculated within each patient, followed by Z-score comparison
# between MPC and NPC groups as in the original analysis.
# ----------------------------
mpc_tf_list <- list()
scenic_cache <- list()

for (p in patients) {
  loom_file <- file.path(scenic_dir, paste0("GSE", p),"out_SCENIC.loom" )

  if (!file.exists(loom_file)) {
    warning("Missing SCENIC loom file for patient: ", p)
    next
  }
  scenic_res <- ReadSCENICloom(loom_file)
  scenic_cache[[p]] <- scenic_res
  meta_sub <- labels[labels[[patient_col]] == p & labels[[label_col]] %in% c("MPC", "NPC"),, drop = FALSE]
  regulon_auc <- scenic_res$regulon.AUC
  common_cells <- intersect(colnames(regulon_auc), rownames(meta_sub))

  if (length(common_cells) < 10) {
    warning("Too few matched cells for patient: ", p)
    next
  }
  auc_sub <- regulon_auc[, common_cells]
  cell_annotation <- meta_sub[common_cells, label_col]

  rss <- calcRSS(AUC = getAUC(auc_sub),cellAnnotation = cell_annotation)
  rss <- na.omit(rss)
  rss_plot <- plotRSS(rss,cluster_columns = FALSE, order_rows = TRUE,thr = 0.1, varName = "Type")
  zscore_df <- rss_plot$df

  mpc_tf <- zscore_df %>%
    group_by(Topic) %>%
    filter(Z == max(Z)) %>%
    filter(Type == "MPC") %>%
    ungroup()

  mpc_tf_list[[p]] <- mpc_tf
}

# ----------------------------
# 3. Summarize TF recurrence across patients
# ----------------------------
tf_presence <- bind_rows(
  lapply(names(mpc_tf_list), function(p) {
    data.frame(
      PatientID = p,
      TF = unique(mpc_tf_list[[p]]$Topic),
      stringsAsFactors = FALSE
    )
  })
)

tf_recurrence <- tf_presence %>%
  distinct(PatientID, TF) %>%
  count(TF, name = "N_patients") %>%
  arrange(desc(N_patients))

write.table( tf_recurrence, file = file.path(output_dir, "MPC_TF_recurrence.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# TFs recurrently activated in >=2 or >=3 patients, respectively.
TF_2 <- tf_recurrence$TF[tf_recurrence$N_patients > 1]
TF_3 <- tf_recurrence$TF[tf_recurrence$N_patients > 2]

saveRDS( list(TF_2 = TF_2, TF_3 = TF_3), file.path(output_dir, "MPC_recurrent_TFs.rds"))

# ----------------------------
# 4. Extract TF-target relationships and importance scores
# ----------------------------
tf_target_all <- list()

for (p in intersect(names(mpc_tf_list), names(scenic_cache))) {
  scenic_res <- scenic_cache[[p]]
  active_tf <- unique(mpc_tf_list[[p]]$Topic)

  target_list <- scenic_res$regulon.list[active_tf]

  tf_target <- bind_rows(
    lapply(names(target_list), function(tf) {
      data.frame(
        TF = gsub("\\(\\+\\)", "", tf),
        Target = target_list[[tf]],
        stringsAsFactors = FALSE
      )
    })
  )

  adj_file <- file.path(scenic_dir, paste0("GSE", p), "adj.sample.csv")

  if (!file.exists(adj_file)) {
    warning("Missing adjacency file for patient: ", p)
    next
  }

  importance_df <- read.csv(adj_file, stringsAsFactors = FALSE)
  colnames(importance_df)[1:3] <- c("TF", "Target", "importance")

  tf_target <- merge( tf_target, importance_df[, c("TF", "Target", "importance")], by = c("TF", "Target"), all = FALSE )
  tf_target$PatientID <- p
  tf_target_all[[p]] <- tf_target
}

all_edges <- bind_rows(tf_target_all)

# ----------------------------
# 5. Build the recurrent MPC regulatory network
# ----------------------------
if (nrow(all_edges) > 0) {
  tf_freq <- all_edges %>%
    distinct(PatientID, TF) %>%
    count(TF, name = "TF_freq")

  target_freq <- all_edges %>%
    distinct(PatientID, Target) %>%
    count(Target, name = "Target_freq")

  edge_freq <- all_edges %>%
    mutate(TF_Target = paste(TF, Target, sep = "-")) %>%
    distinct(PatientID, TF_Target) %>%
    count(TF_Target, name = "TF_Target_freq")

  network <- all_edges %>%
    mutate(TF_Target = paste(TF, Target, sep = "-")) %>%
    left_join(tf_freq, by = "TF") %>%
    left_join(target_freq, by = "Target") %>%
    left_join(edge_freq, by = "TF_Target") %>%
    group_by(TF, Target, TF_freq, Target_freq, TF_Target_freq) %>%
    summarise( importance = mean(importance, na.rm = TRUE),.groups = "drop")

  # Original network-filtering criteria.
recurrent_tf_3 <- unique(gsub("\\(\\+\\)", "", TF_3))
  network_filtered <- network %>%
    filter(
      TF %in% recurrent_tf_3,
      importance > 1
    )
  if (nrow(network_filtered) == 0L) {
    stop("No cross-patient edges remain after manuscript thresholds; check TF naming and importance scale.")
  }

  node_table <- data.frame( node = unique(c(network_filtered$TF, network_filtered$Target)), stringsAsFactors = FALSE)
  node_table$Type <- ifelse( node_table$node %in% network_filtered$TF, "TF", "Target")

  write.table( network_filtered, file = file.path(output_dir, "MPC_regulatory_network.tsv"), sep = "\t",row.names = FALSE, quote = FALSE)

  write.table( node_table, file = file.path(output_dir, "MPC_regulatory_network_nodes.tsv"),sep = "\t", row.names = FALSE,quote = FALSE)
}


# ============================================================
# 6. MPC-NPC differential-edge permutation analysis
#    Methods Section 8 | Run in separate R -> Python -> R stages
# ============================================================
# Inputs:
#   - seuretList: named patient Seurat objects (malignant cells only),
#                 with MPC/NPC labels in metadata column "MPCLable".
#                 Load this object before running Section 6.
#   - network_filtered: TF/Target edges from Section 5.
#   - pySCENIC GRNBoost2: Python implementation embedded in this R file.
#
# Workflow:
#   6.1 R: export observed MPC/NPC expression CSV files
#   6.2 R calls embedded Python: infer observed GRNBoost2 importance TSVs
#   6.3 R: export 1,000 independently permuted MPC/NPC expression CSV pairs
#   6.4 R calls embedded Python: infer permuted GRNBoost2 importance TSVs
#   6.5 R: calculate observed differences, null intervals, empirical P, BH FDR
#
# Each patient is analyzed independently. Full random label permutations
# preserve MPC/NPC group sizes (not the original notebook's 50% label swaps).
# All code is distributed in ONE .R file. Python source is embedded below;
# R writes it to a temporary .py file only when stage 6.2 or 6.4 is enabled.
# pySCENIC runs in Python, not in R. Its internal random state is not
# explicitly fixed by the pySCENIC CLI call; record software versions.
# Run stages sequentially and verify
# each stage finishes before enabling the next stage.
# ============================================================

# ----------------------------
# 6.0 Configuration and helpers
# ----------------------------
# Edit the location of your Seurat list before running this section, e.g.:
# seuretList <- readRDS("path/to/patient_malignant_Seurat_list.rds")
# Alternatively, load an RData file containing the appropriate named list.
#
# Use the original MPCLable metadata labels; do not re-identify MPC/NPC.
edge_group_col <- "MPCLable"
edge_dir <- file.path(output_dir, "differential_edges")
edge_nperm <- 1000L
edge_seed <- 1234L
edge_min_cells_per_gene <- 3L  # As in the original expression-export notebook
dir.create(edge_dir, showWarnings = FALSE, recursive = TRUE)

# The candidate edges come from the existing cross-patient TF-target network.
# If Section 5 did not create network_filtered, run it before Section 6.
stopifnot(exists("network_filtered"),
          all(c("TF", "Target") %in% names(network_filtered)))
edge_candidates <- unique(as.data.frame(network_filtered[, c("TF", "Target")]))
edge_candidates <- subset(edge_candidates,
                          !is.na(TF) & !is.na(Target) & TF != Target)
stopifnot(nrow(edge_candidates) > 0L)
write.csv(edge_candidates, file.path(edge_dir, "candidate_edges.csv"),
          row.names = FALSE)

edge_patient_ids <- function() {
  stopifnot(exists("seuretList"), is.list(seuretList),
            length(seuretList) > 0L, !is.null(names(seuretList)),
            all(nzchar(names(seuretList))))
  names(seuretList)
}

# Keep patient file naming consistent with the original Python analysis.
# Example: "161_3" -> "GSE161_3_MPC.csv".
edge_sample_name <- function(patient, group) {
  paste0(if (startsWith(patient, "GSE")) "" else "GSE",
         patient, "_", group)
}

# Export log-normalized, non-imputed RNA expression (cells x genes).
# To stay close to the original notebook, genes with counts > 0 in fewer
# than three cells are filtered WITHIN each group, then the remaining
# group-specific raw counts are log-normalized (scale factor 10,000),
# using full RNA library sizes to match Seurat LogNormalize.
# This can yield different gene columns across groups/permutations.
edge_export_group <- function(obj, cell_names, genes, filename) {
  stopifnot(length(cell_names) > 0L)
  all_counts <- Seurat::GetAssayData(obj, assay = "RNA", layer = "counts")
  # Compute library size from ALL RNA genes, not just candidate network genes.
  # This reproduces the scale-factor normalization denominator used by
  # Seurat::NormalizeData(normalization.method = "LogNormalize").
  lib <- Matrix::colSums(all_counts[, cell_names, drop = FALSE])
  if (any(lib <= 0)) stop("Zero library size: ", filename)
  keep_genes <- intersect(genes, rownames(all_counts))
  counts <- all_counts[keep_genes, cell_names, drop = FALSE]
  counts <- counts[Matrix::rowSums(counts > 0) >= edge_min_cells_per_gene,
                   , drop = FALSE]
  if (nrow(counts) == 0L) stop("No expressed candidate genes: ", filename)
  log_expr <- log1p(counts %*% Matrix::Diagonal(x = 10000 / lib))
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  write.csv(as.matrix(Matrix::t(log_expr)), filename, row.names = TRUE)
}

edge_prepare_patient <- function(patient) {
  obj <- seuretList[[patient]]
  stopifnot(edge_group_col %in% colnames(obj[[]]))
  labs <- as.character(obj[[]][[edge_group_col]])
  names(labs) <- colnames(obj)
  if (anyNA(labs) || !all(labs %in% c("MPC", "NPC")) ||
      !all(c("MPC", "NPC") %in% labs)) {
    stop("Expected only MPC/NPC labels, both present, for ", patient)
  }
  genes <- intersect(unique(c(edge_candidates$TF, edge_candidates$Target)),
                     rownames(obj[["RNA"]]))
  if (length(genes) == 0L) stop("No network genes in ", patient)
  list(obj = obj, labs = labs, genes = genes)
}

# ----------------------------
# Embedded Python runner (single-file distribution)
# ----------------------------
# LANGUAGE: R. The following lines store Python source code as strings.
# When a Python stage is enabled, R writes the embedded code to a temporary
# .py file and invokes it with system2(). GRNBoost2 still runs in Python.
edge_python_executable <- "python3"  # Python environment with pySCENIC installed
edge_pyscenic_executable <- "pyscenic"
edge_tf_file <- "path/to/allTFs_hg38.txt"  # EDIT before running Python stages
edge_workers <- 4L
edge_python_source <- c(
  "#!/usr/bin/env python3",
  "\"\"\"Run pySCENIC GRNBoost2 for observed or permuted MPC/NPC CSV matrices.",
  "",
  "R exports expression matrices; Python performs network inference; R then",
  "computes edge-wise empirical P values and BH FDR. Use the same TF reference",
  "and pySCENIC version for both stages.",
  "\"\"\"",
  "import argparse",
  "import subprocess",
  "from pathlib import Path",
  "",
  "def main():",
  "    ap = argparse.ArgumentParser()",
  "    ap.add_argument(\"--base-dir\", type=Path, required=True)",
  "    ap.add_argument(\"--tf-file\", type=Path, required=True)",
  "    ap.add_argument(\"--stage\", choices=[\"observed\", \"permuted\"], required=True)",
  "    ap.add_argument(\"--pyscenic\", default=\"pyscenic\")",
  "    ap.add_argument(\"--workers\", type=int, default=4)",
  "    ap.add_argument(\"--overwrite\", action=\"store_true\")",
  "    args = ap.parse_args()",
  "    if not args.tf_file.is_file():",
  "        ap.error(f\"TF reference not found: {args.tf_file}\")",
  "    if args.workers < 1:",
  "        ap.error(\"--workers must be >= 1\")",
  "    if args.stage == \"observed\":",
  "        folders = [args.base_dir / \"TrueLable\"]",
  "    else:",
  "        folders = [args.base_dir / \"perturbed\" / f\"perturbed{i:04d}\"",
  "                   for i in range(1, 1001)]",
  "    total = 0",
  "    for folder in folders:",
  "        input_dir = folder / \"expData\"",
  "        if not input_dir.is_dir():",
  "            raise FileNotFoundError(f\"Missing expression directory: {input_dir}\")",
  "        files = sorted(input_dir.glob(\"*.csv\"))",
  "        if not files:",
  "            raise RuntimeError(f\"No expression CSVs in {input_dir}\")",
  "        output_dir = folder / \"importantScore\"",
  "        output_dir.mkdir(parents=True, exist_ok=True)",
  "        for infile in files:",
  "            outfile = output_dir / f\"{infile.stem}.tsv\"",
  "            if outfile.exists() and outfile.stat().st_size > 0 and not args.overwrite:",
  "                print(f\"SKIP {outfile}\", flush=True)",
  "                continue",
  "            cmd = [args.pyscenic, \"grn\", str(infile), str(args.tf_file),",
  "                   \"-o\", str(outfile), \"--num_workers\", str(args.workers)]",
  "            print(f\"RUN {infile}\", flush=True)",
  "            subprocess.run(cmd, check=True)",
  "            if not outfile.exists() or outfile.stat().st_size == 0:",
  "                raise RuntimeError(f\"Missing/empty GRNBoost2 output: {outfile}\")",
  "            total += 1",
  "    print(f\"Completed {total} GRNBoost2 jobs.\", flush=True)",
  "",
  "if __name__ == \"__main__\":",
  "    main()"
)

edge_run_python_stage <- function(stage) {
  stopifnot(stage %in% c("observed", "permuted"))
  if (!file.exists(edge_tf_file)) {
    stop("Set edge_tf_file to an existing TF reference file.")
  }
  script_file <- tempfile(pattern = "edge_grnboost2_", fileext = ".py")
  on.exit(unlink(script_file), add = TRUE)
  writeLines(edge_python_source, script_file, useBytes = TRUE)
  args <- c(shQuote(script_file),
            "--base-dir", shQuote(normalizePath(edge_dir)),
            "--tf-file", shQuote(normalizePath(edge_tf_file)),
            "--stage", stage, "--pyscenic", shQuote(edge_pyscenic_executable),
            "--workers", as.character(edge_workers))
  status <- system2(edge_python_executable, args = args)
  if (!identical(as.integer(status), 0L)) {
    stop("Python pySCENIC stage failed: ", stage, " (exit status ", status, ")")
  }
}

# ----------------------------
# 6.1 Export observed expression matrices (R)
# ----------------------------
# Run this block first, then run Python stage 6.2.
# The output for each patient consists of *_MPC.csv and *_NPC.csv.
run_edge_stage_6_1 <- FALSE
if (run_edge_stage_6_1) {
  edge_patient_ids()
  for (patient in edge_patient_ids()) {
    x <- edge_prepare_patient(patient)
    for (group in c("MPC", "NPC")) {
      file <- file.path(edge_dir, "TrueLable", "expData",
                        paste0(edge_sample_name(patient, group), ".csv"))
      edge_export_group(x$obj, names(x$labs)[x$labs == group],
                        x$genes, file)
    }
  }
}

# ----------------------------
# 6.2 Infer observed GRNBoost2 networks (Python, launched from R)
# ----------------------------
# LANGUAGE: R orchestration -> Python pySCENIC/GRNBoost2 computation.
# Run after stage 6.1. Requires edge_python_executable, edge_tf_file,
# and edge_pyscenic_executable to point to the correct Python environment.
# Reads TrueLable/expData and writes TrueLable/importantScore.
run_edge_stage_6_2 <- FALSE
if (run_edge_stage_6_2) {
  edge_run_python_stage("observed")
}

# ----------------------------
# 6.3 Export 1,000 label-permuted expression matrices (R)
# ----------------------------
# For each patient, randomly permute its complete label vector while
# preserving the number of MPC and NPC cells. Re-export each group matrix
# using the SAME filtering/normalization steps as the observed labels.
# Each permutation index is a shared iteration across all patients, but
# permutations are independently sampled for each patient.
# This can require substantial disk space.
run_edge_stage_6_3 <- FALSE
if (run_edge_stage_6_3) {
  edge_patient_ids()
  set.seed(edge_seed)
  for (patient in edge_patient_ids()) {
    x <- edge_prepare_patient(patient)
    for (b in seq_len(edge_nperm)) {
      shuffled <- sample(unname(x$labs), length(x$labs), replace = FALSE)
      names(shuffled) <- names(x$labs)
      for (group in c("MPC", "NPC")) {
        file <- file.path(edge_dir, "perturbed",
                          sprintf("perturbed%04d", b), "expData",
                          paste0(edge_sample_name(patient, group), ".csv"))
        edge_export_group(x$obj, names(shuffled)[shuffled == group],
                          x$genes, file)
      }
    }
  }
}

# ----------------------------
# 6.4 Infer permuted GRNBoost2 networks (Python, launched from R)
# ----------------------------
# LANGUAGE: R orchestration -> Python pySCENIC/GRNBoost2 computation.
# Run after stage 6.3. Reads all 1,000 permuted expression folders,
# writing their importantScore TSV outputs. This is computationally costly.
run_edge_stage_6_4 <- FALSE
if (run_edge_stage_6_4) {
  edge_run_python_stage("permuted")
}

# ----------------------------
# 6.5 Differential-edge permutation statistics (R)
# ----------------------------
# For each TF-target edge:
#   Delta_omega_obs = importance_MPC - importance_NPC
#   Delta_omega_perm = permuted_importance_MPC - permuted_importance_NPC
#   95% null interval = quantile(Delta_omega_perm, c(0.025, 0.975))
#   one-sided empirical P = (1 + sum(perm >= obs)) / (1000 + 1)
#   BH FDR is calculated across tested candidate edges WITHIN each patient.
#
# IMPLEMENTATION DETAIL NOT SPECIFIED BY THE MANUSCRIPT:
# GRNBoost2 may omit candidate edges. This script treats an absent edge as
# importance=0, keeping the patient-specific candidate universe fixed.
# This is an explicit modeling assumption, not a manuscript-stated rule.
# The original notebook used intersecting edges; outputs may differ.
# Confirm this choice against the source result tables before publication.
edge_read_weights <- function(file, candidates) {
  if (!file.exists(file) || file.size(file) == 0L) {
    stop("Missing or empty GRNBoost2 result: ", file)
  }
  tab <- read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
  if ("target" %in% names(tab)) names(tab)[names(tab) == "target"] <- "Target"
  stopifnot(all(c("TF", "Target", "importance") %in% names(tab)))
  tab <- tab[, c("TF", "Target", "importance")]
  tab$importance <- as.numeric(tab$importance)
  tab <- tab[is.finite(tab$importance), , drop = FALSE]
  if (nrow(tab)) {
    tab <- aggregate(importance ~ TF + Target, tab, max)
  }
  key <- paste(candidates$TF, candidates$Target, sep = "\r")
  idx <- match(key, paste(tab$TF, tab$Target, sep = "\r"))
  weights <- tab$importance[idx]
  weights[is.na(weights)] <- 0
  weights
}

edge_diff_from_files <- function(folder, patient, candidates) {
  mpc <- file.path(folder, "importantScore",
                   paste0(edge_sample_name(patient, "MPC"), ".tsv"))
  npc <- file.path(folder, "importantScore",
                   paste0(edge_sample_name(patient, "NPC"), ".tsv"))
  edge_read_weights(mpc, candidates) - edge_read_weights(npc, candidates)
}

run_edge_stage_6_5 <- FALSE
if (run_edge_stage_6_5) {
  edge_patient_ids()
  results <- list()
  for (patient in edge_patient_ids()) {
    message("Testing TF-target edges: ", patient)
    # Use only edges with both genes present in this patient's RNA assay.
    obj <- seuretList[[patient]]
    genes <- rownames(obj[["RNA"]])
    candidates <- subset(edge_candidates, TF %in% genes & Target %in% genes)
    stopifnot(nrow(candidates) > 0L)

    observed <- edge_diff_from_files(
      file.path(edge_dir, "TrueLable"), patient, candidates)
    null <- matrix(NA_real_, nrow = nrow(candidates), ncol = edge_nperm)
    for (b in seq_len(edge_nperm)) {
      folder <- file.path(edge_dir, "perturbed",
                          sprintf("perturbed%04d", b))
      null[, b] <- edge_diff_from_files(folder, patient, candidates)
    }

    # Full 1000-permutation null distribution; no missing iterations allowed.
    p_emp <- (1 + rowSums(null >= observed)) / (edge_nperm + 1)
    fdr <- p.adjust(p_emp, method = "BH")
    result <- data.frame(
      PatientID = patient, TF = candidates$TF, Target = candidates$Target,
      effect_size = observed,
      null_lower = apply(null, 1, quantile, probs = 0.025),
      null_upper = apply(null, 1, quantile, probs = 0.975),
      empirical_p = p_emp, FDR = fdr,
      MPC_enriched = observed > 0 & fdr <= 0.05,
      stringsAsFactors = FALSE
    )
    results[[patient]] <- result
  }
  all_edge_results <- dplyr::bind_rows(results)
  write.csv(all_edge_results,
            file.path(edge_dir, "all_patient_edge_permutation_results.csv"),
            row.names = FALSE)
  write.csv(subset(all_edge_results, effect_size > 0 & empirical_p < 0.05),
            file.path(edge_dir, "MPC_enriched_nominal_P_lt_0.05.csv"),
            row.names = FALSE)
  write.csv(subset(all_edge_results, MPC_enriched),
            file.path(edge_dir, "MPC_enriched_BH_FDR_le_0.05.csv"),
            row.names = FALSE)
  message("Nominal MPC-enriched edges: ",
          sum(all_edge_results$effect_size > 0 &
                all_edge_results$empirical_p < 0.05))
  message("FDR-significant MPC-enriched edges: ",
          sum(all_edge_results$MPC_enriched))
}