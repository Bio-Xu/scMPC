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
  network_filtered <- network %>%
    filter(
      TF_freq > 1,
      TF_Target_freq > 1,
      Target_freq > 2,
      importance > 0.1
    )

  node_table <- data.frame( node = unique(c(network_filtered$TF, network_filtered$Target)), stringsAsFactors = FALSE)
  node_table$Type <- ifelse( node_table$node %in% network_filtered$TF, "TF", "Target")

  write.table( network_filtered, file = file.path(output_dir, "MPC_regulatory_network.tsv"), sep = "\t",row.names = FALSE, quote = FALSE)

  write.table( node_table, file = file.path(output_dir, "MPC_regulatory_network_nodes.tsv"),sep = "\t", row.names = FALSE,quote = FALSE)
}

