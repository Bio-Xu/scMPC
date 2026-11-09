# ============================================================
# CellChat analysis of MPC/NPC interactions with the tumor microenvironment
# ============================================================

library(CellChat)
library(Seurat)
library(dplyr)
library(tidyr)
library(future)
library(ggplot2)

# ----------------------------
# User settings
# ----------------------------
seurat_list_rds <- "path/to/PriSeuratObjList.rds"
output_dir <- "results/CellChat"
group_col <- "cellchat"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "patient_objects"), recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# 1. Patient-level CellChat analysis
# Each Seurat object should contain MPC, NPC and TME cell-type labels
# in the metadata column specified by group_col.
# ----------------------------
run_cellchat_patient <- function(seurat_obj, group_col = "cellchat") {
  # Remove normal epithelial cells if present, following the original workflow.
  if (group_col %in% colnames(seurat_obj@meta.data)) {
    keep_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[group_col]] != "Normal_Epithelial"]
    seurat_obj <- subset(seurat_obj, cells = keep_cells)
  }

  cellchat_db <- CellChatDB.human
  exp_data <- GetAssayData(seurat_obj,assay = "RNA",layer = "counts")
  exp_data <- normalizeData(exp_data)
  meta <- seurat_obj@meta.data
  cellchat <- createCellChat( object = exp_data, meta = meta,group.by = group_col)
  cellchat@DB <- cellchat_db
  cellchat <- subsetData(cellchat)
  future::plan("sequential")
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb( cellchat, type = "triMean", population.size = FALSE)
  cellchat <- filterCommunication( cellchat,min.cells = 10 )
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
}
seurat_list <- readRDS(seurat_list_rds)

cellchat_list <- list()

for (p in names(seurat_list)) {
  message("Processing patient: ", p)
  cellchat_obj <- run_cellchat_patient( seurat_list[[p]], group_col = group_col)
  saveRDS(cellchat_obj,file.path(output_dir, "patient_objects", paste0(p, "_CellChat.rds")))
  cellchat_list[[p]] <- cellchat_obj
}
saveRDS(cellchat_list, file.path(output_dir, "CellChat_patient_list.rds"))

# ----------------------------
# 2. Extract MPC/NPC communications across patients
# ----------------------------
all_lr <- list()

for (p in names(cellchat_list)) {
  chat <- cellchat_list[[p]]

  lr_source <- subsetCommunication(chat,sources.use = c("MPC", "NPC"), targets.use = NULL )
  lr_target <- subsetCommunication( chat, sources.use = NULL,targets.use = c("MPC", "NPC"))
  lr <- bind_rows(lr_source, lr_target)

  # Exclude direct MPC/NPC-to-MPC/NPC interactions.
  lr <- lr[ !(lr$source %in% c("MPC", "NPC") & lr$target %in% c("MPC", "NPC")), ,drop = FALSE]

  lr_name_col <- if ("interaction_name_2" %in% colnames(lr)) {
    "interaction_name_2"
  } else {
    "interaction_name"
  }

  lr$LR <- lr[[lr_name_col]]
  lr <- lr %>%
    group_by(LR) %>%
    mutate( LR_unique_in_MPC = ifelse( any(source == "NPC" | target == "NPC"), "NO", "YES") ) %>%
    group_by(pathway_name) %>%
    mutate( pathway_unique_in_MPC = ifelse( any(source == "NPC" | target == "NPC"), "NO", "YES" )) %>%
    ungroup()
  lr$patientID <- p
  lr$Type <- ifelse( lr$source == "NPC" | lr$target == "NPC", "NPC", "MPC" )
  all_lr[[p]] <- lr
}

lr_all <- bind_rows(all_lr)
write.table( lr_all, file = file.path(output_dir, "all_MPC_NPC_communications.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------
# 3. Recurrent MPC-specific ligand-receptor interactions
# ----------------------------
mpc_specific_lr <- lr_all %>%
  filter(LR_unique_in_MPC == "YES") %>%
  group_by(LR) %>%
  summarise( N_patients = n_distinct(patientID), Patients = paste(sort(unique(patientID)), collapse = ","), .groups = "drop") %>%
  arrange(desc(N_patients))
write.table( mpc_specific_lr, file = file.path(output_dir, "MPC_specific_LR_recurrence.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------
# 3.1. Cross-patient recurrence of MPC-specific interactions
# ----------------------------
lr_recurrence <- mpc_specific_lr %>%
  mutate(Group = ifelse(N_patients == 1, "special", "share"))

lr_recurrence_summary <- lr_recurrence %>%
  count(Group, name = "Freq") %>%
  mutate(pro = Freq / sum(Freq), Group_label = paste0(Group, " (", Freq, ")"))

write.table(lr_recurrence, file = file.path(output_dir, "MPC_specific_LR_patient_groups.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(lr_recurrence_summary, file = file.path(output_dir, "MPC_specific_LR_patient_group_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

if (nrow(lr_recurrence_summary) > 0) {
  p <- ggplot(lr_recurrence_summary, aes(x = "", y = pro, fill = Group)) +
    geom_col(width = 0.1, color = "black") +
    coord_polar(theta = "y") +
    theme_void() +
    scale_fill_manual(values = c("special" = "#7F9BD2", "share" = "#C5D99D"),
                      labels = setNames(lr_recurrence_summary$Group_label, lr_recurrence_summary$Group)) +
    labs(title = "MPC-specific LR recurrence", fill = "LR Number") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(output_dir, "MPC_specific_LR_patient_groups.pdf"), p, width = 6, height = 6)
}

# ----------------------------
# 3.2. CellChatDB annotation of MPC-specific interactions
# ----------------------------
annotation_db <- CellChatDB.human$interaction %>%
  select(interaction_name, interaction_name_2, annotation, pathway_name) %>%
  distinct()

annotation_lookup <- annotation_db %>%
  select(LR = interaction_name_2, annotation, pathway_name) %>%
  distinct()

annotated_lr <- lr_all %>%
  filter(
    LR_unique_in_MPC == "YES",
    LR %in% mpc_specific_lr$LR[mpc_specific_lr$N_patients > 1]
  ) %>%
  select(LR, pathway_name, any_of("annotation")) %>%
  distinct()

if ("annotation" %in% names(annotated_lr)) {
  annotated_lr <- annotated_lr %>% select(-annotation)
}
annotated_lr <- annotated_lr %>%
  left_join(annotation_lookup, by = c("LR", "pathway_name")) %>%
  distinct(LR, annotation, pathway_name)

if (anyNA(annotated_lr$annotation)) {
  warning("Some retained interactions have no matching CellChatDB annotation.")
}

annotation_summary <- annotated_lr %>%
  count(annotation, name = "Freq") %>%
  mutate(pro = Freq / sum(Freq), percentage = 100 * pro) %>%
  arrange(desc(Freq))

ecm_summary <- annotation_summary %>%
  summarise(N_total = sum(Freq),
            N_ECM_Receptor = sum(Freq[annotation == "ECM-Receptor"]),
            ECM_Receptor_percent = 100 * N_ECM_Receptor / N_total)

pathway_annotation <- annotated_lr %>%
  group_by(pathway_name) %>%
  summarise(Freq = n(), annotation_list = paste(sort(unique(annotation)), collapse = "; "),
            .groups = "drop") %>%
  arrange(annotation_list, desc(Freq)) %>%
  mutate(pro = 100 * Freq / sum(Freq))

write.table(annotated_lr, file = file.path(output_dir, "MPC_specific_LR_CellChatDB_annotation.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(annotation_summary, file = file.path(output_dir, "MPC_specific_LR_annotation_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ecm_summary, file = file.path(output_dir, "MPC_specific_LR_ECM_receptor_proportion.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pathway_annotation, file = file.path(output_dir, "MPC_specific_LR_pathway_annotation.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

print(annotation_summary)
print(ecm_summary)

if (nrow(pathway_annotation) > 0) {
  base_colors <- c("ECM-Receptor" = "#7F9BD2",
                   "Secreted Signaling" = "#45B69C",
                   "Cell-Cell Contact" = "#C06C84")
  pathway_annotation <- pathway_annotation %>%
    mutate(pathway_name = factor(pathway_name, levels = pathway_name))
  pathway_colors <- unlist(lapply(split(pathway_annotation, pathway_annotation$annotation_list),
                                 function(df) {
                                   base <- base_colors[df$annotation_list[1]]
                                   if (is.na(base)) base <- "#A0A0A0"
                                   setNames(scales::seq_gradient_pal("white", base, "Lab")(
                                     seq(1, 0.1, length.out = nrow(df))),
                                     as.character(df$pathway_name))
                                 }))
  p <- ggplot(pathway_annotation, aes(x = "", y = pro, fill = pathway_name)) +
    geom_col(width = 0.05, color = "black") +
    coord_polar(theta = "y") +
    theme_void() +
    scale_fill_manual(values = pathway_colors) +
    labs(title = "LR Pathway Class", fill = "Pathway") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(output_dir, "MPC_specific_LR_pathway_annotation.pdf"),
         p, width = 8, height = 6)
}

# ----------------------------
# 4. Compare pathway-level information flow for MPC and NPC
# ----------------------------
mpc_flow <- lr_all %>%
  filter(source == "MPC" | target == "MPC") %>%
  group_by(patientID, pathway_name) %>%
  summarise( MPC_Pathway_LRprob = sum(prob, na.rm = TRUE), MPC_Pathway_LRcount = n_distinct(LR),.groups = "drop")
npc_flow <- lr_all %>%
  filter(source == "NPC" | target == "NPC") %>%
  group_by(patientID, pathway_name) %>%
  summarise(  NPC_Pathway_LRprob = sum(prob, na.rm = TRUE), NPC_Pathway_LRcount = n_distinct(LR),.groups = "drop")
information_flow <- full_join( mpc_flow, npc_flow,by = c("patientID", "pathway_name"))

information_flow[is.na(information_flow)] <- 0
denom <- information_flow$MPC_Pathway_LRprob + information_flow$NPC_Pathway_LRprob
information_flow$MPC_Pathway_LRprob_proportion <- ifelse( denom > 0, information_flow$MPC_Pathway_LRprob / denom, NA)
information_flow$NPC_Pathway_LRprob_proportion <- ifelse( denom > 0, information_flow$NPC_Pathway_LRprob / denom,NA)

write.table( information_flow, file = file.path(output_dir, "MPC_NPC_pathway_information_flow.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# ----------------------------
# 5. Receptor-expression and TME-state correlations (Figure 4g)
# Cell-state gene sets were manually curated in the original analysis.
# Candidate LR interactions and representative panels were selected
# manually based on MPC-to-TME communication and biological relevance.
# ----------------------------
library(AUCell)

state_signature_rds <- "path/to/ImmuneStromalState.list.final.RDS"
state_dir <- file.path(output_dir, "receptor_TME_correlations")
dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)
imm_Str <- readRDS(state_signature_rds)

# Patient-specific TME cell types included in the original analysis.
# 每个患者中筛选出状态列表中有的细胞类型进行后续分析--这是seurat对象的命名方法不是状态基因集的命名方式
patient_tme <- list(
    'GSE161529_1' = c('Macrophage','Fibroblast','CD8 Memory T','Memory B','IgA Plasma'),
    'GSE161529_2' = c('Macrophage','CD8 Memory T','Fibroblast','IgG Plasma','IgA Plasma','Memory B','Treg','CD4 T' ),
    'GSE161529_3' = c('Macrophage','Fibroblast'),
    'GSE161529_4' = c('Macrophage','Fibroblast','Memory B','IgG Plasma','IgA Plasma' ),
    'GSE161529_5' = c('Macrophage','IgA Plasma','CD8 Memory T','CD4 T','Fibroblast','Treg','Memory B','IgG Plasma'),
    'GSE161529_6' = c('IgA Plasma','Macrophage','Treg','CD8 Memory T','IgG Plasma','Fibroblast','CD4 T','Memory B'),
    'GSE161529_7' = c('Macrophage'),
    'GSE167036_8' = c('CD8 Memory T','CD4 Helper T','Memory B','IgG Plasma','Monocyte','CD4 Memory T','Treg','Endothelial','CD4 T'),
    'GSE180286_1' = c('IgG Plasma','CD8 Memory T','Macrophage','IgA Plasma','Endothelial','Memory B','Fibroblast','CD4 Helper T','CD4 T'),
    'GSE180286_3' = c('Fibroblast','Endothelial','Macrophage','IgG Plasma','IgA Plasma','CD8 Memory T'),
    'GSE180286_5' = c('Fibroblast','IgG Plasma','Endothelial','Macrophage','IgA Plasma','Memory B','CD8 Memory T','CD4 T' ),
    'GSE195861_1' = c('Treg','Memory B','Macrophage','Monocyte','CD8 Memory T','IgG Plasma','CD4 T'),
    'GSE195861_2' = c('CD8 Memory T','Treg','Macrophage','IgG Plasma','CD4 T','Monocyte','Memory B'),
    'GSE195861_3' = c('Treg','IgG Plasma','CD4 T','Monocyte','Memory B','Macrophage','CD8 Memory T'),
    'GSE195861_5' = c('Macrophage','Monocyte')
)

# Match Seurat cell-type names to the curated signature database.
#因集合的细胞类型名称与我们研究数据中的细胞类型名称不对应,建立seurat对象与sig细胞类型间的对应关系，方便后续取用数据
cell_index <- data.frame(seurat = c('Macrophage','Monocyte','Fibroblast','Endothelial',
                                   'CD8 Memory T','Treg','CD4 T','CD4 Helper T','CD4 Memory T',
                                   'Memory B','IgA Plasma','IgG Plasma'),
                        Sig = c('Macrophage','Monocyte','CAF','Endothelial cell',
                                   'CD8 Memory T cell','CD4 Regulatory T cell','CD4 T cell','CD4 Helper T cell','CD4 Memory T cell',
                                   'Memory B cell','IgA Plasma cell','IgG Plasma cell')
                        )
cell_index

# 5.1. AUCell state scores from SCT counts, separately for each patient and cell type.
# Rankings are computed once per patient/cell type and reused across signatures.
score_list <- list()
for (p in intersect(names(seurat_list), names(patient_tme))) {
  obj <- seurat_list[[p]]
  for (ct in patient_tme[[p]]) {
    sig_type <- cell_index$Sig[match(ct, cell_index$seurat)]
    if (is.na(sig_type)) next
    ids <- rownames(obj@meta.data)[!is.na(obj$CellType) & obj$CellType == ct]
    if (length(ids) == 0) next
    states <- imm_Str$Meta$Name[imm_Str$Meta$CellType == sig_type]
    states <- unique(states[states %in% names(imm_Str$GeneSet)])
    states <- states[vapply(states, function(s) length(imm_Str$GeneSet[[s]]) > 0, logical(1))]
    if (length(states) == 0) next

    counts <- GetAssayData(obj, assay = "SCT", layer = "counts")[, ids, drop = FALSE]
    rankings <- AUCell_buildRankings(as.matrix(counts), plotStats = FALSE)
    scores <- lapply(states, function(s) {
      auc <- AUCell_calcAUC(imm_Str$GeneSet[[s]], rankings)
      as.numeric(getAUC(auc)[1, ids])
    })
    names(scores) <- states
    score_df <- as.data.frame(scores, check.names = FALSE)
    rownames(score_df) <- ids
    if (is.null(score_list[[p]])) score_list[[p]] <- list()
    score_list[[p]][[ct]] <- score_df
  }
}
saveRDS(score_list, file.path(state_dir, "AUCell_TME_state_scores.rds"))

# 5.2. Receptor genes from MPC-specific outgoing LR interactions.
# As in the notebook, receptor complexes are split on underscores.
tme_receptors <- lapply(names(seurat_list), function(p) {
  receptors <- lr_all$receptor[
    lr_all$patientID == p &
      lr_all$source == "MPC" &
      lr_all$LR_unique_in_MPC == "YES"
  ]
  unique(unlist(strsplit(receptors[!is.na(receptors)], "_"), use.names = FALSE))
})
names(tme_receptors) <- names(seurat_list)

# Pearson correlations are computed within each patient and TME cell type.
cor_results <- list()
for (p in names(score_list)) {
  obj <- seurat_list[[p]]
  genes <- tme_receptors[[p]]
  for (ct in names(score_list[[p]])) {
    scores <- score_list[[p]][[ct]]
    ids <- intersect(rownames(scores), colnames(obj))
    genes_use <- intersect(genes, rownames(obj[["SCT"]]))
    if (length(ids) < 3 || length(genes_use) == 0) next
    expr <- GetAssayData(obj, assay = "SCT", layer = "data")[genes_use, ids, drop = FALSE]
    for (g in genes_use) {
      for (s in colnames(scores)) {
        x <- as.numeric(expr[g, ids])
        y <- as.numeric(scores[ids, s])
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) next
        test <- cor.test(x[ok], y[ok], method = "pearson")
        cor_results[[length(cor_results) + 1L]] <- data.frame(
          patientID = p, CellType = ct, receptor = g, state = s,
          N_cells = sum(ok), correlation = unname(test$estimate),
          p_value = test$p.value
        )
      }
    }
  }
}
cor_table <- dplyr::bind_rows(cor_results)
write.table(cor_table, file.path(state_dir, "receptor_TME_state_Pearson_correlations.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 5.3. Plot manually selected significant positive receptor-TME state correlations.
# Select the combinations for visualization after inspecting the correlation table.
# Each row identifies one patient, cell type, receptor, and TME state.
selected_correlations <- data.frame(
  patientID = character(),
  CellType = character(),
  receptor = character(),
  state = character(),
  stringsAsFactors = FALSE
)

plot_list <- list()
plot_stats <- list()
for (i in seq_len(nrow(selected_correlations))) {
  p <- selected_correlations$patientID[i]
  ct <- selected_correlations$CellType[i]
  g <- selected_correlations$receptor[i]
  state <- selected_correlations$state[i]

  if (is.null(score_list[[p]]) || is.null(score_list[[p]][[ct]])) next
  obj <- seurat_list[[p]]
  scores <- score_list[[p]][[ct]]
  if (!state %in% colnames(scores) || !g %in% rownames(obj[["SCT"]])) next

  ids <- intersect(rownames(scores), colnames(obj))
  if (length(ids) < 3) next
  expression <- GetAssayData(obj, assay = "SCT", layer = "data")
  df <- data.frame(
    expression = as.numeric(expression[g, ids]),
    AUCell_score = as.numeric(scores[ids, state])
  )
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 3 || sd(df$expression) == 0 || sd(df$AUCell_score) == 0) next

  test <- cor.test(df$expression, df$AUCell_score, method = "pearson")
  if (unname(test$estimate) <= 0.4 || test$p.value >= 0.05) next

  plot_stats[[length(plot_stats) + 1L]] <- data.frame(
    patientID = p, CellType = ct, receptor = g, state = state,
    N_cells = nrow(df), correlation = unname(test$estimate), p_value = test$p.value
  )
  plot_list[[length(plot_list) + 1L]] <- ggplot(df, aes(expression, AUCell_score)) +
    geom_point(color = "blue", size = 1) +
    geom_smooth(method = "lm", color = "red") +
    labs(
      title = paste0(p, " | ", state, "\\nr = ", round(unname(test$estimate), 2),
                     " (p = ", signif(test$p.value, 3), ")"),
      x = paste(g, "expression"), y = "AUCell score"
    ) +
    theme_minimal(base_size = 7)
}

if (length(plot_list) > 0) {
  ggsave(
    file.path(state_dir, "receptor_TME_state_correlation_plots.pdf"),
    patchwork::wrap_plots(plot_list, ncol = 4),
    width = 12, height = max(4, 2.8 * ceiling(length(plot_list) / 4)),
    limitsize = FALSE
  )
}
if (length(plot_stats) > 0) {
  write.table(
    dplyr::bind_rows(plot_stats),
    file.path(state_dir, "selected_receptor_TME_state_correlations.tsv"),
    sep = "\\t", row.names = FALSE, quote = FALSE
  )
}
