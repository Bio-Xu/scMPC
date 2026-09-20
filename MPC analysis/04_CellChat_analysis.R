# ============================================================
# CellChat analysis of MPC/NPC interactions with the tumor microenvironment
# ============================================================

library(CellChat)
library(Seurat)
library(dplyr)
library(tidyr)
library(future)

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
