# ============================================================
# 07. Identification of the PPI-prioritized MPC signature
# ============================================================
#
# Workflow:
# (1) patient-specific differential expression in R
# (2) RRA integration are implemented in R
# (3) STRING network analysis and Cytoscape/cytoHubba prioritization were performed externally
#####

library(Seurat)
library(RobustRankAggreg)

# Set the metadata column containing the existing MPC/NPC labels.
group_col <- "MPC_group"
mpc_label <- "MPC"
npc_label <- "NPC"


## 1. Patient-specific differential expression

de_tables <- list()
gene_ranks <- list()
tested_genes <- list()

for (p in names(seuretList)) {
  obj <- seuretList[[p]]
  stopifnot(group_col %in% colnames(obj[[]]))
  Idents(obj) <- group_col

  de <- FindMarkers(
    obj, ident.1 = mpc_label, ident.2 = npc_label,
    assay = "SCT", test.use = "wilcox",
    logfc.threshold = 0, min.pct = 0
  )

  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(de))[1]
  stopifnot(!is.na(fc_col))
  tested_genes[[p]] <- rownames(de)

  de$gene <- rownames(de)
  de <- de[!is.na(de$p_val_adj) &
             de$p_val_adj < 0.01 &
             de[[fc_col]] > 2, , drop = FALSE]
  de <- de[order(-de[[fc_col]]), , drop = FALSE]
  de_tables[[p]] <- de
  gene_ranks[[p]] <- de$gene
}

saveRDS(de_tables, "MPC_signature/patient_DE.rds")
stopifnot(all(lengths(gene_ranks) > 0))

## 2. Robust Rank Aggregation
gene_universe <- unique(unlist(tested_genes, use.names = FALSE))
rra <- aggregateRanks(glist = gene_ranks, N = length(gene_universe))
rra <- rra[order(rra$Score), , drop = FALSE]
rra_genes <- subset(rra, Score < 0.05)

write.csv(rra, "MPC_signature/RRA_all.csv", row.names = FALSE)
write.csv(rra_genes, "MPC_signature/RRA_P_lt_0.05.csv",row.names = FALSE)


# 3. External PPI network analysis and cytoHubba prioritization
#
# STRING and cytoHubba analyses were performed externally using
# STRING v12.0 and Cytoscape, rather than directly in R.
#
# (1) PPI network construction
# Genes with RRA Score < 0.05 were submitted to STRING v12.0.
# Interactions with a combined confidence score >= 0.700 were
# retained. Self-loops, duplicate edges, and disconnected nodes
# were removed from the resulting PPI network.
#
# (2) cytoHubba topological analysis
# The PPI network was imported into Cytoscape and analyzed using
# 12 cytoHubba topological algorithms:
# MCC, DMNC, MNC, Degree, EPC, BottleNeck, EcCentricity,
# Closeness, Radiality, Betweenness, Stress, and
# Clustering Coefficient.
# The top 20 ranked genes were selected from each algorithm.
#
# (3) Gene prioritization
# For each gene, the number of algorithms ranking it among the
# top 20 was counted. Gene sets sig_1 through sig_12 were defined
# by requiring support from at least 1 through 12 algorithms,
# respectively.
#
# The sig_4 gene set, comprising genes supported by at least
# four algorithms, was selected as the final PPI-prioritized
# MPC signature, yielding 22 genes.
#