###############################
# TODO: Optimize functional gene sets by integrating different GO terms
#
#' @param GOTermsID: Functional gene set to be optimized (GO term ID).
#' @param AmiGO: GO terms obtained from AmiGO 2 platform (text file).
#' @param output_path: Output directory path for Venn diagram.
#' @param output_filename: Output filename for Venn diagram.
######################
integrationGeneSet <- function(GOTermsID, AmiGO, output_path, output_filename){
  
  # Load required packages
  require(org.Hs.eg.db)
  require(dplyr)
  require(msigdbr)
  require(VennDiagram)
  require(RColorBrewer)
  require(qpdf) 
  
  # 1. Extract GO term gene set from R package
  DNA_geneID <- get(GOTermsID, org.Hs.egGO2ALLEGS)
  GOR <- mget(DNA_geneID, org.Hs.egSYMBOL) %>% unlist() 
  
  # 2. Extract GO term gene set from AmiGO2
  ami <- read.table(AmiGO, sep = "\t") %>% unlist()
  
  # 3. Extract GO term gene set from MsigDB
  Msi <- msigdbr(species = "Homo sapiens", category = "C5")
  M <- Msi[grep(pattern = GOTermsID, Msi$gs_exact_source), ]
  MSIG <- M[, "gene_symbol"] %>% unlist()
 
  # 4. Generate Venn diagram comparing the three sets
  V <- venn.diagram(
    x = list(
      R_GO = GOR,
      AmiGO2 = ami,
      MsigDB = MSIG
    ),
    filename = NULL,
    fill = brewer.pal(7, "Set1")[1:3]
  )
  
  # Save Venn diagram to PDF
  pdf(file.path(output_path, output_filename))
  grid.draw(V)
  dev.off()
  
  # 5. Create union set from all three sources
  G <- union(ami, GOR)
  G1 <- union(G, MSIG)
  
  return(G1)
}