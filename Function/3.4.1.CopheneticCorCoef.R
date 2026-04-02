########################
# TODO: Evaluate Clustering Results Based on Cophenetic Correlation Coefficient
#
#' @description 
#' 1) Calculate the cophenetic correlation coefficient using the `cophcor` function based on the consensus matrix provided by the consensus clustering results;
#' 2) Alternatively, calculate cophenetic distances using the `cophenetic` function from the distance matrix (dist object) and obtain the coefficient;
#' 3) Plot a line chart to show the cophenetic correlation coefficients for various numbers of clusters. The larger the coefficient, the better the clustering performance.
#'
#' @param consensus_result Consensus clustering results, with list names representing labels for different clustering results (e.g., 2...10), and each entry contains the consensus matrix as the first object.
#' @param dist_obj Distance matrix list, with list names representing labels for different distance methods (e.g., various distance methods). Each entry contains a dist object, default is NULL.
#' @param linkage Linkage method used for hierarchical clustering, default is "average", other methods are available in hclust.
#' @param output Output directory (for saving plots).
#' @param plot Plot output format, default is "png", NULL means no plot.
#' 
#' @output A vector of cophenetic correlation coefficients for each number of clusters, with the length equal to the length of `consensus_result`.
########################
CopheneticCorCoef <- function(consensus_result, dist_obj = NULL, linkage = "average", output, plot = "png"){
  library(NMF)
  library(ggplot2)
  
  #####
  # 1 Calculate Based on Consensus Clustering Results
  ####
  if(!is.null(consensus_result)){
    t.ccc <- sapply(consensus_result, function(x){
      cophcor(x[[1]], linkage = linkage)
    })
  }
    
  #####
  # 2 Calculate Based on Hierarchical Clustering Results
  ####
  if(!is.null(dist_obj)){
    t.ccc <- sapply(dist_obj, function(x){
      d1 <- x
      hc <- hclust(d1, linkage)
      d2 <- cophenetic(hc)
      cor(d1, d2)
    })
  }
  
  #####
  # 3 Plot and Output the Results
  ####
  if(!is.null(plot)){
    ccc_df <- tibble::enframe(t.ccc, name = "Label", value = "Coefficient")
    ccc_df$Label <- ordered(ccc_df$Label, levels = ccc_df$Label)
    if(plot == "png"){
      png(file.path(output, "CopheneticCorCoef.png"))
      print(ggplot(ccc_df, mapping = aes(x = Label, y = Coefficient, group = 1)) + ylab("Cophenetic Correlation Coefficient") +
              geom_line(col="steelblue") + geom_point() + theme_bw())
      dev.off()
    }else if(plot == "pdf"){
      pdf(file.path(output, "CopheneticCorCoef.pdf"))
      print(ggplot(ccc_df, mapping = aes(x = Label, y = Coefficient, group = 1)) + ylab("Cophenetic Correlation Coefficient") +
              geom_line(col="steelblue") + geom_point() + theme_bw())
      dev.off()
    }
  }

  return(t.ccc)
}
                                   