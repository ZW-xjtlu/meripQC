#' @title Plot the joint distribution between matrix columns.
#'
#' @description \code{Plot_column_joint} is a function used to non-supervised learn / cluster the relationship between columns. (Insight into the joint distribution between samples)
#' @param M A \code{matrix}.
#' @param METRIC Can be either "euclidean" or "pearson". The later is simply the euclidean distance after rescale the columns.
#' @param VISUAL Can be either "dendrogram" or "MDS".
#' @param HDER The subtitle and the file name of the plot.
#' @param GROUP_LABEL Optional, a vector used for colour labeling or faceting the dendrogram or MDS plot.
#'
#' @details By default, the column names of the matrix M will be used as the sample labels, other wise, it will use V_{1:ncol(M)}.
#'
#' @return A plot for column wised clustering analysis.
#'
#' @examples
#' Matrix_ex <- matrix(rnorm(9000),300,30)
#' Group_lab <- paste0( "V_",rep(1:10,each = 3) )
#' Plot_column_joint( Matrix_ex, "density", "Test1", Group_lab )
#' Plot_column_joint( Matrix_ex, "box", "Test2", Group_lab )
#'
#' @seealso \code{\link{Plot_column_marginal}}
#'
#'
#' @import ggplot2
#' @import ggdendro
#' @importFrom reshape2 melt
#' @export

Plot_column_joint <- function(M, METRIC = "euclidean",VISUAL = "dendrogram", HDER = "", GROUP_LABEL = NULL){
  stopifnot(METRIC %in% c("euclidean","pearson"))
  stopifnot(VISUAL %in% c("dendrogram","MDS"))
  if(is.null(colnames(M))) {colnames(M) = paste0("V_",seq_len(ncol(M)))}

  require("ggdendro")
  rownames(dist_matrix) = Clean_head(rownames(dist_matrix))
  rownames(dist_matrix) = Clean_head(rownames(dist_matrix))

  dend <- as.dendrogram(hclust(as.dist(dist_matrix)))

  dend_data <- dendro_data(dend, type = "rectangle")

  dend_data$labels[["Writer"]] = Extract_Role(dend_data$labels[["label"]])

  p <- ggplot(dend_data$segments) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
    geom_text(data = dend_data$labels, aes(x, y - .02, label = label, colour = Writer),
              hjust = 1, angle = 90, size = 2.3)+
    ylim(-1, 1.2) + theme_classic() + labs(x ="", y = "height", title = "Dendrogram of differential methylations", subtitle = Name) + scale_color_brewer(palette = "Dark2")
  ggsave(paste0(Name,"_dg.pdf"),p,height = 3.2,width = 4.6)
  Dendrogram_Plot(sqrt_dist_matrix,Name = "sqrt_asym_dist")



}





