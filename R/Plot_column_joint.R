#' @title Plot the joint distribution between matrix columns.
#'
#' @description \code{Plot_column_joint} is a function used to non-supervised learn / cluster the relationship between columns. (Insight into the joint distribution between samples)
#' @param M A \code{matrix}.
#' @param METRIC The metric used in clustering, can be one in "euclidean", "pearson", and "binary". The pearson metric is calculated as the euclidean metric after rescaling the columns.
#' @param CLUSTER The clustering method, can be one in "MDS" and  "dendrogram".
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
#' Plot_column_joint( Matrix_ex, "euclidean", "dendrogram", "Test1", Group_lab )
#' Plot_column_joint( Matrix_ex, "pearson", "MDS", "Test2", Group_lab )
#'
#' @seealso \code{\link{Plot_column_marginal}}
#'
#'
#' @import ggplot2
#' @import ggdendro
#' @importFrom reshape2 melt
#' @export

Plot_column_joint <- function(M, METRIC = "euclidean",VISUAL = "MDS", HDER = "", GROUP_LABEL = NULL){
  stopifnot(METRIC %in% c("euclidean","pearson","binary"))
  stopifnot(VISUAL %in% c("dendrogram","MDS"))
  stopifnot(length(HDER) == 1)

  if(is.null(colnames(M))) {colnames(M) = paste0("V_",seq_len(ncol(M)))}

  if(!is.null(GROUP_LABEL)) {stopifnot( is.vector( GROUP_LABEL ) & (length(GROUP_LABEL) ==  ncol(M)) )} else{
    GROUP_LABEL = colnames(M)
  }

  if(METRIC == "pearson"){
  dist_matrix <- dist( t(scale(M)), method = "euclidean")
  } else {
  dist_matrix <- dist( t(M), method = METRIC )
  }

  if (VISUAL == "dendrogram") {

  dend <- hclust(dist_matrix)

  dend_data <- dendro_data(dend, type = "rectangle")

  yrange <- range(dend_data$segments$y)

  ywidth <- (yrange[2]-yrange[1])

  p <- ggplot(dend_data$segments) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = dend_data$labels,
              aes(x, y - .02, label = label, colour = GROUP_LABEL),
              hjust = 1.05, angle = 90, size = 2.3) +
    ylim(-3 -.5 * max(nchar(as.character(colnames(M)))),
         yrange[2]+.25*ywidth) + theme_classic() +
    labs(x ="", y = "height",
         title = paste0("Collumn dendrogram with metric: ", METRIC),
         subtitle = HDER)

  p <- p + theme(plot.margin = margin(1,5,1,5,"cm"),
                   legend.position = "bottom",
                   legend.direction = "horizontal",
                   legend.box = "vertical",
                   legend.justification = "center",
                   axis.text.x = element_text(angle = 310, vjust =.9, hjust = .1))

  figwidth = 5 + 0.08 * max(5,max(nchar(as.character( colnames(M) )))) + .25* length( unique(colnames(M)) )

  figheight = 2.65 + 0.1* max(5,max(nchar(as.character( colnames(M))))) + .1 * length( unique(colnames(M)) )

  ggsave(paste0(HDER,"_dg.pdf"),p, height = figheight,width = figwidth)

  }

  if (VISUAL == "MDS"){

    cmd_scale <- cmdscale(dist_matrix)

    plot_df = data.frame(PC_1 = cmd_scale[,1],
                         PC_2 = cmd_scale[,2],
                         Group = GROUP_LABEL)

    p1 <- ggplot(plot_df,aes(x = PC_1,y=PC_2,colour = Group)) +
      geom_point(size = 2, alpha = .8) +
      theme_classic() +
      theme(plot.margin = margin(t = 2, r = 3, b = 2, l = 3.5, unit = "cm"),
            legend.position = "bottom",legend.box = "vertical") +
      labs(x = "First dimmension",
           y="Second dimmension",
           title = paste0("Collumn 2D plot with metric: ", METRIC),
           subtitle = HDER) +
      ylim(max(plot_df$PC_2) +1 ,min(plot_df$PC_2) - 1)

    figheight =  4.5 + .12* ceiling(length( unique(GROUP_LABEL) )/5)
    ggsave(paste0(HDER,"_mds.pdf"),p1,width = 6.35,height = figheight)

  }

}





