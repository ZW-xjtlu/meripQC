#' @title Plot the column wise data distribution.
#'
#' @description \code{Plot_column_marginal} is a function used to plot the columned wise data distributions. (marginal distribution for each sample)
#' @param M A \code{matrix}.
#' @param TYPE Can be one in "violine", "box", and "density".
#' @param HDER The subtitle and the file name of the plot.
#' @param GROUP_LABEL Optional, a vector used for colour labeling or faceting the box plot and density plot, the length should equal to the column number of M.
#' @param VALUE_LABEL Optional, The label for the entries of M.
#'
#' @details By default, the column names of the matrix M will be used as the sample labels, other wise, it will use V_{1:ncol(M)}.
#'
#' @return A plot for column wised marginal distribution visualization.
#'
#' @examples
#' Matrix_ex <- matrix(rnorm(9000),300,30)
#' Group_lab <- paste0( "V_",rep(1:10,each = 3) )
#'
#' Plot_column_marginal( Matrix_ex, "violine", "Test2", Group_lab )
#' Plot_column_marginal( Matrix_ex, "density", "Test1", Group_lab )
#' Plot_column_marginal( Matrix_ex, "box", "Test2", Group_lab )
#'
#' @seealso \code{\link{Plot_column_joint}}
#'
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export

Plot_column_marginal <- function(M,TYPE = "violine",HDER = "",GROUP_LABEL = NULL,VALUE_LABEL = "value"){

  stopifnot(TYPE %in% c("density","box","violine"))

  if(is.null(colnames(M))) {colnames(M) = paste0("V_",seq_len(ncol(M)))}

  Plot_df <- melt(M)[,2:3]

  colnames(Plot_df) = c("sample","value")

  Plot_df$Group = rep(GROUP_LABEL,each = nrow(M))

  if(TYPE == "density"){

  p1 <- ggplot(Plot_df,aes(x = value, fill = sample)) + geom_density(alpha = .5,linetype = 0) + theme_classic() + labs(x = VALUE_LABEL, title = "stratified density plot", subtitle = HDER)

  figwidth = 5
  figheight = 2.6

  if(!is.null(GROUP_LABEL)) {
    stopifnot( is.vector( GROUP_LABEL ) & (length(GROUP_LABEL) ==  ncol(M)) )

    p1 <- p1 + facet_wrap( ~ Group,
                           ncol = min(length(unique(Plot_df$Group)), 5) )

    p1 <- p1 + theme(plot.margin = margin(1,5,1,5,"cm"),
                     legend.position = "bottom",
                     legend.direction = "horizontal",
                     legend.box = "vertical",
                     legend.justification = "center",
                     axis.text.x = element_text(angle = 310, vjust =.9, hjust = .1))

    figwidth = 8 + 1 * min(5,length(unique(Plot_df$Group)))

    figheight = 3.65 + 1* ceiling(length(unique(Plot_df$Group))/5)

  }

  ggsave(paste0(HDER,"_density.pdf"),p1,width = figwidth,height = figheight)

  }

  if(TYPE == "box"){

  if(is.null(GROUP_LABEL)) {
    Plot_df$label =  Plot_df$sample
  } else{

    stopifnot( is.vector( GROUP_LABEL ) & (length(GROUP_LABEL) ==  ncol(M)) )

    Plot_df$label =  rep(GROUP_LABEL,each = nrow(M))
  }

  p1 <- ggplot(Plot_df,aes(y = value, x = sample, fill = label)) + geom_boxplot(alpha = .5,linetype = 1) + theme_classic() + labs(y = VALUE_LABEL, x = "samples", title = "box plot of sample distribution", subtitle = HDER)

  p1  = p1 + theme(plot.margin = margin(1,5,1,5,"cm"),
                   legend.position = "bottom",
                   legend.direction = "horizontal",
                   legend.box = "vertical",
                   legend.justification = "center",
                   axis.text.x = element_text(angle = 310, vjust =.9, hjust = .1))

  figheight = 3.8 + .3* ceiling(length( unique(colnames(M)) )/4) + .05 * max(nchar(as.character( Plot_df$sample )))
  figwidth = 6.15 + .2 * length( unique(colnames(M)) )
  ggsave(paste0(HDER,"_boxplot.pdf"),p1,width = figwidth,height = figheight)
  }

  if (TYPE == "violine") {
    if(is.null(GROUP_LABEL)) {
      Plot_df$label =  Plot_df$sample
    } else{

      stopifnot( is.vector( GROUP_LABEL ) & (length(GROUP_LABEL) ==  ncol(M)) )

      Plot_df$label =  rep(GROUP_LABEL,each = nrow(M))
    }

    p1 <- ggplot(Plot_df,aes(y = value, x = sample, fill = label)) + geom_violin(alpha = .5,linetype = 0) + theme_classic() + labs(y = VALUE_LABEL, x = "samples", title = "violine plot of sample distribution", subtitle = HDER)

    p1  = p1 + theme(plot.margin = margin(1,5,1,5,"cm"),
                     legend.position = "bottom",
                     legend.direction = "horizontal",
                     legend.box = "vertical",
                     legend.justification = "center",
                     axis.text.x = element_text(angle = 310, vjust =.9, hjust = .1))

    figheight = 3.4 + .3* ceiling(length( unique(colnames(M)) )/4) + .05 * max(nchar(as.character( Plot_df$sample )))
    figwidth = 4 + .4 * length( unique(colnames(M)) )
    ggsave(paste0(HDER,"_violine.pdf"),p1,width = figwidth,height = figheight)
  }
}
