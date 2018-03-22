#' @title Plot the collumn wise data distribution.
#'
#' @description \code{DESeq2_merip} is a function used to plot the collumned wise data distribution.
#' @param M A \code{matrix}.
#' @param TYPE Can be either "density" or "box".
#' @param HDER The subtitle and the file name of the plot.
#'
#' @return A plot of the collumn wised density/boxplot.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export

Plot_collumn_dist <- function(M,TYPE = "box",HDER){
  stopifnot(TYPE %in% c("density","box"))
  Plot_df <- reshape2::melt(M)[,2:3]
  colnames(Plot_df) = c("sample","value")
  if(TYPE == "density"){
  p1 <- ggplot(Plot_df,aes(x = value, fill = sample)) + geom_density(alpha = .5,linetype = 0) + theme_classic() + labs(x = "value", title = "stratified density plot", subtitle = HDER)
  ggsave(paste0(HDER,"_density.pdf"),p1,width = 5,height = 2.6)
  }
  if(TYPE == "box"){
  p1 <- ggplot(Plot_df,aes(y = value, x = sample, fill = sample)) + geom_boxplot(alpha = .5,linetype = 1) + theme_classic() + labs(x = "samples", title = "box plot of sample distribution", subtitle = HDER)
  p1  = p1 + theme(plot.margin = margin(1,5,1,5,"cm"),
                   legend.position = "bottom",
                   legend.direction = "horizontal",
                   legend.box = "vertical",
                   legend.justification = "center",
                   axis.text.x = element_text(angle = 310, vjust =.9, hjust = .1))
  figheight = 4.15 + .2* round(length( unique(colnames(M)) )/4)
  figwidth = 6.15 + .2 * length( unique(colnames(M)) )
  ggsave(paste0(HDER,"_boxplot.pdf"),p1,width = figwidth,height = figheight)
  }
}
