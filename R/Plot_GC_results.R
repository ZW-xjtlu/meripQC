#' @title Plot GC content bias plot based on DESeq2 result.
#'
#' @description \code{DESeq2_merip} is an internal function used to infer methylation and differential methylation given merip datasets.
#' @param DS_RES A \code{DESeqResults} object with additional factor collumn "Decision".
#' @param GC_IDX A numeric vector indicate the GC content of each individual feature.
#' @param HDER Determine the content of the title and the file name of the plot saved.
#'
#' @return A plot of the GC content saved under the current working directory.
#'
#
#' @export

Plot_GC_results <- function(DS_RES,GC_IDX,HDER) {

na_idx <- is.na( DS_RES$log2FoldChange ) | is.na(GC_IDX)

plot_df = data.frame(
  Log2FC = DS_RES$log2FoldChange[!na_idx],
  GC_idx = GC_IDX[!na_idx],
  Label = DS_RES$Decision[!na_idx]
                     )

p1 <- ggplot(plot_df, aes(x =  GC_idx , y = Log2FC )) + stat_density_2d(aes(fill = ..level..), geom = "polygon") + geom_smooth() + theme_classic() + scale_fill_gradient2() +
  labs(x = "GC contents", y = "Log 2 fold changes", title = "GC content bias plot of the inference", subtitle = HDER) + xlim(c(0.25,0.75)) + theme_classic()

p2 <- ggplot(plot_df, aes(x = GC_idx, fill = Label)) + geom_density(linetype = 0, alpha = .4) + theme_classic() + xlim(c(0.25,0.75)) + scale_fill_brewer(palette = "Dark2") + labs(x = "GC contents", title = "GC content distribution", subtitle = HDER)

ggsave(paste0(HDER,"_GC_bias.pdf"),p1,width = 5,height = 2.6)
ggsave(paste0(HDER,"_GC_dist.pdf"),p2,width = 5,height = 2.6)
}
