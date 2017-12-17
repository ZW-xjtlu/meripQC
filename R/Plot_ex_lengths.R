#' @title Plot exon length distribution plot based on the DESeq2 result.
#'
#' @description \code{DESeq2_merip} is an internal function used to infer methylation and differential methylation given merip datasets.
#' @param LIDT_GR A \code{list} of \code{GRanges}.
#' @param TXDB A \code{txdb} object.
#' @param HDER Determine the content of the title and the file name of the plot saved.
#'
#' @return A plot of the exon length distribution saved under the current working directory.
#'
#
#' @export

Plot_ex_lengths <- function(LIST_GR,TXDB,HDER){
ex_txdb <- exons(TXDB)
LST_exlengths <- lapply(LIST_GR, function(x) width( subsetByOverlaps( ex_txdb, x) ) )
plot_df <- data.frame(log2_exon_lengths = log2( unlist(LST_exlengths,use.names = FALSE) ),
           Group = rep(names(LST_exlengths),sapply(LST_exlengths,length)))
p1 <- ggplot(plot_df,aes(x = log2_exon_lengths, fill = Group)) + geom_density(alpha = .5,linetype = 0) + theme_classic() + labs(x = "log2 exon lengths", title = "exon length distribution", subtitle = HDER)
ggsave(paste0(HDER,"_exl_dist.pdf"),p1,width = 5,height = 2.6)
}
