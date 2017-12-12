#' @title Bar Plot of total sequencing depth.
#'
#' @description \code{Plot_Seq_depth} is an internal function used to generate a sequencing depth plot given an input SummarizedExperiment object.
#'
#' @param SE_M A \code{SummarizedExperiment} object with 2 necessary collumns in \code{colData}: c( "IP_input", "SRR_RUN" ).
#' @param HDER Determine the content of the title and the file name of the plot saved.
#'
#' @return A pdf diagram saved under the current working directory.
#'
#
#' @export


Plot_Seq_depth <- function(SE_M,HDER) {

Plot_df = data.frame(Sum_count = colSums(assay(SE_M)),
                     IP_input = colData(SE_M)[,"IP_input"],
                     Bam_names = colData(SE_M)[,"SRR_RUN"])

p1 <- ggplot(Plot_df,aes(x = Bam_names, y = Sum_count,fill = IP_input)) + geom_bar(stat = "identity",width = .55) +
  theme_classic() + labs(x = "Bam files", y = "Total read counts", title = paste0("Sequencing depth of ",HDER)) +
  theme(axis.text.x = element_text(angle = 310, vjust =.9, hjust = .1)) + scale_fill_brewer(palette = "Dark2")

ggsave(paste0(HDER,"_Seq_depth.pdf"),p1,width = 2.3 + ncol(SE_M)/6,height = 3.3)

}
