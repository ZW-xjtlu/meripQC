#' @title Plot GC content distribution of features by columns.
#'
#' @description \code{Plot_Seq_depth} is an internal function used to generate a GC content distribution plot given an input \code{SummarizedExperiment} object.
#'
#' @param SE_M A \code{SummarizedExperiment} object with 2 necessary columns in \code{colData}: c( "IP_input", "SRR_RUN" ).
#' @param GC_IDX A numeric vector indicate the GC content of each individual feature.
#' @param HDER Determine the content of the title and the file name of the plot saved.
#' @param ANNOT Whether annotate the plot with additional batch or treatment information,
#' in this version, it depends on the column of c("Perturbation_detail") in the \code{colData}.
#'
#' @return A pdf diagram saved under the current working directory.
#'
#' @export

Plot_GC_columns <- function(SE_M,GC_IDX,HDER,ANNOT = F) {

Plot_df  <- reshape2::melt(assay(SE_M))
Plot_df$GC_cont = GC_IDX[Plot_df$Var1]
Plot_df$Bam_names = colData(SE_M)[Plot_df$Var2,"SRR_RUN"]
Plot_df$IP_input = colData(SE_M)[Plot_df$Var2,"IP_input"]
if(ANNOT) Plot_df$Treatment = colData(SE_M)[Plot_df$Var2,"Perturbation_detail"]
Plot_df = na.omit(Plot_df)


if(!ANNOT){

p1 <- ggplot(Plot_df,aes(x = GC_cont, y = value, colour = Bam_names, linetype = IP_input))+ geom_smooth(alpha = .09,size = .9) + theme_classic()  + labs(x = "GC content", y = "Reads count per feature", title = "GC content diagnosis plot for samples")
if (ncol(SE_M) <= 11)  p1 = p1 + scale_color_brewer(palette = "Spectral")

}else{

Plot_df = rbind(Plot_df,Plot_df)
p1 <- ggplot(Plot_df,aes(x = GC_cont, y = value, colour = Treatment, linetype = IP_input))+ geom_smooth(alpha = .15,size = .9) + theme_classic() + scale_color_brewer(palette = "Dark2")  + labs(x = "GC content", y = "Reads count per feature", title = "GC content diagnosis plot for treatments")

}

if(!ANNOT){
figheight = 4.15 + .2* round(length( unique( colData(SE_M)[,"SRR_RUN"] ) )/4)
} else {
figheight = 4.15 + .2* round(length( unique( colData(SE_M)[,"Perturbation_detail"] ) )/4)
}

p1 = p1 + theme(plot.margin = margin(1,5,1,5,"cm"),legend.position = "bottom",legend.direction = "horizontal", legend.box = "vertical", legend.justification = "center")
suppressMessages( ggsave(paste0(HDER,"_GC_diagnosis_col.pdf"),p1,width = 7.5,height = figheight) )

}



