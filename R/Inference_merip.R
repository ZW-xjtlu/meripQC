#' @title Analysis MeRIP datasets with DESeq2.
#'
#' @description \code{Inference_merip} is a function to infer methylation and differential methylation given merip datasets.
#' @param SE_M A \code{SummarizedExperiment} object with 1 necessary column in \code{colData} named "IP_input", its content should be a character vector consists of 2 values: c("IP", "input").
#' @param DM_METHOD A character string indicating the statistical method used in differential methylation analysis, could be one in c("DESeq2","QNB").
#' @param PCA Wheather to save the PCA plot after rlog transformation, default is FALSE, the plot will not be generated when \code{DM_METHOD} = "QNB" while the \code{MODE} = "DM".
#' @param HDER What should be the header of the PCA plot, applied when \code{PCA} = TRUE.
#' @param MODE Could be either "Meth" or "DM", the later will conduct differential methylation analysis with the design:
#'
#' log2(Q) = intercept + I(Treated) + I(IP) + I(IP):I(Treated).
#'
#' The result is just the differences in conditioning effect, or the statistics for estimate before the term I(IP):I(Treated).
#'
#' We don't get this by contrast, because the information is already contained in coefficient estimate under the design above, and we don't need to linear combine (or a linear combination by t(c(0,0,0,1))) the estimates to get it.
#'
#' An additional column c("Perturbation") is necessary for this option, the Perturbation column has to include character "C" for control condition.
#'
#' @return A DESeq2 result object for DESeq2 analysis; for other analysis, it will generate a \code{data.frame} object.
#' @import DESeq2
#' @import QNB
#' @import exomePeak
#' @import ggplot2
#
#' @export

Inference_merip <- function(SE_M,
                             MODE = "Meth",
                              DM_METHOD = "DESeq2",
                              PCA = FALSE,
                             HDER = "Unknown",
                            ROW_FILTER = 0) {

SE_M$IPinput = SE_M$IP_input

if( MODE == "DM" & DM_METHOD == "QNB" ) ROW_FILTER = max(1,ROW_FILTER)

Omit_indx <- rowSums( assay(SE_M) ) < ROW_FILTER

SE_M <- SE_M[!Omit_indx,]

if(MODE == "Meth") {

  Cov = ~ IPinput
  dds = suppressMessages( DESeqDataSet(se = SE_M, design = Cov) )
  dds$IPinput <- relevel(dds$IPinput, "input")
  dds <- suppressMessages( DESeq(dds) )

  if(PCA) {
    INTGRP = c("IPinput", if(any(colnames(colData(SE_M)) == "Perturbation_detail")){"Perturbation_detail"}else{NULL})
    rld <- rlog(dds)
    ggsave(paste0(HDER,"_PCA.pdf"), plotPCA(rld,intgroup = INTGRP) + theme_classic() + theme() + scale_color_brewer(palette = "Spectral"), width = 5, height = 5)
  }

  inference_rst <- results(dds)

  } else {

    if (DM_METHOD == "DESeq2") {

  Pert_u = as.character(SE_M$Perturbation)
  Pert_u[Pert_u != "C"] = "Treated"
  SE_M$Perturbation = factor( Pert_u )
  Cov = ~ Perturbation + IPinput + Perturbation:IPinput
  dds = suppressMessages( DESeqDataSet(se = SE_M, design = Cov) )
  dds$IPinput <- relevel(dds$IPinput, "input")
  dds$Perturbation <- relevel(dds$Perturbation, "C")
  dds <- suppressMessages( DESeq(dds) )

  if(PCA) {
    INTGRP = c("IPinput","Perturbation")
    rld <- rlog(dds)
    ggsave(paste0(HDER,"_PCA.pdf"), plotPCA(rld,intgroup = INTGRP) + theme_classic() + theme() + scale_color_brewer(palette = "Spectral"), width = 5, height = 5)
  }

  inference_rst <- results(dds, name="PerturbationTreated.IPinputIP")

  }

  if (DM_METHOD == "QNB") {
    C_indx = SE_M$Perturbation == "C"
    input_indx = SE_M$IPinput == "input"

    C_IP = data.frame( assay( SE_M[,C_indx&!input_indx] ) )
    T_IP = data.frame( assay( SE_M[,!C_indx&!input_indx] ) )
    C_input = data.frame( assay( SE_M[,C_indx&input_indx] ) )
    T_input = data.frame( assay( SE_M[,!C_indx&input_indx] ) )

    #Here, it merges the samples where the copy numbers are redundent.

    Rep_Num <- c(ncol(C_IP),ncol(T_IP),ncol(C_input),ncol(T_input))

   if(any(Rep_Num != min(Rep_Num))) {

     combine_columns_to <- function(DF_x,column_target_num) {
       while(ncol(DF_x) > column_target_num) {
       DF_x[,1] = DF_x[,1] + DF_x[,2]
       DF_x = data.frame( cbind(DF_x[,-2]) )
       }
       return(DF_x)
     }

     C_IP = combine_columns_to(C_IP, min(Rep_Num))

     T_IP = combine_columns_to(C_IP, min(Rep_Num))

     C_input = combine_columns_to(C_input, min(Rep_Num))

     T_input = combine_columns_to(T_input, min(Rep_Num))

   }

    inference_rst = qnbtest(C_IP,T_IP,C_input,T_input,mode="auto")

    colnames(  inference_rst )[4] = "log2FoldChange"

  }
}

inference_rst <- as.data.frame(inference_rst)

inference_return <- data.frame(matrix(NA,ncol = ncol(inference_rst), nrow = length(Omit_indx)))

colnames(inference_return) = colnames(inference_rst)

inference_return[!Omit_indx,] <- inference_rst

return(inference_return)

}
