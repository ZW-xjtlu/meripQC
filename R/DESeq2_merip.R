#' @title Analysis MeRIP datasets with DESeq2.
#'
#' @description \code{DESeq2_merip} is an internal function used to infer methylation and differential methylation given merip datasets.
#' @param SE_M A \code{SummarizedExperiment} object with 1 necessary collumn in \code{colData}: c("IP_input").
#' @param MODE Could be either "Meth" or "DM", the later will conduct differential methylation analysis with the design:
#' @param PCA Wheather to plot a PCA after an rlog transformation, default is FALSE.
#' @param HDER What should be the header of the PCA plot, applied when \code{PCA} = TRUE.
#'
#' log2(Q) = intercept + I(Treated) + I(IP) + I(IP):I(Treated).
#'
#' The result is just the differences in conditioning effect, or the statistics for estimate before the term I(IP):I(Treated).
#'
#' We don't get this by contrast, because the information is already contained in coefficient estimate under the design above, and we don't need to linear combine (or a linear combination by t(c(0,0,0,1))) the estimates to get it.
#'
#' An additional collumn c("Perturbation") is necessary for this option, the Perturbation collumn has to include character "C" for control condition.
#'
#' @return A DESeq2 result object.
#'
#
#' @export

DESeq2_merip <- function(SE_M,MODE = "Meth",PCA = FALSE,HDER = "Unknown"){
SE_M$IPinput = SE_M$IP_input
if(MODE == "Meth") {
  Cov = ~ IPinput
  } else {
  Pert_u = as.character(SE_M$Perturbation)
  Pert_u[Pert_u!="C"] = "Treated"
  SE_M$Perturbation = factor( Pert_u )
  Cov = ~ Perturbation + IPinput + Perturbation:IPinput
  }
dds = suppressMessages( DESeqDataSet(se = SE_M, design = Cov) )
dds$IPinput <- relevel(dds$IPinput, "input")
if(MODE == "DM") {
  dds$Perturbation <- relevel(dds$Perturbation, "C")
}
dds <- suppressMessages( DESeq(dds) )

if(PCA) {
if(MODE == "DM") {INTGRP = c("IPinput","Perturbation")} else {INTGRP = ("IPinput")}
rld <- rlog(dds)
ggsave(paste0(HDER,"_PCA.pdf"), plotPCA(rld,intgroup = INTGRP) + theme_classic() + scale_color_brewer(palette = "Spectral"), width = 5, height = 3)
}

if(MODE == "Meth"){
dds_rst <- results(dds)
} else {
dds_rst <- results(dds, name="PerturbationTreated.IPinputIP")
}
return(dds_rst)
}
