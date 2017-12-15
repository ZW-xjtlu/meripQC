#' @title Generate quality control report of a single MeRIP data site.
#'
#' @description \code{meRIP_mod_QC_report} is used to generate a single quality control report for a summarized experiment object of MeRIP experiment.
#'
#' @details This function is an internal function, and it defines the behavior of a single QC report on a well formated summarized experiment object of count.
#' Under current version, \code{meRIP_mod_QC_report} supports the generation of the following reports.
#'
#' 1. A reads number distribution plot.
#'
#' 2. A GC content diagnosis plot for single collumns of SummarizedExperiment.
#'
#' 3. A methylation profile report in tabular format based on DeSeq2 result.
#'
#' 4. A GC content diagnosis plot for methylation sites.
#'
#' 5. Guitar plot for methylation sites.
#'
#' 6. Exon length distribution for methylation sites.
#'
#'
#' @seealso This function is called by \code{meRIP_QC}
#'
#' @param se_M A \code{SummarizedExperiment} object containing the counts of each modification sites of each bam files. Appropriate \code{colData} and \code{rowRanges} should be available.
#' Specifically, \code{colData} should be a \code{DataFrame} object including the following collumns:
#'
#' \code{SRR_RUN} : a factor variable that uniquely indentify each collumns of the count matrix, could be ID for each bam files.
#'
#' \code{IP_input} : a factor variable indicating whether the collumns belong to IP or input, the levels need to be c("input", "IP").
#'
#' @param txdb \code{TxDb} object of the corresponding \code{rowRanges}, this is either obtained from biocoductor or converted from the user provided GFF files.
#' @param gtcoord A variable containing guitar coordinate, which is defined by the \code{Guitar} package.
#' @param save_dir A character string indicating the directory to save the report, by default it is the current working directory.
#' @param save_title A character string indicating the header of the reports generated.
#' @param DeSeq2_p_threshold A numeric value between 0 to 1, indicating the p value cut off of the Wald test defined by DESeq2, it will be neglected if \code{DeSeq2_fdr_threshold} is not NULL.
#' @param DeSeq2_fdr_threshold A numeric value between 0 to 1, indicating the fdr cut off of the Wald test defined by DESeq2.
#'
#' Specifically, \code{meRIP_mod_QC_report} want to call DESeq2 and infer methylation under the design log2(Q) ~ intercept + I(IP).
#' The Wald test is conducted on the coefficient estimate of the second term I(IP).
#'
#' @param log2FC_cutoff The log2 fold change cutoff of the DESeq2 result, default setting is 0.
#' @param min_num_Mod The minimal number of sites inferred in the Methylation and Control groups, IP bigger than input and vice versa (for control), default is 10000.
#' @param Save_DESeq2_result Whether to save the DESeq2 result, default setting is TRUE.
#' @param GC_idx_feature Optional: The GC content values for each features (rows) of the count matrix.
#' @param DM_analysis Optional: Whether to conduct differential methylation analysis or not, default setting is FALSE.
#' @param Expected_change Optional: could be either "hyper" and "hypo", indicating the expected change of treated condition over input condition,
#' @param PCA_PLOT Whether to plot the PCA plot in DESeq2, default setting is FALSE, it is slow because it requires rlog transformation.
#' this is useful when inferring the target sites of RNA modification writers or erasers from the MeRIP Seq data. Default setting is NULL.
#'
#' @return This function will write files and plots under the directory provided by \code{save_dir}
#'
#' @examples
#' meRIP_mod_QC_report(se_M = se_mm10,
#' txdb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
#' gtcoord = Gtcoord_mm10,
#' min_num_Mod = 1000)
#'
#' @import DESeq2
#' @import ggplot2
#' @import Guitar
#' @import GenomicFeatures
#' @import SummarizedExperiment
#' @export

meRIP_mod_QC_report <-
  function(se_M,
           txdb = NULL,
           save_title = "modX",
           save_dir = getwd(),
           gtcoord = NULL,
           DeSeq2_p_threshold = NULL,
           DeSeq2_fdr_threshold = NULL,
           log2FC_cutoff = 0,
           min_num_Mod = 10000,
           Save_DESeq2_result = TRUE,
           GC_idx_feature = NULL,
           DM_analysis = FALSE,
           Expected_change = NULL,
           PCA_PLOT = FALSE) {
    #0. directory
    setwd(save_dir)

    #1. A reads count bar plot.
    Plot_Seq_depth(se_M,save_title)

    #2. A GC content diagnosis plot for single collumns of SummarizedExperiment.
    if(!is.null(GC_idx_feature)) Plot_GC_collumns(se_M,GC_idx_feature,save_title,DM_analysis)

    #3. A methylation profile report in tabular format based on DeSeq2 result.
    # Run DeSeq2.
    ds_result <- DESeq2_merip(se_M,MODE = ifelse(DM_analysis,"DM","Meth"),PCA = PCA_PLOT,HDER = save_title)

    # Analysis DESeq2 result and generate a decision table:
    Dcs_tb <- Decision_dsresult(ds_result,log2FC_cutoff,DeSeq2_p_threshold,DeSeq2_fdr_threshold,min_num_Mod,DM_analysis,Expected_change,save_title)
    ds_result$Decision = "Insig"

    #Make decisions based on the decision table.
    Control_index =  (ds_result[[Dcs_tb$Cut_By_ctrl]] < Dcs_tb$Cut_Val_ctrl) & (ds_result$log2FoldChange < Dcs_tb$log2FC_cut)
    Expected_index = (ds_result[[Dcs_tb$Cut_By_expected]] < Dcs_tb$Cut_Val_expected) & (ds_result$log2FoldChange > Dcs_tb$log2FC_cut)

    if (Dcs_tb$Expected_dir == "< 0"){
    tmp = Control_index
    Control_index = Expected_index
    Expected_index = tmp
    tmp = NULL
    }

    if(is.null(Expected_change) & DM_analysis){
    ds_result$Decision[Control_index] = "Hypo-Meth"
    ds_result$Decision[Expected_index] = "Hyper-Meth"
    }

    ds_result$Decision[Control_index] = "Control"
    ds_result$Decision[Expected_index] = ifelse(!is.null(Expected_change),"Targeted","Methylated")

    #Save the result of the decided sites.
    write.csv(Dcs_tb, paste0(save_title,"_Dcs_tb.csv"))

    if(Save_DESeq2_result) saveRDS(ds_result, paste0(save_title,"_ds_result.rds"))

    #4. A GC content diagnosis plot for inference.
    if(!is.null(GC_idx_feature)) Plot_GC_results(se_M,ds_result,GC_idx_feature,save_title)

    #5. Guitar plot for methylation sites.
    Plot_ls_Gr <- as.list(split(rowRanges(se_M),ds_result$Decision))
    Plot_ls_Gr = Plot_ls_Gr[ names(Plot_ls_Gr) != "Insig" ]

    if(!(is.null(gtcoord) & is.null(txdb))) {
      if(is.null(gtcoord)) {gtcoord <- Guitar::makeGuitarCoordsFromTxDb(txdb)}
      capture.output( suppressWarnings( Guitar::GuitarPlot(Plot_ls_Gr,gtcoord,saveToPDFprefix = paste0(save_title,"_guitar.rds")) ) )
    }

    #6. Exon lengths distribution
    if(!is.null(txdb)) Plot_ex_lengths(Plot_ls_Gr,txdb,save_title)
}
