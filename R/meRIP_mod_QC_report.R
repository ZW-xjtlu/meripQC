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
#' @param p_threshold A numeric value between 0 to 1, it indicates the p value cut off of the statistical inference, it will be neglected if \code{fdr_threshold} is not NULL.
#' @param fdr_threshold A numeric value between 0 to 1, it indicates the fdr cut off of the statistical inference.
#'
#' By default, \code{meRIP_mod_QC_report} want to call DESeq2 and infer methylation under the design log2(Q) ~ intercept + I(IP).
#' The Wald test is conducted on the coefficient estimate of the second term I(IP).
#'
#' @param log2FC_cutoff The log2 fold change cutoff of the inference result, default setting is 0.
#' @param min_num_Mod The minimal number of sites inferred in the Methylation and Control groups, i.e.IP bigger than input and vice versa (for control), default setting is 10000.
#' @param Save_inference_result Whether to save the result of the inference, default setting is TRUE.
#' @param GC_idx_feature Optional: The GC content values for each features (rows) of the count matrix.
#' @param DM_analysis Optional: Whether to conduct differential methylation analysis or not, default setting is FALSE.
#' @param DM_method Decide the statistical inference method used in differential methylation procedure. The default setting is "DESeq2"; an alternative setting is "QNB", which will use the \pkg{QNB} package to compute the differential methylation statistics.
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
#' #To do:
#' 1. add QNB
#' 2. add cqn (adjust GC content) / probably add GC content adjustment for CHIP-seq (if possible).
#' 3. add plot over-dispersion for both QNB and DESeq2.
#' 4. change the save dir into paste, or record the original dir.
#'
#' @import DESeq2
#' @import ggplot2
#' @import Guitar
#' @import GenomicFeatures
#' @import SummarizedExperiment
#' @import QNB
#' @export

meRIP_mod_QC_report <-
  function(se_M,
           txdb = NULL,
           save_title = "modX",
           save_dir = save_title,
           gtcoord = NULL,
           p_threshold = NULL,
           fdr_threshold = NULL,
           log2FC_cutoff = 0,
           min_num_Mod = 10000,
           Save_inference_result = TRUE,
           GC_idx_feature = NULL,
           DM_analysis = FALSE,
           DM_method = "DESeq2",
           Expected_change = NULL,
           PCA_PLOT = FALSE) {
    #0. directory
    dir_org = getwd()
    if(!dir.exists(save_dir)) dir.create(save_dir)
    setwd(save_dir)

    #1. A reads count bar plot.
    Plot_Seq_depth(se_M,save_title)

    #2. A GC content diagnosis plot for single collumns of SummarizedExperiment.
    if(!is.null(GC_idx_feature)) Plot_GC_collumns(se_M,GC_idx_feature,save_title,DM_analysis)

    #3. A methylation profile report in tabular format based on DeSeq2 result.
    # Run Inference.
    inf_result <- Inference_merip(se_M,MODE = ifelse(DM_analysis,"DM","Meth"),DM_METHOD = DM_method,PCA = PCA_PLOT,HDER = save_title)

    # Analysis Inference result and generate a decision table:
    Dcs_tb <- Decision_infresult(inf_result,log2FC_cutoff,p_threshold,fdr_threshold,min_num_Mod,DM_analysis,Expected_change,save_title)
    inf_result$Decision = "Insig"

    #Make decisions based on the decision table.
    if (Dcs_tb$Expected_dir == "< 0"){
      Control_index =  (inf_result[[as.character(Dcs_tb$Cut_By_ctrl)]] < Dcs_tb$Cut_Val_ctrl) & (inf_result$log2FoldChange > Dcs_tb$log2FC_cut)
      Expected_index = (inf_result[[as.character(Dcs_tb$Cut_By_expected)]] < Dcs_tb$Cut_Val_expected) & (inf_result$log2FoldChange < Dcs_tb$log2FC_cut)
    } else {
      Control_index =  (inf_result[[as.character(Dcs_tb$Cut_By_ctrl)]] < Dcs_tb$Cut_Val_ctrl) & (inf_result$log2FoldChange < Dcs_tb$log2FC_cut)
      Expected_index = (inf_result[[as.character(Dcs_tb$Cut_By_expected)]] < Dcs_tb$Cut_Val_expected) & (inf_result$log2FoldChange > Dcs_tb$log2FC_cut)
    }

    if(is.null(Expected_change) & DM_analysis){
    inf_result$Decision[Control_index] = "Hypo-Meth"
    inf_result$Decision[Expected_index] = "Hyper-Meth"
    }

    inf_result$Decision[Control_index] = "Control"
    inf_result$Decision[Expected_index] = ifelse(!is.null(Expected_change),"Targeted","Methylated")

    #Save the result of the decided sites.
    write.csv(Dcs_tb, paste0(save_title,"_Dcs_tb.csv"))

    if(Save_inference_result) saveRDS(inf_result, paste0(save_title,"_inf_result.rds"))

    #4. A GC content diagnosis plot for inference.
    if(!is.null(GC_idx_feature)) Plot_GC_results(inf_result,GC_idx_feature,save_title)

    #5. Guitar plot for methylation sites.
    Plot_ls_Gr <- as.list(split(rowRanges(se_M),inf_result$Decision))
    Plot_ls_Gr = Plot_ls_Gr[ names(Plot_ls_Gr) != "Insig" ]

    if(!(is.null(gtcoord) & is.null(txdb))) {
      if(is.null(gtcoord)) {gtcoord <- Guitar::makeGuitarCoordsFromTxDb(txdb)}
      capture.output( suppressWarnings( Guitar::GuitarPlot(Plot_ls_Gr,gtcoord,saveToPDFprefix = paste0(save_title,"_guitar.rds")) ) )
    }

    #6. Exon lengths distribution
    if(!is.null(txdb)) Plot_ex_lengths(Plot_ls_Gr,txdb,save_title)

  setwd(dir_org)
}
