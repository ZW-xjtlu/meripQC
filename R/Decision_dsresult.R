#' @title Calculate the decision table for a DESeq2 result.
#'
#' @description \code{Decision_dsresult} is an internal function used to summary the cut-off and the number of positive results used for DESeq2 result..
#'
#' @param DS_RES a \code{DESeqResults} object that for either methylation or differential methylation.
#' @param log2FC_cut The log2 fold change cutoff of the DESeq2 result, default setting is 0.
#' @param P_cut A numeric value between 0 to 1, indicating the p value cut off of the Wald test defined by DESeq2, it will be neglected if \code{Padj_cut} is not NULL.
#' @param Padj_cut A numeric value between 0 to 1, indicating the fdr cut off of the Wald test defined by DESeq2.
#' @param Min_mod Minimum number of features returned, when this is smaller than the cut-off results, additional features are called by the order of p values.
#' @param DM Whether the result is differential methylation analysis, default is FALSE.
#' @param Exp_dir This parameter is filled when DM = TRUE, can be either "hyper" or "hypo".
#' @param HDER Determine the ID collumn of the generated table.
#'
#' @return A \code{data.frame} object indicating the collumn and cut-off value used for desicion, also it includes the number of positive sites in both directions based on the decision.
#
#' @export

Decision_dsresult <- function(DS_RES,log2FC_cut = 0,P_cut = 0.5,Padj_cut = NULL,Min_mod = 10000,DM = FALSE,Exp_dir = NULL,HDER){

DS_RES <- DS_RES[!(is.na(DS_RES$log2FoldChange) | is.na(DS_RES$pvalue)),]

DS_RES <- DS_RES[abs(DS_RES$log2FoldChange) > log2FC_cut,]

result_df <- data.frame(
                        ID = HDER,
                        log2FC_cut = log2FC_cut,
                        Cut_By_ctrl = "pvalue" ,
                        Cut_By_expected = "pvalue" ,
                        Cut_Val_ctrl = 0 ,
                        Cut_Val_expected = 0 ,
                        Discoveries_ctrl = 0 ,
                        Discoveries_expected = 0,
                        Expected_dir = "> 0"
                        )

if(!is.null(Exp_dir)) if(Exp_dir=="hypo") result_df$Expected_dir = "< 0"

#Define index of the expected direction.
 EXP_idx =  eval(parse(text = paste0( "DS_RES$log2FoldChange ",result_df$Expected_dir) ))

#First lt's do regular cut-off

if( !is.null( Padj_cut ) ) {
DS_RES$padj[is.na(DS_RES$padj)] = 1
result_df$Cut_By_ctrl = "padj"
result_df$Cut_By_expected = "padj"
result_df$Cut_Val_ctrl = Padj_cut
result_df$Cut_Val_expected = Padj_cut
result_df$Discoveries_ctrl = sum(!EXP_idx & DS_RES$padj < Padj_cut)
result_df$Discoveries_expected = sum(EXP_idx & DS_RES$padj < Padj_cut)
} else {
result_df$Cut_Val_ctrl = P_cut
result_df$Cut_Val_expected = P_cut
result_df$Discoveries_ctrl = sum(!EXP_idx & DS_RES$pvalue < P_cut)
result_df$Discoveries_expected = sum(EXP_idx & DS_RES$pvalue < P_cut)
}

#Second, we will consider the case while the minimum number is not met.
if(result_df$Discoveries_ctrl < Min_mod) {
  result_df$Cut_By_ctrl = "pvalue"
  result_df$Cut_Val_ctrl = sort( DS_RES$pvalue[!EXP_idx])[Min_mod+1]
  result_df$Discoveries_ctrl = sum(DS_RES$pvalue[!EXP_idx] < result_df$Cut_Val_ctrl)
}

if(result_df$Discoveries_expected < Min_mod) {
  result_df$Cut_By_expected = "pvalue"
  result_df$Cut_Val_expected = sort( DS_RES$pvalue[EXP_idx])[Min_mod+1]
  result_df$Discoveries_expected = sum(DS_RES$pvalue[EXP_idx] < result_df$Cut_Val_expected)
}

return(result_df)
}
