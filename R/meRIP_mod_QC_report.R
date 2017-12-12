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
           txdb,
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
           Expected_change = NULL) {
    #0. directory
    setwd(save_dir)

    #1. A reads count bar plot.
    Plot_Seq_depth(se_M,save_title)

    #2. A GC content diagnosis plot for single collumns of SummarizedExperiment.
    Plot_GC_collumns(se_M,GC_idx_feature,save_title)

    #3. A methylation profile report in tabular format based on DeSeq2 result.
    #Hence we need to write a function of DESeq2 now.

    #2. Run deseq2 to report methylation numbers/ratios
    Meth_count <- se_M[,c(IP_colnames,input_colnames)]

    col.data = data.frame(
      sample.id = c(IP_colnames,input_colnames),
      condition = rep(
        c("IP","Input"),
        c(length(IP_colnames),
          length(input_colnames))
      )
    )

    dds <- DESeqDataSetFromMatrix(Meth_count, col.data, design=~condition)

    dds$type <- factor(col.data$condition)

    dds$type <- relevel(dds$type, "Input") #Control as the denominator of the FC

    design(dds) <- ~ type

    dds <- suppressMessages( DESeq(dds) )

    #ggsave("MA-plot.pdf",plotMA(dds,ylim = c(-6,6)),width = 5, height = 4)

    dds_rst <- results(dds)

    egrid <- expand.grid(c("pvalue","padj"),c(".1",".05",".01",".001"))
    egrid <- egrid[order(egrid$Var1),]

    report_m <- Reduce("rbind",Map(function(x,y) dds_rst_stat(x,y,ds2_rslt = dds_rst),
                                   as.character(egrid$Var1),
                                   as.character(egrid$Var2)))
    rownames(report_m) = NULL

    write.table(report_m,"Report-pow.txt",sep = "\t",row.names = F)

    report_sd <- as.data.frame(rbind(colSums(Meth_count)))
    colnames(report_sd) = rep("input",ncol(report_sd))
    colnames(report_sd)[colnames(Meth_count) == IP_colnames] = "IP"

    write.table(report_sd,"Report-sequencing-depth.txt",sep = "\t",row.names = F)

    ##A code chunk decide the cut-off idx=================================================================================

    if(is.na(DeSeq2_fdr_threshold)) {
      filter_col = "pvalue"
      filter_cut = DeSeq2_p_threshold
    } else {
      filter_col = "padj"
      filter_cut = DeSeq2_fdr_threshold
    }

    idx_sig_M = which(!is.na(dds_rst$log2FoldChange) & dds_rst[[filter_col]] < filter_cut & dds_rst$log2FoldChange > log2FC_cutoff)

    if (length(idx_sig_M) < min_num_Meth) {
      idx_M <- which(dds_rst$log2FoldChange > log2FC_cutoff)
      ranks_stat <- rank(dds_rst[[filter_col]][idx_M])
      names(ranks_stat) = 1:length(ranks_stat)
      idx_sig_M <- idx_M[as.numeric(names(sort(ranks_stat,decreasing = F))[1:min_num_Meth])]
    }

    idx_sig_N = which(!is.na(dds_rst$log2FoldChange) & dds_rst[[filter_col]] < filter_cut & dds_rst$log2FoldChange < -1*log2FC_cutoff)

    if (length(idx_sig_N) < min_num_Nonsence) {
      idx_N <- which(dds_rst$log2FoldChange < -1*log2FC_cutoff)
      ranks_stat <- rank(dds_rst[[filter_col]][idx_N])
      names(ranks_stat) = 1:length(ranks_stat)
      idx_sig_N <- idx_N[as.numeric(names(sort(ranks_stat,decreasing = F))[1:min_num_Nonsence])]
    }

    ##===================================================================================================================

    if(Return_FC_idx){
      return_vec = rep(NA,nrow(dds_rst))
      return_vec[idx_sig_M] = dds_rst$log2FoldChange[idx_sig_M]
      saveRDS(return_vec,"Meth.rds")
    }

    Plot_df2 <- data.frame(
      log2.FC.Deseq2 = c(dds_rst$log2FoldChange[idx_sig_M],
                         dds_rst$log2FoldChange[idx_sig_N],
                         dds_rst$log2FoldChange),

      Group = rep(
        c("Methylated",
          "Nonsense",
          "All"),
        c(length(idx_sig_M),
          length(idx_sig_N),
          length(dds_rst$log2FoldChange))
      )
    )

    Plot_es <- ggplot(Plot_df2,aes(x = log2.FC.Deseq2)) + geom_density(aes(colour = Group,fill = Group), alpha = .4) + labs(title = "distribution of methylation level",subtitle = "in DeSeq2 log2 FC")

    ggsave("DeSeq2_es_plot.pdf",Plot_es,width = 4.8, height = 2.8)


    #2. Exon lengths distribution

    List_GR = list(
      Methylated = gr[idx_sig_M],
      Nonsense = gr[idx_sig_N]
    )

    Overlapped_exonLength_distribution <- function(Meth_list,TXDB,plot_title,plot_subtitle){
      All_exons = exons(TXDB)
      length_lst <- lapply(Meth_list,function(grl_x) {
        fol <- findOverlaps(grl_x,All_exons)
        width(All_exons)[subjectHits(fol)]
      }

      )

      names(length_lst) = paste0(names(length_lst),".")
      vec_length = unlist(length_lst)
      plot_df = data.frame(log2_exon_lengths = log2(vec_length),
                           Class = gsub("\\..*","",names(vec_length)))
      ggplot(plot_df,aes(x = log2_exon_lengths)) + geom_density(aes(color = Class,fill = Class), alpha = .5) + theme_bw() + labs(title = plot_title, subtitle = plot_subtitle)
    }

    Plot_ex <- Overlapped_exonLength_distribution(List_GR,
                                                  TXDB = txdb,
                                                  plot_title = "Distribution of the lengths of overlapped exons",
                                                  plot_subtitle = paste0("M sites DeSeq2 using ",filter_col," < ",filter_cut ))
    ggsave("exons_length_plot.pdf",Plot_ex,width = 4.8, height = 2.8)

    if(is.null(gtcoord)){
      gtcoord = makeGuitarCoordsFromTxDb(txdb)
    }

    #3. Guitar distribution

    capture.output( GuitarPlot(
      List_GR,
      saveToPDFprefix = "GuitarPlot",
      includeNeighborDNA = FALSE,
      GuitarCoordsFromTxDb = gtcoord
    ) )

    #4. 2ndary structure

    Secondary_structure_enrichment <- function(List_GR,Structure_grl,exon_gr,binsize = 51, plot_title, plot_subtitle) {
      Struc_Prop_LS <- lapply(List_GR,function(gr_x){
        x_sub <- suppressWarnings( subsetByOverlaps(gr_x,exon_gr) )
        x_sub <- resize(x_sub,width = binsize,fix = "center")
        dummy_x <- countOverlaps(x_sub, Structure_grl) > 0
        btest <- binom.test(sum(dummy_x),length(dummy_x))
        data.frame(
          Struc_proportion = btest$estimate,
          Conf_Int_lower = btest$conf.int[1],
          Conf_Int_upper = btest$conf.int[2]
        )
      })
      Plot_df <- Reduce("rbind",Struc_Prop_LS)
      Plot_df$Group = names(List_GR)
      ggplot(Plot_df,aes(x = Group,y = Struc_proportion)) + geom_bar(stat = "identity",aes(fill = Group), width = 0.5) + theme_bw() + geom_errorbar(aes(ymin = Conf_Int_lower, ymax = Conf_Int_upper),position = "dodge", width = 0.25) + scale_fill_brewer(palette="Accent")  + labs(title = plot_title, subtitle = plot_subtitle)
    }

    if(is.null(struc_gr)){} else{
      Plot_2nd <- Secondary_structure_enrichment(List_GR,struc_gr,exons(txdb),51,"Proportion of sites mapped to RNA MEA 2ndary Structures","bin size = 51")
      ggsave("2ndary_structure.pdf",Plot_2nd,width = 4, height = 2.8)
    }

    #5. PhastCons score
    if(is.null(pcdb)) {} else {

      PhastCons_distribution <- function(Meth_list,PCDB,bin_size,plot_title,plot_subtitle) {

        Scores_lst <- suppressWarnings( lapply(Meth_list,function(grl_x) {
          gr_x <-  keepStandardChromosomes( unlist(grl_x) )
          gr_x <- resize(gr_x,bin_size,fix = "center")
          Phascores <- scores(PCDB,gr_x)
          return(Phascores)
        }
        )
        )

        names(Scores_lst) = paste0(names(Scores_lst),".")
        vec_Scores = unlist(Scores_lst)
        plot_df = data.frame(PhastCons_scores = vec_Scores,
                             Class = gsub("\\..*","",names(vec_Scores)))
        ggplot(plot_df,aes(x = PhastCons_scores)) + geom_density(aes(color = Class,fill = Class,linetype = Class), alpha = .5) + scale_fill_brewer(palette="Accent") + scale_color_brewer(palette="Accent") + theme_bw() + labs(title = plot_title, subtitle = plot_subtitle)

      }

      Pc_dis <- PhastCons_distribution(List_GR,
                                       PCDB = pcdb,
                                       bin_size = 201,
                                       plot_title = "Distribution of the phastCons scores",
                                       plot_subtitle = "DeSeq2 result, bin size = 201")

      ggsave("PhastCons_plot.pdf",Pc_dis,width = 4.8, height = 2.8)
    }

    #6. Topological clustering effect

    Plot_clustering <- function(List_GR, plot_title, plot_subtitle) {
      Ls_Clustering <- lapply(List_GR,function(x){
        x <- resize(x,201,fix = "center")
        countOverlaps(x,x)
      })

      Plot_df = data.frame(Neighbors_num = Reduce("c",Ls_Clustering),
                           Group = rep(names(List_GR),sapply(List_GR,length)))
      #ggplot(Plot_df, aes(x = Neighbors_num)) + geom_density(aes(colour = Group,fill = Group),alpha = .5)
      ggplot(Plot_df, aes(x = Group, y = Neighbors_num)) + geom_boxplot(aes(colour = Group,fill = Group),alpha = .5) + theme_bw() + labs(title = plot_title, subtitle = plot_subtitle)
    }

    Plot_cluster <- Plot_clustering(List_GR,"Boxplot of the number of neigbors in each site","Bin size = 201")
    ggsave("Clustering.pdf",Plot_cluster,width = 4.8, height = 2.8)
  }
