library(meripQC)
library(SummarizedExperiment)
Gtcoord_mm10 <- readRDS("/Users/zhenwei/Datasets/Gtcoords/Gtcoord_mm10.rds")
se_mm10 <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/B_COUNT_2017_12_5/se_mm10.rds")
Gene_GC_mm10 <- readRDS("/Users/zhenwei/Datasets/GC_content_Genes/Gene_GC_mm10.rds")
fol <- findOverlaps( rowRanges( se_mm10 ), Gene_GC_mm10 )
GC_cont = rep(NA,nrow(se_mm10))
GC_cont[queryHits(fol)] = mcols(Gene_GC_mm10)[subjectHits(fol),]

se_fto <- se_mm10[,grepl("midbrain", colData( se_mm10 )$Publication )]

set.seed(1)

indx <- sample.int(nrow(se_fto),10000,replace = T)


meRIP_QC_report(
                se_M = se_fto[indx,],
                txdb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
                             gtcoord = Gtcoord_mm10,
                             min_num_mod = 1000,
                             save_title = "FTO-3T3L1",
                             DM_analysis = T,
                             expected_change = "hyper",
                             fdr_threshold = .05,
                             GC_idx_feature = GC_cont[indx],
                DM_method = "DESeq2",
                row_minimal_counts = 10
                )

meRIP_QC_report(
  se_M = se_fto[indx,],
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
  gtcoord = Gtcoord_mm10,
  min_num_mod = 1000,
  save_title = "FTO-3T3L1-cqn",
  DM_analysis = T,
  expected_change = "hyper",
  fdr_threshold = .05,
  GC_idx_feature = GC_cont[indx],
  DM_method = "DESeq2",
  row_minimal_counts = 10,
  cqn = T
)


#Check weather the reported RDS is still well linked between granges and decisions.
Inf_3T3L1 <- readRDS("/Users/zhenwei/Documents/GitHub/meripQC/FTO-3T3L1/FTO-3T3L1_inf_result.rds")

Plot_ex_lengths( as.list( split( rowRanges(se_fto)[indx], Inf_3T3L1$Decision) ) ,
                 TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
                 "test_exon" )

meRIP_QC_report(se_M = se_fto[indx,],
                txdb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
                gtcoord = Gtcoord_mm10,
                min_num_mod = 1000,
                save_title = "FTO-3T3L1-QNB",
                DM_analysis = T,
                expected_change = "hyper",
                fdr_threshold = .05,
                GC_idx_feature = GC_cont[indx],
                DM_method = "QNB",
                mod_count_filter = 1
)

#correctly link methylation with cqn
se_fto <- se_mm10[, grepl("midbrain", colData( se_mm10 )$Publication ) & colData( se_mm10 )$Perturbation == "C"]

set.seed(1)

indx <- sample.int(nrow(se_fto),10000,replace = T)


meRIP_QC_report(
  se_M = se_fto,
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
  gtcoord = Gtcoord_mm10,
  min_num_mod = 1000,
  save_title = "NO_CQN",
  DM_analysis = F,
  GC_idx_feature = GC_cont,
  DM_method = "DESeq2",
  row_minimal_counts = 10,
  cqn = F
)


meRIP_QC_report(
  se_M = se_fto,
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
  gtcoord = Gtcoord_mm10,
  min_num_mod = 1000,
  save_title = "WITH_CQN",
  DM_analysis = F,
  GC_idx_feature = GC_cont,
  DM_method = "DESeq2",
  row_minimal_counts = 10,
  cqn = T
)


se_M = se_fto[indx]
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
gtcoord = Gtcoord_mm10
min_num_mod = 1000
save_title = "CQNTEST"
fdr_threshold = .05
GC_idx_feature = GC_cont[indx]
DM_method = "DESeq2"
row_minimal_counts = 10
cqn = T

SE_M = se_M
MODE = ifelse(DM_analysis,"DM","Meth")
DM_METHOD = DM_method
PCA = F
HDER = save_title
ROW_FILTER = row_minimal_counts
CQN = cqn
GC_INDX = GC_idx_feature #Still ~ 10% NAs if I use gene level GC content, some sites failed to mapped to any genes.

