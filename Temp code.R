# Generating GC content idx for each transcript.

# mm10 and hg19.

# test for dataset /Users/zhenwei/Documents/GitHub/TREW-cons/B_COUNT_2017_12_5/se_mm10.rds.

# Try to inherit the arguments in meRIP_mod_QC_report into the MeRIP_QC arguments.

se_mm10 <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/B_COUNT_2017_12_5/se_mm10.rds")

#For developing, testing, and debugging
library(DESeq2)
library(ggplot2)
library(Guitar)
library(GenomicFeatures)
library(SummarizedExperiment)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
se_M <- se_mm10[,grepl("midbrain", colData( se_mm10 )$Publication )]
gtcoord <- readRDS("/Users/zhenwei/Datasets/Gtcoords/Gtcoord_mm10.rds")
p_threshold = NULL
fdr_threshold = NULL
log2FC_cutoff = 0
min_num_Mod = 10000
Save_inference_result = TRUE
GC_idx_feature = NULL

Gene_GC_mm10 <- readRDS("/Users/zhenwei/Datasets/GC_content_Genes/Gene_GC_mm10.rds")
fol <- findOverlaps( rowRanges( se_mm10 ), Gene_GC_mm10 )
GC_cont = rep(NA,nrow(se_mm10))
GC_cont[queryHits(fol)] = mcols(Gene_GC_mm10)[subjectHits(fol),]
mean(is.na(GC_cont)) #Notice that there are ~ 11% sites do not belong to exons.
GC_idx_feature = GC_cont
save_title = "modX"
save_dir = save_title
DM_analysis = T
PCA_PLOT = F
Expected_change = "hyper"

#Extracting exons level GC content (by genes).
Retriev_gene_GC_content <- function(txdb,bsgnm){
exbg <- exonsBy(txdb,by = "gene")
gene_ex_seq <- DNAStringSet( Views(bsgnm,unlist(exbg)) )
GC_cont <- letterFrequency(gene_ex_seq, letters="CG", as.prob = F)
Total_cont <- width(gene_ex_seq)
GC_content_pergene <- tapply(GC_cont,names(gene_ex_seq),sum)/tapply(Total_cont,names(gene_ex_seq),sum)
mcols(exbg) = GC_content_pergene
return(exbg)
}

Gene_GC_mm10 <- Retriev_gene_GC_content(TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
                                        BSgenome.Mmusculus.UCSC.mm10::Mmusculus)

Gene_GC_hg19 <- Retriev_gene_GC_content(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                                        BSgenome.Hsapiens.UCSC.hg19::Hsapiens)


saveRDS(Gene_GC_hg19,"/Users/zhenwei/Datasets/GC_content_Genes/Gene_GC_hg19.rds")
saveRDS(Gene_GC_mm10,"/Users/zhenwei/Datasets/GC_content_Genes/Gene_GC_mm10.rds")

#GC content for each feature (Granges)
Gene_GC_mm10 <- readRDS("/Users/zhenwei/Datasets/GC_content_Genes/Gene_GC_mm10.rds")

fol <- findOverlaps( rowRanges( se_mm10 ), Gene_GC_mm10 )

GC_cont = rep(NA,nrow(se_mm10))

GC_cont[queryHits(fol)] = mcols(Gene_GC_mm10)[subjectHits(fol),]


mean(is.na(GC_cont)) #Notice that there are ~ 11% sites do not belong to exons.

GC_idx_feature = GC_cont


ds_result <- meripQC::DESeq2_merip(se_M,MODE = "DM")
saveRDS(ds_result,"example_dds_DM.rds")
ds_result <- meripQC::DESeq2_merip(se_M,MODE = "Meth")
saveRDS(ds_result,"example_dds_Meth.rds")

drs_DM <- readRDS("example_dds_DM.rds")
se_M <- readRDS("example_dds_Meth.rds")
Decision_dsresult(drs_DM,1,0,0.05,5000,T,"hyper","Fto midbr")

meripQC::meRIP_mod_QC_report(se_M = se_mm10[,1:8],
                    txdb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
                    gtcoord = Gtcoord_mm10,
                    min_num_Mod = 5000,
                    save_title = "FTO-3T3L1",
                    DM_analysis = T,
                    Expected_change = "hyper",
                    fdr_threshold = .05,
                    PCA_PLOT = T,
                    GC_idx_feature = GC_cont)

meripQC::meRIP_mod_QC_report(se_M = se_mm10[,1:8],
                             txdb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
                             gtcoord = Gtcoord_mm10,
                             min_num_Mod = 5000,
                             save_title = "FTO-3T3L1-QNB",
                             DM_analysis = T,
                             DM_method = "QNB",
                             Expected_change = "hyper",
                             fdr_threshold = .05,
                             PCA_PLOT = T,
                             GC_idx_feature = GC_cont)


meripQC::meRIP_mod_QC_report(se_M = se_mm10[,48:59],
                             txdb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
                             gtcoord = Gtcoord_mm10,
                             min_num_Mod = 5000,
                             save_title = "FTO-midbr",
                             DM_analysis = T,
                             Expected_change = "hyper",
                             DeSeq2_fdr_threshold = .05,
                             PCA_PLOT = T,
                             GC_idx_feature = GC_cont)

library(SummarizedExperiment)
se_mm10 <- readRDS("/Users/zhenwei/Documents/GitHub/TREW-cons/B_COUNT_2017_12_5/se_mm10.rds")
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
Gtcoord_mm10 <- readRDS("/Users/zhenwei/Datasets/Gtcoords/Gtcoord_mm10.rds")

library(GenomicRanges)
Gene_GC_mm10 <- readRDS("/Users/zhenwei/Datasets/GC_content_Genes/Gene_GC_mm10.rds")
fol <- findOverlaps( rowRanges( se_mm10 ), Gene_GC_mm10 )
GC_cont = rep(NA,nrow(se_mm10))
GC_cont[queryHits(fol)] = mcols(Gene_GC_mm10)[subjectHits(fol),]

as.data.frame( SummarizedExperiment::colData(se_mm10) )[17:28,]

meripQC::meRIP_mod_QC_report(se_M = se_mm10[,17:28],
                             txdb = TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
                             gtcoord = Gtcoord_mm10,
                             min_num_Mod = 5000,
                             save_title = "Batch-HS",
                             DM_analysis = F,
                             Expected_change = NULL,
                             DeSeq2_fdr_threshold = .05,
                             PCA_PLOT = T,
                             GC_idx_feature = GC_cont)



