##instalations 
#BiocManager::install("vctrs", force=TRUE)
##parameters 
options(timeout = 800) 
register(SerialParam())
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gc()
## libraries 
library(BSgenome.Hsapiens.UCSC.hg38)
library(Rbowtie2)
library(ggplot2)
library(DiffBind)
library(GenomicAlignments)
library(dplyr)
library(BiocParallel)
library(DiffBind)
library(profileplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(profileplyr)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(org.Hs.eg.db)
library(csaw)
library(PCAtools)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GenomicRanges)
library(readxl)  
library(stringr) 
library(pheatmap)
####
####set dir
setwd("D:/MiguelW12/Documents/gsd1a_atac_analysis/peak_data")
output_dir <- "results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
###load samples to dba - xls peak files from MACS2 + bam files
##GSD1a
GSD1a_vs_HC<-dba.peakset(NULL,
                         peaks="762.sorted_no_mt_noDup/762.sorted_no_mt_noDup_peaks_peaks.xls",
                         peak.caller="macs", sampID="GSD1a1",condition = "GSD1a", replicate=1,
                         bamReads = "bam_ready_qc/762.sorted_no_mt_noDup.bam")
##HC                     
GSD1a_vs_HC<-dba.peakset(GSD1a_vs_HC,
                         peaks="17507.sorted_no_mt_noDup/17507.sorted_no_mt_noDup_peaks_peaks.xls",
                         peak.caller="macs", sampID="HC2",condition = "HC", replicate=2,
                         bamReads = "bam_ready_qc/17507.sorted_no_mt_noDup.bam")
### 
## verify sample load
plot(GSD1a_vs_HC)
##set thresh
GSD1a_vs_HC$config$th = 0.1
## counts 
GSD1a_vs_HC_counts <- dba.count(GSD1a_vs_HC, bParallel = TRUE, score=DBA_SCORE_NORMALIZED,
                                bUseSummarizeOverlaps=TRUE, summits=250)
GSD1a_vs_HC_counts$config$th = 0.1
######QC of peaks 
peakdata <- dba.show(GSD1a_vs_HC_counts)$Intervals
peakdata
info <- dba.show(GSD1a_vs_HC_counts)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
libsizes
##
dba.plotVenn(GSD1a_vs_HC_counts,GSD1a_vs_HC_counts$masks$GSD1a, main = "Open chromatic region overlaps in GSD1a")
dba.plotVenn(GSD1a_vs_HC_counts,GSD1a_vs_HC_counts$masks$HC, main = "Open chromatic region overlaps in HC")
##occupancy
HCvsGSD1a_occupancy<-dba.peakset(GSD1a_vs_HC, consensus = DBA_CONDITION,minOverlap = 0.5)
HCvsGSD1a_occupancy
##pca
sset <- dba(HCvsGSD1a_occupancy,bSummarizedExperiment=TRUE)
dataFrame <- as.data.frame(assay(sset))
rownames(dataFrame) <- rowRanges(sset)
p <- pca(dataFrame, metadata = colData(sset), removeVar = 0.5)
biplot(p, showLoadings = FALSE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)

## normalization
GSD1a_vs_HC_counts_n <- dba.normalize(GSD1a_vs_HC_counts, library = DBA_LIBSIZE_PEAKREADS,
                                      method = DBA_DESEQ2)

normlibs <- cbind(FullLibSize=GSD1a_vs_HC_counts_n$norm$DESeq2$lib.sizes, NormFacs=GSD1a_vs_HC_counts_n$norm$DESeq2$norm.facs,
                  NormLibSize=round(GSD1a_vs_HC_counts_n$norm$DESeq2$lib.sizes/GSD1a_vs_HC_counts_n$norm$DESeq2$norm.facs))
rownames(normlibs) <- info$ID
normlibs
#####
##differntial
GSD1a_vs_HC_contrast <- dba.contrast(GSD1a_vs_HC_counts_n,categories=DBA_CONDITION,minMembers = 2)
GSD1a_vs_HC_contrast$config$th = 0.1
dif_GSD1a_vs_HC_counts <- dba.analyze(GSD1a_vs_HC_contrast)
dif_GSD1a_vs_HC_counts
dba.show(dif_GSD1a_vs_HC_counts, bContrasts=TRUE)
res_deseq_hc_vs_gsd <- dba.report(dif_GSD1a_vs_HC_counts, contrast = 1, method=DBA_DESEQ2, th=1)
write.csv(as.data.frame(res_deseq_hc_vs_gsd), file.path(output_dir, "res.csv"), row.names = TRUE)
### annotation
GenomeInfoDb::seqlevels(res_deseq_hc_vs_gsd)
GenomeInfoDb::seqlevels(txdb)
##convert chr names if needed 
newnames <- paste0(c("chr1","chr2","chr3","chr4", "chr5","chr6","chr7", "chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrMT"))
names(newnames) <- paste0(c("NC_000001.11","NC_000002.12","NC_000003.12","NC_000004.12", "NC_000005.10","NC_000006.12","NC_000007.14", "NC_000008.11","NC_000009.12","NC_000010.11","NC_000011.10","NC_000012.12","NC_000013.11","NC_000014.9","NC_000015.10","NC_000016.10","NC_000017.11","NC_000018.10","NC_000019.10","NC_000020.11","NC_000021.9","NC_000022.11","NC_000023.11","NC_000024.10","NC_012920.1"))
##
res_deseq_hc_vs_gsd_level <- renameSeqlevels(res_deseq_hc_vs_gsd,newnames)
GenomeInfoDb::seqlevels(res_deseq_hc_vs_gsd_level)
GenomeInfoDb::seqlevels(txdb)
length(seqlevels(res_deseq_hc_vs_gsd))
length(seqlevels(txdb))
res_deseq_hc_vs_gsd_level2 <- keepSeqlevels(res_deseq_hc_vs_gsd_level, standardChromosomes(res_deseq_hc_vs_gsd_level)[1:22], pruning.mode = "coarse")
seqlevels(res_deseq_hc_vs_gsd_level2)
# Annotate peaks
anno_df <- annotatePeak(res_deseq_hc_vs_gsd_level2, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb='org.Hs.eg.db')
anno_df
p <- as.GRanges(anno_df)
p
write.csv(as.data.frame(anno_df), file.path(output_dir, "res_deseq_DF_withanno.csv"), row.names = TRUE)
write.csv(as.data.frame(p), file.path(output_dir, "res_deseq_peaks_DF_withanno-granges.csv"), row.names = TRUE)
#vis
plotAnnoPie(anno_df)
plotAnnoBar(anno_df)
upsetplot(anno_df)
## vis res
dba.plotMA(dif_GSD1a_vs_HC_counts, method = DBA_DESEQ2, bNormalized = TRUE, sub = "HCvsGSD1a")
dba.plotVolcano(dif_GSD1a_vs_HC_counts, method = DBA_DESEQ2)
##
sum(res_deseq_hc_vs_gsd$Fold<0)
sum(res_deseq_hc_vs_gsd$Fold>0)
pvals <- dba.plotBox(dif_GSD1a_vs_HC_counts, contrast = 1)
pvals
### GSEA
gsea_df <- as.data.frame(anno_df)
original_gene_list <- gsea_df$Fold
names(original_gene_list) <- gsea_df$ENSEMBL
gene_list <- na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gsea_res <- gseGO(geneList=gene_list, 
                  ont ="ALL", 
                  keyType = "ENSEMBL", 
                  minGSSize = 3, 
                  maxGSSize = 1000, 
                  pvalueCutoff = 0.1, 
                  verbose = TRUE, 
                  OrgDb = 'org.Hs.eg.db',
                  pAdjustMethod = "BH"
)

#save file
write.csv(as.data.frame(gsea_res), file.path(output_dir, "gseago-atac.csv"), row.names = TRUE)
## subtract res and vis
id_list <- c("GO:0099503",
             "GO:0005768",
             "GO:1901135",
             "GO:0005764",
             "GO:0030163",
             "GO:0031667",
             "GO:0006914",
             "GO:0003682",
             "GO:0000302",
             "GO:0007005", 
             "GO:0005740", 
             "GO:0005759", 
             "GO:0006325") 

gsea_res@result = gsea_res@result[gsea_res@result$ID %in% id_list, ]
## vis
edo <- setReadable(gsea_res, 'org.Hs.eg.db', 'ENSEMBL')
cnetplot(edo, node_label="all", categorySize="pvalue", cex_label_gene = 0.3, cex_label_category = 1,
         cex_category = 1, cex_gene = 0.3, colorEdge = TRUE, showCategory=10) 

##gviz ( for sample plotting generate granges per grouping+ chr plotting available)

dTrack <- DataTrack(p, name = "uniform")
colnames(mcols(p))
chromosome(dTrack) <- "chr9"
dTrack <- dTrack[1:4, ]
plotTracks(dTrack)
atrack <- AnnotationTrack(anno_range, name = "CpG")
##
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = "chr9")
plotTracks(ideoTrack, from = 127397138, to = 127407855)
##
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = "chr2")
plotTracks(ideoTrack, from = 127397138, to = 127407855)
##HAS1
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = "chr19")
plotTracks(ideoTrack, from = 51711112, to = 51725991)
##arg1
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = "chr6")
plotTracks(ideoTrack, from = 131571226, to = 131586329)
##additional vis options
plotTracks(dtrack, type = "boxplot", showSampleNames = TRUE, 
           cex.sampleNames = 0.6)

plotTracks(dTrack, type = "histogram", groups = rownames(values(dTrack)),
           showSampleNames = TRUE, cex.sampleNames = 0.6, separator = 1)



### motif analysis - performed on linux with homer
## combine with rna sig genes with significant motifs if available 

motif_data_up_sig <- read_excel("ATAC_MOTIFS_UP_GSD1A_sig.xlsx")
motif_data_down_sig <- read_excel("ATAC_MOTIFS_DOWN_GSD1A_sig.xlsx")
rna_sig_genes <- read_excel("gsd1a_rna_2024/rna_sig_0.1.xlsx")  
#Prep data
#columns to uppercase for case-insensitive matching
rna_sig_genes <- rna_sig_genes %>% mutate(SYMBOL = toupper(SYMBOL))
motif_data_up_sig <- motif_data_up_sig %>% mutate(`Motif Name` = toupper(`Motif Name`))
motif_data_down_sig <- motif_data_down_sig %>% mutate(`Motif Name` = toupper(`Motif Name`))
motif_data_up_sig$TF_in_RNA_sig <- NA  
motif_data_down_sig$TF_in_RNA_sig <- NA  

#run through each gene in the rna data and check if it appears in the TF column
for (gene in rna_sig_genes$SYMBOL) {
  
  motif_data_up_sig$TF_in_RNA_sig <- ifelse(str_detect(motif_data_up_sig$`Motif Name`, gene), gene, motif_data_up_sig$TF_in_RNA_sig)
}

for (gene in rna_sig_genes$SYMBOL) {
  
  motif_data_down_sig$TF_in_RNA_sig <- ifelse(str_detect(motif_data_down_sig$`Motif Name`, gene), gene, motif_data_down_sig$TF_in_RNA_sig)
}
#save
write.csv(motif_data_up_sig, "up-motif_with_overlapping_TFs.csv", row.names = FALSE)
write.csv(motif_data_down_sig, "down-motif_with_overlapping_TFs.csv", row.names = FALSE)

## gene expression vis for TFs from motif analysis
## gene table 
gene_df = read.csv("data_integration/star_ut-normatotaldata_40filt.csv")
## choose genes 
gene_df2 <- gene_df[gene_df$X %in% c("ENSG00000075426", "ENSG00000116044", "ENSG00000134954", "ENSG00000105997", "ENSG00000100644",
                                     "ENSG00000177606", "ENSG00000175832", "ENSG00000171872", "ENSG00000111704", "ENSG00000197757"),]
#normalize data for heatmap
normalize_data <- function(df) {
  df[, 2:ncol(df)] <- scale(df[, 2:ncol(df)])
  return(df)
}

#normalize dataframes
n_gene_df2 <- normalize_data(gene_df2)
## choose columns based on groups 
gene_df2$HC_mean <- rowMeans(gene_df2[ , c(2,3,4)], na.rm = TRUE)
gene_df2$GSD_mean <- rowMeans(gene_df2[ , c(5,6,7)], na.rm = TRUE)
rownames(gene_df2) <- gene_df2$X
df_HEATMAP = subset(gene_df2, select = c(8,9))
rownames(df_HEATMAP) <- df_HEATMAP$
  heatmap_motifs <- pheatmap(
    df_HEATMAP,
    cluster_rows = TRUE, 
    cluster_cols = FALSE, 
    show_rownames = TRUE, 
    main = "Non-Annotated Heatmap", color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    scale = ("row")) 



gc()

#################################################################################

