## RNA-seq General data analysis pipeline 
## for install - 
BiocManager::install('')
##options 
options(ggrepel.max.overlaps = Inf)
colramp = colorRampPalette(c(3,"white",2))(20)
orgainsm = 'org.Hs.eg.db'
gc()
## libraries to load 
#libs
library(data.table)  
library(ggplot2)   
library(readr)
library(airway)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(genefilter)
library(PoiClaClu)
library(AnnotationDbi)
library(dplyr)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(dplyr)
library(ggfortify)
library(pcaExplorer)
library(cli)
library(topGO)
library(markdown)
library(PCAtools)
library(DEGreport)
library(radiant.data)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(pathview)
library(cowplot)
library(ReactomePA)
library(DOSE)
library(AnnotationHub)
library(dplyr)
library(tibble)

### 
### full data - QC check
count_dir <- "Z:/Lab stuff/MW_231_assay_protocols_new/omics_analysis/RNA_seq/downstream_counts_proc/example_files/counts_txt"
##
#read and merge count files frOm different samples
merge_count_files <- function(count_dir) {
  #list files in dir
  count_files <- list.files(count_dir, pattern = "\\.txt$", full.names = TRUE)
  print(count_files)
  count_data <- lapply(count_files, function(file) {
    sample_name <- gsub("\\.counts\\.txt", "", basename(file))
    counts <- fread(file, header = FALSE)
    colnames(counts) <- c("GeneID", sample_name)
    return(counts)
    
  })
  
  #merge by GeneID
  print(count_data)
  merged_count_data <- Reduce(function(x, y) merge(x, y, by = "GeneID", all = TRUE), count_data)
  return(merged_count_data)
}

#run
count_data <- merge_count_files(count_dir)
write.csv(as.data.frame(count_data), file="MERGED_count.csv")

##
###basic statistics
summary_stats <- apply(count_data[, -1], 2, summary)
write.csv(as.data.frame(summary_stats), file="count_data_summary_stats.csv")
# count zero counts for each sample
zero_counts <- colMeans(count_data[, -1] == 0, na.rm = TRUE) * 100
print(zero_counts)
### qc met
create_qc_metrics <- function(count_data, output_dir) {
  
  dir.create(output_dir, showWarnings = FALSE)
  
  #data frame for qc metrics
  qc_df <- data.frame(Sample = character(),
                      Total_Genes = numeric(),
                      Num_NAs = numeric(),
                      Num_Zeros = numeric(),
                      stringsAsFactors = FALSE)
  
  for (col_name in colnames(count_data)) {
    counts <- count_data[[col_name]]
    
    total_genes <- length(counts)
    num_nas <- sum(is.na(counts))
    num_zeros <- sum(counts == 0)
    
    qc_df <- rbind(qc_df, data.frame(Sample = col_name,
                                     Total_Genes = total_genes,
                                     Num_NAs = num_nas,
                                     Num_Zeros = num_zeros))
  }
  
  write.csv(qc_df, file.path(output_dir, "qc_metrics.csv"), row.names = FALSE)
}
#run
create_qc_metrics(count_data, "count_stats")
##
## vis qc
# hists
create_count_histograms <- function(count_data, output_dir = "path/to/histograms") {
  dir.create(output_dir, showWarnings = FALSE)
  
  for (col_name in colnames(count_data)) {
    if (col_name == "GeneID") next
    
    counts <- as.numeric(count_data[[col_name]])
    
    if (all(!is.na(counts)) && all(counts >= 0)) {
      counts <- counts[counts > 0]
      
      log_counts <- log10(counts)
      
      png(file.path(output_dir, paste0(col_name, "_histogram.png")))
      hist(log_counts, main = paste("Histogram of Counts -", col_name),
           xlab = "Log10(Counts)", ylab = "Frequency", xlim = c(min(log_counts), max(log_counts)),
           freq = TRUE, breaks = "Sturges")
      dev.off()
      
      missing_values <- sum(is.na(counts))
      if (missing_values > 0) {
        cat(paste("Number of missing values for", col_name, ":", missing_values, "\n"))
      }
      
      cat("histogram saved-", file.path(output_dir, paste0(col_name, "_histogram.png")), "\n")
    } else {
      warning(paste("invalid count data in sample", col_name, "- skipping histogram generation."))
    }
  }
}
##run
create_count_histograms(count_data, output_dir = "count_histograms")
#boxes per sample
create_count_boxplots <- function(count_data, output_dir) {
  dir.create(output_dir, showWarnings = FALSE)
  
  
  for (col_name in colnames(count_data)) {
    if (col_name == "GeneID") next
    
    counts <- na.omit(count_data[[col_name]])
    counts <- counts[counts > 0]  # Filter out zero and negative counts
    
    df <- data.frame(Sample = rep(col_name, length(counts)),
                     Counts = counts)
    
    p <- ggplot(df, aes(x = Sample, y = Counts)) +
      geom_boxplot(fill = "skyblue", color = "black") +
      labs(title = paste("Boxplot of Counts -", col_name),
           x = "Sample", y = "Counts") +
      scale_y_continuous(trans = "log10")
    
    png(paste(output_dir, "/", col_name, "_boxplot.png", sep = ""), width = 800, height = 600)
    print(p)
    dev.off()
  }
}
#run 
create_count_boxplots(count_data, "count_BOXPLOTS")
## boxplot all 
create_count_boxplots2 <- function(count_data, output_dir) {
  dir.create(output_dir, showWarnings = FALSE)
  
  
  
  all_counts <- data.frame(Sample = character(),
                           Counts = numeric(),
                           stringsAsFactors = FALSE)
  
  for (col_name in colnames(count_data)) {
    if (col_name == "GeneID") next
    
    counts <- na.omit(count_data[[col_name]])
    counts <- counts[counts > 0]  
    
    all_counts <- rbind(all_counts, data.frame(Sample = rep(col_name, length(counts)),
                                               Counts = counts,
                                               stringsAsFactors = FALSE))
  }
  
  p <- ggplot(all_counts, aes(x = Sample, y = Counts, fill = Sample)) +
    geom_boxplot(color = "black") +
    labs(title = "Boxplot of Counts",
         x = "Sample", y = "Counts") +
    scale_y_continuous(trans = "log10") +
    theme(legend.position = "top")
  
  output_path <- file.path(output_dir, "all_samples_boxplot.png")
  png(output_path, width = 800, height = 600)
  print(p)
  dev.off()
  
  cat("boxplot saved as", output_path, "\n")
}
## run
create_count_boxplots2(count_data, "count_BOXPLOTS")

### differential analysis prep
## arrange your meta data as shown in the csv example
Count_Data = read.csv('path/to/counts', header=TRUE)
Col_Data = read.csv('path/to/meta', header=TRUE)
rownames(Count_Data) <- Count_Data$GeneID
Count_Data = subset(Count_Data, select = -c(GeneID))
rownames(Col_Data) <- Col_Data$X
Col_Data  = subset(Col_Data, select = -c(X))
# make sure colnames of counts equal rownames of meta
all(colnames(Count_Data) %in% rownames(Col_Data))
######
## desqe build
## choose your design, in this case group 
dds_all <- DESeq2::DESeqDataSetFromMatrix(countData = Count_Data, colData = Col_Data, design = ~ group)
## check model 
dds_all
#add filter for counts - use initial QC to decide
# run deseq
keep <- rowMeans(counts(dds_all)) >= 15
dds_all <- dds_all[keep,]
## make sure filter did not remove to much data
dds_all
## run deseq
dds_all_run <- DESeq2::DESeq(dds_all)
resultsNames(dds_all_run)
## optional to use LRT
dds_lrt <- DESeq(dds_all_2vs2, test="LRT", reduced = ~ 1)
res_LRT <- results(dds_lrt, contrast = c("divided","72_GS_UT", "72_HC_UT"))
### Any comparison is optiona, in this case - hcvsgsd
res_all <- DESeq2::results(dds_all_run, contrast = c("group","HC_UT", "GSD_UT"))
#### visualise results
hist(res_all$pvalue)
hist(res_all$padj)
hist(res_all$log2FoldChange)
summary(res_all)
# save results
write.csv(as.data.frame(res_all), file="diff_analysis_res_hcvsgsd_ut.csv")
###add annotation 
## choose organism- human in this case
orgainsm = 'org.Hs.eg.db'
head(res_all)
rownames(res_all)
## choose relevant annotations... keytype should be your current annotation of genes. 
anno <- AnnotationDbi::select(org.Hs.eg.db, rownames(res_all), columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME", "CHR"),keytype = "ENSEMBL")
keytypes(org.Hs.eg.db)
keys(org.Hs.eg.db)
head(anno_2)
res_anno = cbind( ENSEMBL = rownames( res_all), res_all)
outtable <- left_join(as.data.frame(res_anno), anno )
outtable <- outtable[order(outtable$padj),]
#save
write.csv(as.data.frame(outtable), file="diff_analysis_annotated_res.csv")
#Enhanced Volcano
## change title etc according to your data
EnhancedVolcano::EnhancedVolcano(as.data.frame(outtable), lab = outtable$SYMBOL, x = 'log2FoldChange', y = 'padj', xlim = c(-8,8), title = 'HC vs Gsd', pCutoff = 0.1, FCcutoff = 1, pointSize = 2.0, labSize = 4.0)
##export normlized data 
norm_data_star <- counts(dds_all_run, normalized=TRUE)
write.csv(as.data.frame(norm_data_star ), file="normal_data_deseq.csv")
### sig genes 
resSig = res_all[ which(res_all$padj < 0.1 & res_all$log2FoldChange > 1), ]
resSig2 = res_all[ which(res_all$padj < 0.1 & res_all$log2FoldChange < -1), ]
resSig_all = rbind(resSig, resSig2)
write.csv(as.data.frame(resSig_all), file="resSig_all.csv")

## plot chosen main genes following examining lists
## load dfs 
Imp_genes_fold = read.csv('path/to/your/fold/genes', header=TRUE)
Imp_genes_class = read.csv('path/to/genes/based/on/class', header=TRUE)
sorted_imp_genes_df <- merge(Imp_genes_fold, Imp_genes_class, by = "Gene", all.x = TRUE)
sorted_imp_genes_df$Gene <- factor(sorted_imp_genes_df$Gene, levels = sorted_imp_genes_df$Gene[order(sorted_imp_genes_df$Class)])
# colors for each class
class_colors <- c("Hox_genes" = "skyblue", "Metabolic_genes" = "orange", "Fibroblast_genes" = "yellow", "Wnt_pathway_genes" = "purple", "Lysosomal/Mitochondrial_genes" = "red")  
#plot
level_order <- c("Hox_genes", "Metabolic_genes", "Fibroblast_genes", "Wnt_pathway_genes", "Lysosomal/Mitochondrial_genes")
ggplot(sorted_imp_genes_df, aes(x = Gene, y = L2FC, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = class_colors) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5)) +
  labs(title = "Genes", x = "Gene", y = "L2FC(HC/GSD1a")

## other simple ploting options 
#PLOT SPECFIC GENES 
plotCounts(dds_all, gene='ENSG00000181019', intgroup="divided", col())
plotCounts(dds_all, gene=which.min(res$padj), intgroup="divided")
d <- plotCounts(dds_all, gene='ENSG00000132170', intgroup="divided", 
                returnData=TRUE)

e <- ggplot(d, aes(x = divided, y = count, color = divided)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("peroxisome proliferator activated receptor gamma") +
  theme(plot.title = element_text(hjust = 0.5))
e + theme(text = element_text(size = 20))   


##heatmap 
## simple vs
vst_mean <- vst(dds_all)
mat <- assay(vst_mean)
ids <- resSig_all$ENSEMBL
DEGENES <- mat[ids, ]
annotation <- as.data.frame(colData(vst_mean)[, c("divided")])
rownames(annotation) <- colnames(vst_mean)
pheatmap(DEGENES, scale = "row", show_rownames = FALSE, clustering_distance_rows = "correlation", annotation_col = annotation, main ="Differentially Expressed Genes")
# heatmap( by variance and mean of groups)
dds_all_rlog  <- rlog(dds_all)
variances <- apply(assay(dds_all_rlog), 1, var)
upper_var <- quantile(variances, 0.75)
df_by_var <- data.frame(assay(dds_all_rlog )) %>%  dplyr::filter(variances > upper_var)
## choose samples based on your data...
df_by_var$HC_48_MEAN <- rowMeans(df_by_var[ , c(5,6)], na.rm = TRUE)
df_by_var$HC_72_MEAN <- rowMeans(df_by_var[ , c(3,4)], na.rm = TRUE)
df_by_var$GSD_72_MEAN <- rowMeans(df_by_var[ , c(1,2)], na.rm = TRUE)
df_by_var$GSD_48_MEAN <- rowMeans(df_by_var[ , c(7,8)], na.rm = TRUE)
df_by_var_HEATMAP = subset(df_by_var, select = c(9,10,11,12))
heatmap_all_t <- pheatmap(
  df_by_var_HEATMAP,
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  show_rownames = FALSE, 
  main = "Non-Annotated Heatmap", color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  scale = ("row")) 
### PCA 
## SIMPLE VS
#pca's
rld_mean <- rlog(dds_all)
plotPCA(rld_mean, intgroup = c("divided"))
vst_mean <- vst(dds_all)
plotPCA(vst_mean, intgroup = c("divided"))
## better vs 
## you can use vst to scale as well. 
rdat <- assay(rlog(dds_all))
p <- pca(rdat, metadata = colData(dds_all), removeVar = 0.5)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
elbow <- findElbowPoint(p$variance)
elbow
horn <- parallelPCA(rdat)
horn$n
biplot(p,
       colby = 'group',
       colLegendTitle = 'pca_rnaseq',
       # encircle config
       encircle = TRUE,
       encircleFill = TRUE,
       ellipseLevel = 0.95,
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)

##deg reports cluster analysis 
res_clust <- results(dds_all)
## choose sig genes based on thresh
sig_res <- res_clust %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.1)

# Get sig gene lists
sig_genes <- sig_res %>% 
  pull(gene)
length(sig_genes)
## choose amount of sig genes
clustering_sig_genes <- sig_res %>%
  arrange(padj) %>%
  head(n=400)
## norm
rld <- rlog(dds_all)
rld_mat <- assay(rld)

cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]

clusters_final <- degPatterns(cluster_rlog, metadata = Col_Data,
                           time = "CELL.LINE", col="divided")
## vis
color = colorRampPalette(c("red","orange","blue"))
ggplot(clusters_final[["normalized"]],
       mapping = aes(divided, value, fill = time)) 
+
  geom_boxplot(colour = colorRampPalette(c("red","orange","blue" ))) 
##cluster 1 extraction
cluster_groups <- clusters_final$df
group1 <- clusters_final$df %>%
  filter(cluster == 1)
## annotate cluster 
anno_2 <- AnnotationDbi::select(org.Hs.eg.db, rownames(group1), columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"),keytype = "ENSEMBL")
head(anno_2)
group1 = cbind( ENSEMBL = rownames( group1), group1)
outtable_clusterg1 <- left_join(as.data.frame(group1), anno_2 )

write.csv(as.data.frame(outtable_clusterg1), file="group1-clusterano.csv")
####



#gene_set_enrichment_analysis
##Gene ontology 
new_df = read.csv('path/to/full/data', header=TRUE)
original_gene_list <- new_df$log2FoldChange
names(original_gene_list) <- new_df$X
gene_list <- na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
## choose parameters as needed
gse_ut <- gseGO(geneList=gene_list, 
                   ont ="ALL", 
                   keyType = "ENSEMBL", 
                   minGSSize = 20, 
                   maxGSSize = 10000, 
                   pvalueCutoff = 0.1,
                   nPerm = 10000,
                   verbose = TRUE, 
                   OrgDb = 'org.Hs.eg.db',
                   pAdjustMethod = "fdr")

#save results 
write.csv(as.data.frame(gse_ut), file="gseago.csv")
## select pathways for visualzation 

id_list_gsea_plot <- c("GO:1902493",
                       "GO:0004402",
                       "GO:0000123",
                       "GO:0061733",
                       "GO:0010628",
                       "GO:0031248",
                       "GO:0008276",
                       "GO:0008135",
                       "GO:0090079",
                       "GO:0006418") 

## vis 
gse_ut@result = gse_ut@result[gse_ut@result$ID %in% id_list_gsea_plot, ]
dotplot(gse_ut, showCategory=30, font.size = 8,label_format = 60, color = "p.adjust", title = "dotplot_GSEA", split=".sign") + facet_grid(.~.sign)
ema <- pairwise_termsim(gse_ut, method="JC", semData = NULL, showCategory = 30)
emapplot(ema, cex_category=1, cex_label_category=0.8, cex_line=0.8, color="pvalue", layout="kk" )
treeplot(ema, showCategory = 30, nWords = 4, fontsize = 2, nCluster = 5)
edo <- setReadable(gse_ut, 'org.Hs.eg.db', 'ENSEMBL')
heatplot(edo,  foldChange=gene_list, label_format = 8, showCategory=5)
p1 <- cnetplot(edo, foldChange=gene_list)
p2 <- cnetplot(edo, categorySize="pvalue", foldChange=gene_list)
p3 <- cnetplot(edo, foldChange=gene_list, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
cnetplot(edo, node_label="all", categorySize="pvalue", cex_label_gene = 0.6, cex_label_category = 1,
         cex_category = 1, cex_gene = 0.5, colorEdge = TRUE, showCategory=10) 
ridgeplot(gse_ut, showCategory = 37, label_format = 60) + labs(x = "enrichment distribution")
upsetplot(gse_ut) 
gseaplot2(edo, geneSetID = 20:23, pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")

#kegg 
View(original_gene_list)
ids <- bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
df_kegg = df_kegg[df_kegg$X %in% df_kegg$ENSEMBL,]
df_kegg$Y = dedup_ids$ENTREZID
kegg_gene_list <- df_kegg$log2FoldChange
names(kegg_gene_list) <- df_kegg$Y
kegg_gene_list <- na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "hsa"
kegg_gse <- gseKEGG(geneList     = kegg_gene_list,
                             organism     = kegg_organism,
                             minGSSize    = 3,
                             maxGSSize    = 1000,
                             pvalueCutoff = 0.25,
                             pAdjustMethod = "fdr",
                             keyType       = "kegg")
#save file
write.csv(as.data.frame(kegg_gse), file = "kegg_gse.csv")
#vis
## filter pathways
kegg_id_list <- c("hsa04630", "hsa04110", "hsa04621", "hsa00350", "hsa04020", "hsa04060", "hsa04668", "hsa00010")

kegg_gse@result = kegg_gse@result[kegg_gse@result$ID %in% kegg_id_list, ]

ridgeplot(kegg_gse, showCategory = 37, label_format = 60) + labs(x = "kegg gene set enrichment distribution")
## choose pathway for kegg vis 
dme_try_meanf10_ut_72 <- pathview(gene.data=kegg_gene_list, pathway.id="hsa00520", species = kegg_organism)
knitr::include_graphics("hsa00520.pathview.png")

#Wikipathway
GSWP <- gseWP(kegg_gene_list, nPerm = 2000,
              organism = "Homo sapiens", pvalueCutoff = 1, pAdjustMethod = 'BH')

write.csv(as.data.frame(GSWP), file = "gse_wiki.csv")

# REACTOME + DO - CREATE DE' LIST WITH ENTREZ ID INSTEAD OF ENSEMBL 
gene_list_entrez = new_df$X
eg = bitr(gene_list_entrez, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
names(original_gene_list) <- eg$ENTREZID
original_gene_list = sort(original_gene_list, decreasing = TRUE)
y <- gsePathway(original_gene_list, 
                pvalueCutoff = 0.1,
                pAdjustMethod = "fdr",
                verbose = FALSE)
y <- setReadable(y, 'org.Hs.eg.db')
write.csv(as.data.frame(y), file = "REACTOME_gse.csv")
# vis
viewPathway("Activation of AMPK downstream of NMDARs",
            readable = FALSE,
            layout = "nicely",
            foldChange = original_gene_list)
cnetplot(y, node_label="all", categorySize="p.adjust", cex_label_gene = 0.6, cex_label_category = 0.8,
         cex_category = 1, cex_gene = 0.5, colorEdge = TRUE, showCategory=20) 
# DO 
y2 <- gseDO(original_gene_list,
           nPerm = 2000, 
           minGSSize     = 1,
           pvalueCutoff  = 0.25,
           pAdjustMethod = "fdr",
           verbose       = FALSE)
write.csv(as.data.frame(y2), file = "gse_DO.csv")
# select pathways
do_id_list <- c("DOID:4195", "DOID:10603", "DOID:0080014")
y2@result = y2@result[y2@result$ID %in% do_id_list, ]

do_vis <- setReadable(y2, 'org.Hs.eg.db', 'ENTREZID')
heatplot(do_vis,  foldChange=original_gene_list, label_format = 8, showCategory=5)
## ncg
ncg <- gseNCG(original_gene_list,
              pvalueCutoff  = 0.5,
              pAdjustMethod = "none",
              verbose       = FALSE)
ncg <- setReadable(ncg, 'org.Hs.eg.db')
write.csv(as.data.frame(ncg), file = "gse_DO_ncg.csv")

## dgn 
dgn <- gseDGN(original_gene_list,
              pvalueCutoff  = 0.01,
              pAdjustMethod = "none",
              verbose       = FALSE)
write.csv(as.data.frame(dgn), file = "gse_DO_dgn.csv")
## enrichemnt analysis ORA 
#GO
original_gene_list <- new_df$log2FoldChange
names(original_gene_list) <- new_df$X
gene_list <-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
## choose thresh
sig_genes = subset(new_df, padj < 0.1)
genes <- sig_genes$log2FoldChange
names(genes) <- sig_genes$X
genes <- na.omit(genes)
genes_up <- names(genes)[(genes) > 1]
genes_down <- names(genes)[(genes) < -1]

## try spliting to down reg and up reg 
go_enrich <- enrichGO(gene = genes_up,
                               universe = names(gene_list),
                               OrgDb = org.Hs.eg.db, 
                               keyType = 'ENSEMBL',
                               readable = T,
                               ont = "ALL",
                               pvalueCutoff = 0.1, 
                               pAdjustMethod = "BH")
#
write.csv(as.data.frame(go_enrich), file="go_enrich.csv")
## filter pathways
id_up <- c("GO:0030178", "GO:0045725",
                                           "GO:0010906",
                                           "GO:0042800",
                                           "GO:0008276",
                                           "GO:0098531",
                                           "GO:0018024",
                                           "GO:0051568")

go_enrich@result = go_enrich@result[go_enrich@result$ID %in% id_up, ]
## vis
goplot(go_enrich,
       showCategory = 170,
       color = "p.adjust",
       layout = "sugiyama",
       geom = "label")


#########kegg
ids <- bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
df = new_df[new_df$X %in% dedup_ids$ENSEMBL,]
df$Y = dedup_ids$ENTREZID
kegg_gene_list <- df$log2FoldChange
names(kegg_gene_list) <- df$Y
kegg_gene_list <- na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_sig_genes = subset(df, padj < 0.1)
kegg_genes <- kegg_sig_genes$log2FoldChange
names(kegg_genes) <- kegg_sig_genes$Y
kegg_genes <- na.omit(kegg_genes)
kegg_genes_up <- names(kegg_genes)[abs(kegg_genes) > 1]
kegg_genes_down <- names(kegg_genes)[(kegg_genes) < -1]

kegg_organism = "hsa"
kegg_enrich <- enrichKEGG(gene=kegg_genes_up, universe=names(kegg_gene_list),organism=kegg_organism,
                        pvalueCutoff = 0.1, keyType = "ncbi-geneid", pAdjustMethod = "BH")
write.csv(as.data.frame(kegg_enrich), file="kegg_enrich.csv")

# filter pathways 
id_list <- c("hsa05225",
                   "hsa04932",
                   "hsa03050",
                   "hsa05415",
                   "hsa00190",
                   "hsa05022",
                   "hsa05010"
)

kegg_enrich@result = kegg_enrich@result[kegg_enrich@result$ID %in% id_list, ]
## vis
barplot(kegg_enrich, 
        showCategory = 10, 
        title = "kegg_enriched_pathways-",
        font.size = 10)

#Wiki 
wp <- enrichWP(kegg_genes_up, pvalueCutoff = 0.1, organism = "Homo sapiens", pAdjustMethod = "BH")
write.csv(as.data.frame(wp), file="wp_up.csv")

## Reactome + DO 
sig_RA <- subset(new_df, padj < 0.01)
gene_list_entrez = sig_RA$X
original_gene_list_ENTREZ <- sig_RA$log2FoldChange
eg = bitr(gene_list_entrez, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dedup_list = eg[!duplicated(eg[c("ENTREZID")]),]
names(original_gene_list_ENTREZ) <- eg$ENTREZID
de_up <- names(original_gene_list_ENTREZ)[(original_gene_list_ENTREZ) > 1]
de_down <- names(original_gene_list_ENTREZ)[(original_gene_list_ENTREZ) < -1]

x <- enrichPathway(gene=de_d, pvalueCutoff = 0.1, readable=TRUE, pAdjustMethod = "BH")
x <- setReadable(x, 'org.Hs.eg.db')
write.csv(as.data.frame(x), file = "down_RA_enrich.csv")

##DO
d <- enrichDO(gene          = de_down,
              ont           = "DO",
              pvalueCutoff  = 0.1,
              pAdjustMethod = "BH",
              universe      = names(original_gene_list_72_ut_mean_f10_ENTREZ),
              minGSSize     = 3,
              maxGSSize     = 1000,
              readable      = TRUE)
write.csv(as.data.frame(d), file = "down_DO_enrich.csv")







##############################################################
