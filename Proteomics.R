###
#auto installation of all required packages
required_packages <- c("ggplot2", "dplyr", "readr", "plotly", "data.table", "magrittr",
                       "tidyr", "ggthemes", "tidyverse", 
                       "devtools", "corrplot", "cli")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("Package", pkg, "could not be installed. Check library paths and permissions."))
    }
  }
}
## installing using biocmanager 
BiocManager::install('patchwork')
# loading all required libraries
library("DEP")
library(ggrepel)
library(ggVennDiagram)
library(tidyr)
library(scales)   
library(patchwork)
library(DEqMS)
library(matrixStats)
library(MSnbase)
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(Biobase)
library(limma)
library(QFeatures)
library(msqrob2)
library(plotly)
library(gridExtra)
library(ggrepel)
library(devtools)
library(dplyr)
library(readr)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(pathview)
library(cowplot)
library(ReactomePA)
library(DOSE)
library(AnnotationHub)
gc()
####
run_app("LFQ")
## create results dir 
results_dir <- file.path("D:/MiguelW12/Documents/prot_test/imp_test_folder")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
##load data based on specfied folder with data + meta data
t_data = read.csv('D:/MiguelW12/Documents/prot_test/dat.csv', header=TRUE)
experimental_design = read.csv('D:/MiguelW12/Documents/prot_test/meta.csv', header=TRUE)

## clean uniqe peptides based on filter 
filter_rows_u <- function(df) {
  # define threshold
  threshold <- as.numeric(readline(prompt = "Enter threshold value for unique peptides: "))
  
  # find all columns containing 'Unique.peptides'
  peptide_cols <- grep("Unique\\.peptides", names(df), value = TRUE)
  if (length(peptide_cols) == 0) {
    stop("No columns containing 'Unique.peptides' found.")
  }
  
  # remove if ANY column is below threshold
  row_to_delete_any <- apply(df[, peptide_cols], 1, function(x) any(x < threshold))
  scenario1_df <- df[!row_to_delete_any, ]
  
  #remove if ≥80% of columns are below threshold
  row_to_delete_80 <- apply(df[, peptide_cols], 1, function(x) {
    mean(x < threshold) >= 0.8
  })
  scenario2_df <- df[!row_to_delete_80, ]
  
  # print outcomes
  cat("\nScenario 1: Remove rows if ANY column is below threshold")
  cat("\nRows removed:", sum(row_to_delete_any), " | Remaining rows:", nrow(scenario1_df), "\n")
  
  cat("\nScenario 2: Remove rows if ≥80% of columns are below threshold")
  cat("\nRows removed:", sum(row_to_delete_80), " | Remaining rows:", nrow(scenario2_df), "\n")
  
  # which scenario to use ? 
  choice <- readline(prompt = "Choose scenario to proceed with (1 = any below threshold, 2 = ≥80% below threshold, 0 = cancel): ")
  
  if (choice == "1") {
    return(scenario1_df)
  } else if (choice == "2") {
    return(scenario2_df)
  } else {
    cat("No filtering applied. Returning original data.\n")
    return(df)
  }
}
data_uniqe <- filter_rows_u(t_data)
## clean lfq based on filter 
filter_rows_l <- function(df) {
  # define threshold
  threshold <- as.numeric(readline(prompt = "Enter threshold value for LFQ intensity: "))
  
  # find all columns containing 'LFQ.intensity'
  peptide_cols <- grep("LFQ\\.intensity", names(df), value = TRUE)
  if (length(peptide_cols) == 0) {
    stop("No columns containing 'LFQ.intensity' found.")
  }
  
  # remove if ANY column is below threshold
  row_to_delete_any <- apply(df[, peptide_cols], 1, function(x) any(x < threshold))
  scenario1_df <- df[!row_to_delete_any, ]
  
  #remove if ≥80% of columns are below threshold
  row_to_delete_80 <- apply(df[, peptide_cols], 1, function(x) {
    mean(x < threshold) >= 0.8
  })
  scenario2_df <- df[!row_to_delete_80, ]
  
  # print outcomes
  cat("\nScenario 1: Remove rows if ANY column is below threshold")
  cat("\nRows removed:", sum(row_to_delete_any), " | Remaining rows:", nrow(scenario1_df), "\n")
  
  cat("\nScenario 2: Remove rows if ≥80% of columns are below threshold")
  cat("\nRows removed:", sum(row_to_delete_80), " | Remaining rows:", nrow(scenario2_df), "\n")
  
  # which scenario to use ? 
  choice <- readline(prompt = "Choose scenario to proceed with (1 = any below threshold, 2 = ≥80% below threshold, 0 = cancel): ")
  
  if (choice == "1") {
    return(scenario1_df)
  } else if (choice == "2") {
    return(scenario2_df)
  } else {
    cat("No filtering applied. Returning original data.\n")
    return(df)
  }
}
data_uniqe2 <- filter_rows_l(data_uniqe)

##save intermediate filtered data and check 
write_csv(as.data.frame(data_uniqe2), file.path(results_dir, "step1_prot_proc.csv"))

# remove NA values
data_final <- drop_na(data_uniqe2)

## remove duplicated gene names and make sure they do not exist- should see FALSE printed

data_final$Gene.names %>% duplicated() %>% any()
data_final %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
data_final2 <- make_unique(data_final, "Gene.names", "Protein.IDs", delim = ";")
data_final2$name %>% duplicated() %>% any()

##prepare the data design for proteomics analysis using the DEP package
## run the below line if you want information on the package.
??DEP
##
# get LFQ column numbers
LFQ_columns <- grep("LFQ.", colnames(data_final2)) 
data_se <- make_se(data_final2, LFQ_columns, experimental_design)
##view data design
data_se
qc_dir <- file.path(results_dir, "qc_res")
if (!dir.exists(qc_dir)) {
  dir.create(qc_dir)
}
# pre analysis processing and plots 
# this will show a barplot of proteins overlap among all samples
pdf(file = file.path(qc_dir, "freq_plot.pdf"))
plot_frequency(data_se)
dev.off()
## here we filter proteins that are identified in all replicates of at least one condition 
data_filt2 <- filter_missval(data_se, thr = 0)
## or less stringent,  proteins that are identified in 2 out of 3 replicates of at least one condition
data_filt2 <- filter_missval(data_se, thr = 1)
### barplot of identified proteins per samples
pdf(file = file.path(qc_dir, "proteins_per_sample_plot.pdf"))
plot_numbers(data_filt2)
dev.off()
#### barplot of protein overlap between samples
pdf(file = file.path(qc_dir, "proteins_overlap_per_sample_plot.pdf"))
plot_coverage(data_filt2)
dev.off()

##### normalize data uisng variance stabilizing transformation 
data_norm <- normalize_vsn(data_filt2)
##### examine normalization by boxplots before and after normalization
pdf(file = file.path(qc_dir, "before_and_after_normdata_per_sample.pdf"))
plot_normalization(data_filt2, data_norm)
dev.off()

###### prior to imputations for missing values this will create a heatmap of proteins with missing values
pdf(file = file.path(qc_dir, "heatmap_missing_values_per_sample.pdf"))
plot_missval(data_filt2)
dev.off()
###### intensity distributions and cumulative fraction of proteins with and without missing values
pdf(file = file.path(qc_dir, "intensity_distributions_per_sample.pdf"))
plot_detect(data_filt2)
dev.off()
## imputation & statistics directory 
stats_dir <- file.path(results_dir, "stats_res")
if (!dir.exists(stats_dir)) {
  dir.create(stats_dir)
}

##imputations
## different strategies
#for data where proteins are not quantified in specific conditions (e.g. in the control samples)-  missing not at random (MNAR)
##MNMR
#impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- DEP::impute(data_norm, fun = "MinProb", q = 0.01)
plot_imputation(data_filt2, data_imp)
#impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- DEP::impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)
plot_imputation(data_norm, data_imp_man)
# for data where proteins are quantified in some replicates but not in others. imputing by missing at random (MAR) is preferable.  
#impute missing data using the k-nearest neighbour approach (for MAR)
#MAR
data_imp_knn <- DEP::impute(data_norm, fun = "knn", rowmax = 0.9)
plot_imputation(data_norm, data_imp_knn)
data_imp <- data_imp_knn
## if you want to save imputed plot 
pdf(file = file.path(stats_dir, "imputation_plot.pdf"))
plot_imputation(data_norm, data_imp)
dev.off()

##statistical analysis+plots
# comparisions are based on protein-wise linear models combined with empirical Bayes statistics 
## testing every sample against control
data_diff <- test_diff(data_imp , type = "control", control = "ctrl")
#test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_imp, type = "all")
#test manually defined comparisons
data_diff_manual <- test_diff(data_imp, type = "manual", 
                              test = c("sample1_vs_ctrl", "sample2_vs_ctrl"))
## or 
data_diff_manual <- test_diff(data_imp, type = "manual", 
                              test = c("ctrl_vs_sample1"))
## add rejections is optional as needed
dep <- add_rejections(data_diff_manual, alpha = 0.1, lfc = log2(1))
## no rejecetion threshold
dep <- add_rejections(data_diff_all_contrasts)

## plot a frequency plot of significant proteins for the different conditions
pdf(file = file.path(stats_dir, "significant_overlaps_plot.pdf"))
plot_cond(dep)
dev.off()

stats_dir_main_file <- file.path(stats_dir, "main_file")
if (!dir.exists(stats_dir_main_file)) {
  dir.create(stats_dir_main_file)
}

##results TABLES
data_results <- get_results(dep)
##output_dataframe general
write_csv(as.data.frame(data_results), file.path(stats_dir_main_file, "data_results_all_comparisions.csv"))
## long and wide structures
df_long <- get_df_long(dep)
df_wide <- get_df_wide(dep)
write_csv(as.data.frame(df_long), file.path(stats_dir, "df_long_all_main_comparisions.csv"))
write_csv(as.data.frame(df_wide), file.path(stats_dir, "df_wide_all_main_comparisions.csv"))


stats_plot_dir <- file.path(results_dir, "stats_plots")
if (!dir.exists(stats_plot_dir)) {
  dir.create(stats_plot_dir)
}
##plot specific protein/proteins
pdf(file = file.path(stats_plot_dir, "xx_proteins_plot.pdf"))
plot_single(data_diff_all_contrasts, proteins = c("CYP51A1", "COPE"), type = "centered")
dev.off()

# volcanos for each group comparision 
venn_plot_dir <- file.path(stats_plot_dir, "venns")
if (!dir.exists(venn_plot_dir)) {
  dir.create(venn_plot_dir)
}
create_volcano_plots <- function(input_folder, output_folder) {
  
  # choose column suffix and cutoffs
  p_suffix <- readline(prompt = "Choose p-value column suffix (_p.val or _p.adj): ")
  p_suffix <- trimws(p_suffix)
  
  p_cutoff <- as.numeric(readline(prompt = "enter p-value cutoff: "))
  fc_cutoff <- as.numeric(readline(prompt = "enter log2 fold-change cutoff: "))
  
  # list csvs in the folder
  csv_files <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
  
  for (csv_file in csv_files) {
    #read each csv
    data <- read_csv(csv_file, show_col_types = FALSE)
    
    #all comparisons by finding *_ratio columns
    ratio_cols <- grep("_ratio$", names(data), value = TRUE)
    comparisons <- sub("_ratio$", "", ratio_cols)
    
    for (comp in comparisons) {
      ratio_col <- paste0(comp, "_ratio")
      pval_col <- paste0(comp, p_suffix)
      
      # protrct for columns that are missing
      if (!all(c(ratio_col, pval_col) %in% names(data))) {
        next
      }
      
      # prep df for plotting
      df_plot <- data[, c("name", ratio_col, pval_col)]
      names(df_plot) <- c("name", "ratio", "pval")
      
      #points based on thresholds
      df_plot$color <- ifelse(df_plot$pval < p_cutoff & df_plot$ratio > fc_cutoff, "Upregulated",
                              ifelse(df_plot$pval < p_cutoff & df_plot$ratio < -fc_cutoff, "Downregulated", "Not Significant"))
      
      # Volcano
      p <- ggplot(df_plot, aes(x = ratio, y = -log10(pval), color = color)) +
        geom_point() +
        geom_label_repel(
          data = subset(df_plot, color %in% c("Upregulated", "Downregulated")),
          aes(label = name),
          box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50',
          size = 4, hjust = 0, nudge_y = 0.1, max.overlaps = 1000
        ) +
        theme_minimal() +
        labs(
          title = comp,
          x = "log2 Fold Change",
          y = paste0("-log10(", p_suffix, ")"),
          color = "Legend"
        ) +
        theme(
          plot.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)
        ) +
        guides(color = guide_legend(override.aes = aes(label = "O")))
      
      # Save 
      output_file <- file.path(output_folder, paste0(comp, ".pdf"))
      ggsave(output_file, plot = p, width = 12, height = 8, dpi = 600)
    }
  }
}
input_folder_path <- stats_dir_main_file
output_folder_path <- venn_plot_dir
create_volcano_plots(input_folder_path, output_folder_path)

## venn diagrams
## here you should based on pval/adjpval, which group comparisions to show, 
plot_sig_venn <- function(results_dir, data_results) {
  #suffix and thresh prompts
  suffix <- readline(prompt = "Enter column suffix ( _p.val or _p.adj'): ")
  threshold <- as.numeric(readline(prompt = "Enter significance threshold (e.g. 0.05): "))
  
  #columns ending with suffix
  cols <- grep(paste0(suffix, "$"), names(data_results), value = TRUE)
  if (length(cols) == 0) {
    stop(paste("No columns ending with", suffix, "found."))
  }
  
  # list of significant proteins per comparison
  sig_lists <- list()
  for (col in cols) {
    comp_name <- sub(paste0("_", suffix, "$"), "", col)  # remove suffix for title
    sig_prots <- data_results[data_results[[col]] < threshold, "name"]
    sig_lists[[comp_name]] <- sig_prots
  }
  
  #random colors
  random_colors <- grDevices::colors()[sample(1:657, length(sig_lists))]
  
  # venn
  venn_plot <- ggVennDiagram(sig_lists, label_alpha = 0.1) +
    scale_fill_gradient(low = sample(random_colors, 1), high = sample(random_colors, 1)) +
    theme(legend.position = "none")
  
  # Save
  filename <- paste0("venn_sig_", suffix, "_", threshold, ".pdf")
  filepath <- file.path(results_dir, filename)
  
  pdf(file = filepath)
  print(venn_plot)
  dev.off()
  
  cat("venn diagram saved to:", filepath, "\n")
  
  print(venn_plot)
  return(list(sig_lists = sig_lists, venn_plot = venn_plot))
}
plot_sig_venn(output_folder_path, data_results)

??rda
general_plots_dir <- file.path(results_dir, "general_plots")
if (!dir.exists(general_plots_dir)) {
  dir.create(general_plots_dir)
}
## some general plotting n- equals number of proteins point size is used for visual clarity
pdf(file = file.path(general_plots_dir, "pca_plot.pdf"))
plot_pca(dep, x = 1, y = 2, n = 92, point_size = 4)
dev.off()
## correlation plot between samples plotting use signficance to focus on significant proteins
pdf(file = file.path(general_plots_dir, "cor_plot.pdf"))
plot_cor(dep, significant = FALSE, lower = 0, upper = 1, pal = "Reds")
dev.off()
## heatmap plot between samples plotting, k means clustering 
pdf(file = file.path(general_plots_dir, "heatmap_plot1.pdf"))
heat <- DEP::plot_heatmap(dep, type = "centered", kmeans = TRUE, 
                          k = 3, col_limit = 3, show_row_names = TRUE,
                          indicate = c("condition", "replicate"))
dev.off()
## another options for heatmap 
#heatmap of all significant proteins (rows) and the tested contrasts (columns)
pdf(file = file.path(general_plots_dir, "heatmap_plot2.pdf"))
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = FALSE)
dev.off()



enrich_dir <- file.path(results_dir, "enrichment")
if (!dir.exists(enrich_dir)) {
  dir.create(enrich_dir)
}
##enrichment and gene set enrihmcnet 
## GSEA
## choose the groups to compare on 
original_gene_list <- data_results$sample2_vs_ctrl_ratio
#process
names(original_gene_list) <- data_results$name
gene_list <- na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
head(gene_list)
gene_list = na.omit(gene_list)

## run - yuo can choose ont- as CC,MF,BP or ALL of them, you can also swap the padj method.
gsea <- gseGO(geneList=gene_list, 
                 ont ="ALL", 
                 keyType = "SYMBOL", 
                 minGSSize = 5, 
                 maxGSSize = 1000, 
                 pvalueCutoff = 1, 
                 verbose = TRUE, 
                 OrgDb = 'org.Hs.eg.db',
                 pAdjustMethod = "fdr")

#save file
write.csv(as.data.frame(gsea), file.path(enrich_dir, "gseago_ctrlvssamp2.csv"))

##visulazation options
##1
pdf(file = file.path(enrich_dir, "dotplot_gsea.pdf"))
dotplot(gsea, showCategory=30, font.size = 8,label_format = 60, color = "p.adjust", title = "dotplot_ctrlvs2", split=".sign") + facet_grid(.~.sign)
dev.off()
##2
ema <- pairwise_termsim(gsea, method="JC", semData = NULL, showCategory = 30)
pdf(file = file.path(enrich_dir, "emapplot_gsea.pdf"))
emapplot(ema)
emapplot(ema, cex_category=1, cex_label_category=0.8, cex_line=0.8, color="pvalue", layout="kk" )
dev.off()

##3
pdf(file = file.path(enrich_dir, "treeplot_gsea.pdf"))
treeplot(ema, showCategory = 30, nWords = 4, fontsize = 2, nCluster = 5)
dev.off()
##4
pdf(file = file.path(enrich_dir, "heat_plot_gsea.pdf"))
edo <- setReadable(gse_glioma, 'org.Hs.eg.db', 'ENSEMBL')
dev.off()
heatplot(edo,  foldChange=gene_list, label_format = 8, showCategory=5)
##5

p1 <- cnetplot(edo, foldChange=gene_list)
p2 <- cnetplot(edo, categorySize="p.adjust", foldChange=gene_list)
p3 <- cnetplot(edo, foldChange=gene_list, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
pdf(file = file.path(enrich_dir, "net_plot_gsea.pdf"))
cnetplot(edo, node_label="all", categorySize="p.adjust", cex_label_gene = 0.6, cex_label_category = 1,
         cex_category = 1, cex_gene = 0.5, colorEdge = TRUE, showCategory=30) 

dev.off()

###Enrichment 
#again, choose comparisons
data_results$sample1_vs_ctrl_p.adj
#process
gene_list <- data_results$sample2_vs_ctrl_ratio
names(gene_list) <- data_results$name
gene_list_2<-na.omit(gene_list)
gene_list_2 = sort(gene_list_2, decreasing = TRUE)
## adjust p value threshold as needed
sig_genes = subset(data_results, sample2_vs_ctrl_p.adj < 0.1)
genes <- sig_genes$sample2_vs_ctrl_ratio
names(genes) <- sig_genes$name
genes <- na.omit(genes)
head(names(gene_list_2))
head(genes)
dim(sig_genes)
names(gene_list_2)
dim(gene_list_2)
## adjust log 2 foldchange as needed
genes <- names(genes)[abs(genes) > 0.3]
dim(genes)
##option to split for up/down regulated proteins
genes_up <- names(genes)[abs(genes) > 0.3]
genes_down <- names(genes)[abs(genes) < -0.3]

## choose ont, p value adjustment method etc. 
go_enrich <- enrichGO(gene = genes,
                        universe = names(gene_list_2),
                        OrgDb = 'org.Hs.eg.db', 
                        keyType = 'SYMBOL',
                        readable = T,
                        ont = "ALL",
                        pvalueCutoff = 1, 
                        pAdjustMethod = "none")
#
write.csv(as.data.frame(go_enrich), file.path(enrich_dir, "go_enrich_ctrlvs2.csv"))
##visulazation options
#1
pdf(file = file.path(enrich_dir, "heat_plot_enrich.pdf"))
edo <- setReadable(go_enrich_t, 'org.Hs.eg.db', 'ENSEMBL')
heatplot(edo,  foldChange=gene_list, label_format = 8, showCategory=20)
dev.off()
#2
pdf(file = file.path(enrich_dir, "upset_plot_enrich.pdf"))
upsetplot(go_enrich_t) 
dev.off()
#3
pdf(file = file.path(enrich_dir, "tree_plot_enrich.pdf"))
treeplot(go_enrich_t, showCategory = 30, nWords = 4, fontsize = 2, nCluster = 5)
dev.off()
