## Analysis for EPICV2 array
##libraries
library(Rtools)
library(ggplot2)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(minfi)
library(minfiData)
library(sesame)
library(ExperimentHub)
library(parallel)
library(pals)
library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(missMethyl)
library(Gviz)
library(DMRcate)
library(stringr)
library(devtools)
library(minfiDataEPIC)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(quantro)
library(PCAtools)
library(methylGSA)
library(pheatmap)
library(creditmodel)
library(ggplot2)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
#####
gc()
##create data dir - uri edit for example
dataDirectory <- "D:/MiguelW12/Documents/GSD1a_met_minfi -t_out"
### create annotation data
annEPICv2 = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
head(annEPICv2)
#list files in dir
list.files(dataDirectory, recursive=TRUE)
#meta+data prep - uri edit for example
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet.csv")
targets
rgSet_data <- read.metharray.exp(targets=targets)
rgSet_data
##add EPICV2 annotation
annotation(rgSet_data)
rgSet_data@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
rgSet_data
##rename samples
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet_data) <- targets$ID
rgSet_data
###proc+qc
### in addition to Genome_Studio QC and BACR - see link at the beggining of script -  - uri edit for example
##outlier samples?
detP_data <- detectionP(rgSet_data)
head(detP_data)
dim(detP_data)
dim(rgSet_data)
##plot
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP_data), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

barplot(colMeans(detP_data), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
       bg="white")

## generate report
qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="qcReport-gsd1a-real.pdf")

## any samples need to be excluded? 
keep1 <- colMeans(detP_data) < 0.05
rgSet_data <- rgSet_data[,keep1]
rgSet_data
targets <- targets[keep1,]
detP_data <- detP_data[,keep1]
dim(detP_data)
dim(rgSet_data)

## normalization- qunatro analysis
p <- getBeta(rgSet, offset = 100)
pd <- pData(rgSet)
dim(p)
head(pd)

##plots
matdensity(p, groupFactor = pd$Sample_Group, xlab = "Beta", ylab = "density",
           main = "Beta Values", brewer.n = 8, brewer.name = "Dark2")
legend('top', c("HC", "GSD1A"), col = c(1, 2), lty = 1, lwd = 3)

matboxplot(p, groupFactor = pd$Sample_Group, xaxt = "n", main = "Beta Values")
legend('top', c("HC", "GSD1A"), col = c(1, 2), lty = 1, lwd = 3)

## testing
qtest <- quantro(object = p, groupFactor = pd$Sample_Group)
qtest
summary(qtest)
##
anova(qtest)
quantroStat(qtest)
## add perm
qtest_p <- quantro(object = p, groupFactor = pd$Sample_Group, B = 1000)
qtest_p
quantroPlot(qtest_p)
## Implement normalization
NSetSqt_data <- preprocessFunnorm(rgSet_data)
##no norm to compare
RSetsq_data <- preprocessRaw(rgSet_data)
## examine data
dim(rgSet_data)
dim(NSetSqt_data)
dim(RSetsq_data)
dim(detP_data)
##sex idenitification
NSetSqt_data$predictedSex
plotSex(NSetSqt_data)
## plot norm data - bef and aft
par(mfrow=c(1,2))
densityPlot(getBeta(RSetsq_data), sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))

densityPlot(getBeta(NSetSqt_data), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))

# MDS plots to look at largest sources of variation
par(mfrow=c(1,3))
plotMDS(getM(NSetSqt_data), top=900000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

##other components
par(mfrow=c(1,3))
plotMDS(getM(NSetSqt_data), top=900000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(NSetSqt_data), top=900000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(NSetSqt_data), top=900000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(5,4))
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")


## filtering low quality probes 
# ensure probes are in the same order in the  objects
dim(detP_data)
dim(NSetSqt_data)
length(featureNames(NSetSqt_data))
detP_data <- detP_data[match(featureNames(NSetSqt_data),rownames(detP_data)),]
dim(detP_data)
dim(NSetSqt_data)
# remove any probes that have failed in one or more samples
keep2 <- rowSums(detP_data < 0.01) == ncol(NSetSqt_data)
table(keep2)
NSetSqFlt_data <- NSetSqt_data[keep2,]
NSetSqFlt_data
dim(NSetSqFlt_data)

#remove probes on the sex chromosomes
keep3 <- !(featureNames(NSetSqFlt_data) %in% annEPICv2$Name[annEPICv2$chr %in%
                                                        c("chrX","chrY")])
table(keep3)
NSetSqFlt_data <- NSetSqFlt_data[keep3,]
dim(NSetSqFlt_data)
##SNPs removal
NSetSqFlt_data <- dropLociWithSnps(NSetSqFlt_data)
NSetSqFlt_data
dim(NSetSqFlt_data)

#post filter MDS
par(mfrow=c(1,3))
plotMDS(getM(NSetSqFlt_data), top=900000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], cex=0.8)
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.65, bg="white")


par(mfrow=c(1,3))
#other components
plotMDS(getM(NSetSqFlt_data), top=900000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(1,3))
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(NSetSqFlt_data), top=900000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(NSetSqFlt_data), top=900000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(4,5))
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")


# PCA
mVals <- getM(NSetSqFlt_data)
head(mVals[,1:5])
p <- pca(mVals, metadata = colData(NSetSqFlt_data), removeVar = 0.1)
plotloadings(p, labSize = 3)
sam_pc <- p$rotated
write.csv(as.data.frame(sam_pc), file="samo_pc.csv")
load_pc <- p$loadings
write.csv(as.data.frame(load_pc), file="load_pc.csv")

horn <- parallelPCA(bVals_data)
horn$n
h <- parallelPCA(mVals_data)
h$n
elbow <- findElbowPoint(p$variance)
elbow
screeplot(p,
          components = getComponents(p, 1:20),
          vline = c(horn$n, elbow)) +
  
  geom_label(aes(x = horn$n + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))

which(cumsum(p$variance) > 80)[1]

biplot(p,
       colby = 'Sample_Group', colkey = c('HC' = 'forestgreen', 'GSD1A' = 'purple'),
       colLegendTitle = 'Group',
       # ellipse config
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.95,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 0,
       ellipseFillKey = c('HC' = 'yellow', 'GSD1A' = 'pink'),
       xlim = c(-2000,2000), ylim = c(-2000, 2000),
       hline = 0, vline = c(-25, 0, 25), pointSize = 1,
       legendPosition = 'left', legendLabSize = 16, legendIconSize = 8.0,
       shape = 'Sample_Group', shapekey = c('HC'=15, 'GSD1A'=17),
       title = 'PCA bi-plot',
       subtitle = 'PC1 versus PC2',
       caption = '5 PCs ≈ 80%,
       circle-95% CI')


######Extract beta + m values 
mVals_data <- getM(NSetSqFlt_data)
head(mVals_data[,1:5])
dim(mVals_data)

bVals_data <- getBeta(NSetSqFlt_data)
head(bVals_data[,1:5])

##plot
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values",
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values",
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))


######DMP's
#factor of interest
targets$Sample_Group
cellType <- factor(targets$Sample_Group)
# individual effect( not needed) 
#individual <- factor(targets$)
##design
##option1
design <- model.matrix(~0+cellType, data=targets)
colnames(design) <- c(levels(cellType))
design
##option2
Group <- factor(targets$Sample_Group, levels=c("HC","GSD1A"))
des <- model.matrix(~Group)
colnames(des) <- c("HC","GSD1AvsHC")
des

# fit the linear model
fit1 <- lmFit(mVals_data, design)
# contrast matrix for specific comparisons
contMatrix <- makeContrasts(HC-GSD1A,
                            levels=design)
contMatrix
fit2 <- contrasts.fit(fit1, contMatrix)
fit3 <- eBayes(fit2)

#numbers of DM CpGs at FDR < 0.1 -BH correction 
summary(decideTests(fit3, p.value = 0.1, adjust.method = "BH"))

#table of results for the contrast with annotation
annEPICv2kSub <- annEPICv2[match(rownames(mVals_data),annEPICv2$Name),
                             c(1:4,12:19,24:ncol(annEPICv2))]
DMP_data <- topTable(fit3,  num=Inf, coef=1, genelist=annEPICv2kSub)
head(DMP_data)
write.table(DMP_data, file="DMP_data_GSD1A.csv", sep=",", row.names=FALSE)
# plot significantly differentially methylated CpGs
par(mfrow=c(2,2))
sapply(rownames(DMP_data)[1:4], function(cpg){
  plotCpg(bVals_data, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")
})
##plot specfic cpg's by location in frame 
cpg = rownames(DMP_data)[30893:30894]
plotCpg(bVals_data, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")

## volcano
#Function to create a volcano plot- res_df
create_volcano_plot <- function(results_df, p_value_cutoff = 0.1, logfc_cutoff = 1.5) {
  #results based on p-value and logFC 
  significant_results <- subset(results_df, abs(logFC) > logfc_cutoff & P.Value < p_value_cutoff)
  
  volcano_plot <- ggplot(results_df, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(color = "black", alpha = 0.6) +
    geom_point(data = significant_results, color = "forestgreen", alpha = 0.8) +  
    geom_hline(yintercept = -log10(p_value_cutoff), linetype = "dashed", color = "blue") +  
    geom_vline(xintercept = c(logfc_cutoff, -logfc_cutoff), linetype = "dashed", color = "blue") +  
    labs(title = "Volcano Plot",
         x = "logFC",
         y = "-log10(P.Value)",
          ) +
    scale_color_manual(name='Legend',
                       breaks=c('Significant', 'Not Significant'),
                       values=c('Significant'='forestgreen', 'Not Significant'='black'))+
    theme(legend.title=element_text(size=20),
                  legend.text=element_text(size=14),
                  legend.position = "right")
  
  return(volcano_plot)
}
# run 
volcano_plot <- create_volcano_plot(DMP_data, p_value_cutoff = 0.1, logfc_cutoff = 1.5)
print(volcano_plot)
## clear space
gc()
#### DMR analysis
myAnnotation <- cpg.annotate(object = NSetSqFlt_data, datatype = "array", what = "M",
                             analysis.type= "differential", design = design,
                             contrasts = TRUE, cont.matrix = contMatrix, fdr = 0.1,
                             coef = "HC - GSD1A")

DMRs <- dmrcate(myAnnotation, lambda=1000, C=2, pcutoff = 0.1)
results.ranges <- extractRanges(DMRs, genome = "hg38")
head(results.ranges)
##save res
write.table(results.ranges, file="results.ranges_true.csv", sep=",", row.names=FALSE)
##plot
groups <- c(GSD1A="firebrick1", HC="dodgerblue2")
type <- factor(targets$Sample_Group)
cols <- groups[as.character(type)]
cols
DMR.plot(ranges=results.ranges, dmr=1, CpGs=NSetSqFlt_data, what="Beta", pch=12, toscale=TRUE, plotmedians=TRUE,
         phen.col=cols, genome="hg38")

##methylRRA gene set enrichment analysis

#remove characters after "_"- designed to adjust data to EPICV2 structure
remove_after_underscore <- function(text) {
  gsub("_.*", "", text)
}
# apply on original column to create modified column for EPIVv2
DMP_data$N_name <- remove_after_underscore(DMP_data$Name)
dim(DMP_data)
head(DMP_data$N_name)
###
gt <- unlist(split(as.numeric(DMP_data$adj.P.Val), as.character(DMP_data$N_name)))

######### methylRRA gsea 
##gsea
res = methylRRA(cpg.pval = gt, method = "GSEA", GS.type = "GO",
                 minsize = 2, maxsize = 1000)
write.csv(as.data.frame(res), file="res-gsd1a_met_GSEA.csv")
####
dim(res)
type(res)
##vis
##GSEA PLOT
res_gsea <- subset(res, ID %in% c('GO:0005667',
                               'GO:0000779',
                               'GO:0016055',
                               'GO:0031902',
                               'GO:0042826',
                               'GO:0006638',
                               'GO:0003333',
                               'GO:0030203',
                               'GO:0071837',
                               'GO:0006639',
                               'GO:0010008',
                               'GO:0005770',
                               'GO:0016197',
                               'GO:0006641',
                               'GO:0032418',
                               'GO:0045913',
                               'GO:0035065',
                               'GO:0005741',
                               'GO:0005776',
                               'GO:0010634',
                               'GO:0008170',
                               'GO:0016667',
                               'GO:1902686'))


##add any other prefered visulazitions from clusterprofiler 
e <- barplot(res_gsea, num = 23, colorby = "padj", xaxis = , 
              title = "enrichment-Gsea")
e

#### clustering 
df <- as.data.frame(bVals_data)
df <- t(df)
scaled_data <- scale(df)
clust <- pheatmap(scaled_data, cluster_cols = F, cluster_rows = T, main = "Cluster_Heatmap-EPIC", kmeans_k = 2,
                  fontsize = 8, cutree_rows = 2,
                  show_colnames = FALSE)
clust
t <- clust$kmeans$cluster
write.csv(as.data.frame(t), file.path(output_dir, "clustered_data$kmeans.csv"), row.names = FALSE)


##### scatter_regression analysis -run per two samps
##
df <- as.data.frame(bVals_data)
write.csv(as.data.frame(bVals_data), file="EPIC-bVals_data.csv")
table(df$HC.1, df$HC.2, df$HC.3, df$HC.4)
table(df$GSD1A.5, df$GSD1A.6, df$GSD1A.7, df$GSD1A.8)
## for each pair of samples
df_plot <- df %>% select(GSD1A.7, GSD1A.8)
df_plot_adj <- df_plot * 100
sample1 <- "GSD1A.7"
sample2 <- "GSD1A.8"

ggplot(df_plot_adj,aes(x = GSD1A.7, y = GSD1A.8)) + 
  geom_point(size=1) + 
  labs(title = paste("Feature Comparison:", sample1, "vs", sample2),
       x = paste("Features of", sample1),
       y = paste("Features of", sample2)) +
  geom_smooth(method = "lm", se=TRUE, ) +
  stat_regline_equation(label.y = 105, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 100, aes(label = ..rr.label..))

gc()




####end 
