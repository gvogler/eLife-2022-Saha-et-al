
# 
# ## ----Libraries, echo=TRUE, message=FALSE, results=FALSE-------------------------
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")
# BiocManager::install("tximportData")
# BiocManager::install("sparseMatrixStats")
# BiocManager::install("DEGreport")
# BiocManager::install("apeglm")
# BiocManager::install("biomaRt")
# BiocManager::install("GO.db")
# BiocManager::install("GOstats")
# BiocManager::install("compEpiTools")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("pathview")
# BiocManager::install("vsn")


# Validate Bioconductor:
# BiocManager::valid()
library(ggrepel)
library(cowplot)
library(ggplotify)
library(gridExtra)
library(RColorBrewer)
library(pheatmap)
library(pathview)
library(KEGGREST)
library(org.Hs.eg.db)
library(compEpiTools)
library(GOplot)
library(GO.db)
library(GOstats)
library(biomaRt)
library(ggpubr)
library(ggplot2)
library(DESeq2)
library(tximport)
library(tximportData)
library(sparseMatrixStats)
library(DEGreport)
library(AnnotationDbi)
library(gprofiler2)
library(tidyverse)


## ---- message=FALSE, cache = TRUE-----------------------------------------------
salmon_merged.gene_counts <- readRDS("RNAseq VCM Perrin Colas/salmon.merged.gene_counts.rds")
head(salmon_merged.gene_counts) # only first four experiments shown


## ---- message=FALSE-------------------------------------------------------------
raw_counts <- assays(salmon_merged.gene_counts)$counts
head(raw_counts)[,1:6]


## ---- message=FALSE-------------------------------------------------------------
subset_cols <- names(raw_counts)[c(1:3, 7:9)]
subset_data <- raw_counts[,subset_cols]
head(subset_data)


## -------------------------------------------------------------------------------
subset_data <- round(subset_data)
head(subset_data)


## ---- message=FALSE-------------------------------------------------------------
ExpDesign <- data.frame(row.names=colnames(subset_data), condition=c("ctrl_VCM", "ctrl_VCM", "ctrl_VCM", "PAX9_VCM", "PAX9_VCM", "PAX9_VCM"), libType = c("paired-end","paired-end","paired-end","paired-end","paired-end","paired-end"))
ExpDesign


## ---- collapse = TRUE-----------------------------------------------------------
dds <- DESeq2::DESeqDataSetFromMatrix(countData = subset_data, colData = ExpDesign, design = ~ condition)
dds


## ---- collapse = TRUE-----------------------------------------------------------
dds <- DESeq2::DESeq(dds)

## ---- echo=FALSE, fig.cap="Figure 1. Gene counts before and after normalization."----
log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1)
ggplot(data = data.frame(log2(counts(dds))[,1:2]), aes(x = CTRL_R1, y = CTRL_R2)) +
  geom_point(cex = 0.1, color = "red") +
  geom_abline(intercept = 0, slope = 1, color = "green", cex = 0.5) +
  geom_point(data = data.frame(log.norm.counts[,1:2]), aes(x = CTRL_R1, y = CTRL_R2), color = "blue", cex = 0.1, alpha = 0.2) 


## ---- fig.cap="Figure 2. Dispersion Estimates."---------------------------------
plotDispEsts(dds)


## ----PCA------------------------------------------------------------------------
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "libType"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
plot1 <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=libType)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


## ----Heatmap Clustering, results='hold'-----------------------------------------
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
plot2 <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, silent = TRUE)


## ----PCA and Heatmap Figure, echo=FALSE, fig.cap="Figure 4. PCA and Heatmap Cluster"----
plot_grid(plot1, as.ggplot(plot2), align = "h", nrow = 1, rel_widths = c(1.5,1))


## -------------------------------------------------------------------------------
#### Run these lines if the Rdata file is not present
# ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
# geneset_ensembl <- getBM(attributes = c('entrezgene_id', 'external_gene_name', 'ensembl_gene_id'),  filters = 'ensembl_gene_id', values = row.names(res), mart = ensembl)
# saveRDS(geneset_ensembl, file = "geneset_ensembl.Rdata")
geneset_ensembl <- readRDS(file = "geneset_ensembl.Rdata")



## ---- collapse = TRUE-----------------------------------------------------------
res <- DESeq2::results(dds)


## ---- collapse=TRUE-------------------------------------------------------------
res1 = results(dds, contrast=c("condition","ctrl_VCM", "PAX9_VCM"))



## ---- collapse = TRUE-----------------------------------------------------------
ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
resLRT <- results(ddsLRT)

# Add Entrez/Symbol annotation to the results
counts_table <- as.data.frame(counts(dds))
counts_table$Symbol <- geneset_ensembl$external_gene_name[match(row.names(counts_table), geneset_ensembl$ensembl_gene_id)]


res$Symbol <- geneset_ensembl$external_gene_name[match(row.names(res), geneset_ensembl$ensembl_gene_id)]
res$entrez <- geneset_ensembl$entrezgene_id[match(row.names(res), geneset_ensembl$ensembl_gene_id)]

resLRT$Symbol <- geneset_ensembl$external_gene_name[match(row.names(resLRT), geneset_ensembl$ensembl_gene_id)]
resLRT$entrez <- geneset_ensembl$entrezgene_id[match(row.names(resLRT), geneset_ensembl$ensembl_gene_id)]

head(resLRT)

## ---- echo = TRUE, message=FALSE, cache = TRUE----------------------------------
res1LFC <- lfcShrink(dds, coef="condition_PAX9_VCM_vs_ctrl_VCM", type="apeglm")
dds$condition <- relevel(dds$condition, ref = "PAX9_VCM")
dds <- nbinomWaldTest(dds) # <- to updated conditions



## ---- include=FALSE-------------------------------------------------------------
ggma_plot_res <- function(x)
  {
  z <- x@elementMetadata[4,2]
  heading_plot <- substr(z, gregexpr("value: condition", z)[[1]][1]+17, nchar(z))
  ggmaplot(x, main = heading_plot,
         fdr = 0.05, fc = 2, size = 1, alpha = 0.5,
         palette = c("#B31B21", "#1465AC", "gray80"),
         genenames = as.vector(x$name),
         legend = "top", top = 0,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal()) + theme(plot.title = element_text(size = rel(1.5), hjust = 0.5, face = "bold"), text = element_text(size=12))
}
  
volcano_plot <- function(x)
  {
  z <- x@elementMetadata[4,2]
  heading_plot <- substr(z, gregexpr("value: condition", z)[[1]][1]+17, nchar(z))
  x$symbol <- geneset_ensembl$external_gene_name[match(row.names(x), geneset_ensembl$ensembl_gene_id)]
  
  res_tableOE_tb <- as_tibble(x, rownames = "gene") %>% 
    mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 1) # 1.7-fold change
  res_tableOE_tb <- res_tableOE_tb %>% arrange(padj) %>% mutate(genelabels = "")
  res_tableOE_tb$genelabels[1:10] <- res_tableOE_tb$symbol[1:10]
  
  ggplot(res_tableOE_tb, aes(x = log2FoldChange, y = -log10(padj))) +
    scale_color_manual(values=c("#25B025", "#9D1856")) +
    geom_point(aes(colour = threshold_OE), size = 0.5) +
    geom_text_repel(aes(label = genelabels), point.size = NA) +
    ggtitle(heading_plot) +
    xlab("log2 fold change") + 
    ylab("-log10 adjusted p-value") +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5, face = "bold"),
          axis.title = element_text(size = rel(1.25))) 
  
  } 


## ---- echo=FALSE, warning=FALSE, fig.asp = .62, figures-side1, fig.show="hold", out.width="50%", fig.cap="Figure 4A. Log2-fold changes before (left) and after log fold change shrinkage (right)."----

plot1 <- ggma_plot_res(res1)
plot2 <- volcano_plot(res1)
plot3 <- ggma_plot_res(res1LFC)
plot4 <- volcano_plot(res1LFC)
par(mfrow=c(1,4), mai = c(1, 0.1, 0.1, 0.1))
plot1
plot2
plot3
plot4


## ---- echo=FALSE, warning=FALSE, figures-side2, fig.show="hold", out.width="50%", fig.cap="Figure 4B. Log2-fold changes before (left) and after log fold change shrinkage (right)."----

plot1 <- ggma_plot_res(res2)
plot2 <- volcano_plot(res2)
plot3 <- ggma_plot_res(res2LFC)
plot4 <- volcano_plot(res2LFC)
par(mar = c(4, 4, 4, 4))
plot1
plot2
plot3
plot4


## ---- echo=FALSE, warning=FALSE, figures-side3, fig.show="hold", out.width="50%", fig.cap="Figure 4C. Log2-fold changes before (left) and after log fold change shrinkage (right)."----

plot1 <- ggma_plot_res(res3)
plot2 <- volcano_plot(res3)
plot3 <- ggma_plot_res(res3LFC)
plot4 <- volcano_plot(res3LFC)
par(mar = c(4, 4, 4, 4))
plot1
plot2
plot3
plot4

