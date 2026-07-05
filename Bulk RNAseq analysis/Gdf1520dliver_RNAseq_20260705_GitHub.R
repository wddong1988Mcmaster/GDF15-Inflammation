########################### GDF15 20d liver RNAseq DESeq2-salmon

rm(list = ls())
options(stringsAsFactor = F)

####Step: import salmon data

BiocManager::install("tximport")
install.packages("tidyverse")
BiocManager::install("tximportData")
BiocManager::install("ggplot2")
install.packages("ggplot2")
install.packages("htmltools")

library(GEOquery)
library(dplyr)
library(pathview)
library("tximport")
library("readr")
library("tximportData")

library("DESeq2")
library("export")


dir <- "C:/BaiduSyncdisk/R analysis/Gdf1520dLiver" # change the folder where you doposit your data
samples <- read.table("sample_gdf15liver.txt", header=TRUE)
head (samples)
# construct condition
#samples$condition <- factor(rep(c("A","B"),each=3))
rownames(samples) <- samples$group
#samples[,c("pop","center","run","condition")]

###Next we specify the path to the files using the appropriate columns of samples,
###and we read in a table that links transcripts to genes for this dataset.

# find the files for each sample, and change the colume name to sample name 
files <- file.path(dir,"salmon", samples$run, "quant.sf")
head (files)
names(files) <- samples$group


######I am using this one******
tx2gene <- read.csv("tx2symbol.csv", header=TRUE)
head (tx2gene)

txi <- tximport(files, type="salmon", tx2gene=tx2gene)
names(txi)
head(txi)
head(txi$abundance)

write.csv(as.data.frame(txi$abundance), "1-TPM.csv")



### Finally, we can construct a DESeqDataSet from the txi object 
### and sample information in samples.
###The ddsTxi object here can then be used as dds in the following analysis steps.
BiocManager::install("DESeq2")
library("DESeq2")

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
dds<- ddsTxi
dds



######## annotation for this dataset "dds"
library ("org.Mm.eg.db")
genenames_dds <- mapIds(org.Mm.eg.db,keys = rownames(dds),column = "SYMBOL",keytype="ENTREZID")
annotation_dds <- data.frame(gene_name = genenames_dds,
                                row.names = rownames(dds),
                                stringsAsFactors = FALSE)
head(annotation_dds)
write.csv(annotation_dds, 'annotation_dds_SYMBOL.csv')




### Spre-filter low count genes before running the DESeq2 functions
### Here we perform a minimal pre-filtering to keep only rows that have at least 10 reads total.


## Metabolism Class showing for
oridata = counts(dds)
dim(oridata)

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds

write.csv(as.data.frame(counts(dds)), "1-counts.csv")


###PCA
##First we need to transform the raw count data
## Data transformations and visualization
##for visualization or clustering – it might be useful to work with transformed versions of the count data.
## rlog function might take too long, variance stabilizing transformations (VST) function will be a faster choice.
## By setting blind to FALSE
## vst function will perform variance stabilizing transformation
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)

vsd1=assay(vsd)

write.csv(as.data.frame(vsd1), "1-countsNorm.csv")
#using the DESEQ2 plotPCA fxn we can
plotPCA(vsd, intgroup="condition")

library(export)

graph2eps(x = NULL, file='PCA', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 3, bg = "transparent")

graph2svg(x = NULL, file='PCA', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 3.5, bg = "transparent")

graph2pdf(x = NULL, file='PCA', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 3.5, bg = "transparent")



###PCA2
BiocManager::install("pcaExplorer")
library("pcaExplorer")
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
rld
pcaplot(rld,intgroup = "condition",
        pcX = 1, pcY = 2, title = "",
        ellipse = T, text_labels = F)

?pcaplot

install.packages("export")
library(export)

graph2eps(x = NULL, file='PCA_1', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 4, bg = "transparent")

graph2svg(x = NULL, file='PCA_nocir', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 4, bg = "transparent")

graph2pdf(x = NULL, file='PCA', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 4, bg = "transparent")




####### Step: Differential expression analysis
### specifying the reference level
#dds$condition <-factor(dds$condition, levels = c("Vehicle","GDF15","Pair-fed"))
dds$condition <- relevel(dds$condition, ref = "Vehicle")
dds$condition
#Speed-up
library("BiocParallel")
register(MulticoreParam(4))
#Levels: Vehicle GDf15 Pair-fed, Vehcile should be in the former position
dds <- DESeq(dds)
res <- results(dds)
res
write.csv(as.data.frame(res), "DEG_results.csv")
# make sure that condition treated vs untreated


#### Contrasts: containing a factor with three levels, say A, B and C
## GDF15 vs Vehicle
res_gvv=results(dds, contrast=c("condition","GDF15","Vehicle"))
res_gvv
write.csv(as.data.frame(res_gvv), "resGvV_results.csv")

## PairFed vs Vehicle
res_pvv=results(dds, contrast=c("condition","PairFed","Vehicle"))
res_pvv
write.csv(as.data.frame(res_pvv), "resPvV_results.csv")

## GDF5 vs PairFed
res_gvp=results(dds, contrast=c("condition","GDF15", "PairFed"))
res_gvp
write.csv(as.data.frame(res_gvp), "resGvP_results.csv")





#### Log fold change shrinkage for visualization and ranking
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm")
library("apeglm")

resultsNames(dds)
## [1] "Intercept"     "condition_treated_vs_untreated"
resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")
resLFC

## GDF15 vs Vehicle
resLFC_gvv <- lfcShrink(dds, coef="condition_GDF15_vs_Vehicle", type="apeglm")
resLFC_gvv
write.csv(as.data.frame(resLFC_gvv), "resLFC_gvv_results.csv")

## PairFed vs Vehicle
resLFC_pvv <- lfcShrink(dds, coef="condition_PairFed_vs_Vehicle", type="apeglm")
resLFC_pvv
write.csv(as.data.frame(resLFC_pvv), "resLFC_pvv_results.csv")



### p-values and adjusted p-values
res_gvvOrdered <- res_gvv[order(res$pvalue),]
res_pvvOrdered <- res_pvv[order(res$pvalue),]
res_gvpOrdered <- res_gvp[order(res$pvalue),]


#summarize some basic tallies using the summary function
summary(res_gvv, res_gvv$padj < 0.1 & abs(res_gvv$log2FoldChange)>0.5)
summary(res_pvv, res_pvv$padj < 0.1 & abs(res_pvv$log2FoldChange)>0.5)
summary(res_gvp, res_gvp$padj < 0.1 & abs(res_gvp$log2FoldChange)>0.5)

summary(res_gvp, res_gvp$padj < 0.12 & abs(res_gvp$log2FoldChange)>0.25)


#How many adjusted p-values were less than 0.1?
sum(res_gvv$padj < 0.1, na.rm=TRUE)
sum(res_gvv$log2FoldChange>0.8, na.rm=TRUE)

##Need to filter on adjusted p-values, not p-values, to obtain FDR control. 10% FDR (alpha) is common
##because RNA-seq experiments are often exploratory and having 90% true positives in the gene set is ok
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)





############## Volcano Plot with label gene names.


install.packages("ggpubr")
install.packages("ggthemes")

library(ggpubr)
library(ggthemes)

###GDF15 vs Vehicle
#res_gvv1=as.data.frame(res_gvv)
#res_gvv1=as.data.frame(res_pvv)
res_gvv1=as.data.frame(res_gvp)

deg.data <-res_gvv1

head(deg.data)
deg.data$logQ <- -log10(deg.data$padj)
deg.data$log2FC=deg.data$log2FoldChange
deg.data$FDR=deg.data$padj
deg.data$ID=rownames(deg.data)

head(deg.data)

# plot basic volcano
ggscatter(deg.data, x = "log2FC", y = "logQ") + theme_base()

# set a subcolumn and set up and down regulated genes
deg.data$Group = "normal"

deg.data$Group[which( (deg.data$FDR < 0.1) & (deg.data$log2FC > 0.5) )] = "up"
deg.data$Group[which( (deg.data$FDR < 0.1) & (deg.data$log2FC < -0.5) )] = "down"

# see how many genes up and down
table(deg.data$Group)

# plot new volcano
ggscatter(deg.data, x = "log2FC", y = "logQ",
          color = "Group") + theme_base()

# change colors (palette) and plot size(size)
ggscatter(deg.data, x = "log2FC", y = "logQ",
          color = "Group", 
          palette = c("green", "gray", "red"),
          size = 1.5) + theme_base()

# 为火山图添加logP分界线条(geom_hline)和logFC分界线条(geom_vline)
ggscatter(deg.data, x = "log2FC", y = "logQ",
          color = "Group", 
          palette = c("green", "gray", "red"),
          size = 1) + theme_base() + 
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_vline(xintercept = c(-0.5,0.5), linetype="dashed")


# add a new collumn "Label"
deg.data$Label = ""

# get top 20 DEGs according to FC
deg.data <- deg.data[order(abs(deg.data$log2FC),decreasing = T), ]
log2FC.genes <- head(deg.data$ID, 20)

head(deg.data)
# get top 20 DEGs according to adjp value
# deg.data <- deg.data[order(deg.data$logQ)]
# in high expression genes，select the lowest FDR 20
deg.data <- deg.data[order(abs(deg.data$logQ),decreasing = T), ]
fdr.genes <- head(deg.data$ID, 20)

# combine log2FC.genes and fdr.genes，and put them in 'Label'
deg.top20.genes <- c(as.character(log2FC.genes), as.character(fdr.genes))

write.csv(as.data.frame(deg.top20.genes), "top20gene_gvp.csv")

deg.data$Label[match(deg.top20.genes, deg.data$ID)] <- deg.top20.genes

#print (deg.data$Label)
table (deg.data$Label)


# change color and axis labelling to make figure looking better
ggscatter(deg.data, x = "log2FC", y = "logQ",
          color = "Group", 
          palette = c("#00BA38", "#BBBBBB", "#F8766D"),
          size = 2,
          label = deg.data$Label, 
          font.label = 8, 
          repel = T,
          xlab = "log2FC", 
          ylab = "-log10(FDR)") + 
  theme_base() + 
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_vline(xintercept = c(-0.5,0.5), linetype="dashed")

library(export)
graph2svg(x = NULL, file='Volcano_gvp_label', font = "Arial", cairo = TRUE,   
          width = 6, height = 3.5, bg = "transparent")





########## Venn plotting and dotpot for DEGs ####

res_gvv=read.csv('resGvV_results.csv', header = T)
res_pvv=read.csv('resPvV_results.csv', header = T)
row.names(res_gvv) =res_gvv$X
res_gvv=res_gvv[,-1]
head(res_gvv)

row.names(res_pvv) =res_pvv$X
res_pvv=res_pvv[,-1]
head(res_pvv)


# Apply cutoff to each result set
sig_gvv <- subset(res_gvv, abs(log2FoldChange) > 0.5 & padj < 0.1)
sig_pvv <- subset(res_pvv, abs(log2FoldChange) > 0.5 & padj < 0.1)

head(sig_gvv)
head(sig_pvv)

## Venn figure
# 1) Define gene sets
set_gvv <- rownames(sig_gvv)
set_pvv <- rownames(sig_pvv)
sets <- list(GVV = set_gvv, PVV = set_pvv)

# 2) Quick stats
n_gvv <- length(set_gvv)
n_pvv <- length(set_pvv)
n_inter <- length(intersect(set_gvv, set_pvv))
cat("# sig_gvv:", n_gvv, "\n# sig_pvv:", n_pvv, "\n# overlap:", n_inter, "\n")

# Save the common/unique gene lists if you like
common_genes <- intersect(set_gvv, set_pvv)
unique_gvv   <- setdiff(set_gvv, set_pvv)
unique_pvv   <- setdiff(set_pvv, set_gvv)

# 3) Classic Venn
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
library(VennDiagram)
library(grid)

grid.newpage()
vp <- viewport(width = 0.8, height = 0.6)  # squish vertically → ellipses
pushViewport(vp)

# draw the Venn WITHOUT category labels
draw.pairwise.venn(
  area1 = length(set_gvv),
  area2 = length(set_pvv),
  cross.area = length(intersect(set_gvv, set_pvv)),
  category = c("", ""),                 # <- hide labels here
  fill = c("#B7E2FC", "#FFC6AF"),
  alpha = c(0.70, 0.70),
  lty = "blank",
  cex = 1.4
)

# ---- add custom labels ABOVE each ellipse ----
# tweak x positions slightly if needed (0..1 in the current viewport)
grid.text("DEGs_GDF15 vs Vehicle",
          x = 0.29, y = 0.90,            # above the blue ellipse
          gp = gpar(fontsize = 10, fontface = "bold", col = "#2B6E9E"))

grid.text("DEGs_CR vs Vehicle",
          x = 0.78, y = 0.90,            # above the orange ellipse
          gp = gpar(fontsize = 10, fontface = "bold", col = "#C55518"))

graph2png(x = NULL, file='VennPlot_commonGene', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 4.5, bg = "transparent")

graph2svg(x = NULL, file='VennPlot_commonGene', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 4.5, bg = "transparent")











########## Show gene expression related with these pathways using dotplot ######
## GDF15 treatment, 
# Macropahge activation: Ccl2, Ccl8, Icam1, Gbp2, Gbp4, Gbp6, Gbp9, Ass1, Syk, Fgr, Ubd, Cd274
# Antigen processing: Cd74, H2-Aa, H2-Ab1, H2-DMa, H2-DMb1, H2-Ob
# T cell activation: Cd3g, Lck, Zap70, Cd8a

## Pair-fed:
# Lipid Metabolism: Cpt1a, Acot1, Acot3, Acot4, Pnpla2, Lpl, Pla2g12a, Ppard


# Tpm=read.csv('1-TPM.csv', header = T)
# head(res_gvv)

res_gvv=read.csv('resGvV_results.csv', header = T)
res_pvv=read.csv('resPvV_results.csv', header = T)
row.names(res_gvv) =res_gvv$X
res_gvv=res_gvv[,-1]
head(res_gvv)

row.names(res_pvv) =res_pvv$X
res_pvv=res_pvv[,-1]
head(res_pvv)


# packages
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(forcats)

# ---------------------------
# 1) Define gene classes (updated)
# ---------------------------
gene_sets <- list(
  `Macrophage activation` = c("Ccl2","Ccl8","Icam1","Gbp2","Gbp4","Gbp6","Gbp9","Ass1","Syk","Fgr","Ubd","Cd274"),
  `Antigen processing`    = c("Cd74","H2-Aa","H2-Ab1","H2-DMa","H2-DMb1","H2-Ob"),
  `T cell activation`     = c("Cd3g","Lck","Zap70","Cd8a","Adcy7")
  
)
genes_all <- unique(unlist(gene_sets))

# ---------------------------
# 2) Helper to extract from a DESeq2 result (res_pvv/res_gvv)
# ---------------------------
prep_df <- function(res, dataset_label){
  as_tibble(res, rownames = "gene") %>%
    filter(gene %in% genes_all) %>%
    mutate(
      dataset = dataset_label,
      padj = ifelse(is.na(padj), 1, pmin(padj, 1)),
      neglog10padj = -log10(pmax(padj, 1e-300)),
      category = case_when(
        gene %in% gene_sets$`Macrophage activation` ~ "Macrophage activation",
        gene %in% gene_sets$`Antigen processing`    ~ "Antigen processing",
        gene %in% gene_sets$`T cell activation`     ~ "T cell activation",
        
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::select(gene, category, dataset, log2FoldChange, padj, neglog10padj)
}

df_left  <- prep_df(res_pvv, "CR")
df_right <- prep_df(res_gvv, "GDF15")
df <- bind_rows(df_left, df_right)

# Optional: report missing genes per dataset
missing_left  <- setdiff(genes_all, df_left$gene)
missing_right <- setdiff(genes_all, df_right$gene)
if (length(missing_left))  message("Missing in CR: ", paste(missing_left, collapse = ", "))
if (length(missing_right)) message("Missing in GDF15: ",   paste(missing_right, collapse = ", "))

# ---------------------------
# 3) Aesthetics with thresholds
#    - Size maps to –log10(padj), but clamp small values (<1) to the minimum size bucket
#    - Color maps to log2FC, but if |log2FC| < 0.5 => gray
# ---------------------------
df <- df %>%
  mutate(
    # keep dataset order: Pair-fed (left) | GDF15 (right)
    dataset = factor(dataset, levels = c("CR","GDF15")),
    # size mapping: enforce smallest size when -log10(padj) < 1
    size_mapped = ifelse(neglog10padj < 1, 1, neglog10padj),
    # color mapping: gray when |log2FC| < 0.5
    lfc_for_color = ifelse(abs(log2FoldChange) < 0.5, NA, log2FoldChange)
  )

# ---------------------------
# 4) Order categories/genes for display
# ---------------------------
cat_order  <- c("Macrophage activation","Antigen processing","T cell activation")
gene_order <- df %>%
  distinct(category, gene) %>%
  mutate(category = factor(category, levels = cat_order)) %>%
  arrange(category, gene) %>%
  pull(gene)

df <- df %>%
  mutate(
    category = factor(category, levels = cat_order),
    gene     = factor(gene, levels = rev(unique(gene_order))) # reverse so top-to-bottom within facets
  )

# ---------------------------
# 5) Plot
#    - size = size_mapped (range controls absolute pixel size)
#    - color = lfc_for_color (NA -> gray via na.value)
#    - gradient2 for smooth blue→white→red; midpoint at 0
# ---------------------------
p <- ggplot(df, aes(x = dataset, y = gene)) +
  geom_point(aes(size = size_mapped, color = lfc_for_color)) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  scale_size(
    range  = c(1.8, 9),
    name   = expression(-log[10]("adj P")),
    breaks = c(1, 2, 3, 5, 10),
    labels = c("<1 (min)", "2", "3", "5", "10")
  ) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    na.value = "grey70",  # used when |log2FC| < 0.5
    name = "log2FC"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Selected immune & lipid metabolism genes\n(Left: CR | Right: GDF15)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text.y     = element_text(face = "bold"),
    axis.text.x      = element_text(face = "bold"),
    panel.grid.major = element_line(size = 0.2, color = "grey90")
  )

p
graph2svg(x = p, file='Dotplot_GDF15Immune_PFMet', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 6, bg = "transparent")

print(p)
ggsave("Dotplot_GDF15Immune_PFMet.png", p, width = 3.5, height = 6, dpi = 300)
ggsave("Dotplot_GDF15Immune_PFMet_1.svg", p, width = 3.5, height = 6, dpi = 300)






#### Enrichment analysis using common and different DEGs between gvv and pvv

## ---- deps ----
if (!requireNamespace("clusterProfiler", quietly=TRUE)) install.packages("clusterProfiler")
if (!requireNamespace("org.Mm.eg.db", quietly=TRUE)) BiocManager::install("org.Mm.eg.db")
if (!requireNamespace("enrichplot", quietly=TRUE)) install.packages("enrichplot")

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

## ---- 1) define sets with your cutoff ----
cutoff <- function(x) {
  subset(x, !is.na(padj) & abs(log2FoldChange) > 0.5 & padj < 0.1)
}

sig_gvv <- cutoff(res_gvv)
sig_pvv <- cutoff(res_pvv)

genes_gvv <- rownames(sig_gvv)
genes_pvv <- rownames(sig_pvv)

## basic sets for Venn logic
common_all   <- intersect(genes_gvv, genes_pvv)
gvv_only     <- setdiff(genes_gvv, genes_pvv)
pvv_only     <- setdiff(genes_pvv, genes_gvv)

head(sig_gvv)
head(gvv_only)
head(pvv_only)



## (optional but useful) split by direction
dir <- function(df) {
  list(
    up   = rownames(df[df$log2FoldChange >  0.5, , drop=FALSE]),
    down = rownames(df[df$log2FoldChange <  0.5, , drop=FALSE])
  )
}
d_gvv <- dir(sig_gvv)
d_pvv <- dir(sig_pvv)

common_up_same    <- intersect(d_gvv$up,   d_pvv$up)
common_down_same  <- intersect(d_gvv$down, d_pvv$down)
common_opposite   <- c(intersect(d_gvv$up, d_pvv$down),
                       intersect(d_gvv$down, d_pvv$up))  # discordant direction

gvv_only_up       <- setdiff(d_gvv$up,   genes_pvv)
gvv_only_down     <- setdiff(d_gvv$down, genes_pvv)
pvv_only_up       <- setdiff(d_pvv$up,   genes_gvv)
pvv_only_down     <- setdiff(d_pvv$down, genes_gvv)

## ---- 2) gene universe (recommended for GO) ----
universe <- unique(c(rownames(res_gvv), rownames(res_pvv)))

## ---- 3) helper: run GO:BP enrichment on a gene vector of symbols ----
run_go <- function(genes_symbol, universe_symbol, q=0.05, p=0.05) {
  if (length(genes_symbol) < 5) return(NULL)  # too few for enrichment
  enrichGO(
    gene          = genes_symbol,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    universe      = universe_symbol,
    ont           = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff  = p,
    qvalueCutoff  = q,
    readable      = TRUE
  )
}

## ---- 4) run GO for sets you care about ----
ego_common_all      <- run_go(common_all,     universe)
head(ego_common_all)

# dotplot
dotplot(ego_common_all, showCategory=15)+ggtitle("GO enrichment analysis")
# dotplot1=dotplot(ego_all, showCategory=30)+ggtitle("GO enrichment analysis")
# write_fig(dotplot1, "GO enrichment analysis_gvp_2.tif", show = FALSE, devices = "tif", width = 12,
#           height = 17)
dev.off()


# run plot for gvv
ego_gvv_only        <- run_go(gvv_only,       universe)
head(ego_gvv_only)
write.csv(ego_gvv_only,'ego_gvv_only.csv')

# dotplot
dotplot(ego_gvv_only, showCategory=15,
        x = "GeneRatio", orderBy = "Count")+ggtitle("GO enrichment analysis") #defalt: orderBy = "p.adjust"

# dotplot1=dotplot(ego_all, showCategory=30)+ggtitle("GO enrichment analysis")
# write_fig(dotplot1, "GO enrichment analysis_gvp_2.tif", show = FALSE, devices = "tif", width = 12,
#           height = 17)
dev.off()

# run plot for pvv
ego_pvv_only        <- run_go(pvv_only,       universe)
head(ego_pvv_only)
write.csv(ego_pvv_only,'ego_pvv_only.csv')

# dotplot
dotplot(ego_pvv_only, showCategory=15,
        x = "GeneRatio", orderBy = "Count")+ggtitle("GO enrichment analysis") #defalt: orderBy = "p.adjust"




## (direction-aware, often more interpretable)
ego_common_up       <- run_go(common_up_same,   universe)
ego_common_down     <- run_go(common_down_same, universe)
ego_gvv_only_up     <- run_go(gvv_only_up,      universe)
ego_gvv_only_down   <- run_go(gvv_only_down,    universe)
ego_pvv_only_up     <- run_go(pvv_only_up,      universe)
ego_pvv_only_down   <- run_go(pvv_only_down,    universe)

## ---- 5) quick summaries / plots ----
if (!is.null(ego_common_all))   print(head(as.data.frame(ego_common_all), 10))
if (!is.null(ego_common_up))    print(head(as.data.frame(ego_common_up), 10))
if (!is.null(ego_common_down))  print(head(as.data.frame(ego_common_down), 10))

# Example dotplots (show top terms)
if (!is.null(ego_common_up))    print(dotplot(ego_common_up,  showCategory = 15, title = "Common UP (GVV ∩ PVV)"))
if (!is.null(ego_common_down))  print(dotplot(ego_common_down, showCategory = 15, title = "Common DOWN (GVV ∩ PVV)"))
if (!is.null(ego_gvv_only))     print(dotplot(ego_gvv_only,    showCategory = 15, title = "GVV-only"))
if (!is.null(ego_pvv_only))     print(dotplot(ego_pvv_only,    showCategory = 15, title = "PairFed-only"))

## ---- 6) (optional) save gene lists for records ----
writeLines(common_all,          "genes_common_all.txt")
writeLines(gvv_only,            "genes_gvv_only.txt")
writeLines(pvv_only,            "genes_pvv_only.txt")
writeLines(common_up_same,      "genes_common_up.txt")
writeLines(common_down_same,    "genes_common_down.txt")
writeLines(common_opposite,     "genes_common_opposite.txt")
writeLines(gvv_only_up,         "genes_gvv_only_up.txt")
writeLines(gvv_only_down,       "genes_gvv_only_down.txt")
writeLines(pvv_only_up,         "genes_pvv_only_up.txt")
writeLines(pvv_only_down,       "genes_pvv_only_down.txt")








####======= Immune cell Composition Deconvolution using ImmuCC ======
## Build mouse liver Reference from Tabula Muris

## ============ install ============
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs <- c("TabulaMurisSenisData","SingleCellExperiment","scuttle",
          "AnnotationDbi","org.Mm.eg.db","S4Vectors")
for (p in pkgs) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE)
invisible(lapply(pkgs, library, character.only=TRUE))

## ============ check tissue name（confirm 'Liver' writing） ============
dro_tissues <- listTabulaMurisSenisTissues("Droplet")
facs_tissues <- listTabulaMurisSenisTissues("FACS")
cat("Droplet tissues (head):", head(dro_tissues), "\n")
cat("FACS tissues   (head):", head(facs_tissues), "\n")

## ============ get Liver（Droplet & FACS；If a modality is missing, it will be skipped） ============
get_liver <- function(mode=c("Droplet","FACS"), processedCounts=FALSE){
  mode <- match.arg(mode)
  if (mode == "Droplet") {
    lst <- TabulaMurisSenisDroplet(tissues="Liver",
                                   processedCounts=processedCounts,
                                   reducedDims=FALSE)
  } else {
    lst <- TabulaMurisSenisFACS(tissues="Liver",
                                processedCounts=processedCounts,
                                reducedDims=FALSE)
  }
  if (!("Liver" %in% names(lst))) stop(mode, ": No return: Liver；firstly use listTabulaMurisSenisTissues(\"", mode, "\") to see details")
  sce <- lst[["Liver"]]
  if (!is(sce, "SingleCellExperiment")) sce <- as(sce, "SingleCellExperiment")
  sce
}

sce_d <- tryCatch(get_liver("Droplet", processedCounts=FALSE), error=function(e) NULL)
sce_f <- tryCatch(get_liver("FACS",    processedCounts=FALSE), error=function(e) NULL)

if (is.null(sce_d) && is.null(sce_f)) stop("Droplet 与 FACS 都未获取到 Liver。")

## ============ Metadata standardization（celltype / sample） ============
standardize_meta <- function(sce, prefix=NULL){
  cd <- as.data.frame(colData(sce))
  pick1 <- function(cands) { hit <- cands[cands %in% names(cd)]; if (length(hit)) cd[[hit[1]]] else NULL }
  celltype <- pick1(c("cell_ontology_class","free_annotation","celltype","celltype1","major_cell_type","label"))
  sampleid <- pick1(c("mouse.id","donor","subject","individual","mouse","sample","Sample","sample_id","batch"))
  if (is.null(celltype)) celltype <- rep("unknown", ncol(sce))
  if (is.null(sampleid)) sampleid <- paste0("pseudo_", seq_len(ncol(sce)))
  if (!is.null(prefix)) celltype <- paste(prefix, celltype, sep=":")
  colData(sce)$celltype <- as.character(celltype)
  colData(sce)$sample   <- as.character(sampleid)
  sce
}

if (!is.null(sce_d)) sce_d <- standardize_meta(sce_d, "Droplet")
if (!is.null(sce_f)) sce_f <- standardize_meta(sce_f, "FACS")

## ============ Convert gene names to mouse gene symbols and merge duplicated genes ============
emap_ens_to_symbol <- function(ens_vec){
  ens <- sub("\\.\\d+$","", ens_vec)
  AnnotationDbi::mapIds(org.Mm.eg.db, keys=ens, keytype="ENSEMBL", column="SYMBOL", multiVals="first")
}

to_symbol_and_collapse <- function(sce){
  rn  <- rownames(sce)
  sym <- ifelse(grepl("^ENSMUSG", rn, ignore.case=TRUE), emap_ens_to_symbol(rn), rn)
  # Fallback: keep the original gene name if no mapping is found
  sym[is.na(sym) | !nzchar(sym)] <- rn[is.na(sym) | !nzchar(sym)]
  keep <- !is.na(sym) & nzchar(sym)
  sce  <- sce[keep, ]
  sym  <- sym[keep]
  if (!"counts" %in% assayNames(sce)) {
    # Some objects contain only logcounts; an approximate inverse transformation is applied.
    if ("logcounts" %in% assayNames(sce)) {
      assay(sce, "counts") <- round(pmax(0, 2^assay(sce,"logcounts") - 1))
    } else stop("No counts/logcounts，cannot continue")
  }
  # merge duplicated symbol（sum）
  cts <- as.matrix(assay(sce,"counts"))
  agg <- rowsum(cts, group=sym, reorder=FALSE)
  sce2 <- SingleCellExperiment(list(counts=agg), colData=colData(sce))
  scuttle::logNormCounts(sce2)  # generate logcounts
}

if (!is.null(sce_d)) sce_d <- to_symbol_and_collapse(sce_d)
if (!is.null(sce_f)) sce_f <- to_symbol_and_collapse(sce_f)

## ============ Merge Droplet + FACS ============
## 1) Align genes by rows
common <- intersect(rownames(sce_d), rownames(sce_f))
sce_d2 <- sce_d[common, ]
sce_f2 <- sce_f[common, ]

## 2) Optional: add a prefix to avoid duplicate column names
colnames(sce_d2) <- paste0("Droplet_", colnames(sce_d2))
colnames(sce_f2) <- paste0("FACS_",    colnames(sce_f2))

## 3) Merge counts
cnt <- cbind(assay(sce_d2, "counts"), assay(sce_f2, "counts"))

## 4) Standardize the colData columns
cd_d <- as.data.frame(colData(sce_d2))
cd_f <- as.data.frame(colData(sce_f2))
all_cols <- union(names(cd_d), names(cd_f))

# Fill missing columns with NA
for (nm in setdiff(all_cols, names(cd_d))) cd_d[[nm]] <- NA
for (nm in setdiff(all_cols, names(cd_f))) cd_f[[nm]] <- NA

# Standardize column order; convert factors to characters to avoid inconsistent levels.
cd_d <- lapply(cd_d[all_cols], function(x) if (is.factor(x)) as.character(x) else x)
cd_f <- lapply(cd_f[all_cols], function(x) if (is.factor(x)) as.character(x) else x)

cd <- S4Vectors::DataFrame( rbind( as.data.frame(cd_d), as.data.frame(cd_f) ) )

## 5) Generate the reference and normalize it.
ref <- SingleCellExperiment(list(counts = cnt), colData = cd)
ref <- scuttle::logNormCounts(ref)

## 6) Self-check
stopifnot(all(c("celltype","sample") %in% colnames(colData(ref))))

## ============ save and check ============
stopifnot(all(c("celltype","sample") %in% colnames(colData(ref))))
print(head(sort(table(ref$celltype), decreasing=TRUE), 20))
cat("Unique samples:", length(unique(ref$sample)), "\n")
saveRDS(ref, "ref_MuSic.rds")
cat("Saved:", normalizePath("ref.rds"), "\n")


## ================== 6. (Optional) Retain only immune-related cells to reduce noise ==================
# You may also retain all liver cell types (e.g., hepatocytes, cholangiocytes, HSCs),
# as MuSiC can work with the complete dataset.
# The example below retains only immune-related cells (uncomment if needed):
keep_ct <- c("Kupffer cell","myeloid leukocyte","NK T cell", "macrophage","monocyte","B cell","plasma cell",
             "T cell","NK cell","neutrophil","dendritic cell", "hepatic stellate cell")
mask <- sapply(strsplit(ref$celltype, ":"), function(x) paste(x[-1], collapse=":"))  # 去掉前缀
ref   <- ref[, tolower(mask) %in% tolower(keep_ct)]
saveRDS(ref, "ref_MuSic_Immune.rds")



###==== Run MuSic ========
# Install MuSiC and common dependencies
BiocManager::install("TOAST")
library(TOAST)
devtools::install_github('xuranw/MuSiC')

library(MuSiC)
packageVersion("MuSiC")


# Load the prepared Tabula Muris reference
ref <- readRDS("ref_MuSic_Immune.rds")
table(ref$celltype)


# Load the bulk expression matrix (genes × samples; TPM/CPM recommended; **do not use log-transformed data**)
expr <- read.csv("1-TPM.csv", header = T)
head(expr)
## ===  Arrange expr ===
stopifnot("X" %in% colnames(expr))    
rn <- make.unique(as.character(expr$X))
expr$X <- NULL
# Set the first column as row names and convert to a numeric matrix
expr <- as.data.frame(expr, check.names = FALSE)
expr[] <- lapply(expr, function(x) as.numeric(as.character(x)))  # Prevent character/factor issues
expr <- as.matrix(expr)
rownames(expr) <- rn
# If duplicate gene names are present, merge them by summing their expression values
if (any(duplicated(rownames(expr)))) {
  expr <- rowsum(expr, group = rownames(expr), reorder = FALSE)
}
## Quick check whether the expression matrix is on a linear scale (not log-transformed)
summary(as.numeric(expr))
quantile(expr, c(.5,.9,.99), na.rm=TRUE)
head(expr)


# Ensure gene names are aligned
common <- intersect(rownames(expr), rownames(ref))
expr   <- as.matrix(expr[common, , drop = FALSE])
ref    <- ref[common, ]

# Deconvolution
res <- music_prop(bulk.mtx = expr,
                  sc.sce   = ref,
                  clusters = "celltype",
                  samples  = "sample",
                  verbose  = TRUE)

head(res$Est.prop.weighted)     
cell_fraction <- res$Est.prop.weighted
head(cell_fraction)


## Visualization
library(reshape2); library(ggplot2)
# Bar figure
df <- as.data.frame(t(cell_fraction))
df$celltype <- rownames(df)
m  <- melt(df, id.vars="celltype", variable.name="group", value.name="fraction")
head(m)
table(m$celltype)
ggplot(m, aes(group, fraction, fill=celltype)) +
  geom_bar(stat="identity") + theme_bw() + coord_flip()


## Clean group information
library(dplyr)
library(stringr)
# Extract common cell type names (remove the "Droplet:" and "FACS:" prefixes)
m <- m %>%
  mutate(celltype_clean = str_remove(celltype, "^(Droplet:|FACS:)"))

# Summarize fractions by group and celltype_clean
m_sum <- m %>%
  group_by(group, celltype_clean) %>%
  summarise(fraction = sum(fraction, na.rm = TRUE), .groups = "drop")


# View the results
print(m_sum)
m_sum$celltype = m_sum$celltype_clean
m_sum=m_sum[,-2]

tail(m_sum)
library(dplyr)
m_sum <- m_sum %>%
  mutate(group = sub("PF", "CR", group))

write.csv(m_sum, "CelltypeDeconv.csv")

# Bar figure

head(m_sum)
table(m_sum$celltype)
ggplot(m_sum, aes(group, fraction, fill=celltype)) +
  geom_bar(stat="identity") + theme_bw() + coord_flip()

graph2tif(x = NULL, file='1_CelltypeDeco_pct', font = "Arial", cairo = TRUE,   
          width = 6.5, height = 3, bg = "transparent")

## get better data frame
library(dplyr)  
m2 <- m_sum %>%
  mutate(condition = str_replace(group, "\\..*$", "")) %>%
  dplyr::select(condition, group, celltype, fraction)
head(m2)
write.csv(m2, "CelltypePencFigMaking.csv")



## Bar figure for celltype percentage.
library(dplyr)
library(ggplot2)
library(scales)
# Assume your data frame is named `df` and contains the following columns:
# `condition`, `group`, `celltype`, and `fraction`.
# If your data frame has a different name, first assign it using:
# df <- <your_data_frame>

# Retain only the three conditions of interest and set their display order
cond_order <- c("Veh", "GDF15", "CR")

plot_df <- m2 %>%
  filter(condition %in% cond_order) %>%
  mutate(condition = factor(condition, levels = cond_order)) %>%
  group_by(condition, celltype) %>%
  summarise(fraction = mean(fraction, na.rm = TRUE), .groups = "drop") %>%
  # Normalize within each condition so that each stacked bar sums to 1
  group_by(condition) %>%
  mutate(fraction = fraction / sum(fraction)) %>%
  ungroup()

# Figure）
fs <- 12  
p <- ggplot(plot_df, aes(x = condition, y = fraction, fill = celltype)) +
  geom_col(position = "fill", width = 0.7) +  
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_discrete(labels = c(Veh = "Vehicle", GDF15 = "GDF15", CR = "CR")) +
  labs(x = "group", y = "fraction", fill = "celltype") +
  theme_minimal(base_size = 12) +
  theme_minimal(base_size = fs) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
    
    
    text        = element_text(size = fs, face = "plain"),
    axis.title  = element_text(size = fs, face = "plain"),
    axis.text   = element_text(size = fs, face = "plain"),
    legend.title= element_text(size = fs, face = "plain"),
    legend.text = element_text(size = fs, face = "plain"),
    strip.text  = element_text(size = fs, face = "plain"),
    plot.title  = element_text(size = fs, face = "plain"),
    plot.subtitle = element_text(size = fs, face = "plain"),
    plot.caption  = element_text(size = fs, face = "plain")
  )

print(p)

graph2tif(x = NULL, file='1_CelltypeDeco_pct_1', font = "Arial", cairo = TRUE,   
          width = 6.5, height = 2.5, bg = "transparent")

graph2svg(x = NULL, file='1_CelltypeDeco_pct_1', font = "Arial", cairo = TRUE,   
          width = 6.5, height = 2.5, bg = "transparent")
# ggsave("Dotplot_GDF15Immune_PFMet_1.svg", p, width = 3.5, height = 6, dpi = 300)





####======= Immune cell Functional/Pathway Enrichment ======
library(GSVA)
library(GSEABase)

## load GMT file
library(clusterProfiler)
gene_sets = read.gmt("immune_signatures.gmt")
head(gene_sets)
gene_sets_list <- split(gene_sets$gene, gene_sets$term) #Transfer data.frame(term, gene) into named list


## load expression data: log2 TPM matrix
expr <- read.csv("1-TPM.csv", header = T)
head(expr)
## ===  Arrange expr ===
stopifnot("X" %in% colnames(expr))    
rn <- make.unique(as.character(expr$X))
expr$X <- NULL

expr <- as.data.frame(expr, check.names = FALSE)
expr[] <- lapply(expr, function(x) as.numeric(as.character(x)))  
expr <- as.matrix(expr)
rownames(expr) <- rn

if (any(duplicated(rownames(expr)))) {
  expr <- rowsum(expr, group = rownames(expr), reorder = FALSE)
}

summary(as.numeric(expr))
quantile(expr, c(.5,.9,.99), na.rm=TRUE)
head(expr)
expr_mat <- as.matrix(expr)
colnames(expr_mat) <- sub("^PF\\.", "CR.", colnames(expr_mat)) #把列名里的 PF.* 统一改成 CR.*，与前面分组一致
head(expr_mat)
write.csv(expr_mat, file = "expr_mat.csv", row.names = TRUE)


## GSVA
# Install GSVA
if (!requireNamespace("GSVA", quietly = TRUE)) {
  BiocManager::install("GSVA")
}
library(GSVA)
library(BiocParallel)

stopifnot(is.matrix(expr_mat))
mode(expr_mat) <- "numeric"  # 避免字符/因子
stopifnot(!is.null(rownames(expr_mat)))
expr_log <- log2(expr_mat + 1)

# Optional: remove gene sets with no overlap with the expression matrix
gene_sets_list <- lapply(gene_sets_list, function(x) intersect(x, rownames(expr_mat)))
gene_sets_list <- gene_sets_list[lengths(gene_sets_list) >= 2]
# Create the ssGSEA parameter object
param <- GSVA::ssgseaParam(
  exprData   = expr_log,
  geneSets   = gene_sets_list,
  minSize    = 2,
  maxSize    = 5000,
  normalize  = TRUE      
  
)


gsva_res <- GSVA::gsva(
  param,
  BPPARAM = BiocParallel::SnowParam(workers = 1),  
  verbose = FALSE
)

dim(gsva_res)
head(gsva_res[, 1:5])

# Save
mat <- as.matrix(gsva_res)
write.csv(mat, file = "gsva_scores_ssgsea_wide.csv", row.names = TRUE)



## $$$$$$$$$ Make figures $$$$$$$$$$$$$$$
## ========= 0) packages =========
pkgs <- c("dplyr","tidyr","stringr","forcats","ggplot2","rstatix","pheatmap")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

## ========= 1) # Prepare the data: convert to long format and assign groups =========
stopifnot(is.matrix(gsva_res))
scores_long <- as.data.frame(gsva_res) |>
  tibble::rownames_to_column("pathway") |>
  tidyr::pivot_longer(-pathway, names_to = "sample", values_to = "score") |>
  mutate(group = dplyr::case_when(
    str_starts(sample, "Veh\\.")   ~ "Veh",
    str_starts(sample, "CR\\.")    ~ "CR",
    str_starts(sample, "GDF15\\.") ~ "GDF15",
    TRUE ~ "Other"
  )) |>
  filter(group != "Other") |>
  mutate(group = factor(group, levels = c("Veh","CR","GDF15")))

## ========= # 2) Statistical analysis: Kruskal–Wallis test (overall difference among three groups) + Dunn's post hoc pairwise comparisons =========
kw <- scores_long |>
  group_by(pathway) |>
  rstatix::kruskal_test(score ~ group) |>
  rstatix::adjust_pvalue(method = "BH") |>
  mutate(sig = case_when(
    p.adj < 0.001 ~ "***",
    p.adj < 0.01  ~ "**",
    p.adj < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

## Dunn 
pairs <- scores_long |>
  group_by(pathway) |>
  rstatix::dunn_test(score ~ group, p.adjust.method = "BH") |>
  ungroup() |>
  mutate(comparison = paste(group2, "vs", group1),
         neglog10FDR = -log10(p.adj))

# Calculate the median for each group to determine the direction of the effect (positive or negative difference)
meds <- scores_long |>
  group_by(pathway, group) |>
  summarise(med = median(score, na.rm = TRUE), .groups="drop") |>
  tidyr::pivot_wider(names_from = group, values_from = med)

pairs <- pairs |>
  left_join(meds, by = "pathway") |>
  mutate(diff = case_when(
    comparison == "CR vs Veh"     ~ CR - Veh,
    comparison == "GDF15 vs Veh"  ~ GDF15 - Veh,
    comparison == "GDF15 vs CR"   ~ GDF15 - CR,
    TRUE ~ NA_real_
  ))

## ========= 3) Visualization A: Group Heatmap (Group Medians, Row-Scaled) =========
# Calculate the median of each pathway across the three groups,
# perform row-wise Z-score normalization, and plot the heatmap.
# Pathways are ordered by Kruskal–Wallis significance.
mat_group <- meds |>
  arrange(match(pathway, kw |> arrange(p.adj) |> pull(pathway)))  
mat <- as.matrix(mat_group[, c("Veh","CR","GDF15")])
rownames(mat) <- mat_group$pathway

# Row-wise scaling (facilitates comparison of relative enrichment)
mat_scaled <- t(scale(t(mat)))
rownames(mat_scaled) <- rownames(mat)

N <- 40
top_idx <- seq_len(min(N, nrow(mat_scaled)))
anno_row <- kw |> select(pathway, `KW_FDR` = p.adj, sig) |> as.data.frame()
rownames(anno_row) <- anno_row$pathway
anno_row <- anno_row[rownames(mat_scaled)[top_idx], c("KW_FDR","sig"), drop=FALSE]

## Heatmap
mat_plot <- mat_scaled
mat_plot[mat_plot < -1] <- -1
mat_plot[mat_plot >  1] <-  1

pal <- colorRampPalette(c("#CFE8FF",  "#FFCCCC"))(100)
brks <- seq(-1, 1, length.out = length(pal) + 1)

# figure
fs <- 9  
pheatmap::pheatmap(
  mat_plot[top_idx, ],
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color  = pal,
  breaks = brks,
  main   = "ssGSEA group medians (row-scaled)",
  fontsize = fs,           
  fontsize_row = fs,       
  fontsize_col = fs        
)

graph2tif(x = NULL, file='1_Immunefunction_heatmap', font = "Arial", cairo = TRUE,   
          width = 5, height = 2.5, bg = "transparent")
graph2svg(x = NULL, file='1_Immunefunction_heatmap', font = "Arial", cairo = TRUE,   
          width = 5, height = 2.5, bg = "transparent")


## ========= 4) Visualization B: Faceted Violin/Box Plots 
##            (Three-Group Distributions + KW FDR Labels) =========
# Select the top significantly altered pathways for faceted plotting
top_pathways <- kw |> arrange(p.adj) |> slice(1:12) |> pull(pathway)

ymax <- scores_long |>
  filter(pathway %in% top_pathways) |>
  group_by(pathway) |>
  summarise(y = max(score, na.rm = TRUE), .groups="drop")

lab <- kw |>
  filter(pathway %in% top_pathways) |>
  transmute(pathway, label = paste0("KW FDR = ", signif(p.adj, 3))) |>
  left_join(ymax, by = "pathway") |>
  mutate(x = 2)  

ggplot(scores_long |> filter(pathway %in% top_pathways),
       aes(group, score, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.9) +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  stat_summary(fun = median, geom = "point", size = 1.2, color = "white") +
  facet_wrap(~ pathway, scales = "free_y", ncol = 4) +
  geom_text(data = lab, aes(x = x, y = y * 1.05, label = label),
            inherit.aes = FALSE, size = 3.2) +
  labs(title = "ssGSEA scores by group (Violin + Box)",
       x = "Group", y = "ssGSEA score") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(size = 9))

## ========= 5) Visualization C: Pairwise Comparison Bubble Plot Summary
##            (Three Pairwise Comparisons Among the Three Groups) =========
pairs$neglog10FDR[!is.finite(pairs$neglog10FDR)] <- 0
library(ggplot2)
library(forcats)
fs <- 12  
p <- ggplot(pairs, aes(
  x = forcats::fct_inorder(comparison),
  y = forcats::fct_reorder(pathway, -abs(diff), .fun = max, .na_rm = TRUE),
  size = neglog10FDR, fill = diff
)) +
  geom_point(shape = 21, color = "black", alpha = 0.9) +
  scale_size_continuous(name = expression(-log[10]("FDR")), range = c(1.5, 8)) +
  scale_fill_gradient2(name = "Median diff", low = "#3B8EEA", mid = "white",
                       high = "#E64B35", midpoint = 0) +
  labs(x = "Pairwise comparison", y = "Pathway",
       title = "Pairwise Dunn tests (3 groups)") +
  theme_minimal(base_size = fs) +
  theme(
     
    text         = element_text(size = fs, face = "plain"),
    axis.title   = element_text(size = fs, face = "plain"),
    axis.text    = element_text(size = fs, face = "plain"),
    legend.title = element_text(size = fs, face = "plain"),
    legend.text  = element_text(size = fs, face = "plain"),
    plot.title   = element_text(size = fs, face = "plain", hjust = 0.5),
    axis.text.x  = element_text(angle = 30, hjust = 1),
    
    
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.9)
    
  )

print(p)

graph2tif(x = NULL, file='1_Immunefunction_dotplot', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.2, bg = "transparent")
graph2svg(x = NULL, file='1_Immunefunction_dotplot', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.2, bg = "transparent")
ggsave("1_Immunefunction_dotplot_1.svg", p, width = 6, height = 4.2, dpi = 300)


## ========= 6) Visualization D: Sample-Level PCA
##            (Overall Separation Among the Three Groups) =========
pc <- prcomp(t(gsva_res), scale. = TRUE)
pc_df <- data.frame(pc$x[, 1:2], sample = colnames(gsva_res)) |>
  mutate(group = case_when(
    str_starts(sample, "Veh\\.")   ~ "Veh",
    str_starts(sample, "CR\\.")    ~ "CR",
    str_starts(sample, "GDF15\\.") ~ "GDF15",
    TRUE ~ "Other"
  ),
  group = factor(group, levels = c("Veh","CR","GDF15")))

ggplot(pc_df, aes(PC1, PC2, shape = group, color = group, label = sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(show.legend = FALSE, max.overlaps = 20) +
  labs(title = "PCA of ssGSEA scores", x = "PC1", y = "PC2") +
  theme_minimal(base_size = 12)






