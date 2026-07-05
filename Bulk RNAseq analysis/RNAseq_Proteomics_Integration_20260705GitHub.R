########################### Liver RNAseq and serum proteomics integration analysis ##################

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


dir <- "C:/BaiduSyncdisk/R analysis/Gdf1520dLiver"
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
dds$condition <- relevel(dds$condition, ref = "Vehicle")

dds <- DESeq(dds)
resultsNames(dds)

# GDF15 vs Vehicle
rna_GDF15_vs_Veh <- results(dds, name = "condition_GDF15_vs_Vehicle")

# PairFed vs Vehicle
rna_CR_vs_Veh <- results(dds, name = "condition_PairFed_vs_Vehicle")

# GDF15 vs PairFed (use contrast)
rna_GDF15_vs_CR <- results(dds, contrast = c("condition", "GDF15", "PairFed"))




#########################################

## Proteomics data input and analysis

############################################

library(GEOquery)
library(dplyr)
library(pathview)
library("tximport")
library("readr")
library("tximportData")

library("DESeq2")
library("export")

library(ggplot2)
library(ggrepel)
library(dplyr)

# load data

## integration analysis
library(readr)
df <- read_csv("PlasmaProt_Mouse_GDF15CR_clean.csv")
head(df)

# Proteomics data pre-process
library(limma)
library(edgeR)
library(tidyverse)

# 1. Extract sample columns
sample_cols <- grep("^(Veh|GDF15|CR)\\.", colnames(df), value = TRUE)

prot_mat <- df %>%
  select(PG.Genes, Protein, all_of(sample_cols)) %>%
  mutate(Gene = PG.Genes) %>%
  filter(!is.na(Gene), Gene != "") %>%
  mutate(Gene = str_split(Gene, ";")) %>%
  unnest(Gene) %>%
  mutate(Gene = str_trim(Gene)) %>%
  filter(!is.na(Gene), Gene != "") %>%
  select(Gene, all_of(sample_cols))

# 2. Merge replicated gene
prot_gene <- prot_mat %>%
  group_by(Gene) %>%
  summarise(across(all_of(sample_cols), ~ median(.x, na.rm = TRUE))) %>%
  ungroup()

prot_matrix <- prot_gene %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# 3. Log2 transformation
prot_log2 <- log2(prot_matrix + 1)

# Filter missing values
# Retain proteins detected in at least 50% of samples within at least one group.
group <- case_when(
  grepl("^Veh", sample_cols) ~ "Veh",
  grepl("^GDF15", sample_cols) ~ "GDF15",
  grepl("^CR", sample_cols) ~ "CR"
)

keep <- apply(prot_log2, 1, function(x) {
  any(tapply(!is.na(x), group, mean) >= 0.5)
})

prot_filt <- prot_log2[keep, ]

# 5. Missing value imputation
set.seed(123)

impute_low <- function(x) {
  if (all(is.na(x))) return(x)
  min_val <- min(x, na.rm = TRUE)
  x[is.na(x)] <- rnorm(sum(is.na(x)), mean = min_val - 1.8, sd = 0.3)
  x
}

prot_imp <- t(apply(prot_filt, 1, impute_low))

# 6. Differential protein analysis
metadata <- data.frame(
  sample = sample_cols,
  group = factor(group, levels = c("Veh", "CR", "GDF15"))
)

design <- model.matrix(~ 0 + group, metadata)
colnames(design) <- levels(metadata$group)

fit <- lmFit(prot_imp, design)

contrast.matrix <- makeContrasts(
  GDF15_vs_Veh = GDF15 - Veh,
  CR_vs_Veh    = CR - Veh,
  GDF15_vs_CR  = GDF15 - CR,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

prot_GDF15_vs_Veh <- topTable(fit2, coef = "GDF15_vs_Veh", number = Inf)
prot_CR_vs_Veh    <- topTable(fit2, coef = "CR_vs_Veh", number = Inf)
prot_GDF15_vs_CR  <- topTable(fit2, coef = "GDF15_vs_CR", number = Inf)
head(prot_GDF15_vs_Veh)
head(prot_CR_vs_Veh)
head(prot_GDF15_vs_CR)

# 7. integration
rna_GDF15_vs_Veh
rna_CR_vs_Veh
rna_GDF15_vs_CR

library(tibble)
library(dplyr)
rna_CR_vs_Veh_df <- as.data.frame(rna_CR_vs_Veh) %>%
  rownames_to_column("Gene") %>%
  filter(!is.na(padj))

rna_GDF15_vs_Veh_df <- as.data.frame(rna_GDF15_vs_Veh) %>%
  rownames_to_column("Gene") %>%
  filter(!is.na(padj))

rna_GDF15_vs_CR_df <- as.data.frame(rna_GDF15_vs_CR) %>%
  rownames_to_column("Gene") %>%
  filter(!is.na(padj))


prot_GDF15_vs_Veh$Gene <- rownames(prot_GDF15_vs_Veh)
prot_CR_vs_Veh$Gene <- rownames(prot_CR_vs_Veh)
prot_GDF15_vs_CR$Gene <- rownames(prot_GDF15_vs_CR)

#################
# integration GDF15_vs_Veh data
######################
integration_GDF15 <- merge(
  rna_GDF15_vs_Veh_df,
  prot_GDF15_vs_Veh,
  by = "Gene",
  suffixes = c("_RNA", "_Protein")
)

colnames(integration_GDF15)

# Check the consistency of the effect direction
integration_GDF15 <- integration_GDF15 %>%
  dplyr::rename(
    log2FC_RNA = log2FoldChange,
    padj_RNA = padj,
    stat_RNA = stat,
    logFC_Protein = logFC,
    padj_Protein = adj.P.Val,
    t_Protein = t
  )

integration_GDF15 <- integration_GDF15 %>%
  mutate(
    direction = case_when(
      log2FC_RNA > 0 & logFC_Protein > 0 ~ "Both up",
      log2FC_RNA < 0 & logFC_Protein < 0 ~ "Both down",
      log2FC_RNA > 0 & logFC_Protein < 0 ~ "RNA up, protein down",
      log2FC_RNA < 0 & logFC_Protein > 0 ~ "RNA down, protein up",
      TRUE ~ "No clear pattern"
    )
  )
colnames(integration_GDF15)

# A simple fig
ggplot(integration_GDF15,
       aes(x = log2FC_RNA, y = logFC_Protein, color = direction)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()


# make RNA-protein concordance plot
integration_GDF15 <- integration_GDF15 %>%
  mutate(
    direction = case_when(
      log2FC_RNA > 0 & logFC_Protein > 0 ~ "Both up",
      log2FC_RNA < 0 & logFC_Protein < 0 ~ "Both down",
      log2FC_RNA > 0 & logFC_Protein < 0 ~ "RNA up, protein down",
      log2FC_RNA < 0 & logFC_Protein > 0 ~ "RNA down, protein up",
      TRUE ~ "No clear pattern"
    ),
    significance = case_when(
      padj_RNA < 0.1 & padj_Protein < 0.1 ~ "Both significant",
      padj_RNA < 0.1 ~ "RNA only",
      padj_Protein < 0.1 ~ "Protein only",
      TRUE ~ "Not significant"
    )
  )

library(ggplot2)
library(ggrepel)

top_genes <- integration_GDF15 %>%
  filter(padj_RNA < 0.1 & padj_Protein < 0.1) %>%
  arrange(padj_RNA + padj_Protein) %>%
  head(30)

ggplot(integration_GDF15,
       aes(x = log2FC_RNA, y = logFC_Protein)) +
  
  geom_point(color = "grey80", alpha = 0.5, size = 1.8) +
  
  geom_point(
    data = subset(integration_GDF15, significance == "Both significant"),
    aes(color = direction),
    size = 2.8,
    alpha = 0.9
  ) +
  
  geom_text_repel(
    data = top_genes,
    aes(label = Gene),
    size = 4   
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_classic() +
  
  theme(
    text = element_text(size = 14),          
    axis.title = element_text(size = 14),   
    axis.text = element_text(size = 14),    
    plot.title = element_text(size = 14, face = "bold"), 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  
  labs(
    x = "Liver RNA-seq log2FC",
    y = "Serum proteomics log2FC",
    title = "GDF15 vs Vehicle",
    color = "Direction"
  )

library(export)
graph2tif(x = NULL, file='Integrat_GDF15_Veh', font = "Arial", cairo = TRUE,   
          width = 7.5, height = 6, bg = "transparent")

graph2svg(x = NULL, file='Integrat_GDF15_Veh', font = "Arial", cairo = TRUE,   
          width = 7.5, height = 6, bg = "transparent")

sig_both_GDF15 <- integration_GDF15 %>%
  filter(padj_RNA < 0.1 & padj_Protein < 0.1) %>%
  arrange(direction, padj_RNA + padj_Protein)

View(sig_both_GDF15)



#################
# integration CR_vs_Veh data
######################

integration_GDF15 <- merge(
  rna_CR_vs_Veh_df,
  prot_CR_vs_Veh,
  by = "Gene",
  suffixes = c("_RNA", "_Protein")
)

colnames(integration_GDF15)

# Check the consistency of the effect direction
integration_GDF15 <- integration_GDF15 %>%
  dplyr::rename(
    log2FC_RNA = log2FoldChange,
    padj_RNA = padj,
    stat_RNA = stat,
    logFC_Protein = logFC,
    padj_Protein = adj.P.Val,
    t_Protein = t
  )

integration_GDF15 <- integration_GDF15 %>%
  mutate(
    direction = case_when(
      log2FC_RNA > 0 & logFC_Protein > 0 ~ "Both up",
      log2FC_RNA < 0 & logFC_Protein < 0 ~ "Both down",
      log2FC_RNA > 0 & logFC_Protein < 0 ~ "RNA up, protein down",
      log2FC_RNA < 0 & logFC_Protein > 0 ~ "RNA down, protein up",
      TRUE ~ "No clear pattern"
    )
  )
colnames(integration_GDF15)

# A simple fig
ggplot(integration_GDF15,
       aes(x = log2FC_RNA, y = logFC_Protein, color = direction)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic()


# 9. make RNA-protein concordance plot
integration_GDF15 <- integration_GDF15 %>%
  mutate(
    direction = case_when(
      log2FC_RNA > 0 & logFC_Protein > 0 ~ "Both up",
      log2FC_RNA < 0 & logFC_Protein < 0 ~ "Both down",
      log2FC_RNA > 0 & logFC_Protein < 0 ~ "RNA up, protein down",
      log2FC_RNA < 0 & logFC_Protein > 0 ~ "RNA down, protein up",
      TRUE ~ "No clear pattern"
    ),
    significance = case_when(
      padj_RNA < 0.1 & padj_Protein < 0.1 ~ "Both significant",
      padj_RNA < 0.1 ~ "RNA only",
      padj_Protein < 0.1 ~ "Protein only",
      TRUE ~ "Not significant"
    )
  )

library(ggplot2)
library(ggrepel)

top_genes <- integration_GDF15 %>%
  filter(padj_RNA < 0.1 & padj_Protein < 0.1) %>%
  arrange(padj_RNA + padj_Protein) %>%
  head(30)

ggplot(integration_GDF15,
       aes(x = log2FC_RNA, y = logFC_Protein)) +
  
  geom_point(color = "grey80", alpha = 0.5, size = 1.8) +
  
  geom_point(
    data = subset(integration_GDF15, significance == "Both significant"),
    aes(color = direction),
    size = 2.8,
    alpha = 0.9
  ) +
  
  geom_text_repel(
    data = top_genes,
    aes(label = Gene),
    size = 4   
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_classic() +
  
  theme(
    text = element_text(size = 14),          
    axis.title = element_text(size = 14),    
    axis.text = element_text(size = 14),     
    plot.title = element_text(size = 14, face = "bold"),  
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  
  labs(
    x = "Liver RNA-seq log2FC",
    y = "Serum proteomics log2FC",
    title = "CR vs Vehicle",
    color = "Direction"
  )

library(export)
graph2tif(x = NULL, file='Integrat_CR_Veh', font = "Arial", cairo = TRUE,   
          width = 7.5, height = 6, bg = "transparent")

graph2svg(x = NULL, file='Integrat_CR_Veh', font = "Arial", cairo = TRUE,   
          width = 7.5, height = 6, bg = "transparent")





# Recommended pathway-level integration
# Because the RNA-seq and proteomics datasets were generated from different mouse experiments,
# pathway-level integration is the most robust and reliable approach.
library(fgsea)

## For GDF15 and Vehicle

rna_GDF15_vs_Veh_df
prot_GDF15_vs_Veh

library(dplyr)
library(fgsea)
library(ggplot2)

## Pathway loading
library(GSEABase)
library(clusterProfiler)
geneset = read.gmt("MetabolicPathway3.gmt")
head(geneset)


library(dplyr)
library(tibble)
library(fgsea)
library(ggplot2)
library(ggrepel)

# 1. Convert the gene sets into the format required by fgsea
geneset_clean <- geneset %>%
  filter(!is.na(term), !is.na(gene)) %>%
  distinct(term, gene)

pathway_list <- split(geneset_clean$gene, geneset_clean$term)



# 2. RNA ranking
rna_rank <- rna_GDF15_vs_Veh_df %>%
  filter(!is.na(Gene), !is.na(stat)) %>%
  group_by(Gene) %>%
  summarise(stat = mean(stat, na.rm = TRUE), .groups = "drop") %>%
  deframe()

rna_rank <- sort(rna_rank, decreasing = TRUE)


# 3. Protein ranking
prot_rank <- prot_GDF15_vs_Veh %>%
  filter(!is.na(Gene), !is.na(t)) %>%
  group_by(Gene) %>%
  summarise(t = mean(t, na.rm = TRUE), .groups = "drop") %>%
  deframe()

prot_rank <- sort(prot_rank, decreasing = TRUE)


# 4. fgsea
fgsea_rna <- fgsea(
  pathways = pathway_list,
  stats = rna_rank,
  minSize = 5,
  maxSize = 500
)

fgsea_prot <- fgsea(
  pathways = pathway_list,
  stats = prot_rank,
  minSize = 5,
  maxSize = 500
)


# 5. pathway integration
pathway_integration <- merge(
  fgsea_rna[, c("pathway", "NES", "padj")],
  fgsea_prot[, c("pathway", "NES", "padj")],
  by = "pathway",
  suffixes = c("_RNA", "_Protein")
)

pathway_integration <- pathway_integration %>%
  mutate(
    direction = case_when(
      NES_RNA > 0 & NES_Protein > 0 ~ "Both up",
      NES_RNA < 0 & NES_Protein < 0 ~ "Both down",
      NES_RNA > 0 & NES_Protein < 0 ~ "RNA up, protein down",
      NES_RNA < 0 & NES_Protein > 0 ~ "RNA down, protein up"
    ),
    sig = case_when(
      padj_RNA < 0.15 & padj_Protein < 0.15 ~ "Both significant",
      padj_RNA < 0.15 ~ "RNA only",
      padj_Protein < 0.15 ~ "Protein only",
      TRUE ~ "Not significant"
    )
  )


# 6. Label only significantly enriched pathways
top_pathways <- pathway_integration %>%
  filter(sig != "Not significant") %>%
  arrange(padj_RNA + padj_Protein) %>%
  head(30)

# 7. Add quadrant annotations
pathway_integration <- pathway_integration %>%
  mutate(
    quadrant = case_when(
      NES_RNA > 0 & NES_Protein > 0 ~ "RNA↑ Protein↑",
      NES_RNA < 0 & NES_Protein < 0 ~ "RNA↓ Protein↓",
      NES_RNA > 0 & NES_Protein < 0 ~ "RNA↑ Protein↓",
      NES_RNA < 0 & NES_Protein > 0 ~ "RNA↓ Protein↑",
      TRUE ~ "Other"
    )
  )

quadrant_colors <- c(
  "RNA↑ Protein↑" = "#D55E00",
  "RNA↓ Protein↓" = "#0072B2",
  "RNA↑ Protein↓" = "#CC79A7",
  "RNA↓ Protein↑" = "#009E73"
)

# 8. Plot: color only significant pathways
ggplot(pathway_integration, aes(x = NES_RNA, y = NES_Protein)) +
  
  geom_point(
    data = subset(pathway_integration, sig == "Not significant"),
    color = "grey80",
    alpha = 0.5,
    size = 2
  ) +
  
  geom_point(
    data = subset(pathway_integration, sig != "Not significant" & quadrant != "Other"),
    aes(color = quadrant),
    size = 3,
    alpha = 0.9
  ) +
  
  geom_text_repel(
    data = top_pathways,
    aes(label = pathway),
    size = 4,
    max.overlaps = 30
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_classic(base_size = 14) +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  
  scale_color_manual(values = quadrant_colors) +
  
  labs(
    x = "Liver RNA-seq pathway NES",
    y = "Serum proteomics pathway NES",
    title = "GDF15 vs Vehicle pathway-level integration",
    color = "Quadrant"
  )



library(export)
graph2tif(x = NULL, file='Integrat_GDF15_Veh_pathway', font = "Arial", cairo = TRUE,   
          width = 7, height = 5, bg = "transparent")

graph2svg(x = NULL, file='Integrat_GDF15_Veh_pathway', font = "Arial", cairo = TRUE,   
          width = 7, height = 5, bg = "transparent")







## For CR and Vehicle

rna_CR_vs_Veh_df
prot_CR_vs_Veh

library(dplyr)
library(fgsea)
library(ggplot2)

## Pathway loading
library(GSEABase)
library(clusterProfiler)
geneset = read.gmt("MetabolicPathway3.gmt")
head(geneset)


library(dplyr)
library(tibble)
library(fgsea)
library(ggplot2)
library(ggrepel)

# 1. Convert gene sets to the list format required by fgsea
geneset_clean <- geneset %>%
  filter(!is.na(term), !is.na(gene)) %>%
  distinct(term, gene)

pathway_list <- split(geneset_clean$gene, geneset_clean$term)



# 2. RNA ranking
rna_rank <- rna_CR_vs_Veh_df %>%
  filter(!is.na(Gene), !is.na(stat)) %>%
  group_by(Gene) %>%
  summarise(stat = mean(stat, na.rm = TRUE), .groups = "drop") %>%
  deframe()

rna_rank <- sort(rna_rank, decreasing = TRUE)


# 3. Protein ranking
prot_rank <- prot_CR_vs_Veh %>%
  filter(!is.na(Gene), !is.na(t)) %>%
  group_by(Gene) %>%
  summarise(t = mean(t, na.rm = TRUE), .groups = "drop") %>%
  deframe()

prot_rank <- sort(prot_rank, decreasing = TRUE)


# 4. fgsea
fgsea_rna <- fgsea(
  pathways = pathway_list,
  stats = rna_rank,
  minSize = 5,
  maxSize = 500
)

fgsea_prot <- fgsea(
  pathways = pathway_list,
  stats = prot_rank,
  minSize = 5,
  maxSize = 500
)


# 5. pathway integration
pathway_integration <- merge(
  fgsea_rna[, c("pathway", "NES", "padj")],
  fgsea_prot[, c("pathway", "NES", "padj")],
  by = "pathway",
  suffixes = c("_RNA", "_Protein")
)

pathway_integration <- pathway_integration %>%
  mutate(
    direction = case_when(
      NES_RNA > 0 & NES_Protein > 0 ~ "Both up",
      NES_RNA < 0 & NES_Protein < 0 ~ "Both down",
      NES_RNA > 0 & NES_Protein < 0 ~ "RNA up, protein down",
      NES_RNA < 0 & NES_Protein > 0 ~ "RNA down, protein up"
    ),
    sig = case_when(
      padj_RNA < 0.15 & padj_Protein < 0.15 ~ "Both significant",
      padj_RNA < 0.15 ~ "RNA only",
      padj_Protein < 0.15 ~ "Protein only",
      TRUE ~ "Not significant"
    )
  )


# 6. mark pathway
top_pathways <- pathway_integration %>%
  filter(sig != "Not significant") %>%
  arrange(padj_RNA + padj_Protein) %>%
  head(30)

# 7. Add quadrant annotations
pathway_integration <- pathway_integration %>%
  mutate(
    quadrant = case_when(
      NES_RNA > 0 & NES_Protein > 0 ~ "RNA↑ Protein↑",
      NES_RNA < 0 & NES_Protein < 0 ~ "RNA↓ Protein↓",
      NES_RNA > 0 & NES_Protein < 0 ~ "RNA↑ Protein↓",
      NES_RNA < 0 & NES_Protein > 0 ~ "RNA↓ Protein↑",
      TRUE ~ "Other"
    )
  )

quadrant_colors <- c(
  "RNA↑ Protein↑" = "#D55E00",
  "RNA↓ Protein↓" = "#0072B2",
  "RNA↑ Protein↓" = "#CC79A7",
  "RNA↓ Protein↑" = "#009E73"
)

# 8. Plot: color only significant pathways
ggplot(pathway_integration, aes(x = NES_RNA, y = NES_Protein)) +
  
  geom_point(
    data = subset(pathway_integration, sig == "Not significant"),
    color = "grey80",
    alpha = 0.5,
    size = 2
  ) +
  
  geom_point(
    data = subset(pathway_integration, sig != "Not significant" & quadrant != "Other"),
    aes(color = quadrant),
    size = 3,
    alpha = 0.9
  ) +
  
  geom_text_repel(
    data = top_pathways,
    aes(label = pathway),
    size = 4,
    max.overlaps = 30
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme_classic(base_size = 14) +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  
  scale_color_manual(values = quadrant_colors) +
  
  labs(
    x = "Liver RNA-seq pathway NES",
    y = "Serum proteomics pathway NES",
    title = "CR vs Vehicle pathway-level integration",
    color = "Quadrant"
  )


library(export)
graph2tif(x = NULL, file='Integrat_CR_Veh_pathway', font = "Arial", cairo = TRUE,   
          width = 7, height = 5, bg = "transparent")

graph2svg(x = NULL, file='Integrat_CR_Veh_pathway', font = "Arial", cairo = TRUE,   
          width = 7, height = 5, bg = "transparent")


