########################### Serum Proteomics data analysis for GDF15, pair-feeding and vehicle

rm(list = ls())
options(stringsAsFactor = F)


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



######=========== Run again using raw data ==============######
# Load libraries
library(ggplot2)
library(ggrepel)
library(dplyr)

# Load datasets
Prot <- read.csv("PlasmaProt_Mouse_GDF15CR.csv", header=TRUE)
head(Prot)

SamInfo <- read.csv("PlasmaProt_Mouse_GDF15CRSamInfo.csv", header=TRUE)
head(SamInfo)


# Clean dataframe
df = Prot
clean_pg <- function(v) {
  sapply(v, function(x) {
    if (is.na(x) || trimws(x) == "" || tolower(x) == "nan") return(NA_character_)
    toks <- unlist(strsplit(x, ";", fixed = TRUE))
    toks <- trimws(toks)
    toks <- toks[toks != "" & !is.na(toks)]
    toks <- toupper(toks)            
    toks <- toks[!duplicated(toks)]  
    if (length(toks) == 0) return(NA_character_)
    paste(toks, collapse = ";")
  }, USE.NAMES = FALSE)
}

df$Protein <- clean_pg(df$PG.Genes)
head(df)
write.csv(df, "PlasmaProt_Mouse_GDF15CR_clean.csv", row.names = FALSE, na = "")
Prot1=df
head(Prot1)




### PCA analysis
# apply a log2 transformation,
# perform left-censored (MNAR) imputation for missing values,
# and carry out median normalization across samples;
# then plot a PC1–PC2 scatter with a 68% confidence ellipse and
# a variance-explained plot for the top 10 principal components.

## ===== 0) Packages =====
pkgs <- c("dplyr","ggplot2","ggrepel")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

## ===== 1) Extract the expression matrix and group labels (Protein as row names)=====
prot <- Prot1
stopifnot("Protein" %in% names(prot))

# Sample column：Veh./CR./GDF15.
sample_cols <- grep("^(Veh|CR|GDF15)\\.", names(prot), value = TRUE)

# Use Protein as row names; make duplicate row names unique if necessary
rn <- make.unique(prot$Protein)
X  <- as.matrix(prot[, sample_cols])
mode(X) <- "numeric"
X[is.nan(X)] <- NA_real_
rownames(X) <- rn

samples <- colnames(X)
group   <- factor(sub("\\..*$", "", samples), levels = c("Veh","CR","GDF15"))

## ## ===== 2) Filter by detection rate (quantified in ≥3 samples in at least two groups) =====
min_n <- 3
by_grp_ok <- sapply(levels(group), function(g){
  rowSums(!is.na(X[, group==g, drop=FALSE])) >= min_n
})
keep <- rowSums(by_grp_ok) >= 2
X <- X[keep, , drop=FALSE]

## ===== 3) Log2 Transformation + Left-Censored Imputation + Median Normalization Across Samples =====
pos_vals <- X[X > 0 & is.finite(X)]
pseudo   <- as.numeric(quantile(pos_vals, 0.01, na.rm = TRUE))/2
L <- log2(X + pseudo)

# Save the pre-imputation matrix for quality control (QC)
L_preimp <- L

# Left-censored imputation (1st percentile of each column − 1)
for (j in seq_len(ncol(L))) {
  colj  <- L[, j]
  q1    <- suppressWarnings(quantile(colj, 0.01, na.rm = TRUE))
  fillv <- if (is.finite(q1)) q1 - 1 else min(colj, na.rm = TRUE) - 1
  colj[!is.finite(colj) | is.na(colj)] <- fillv
  L[, j] <- colj
}

Lnorm <- sweep(L, 2, apply(L, 2, median, na.rm = TRUE), FUN = "-")

nzv <- apply(Lnorm, 1, sd) > 0
Lnorm <- Lnorm[nzv, , drop=FALSE]

## ===== 3.5) QC Metrics (Using the Pre-Imputation Matrix: L_preimp) =====
qc_tab <- data.frame(
  Sample       = samples,
  Group        = group,
  missing_frac = colMeans(!is.finite(L_preimp) | is.na(L_preimp)),
  median_int   = apply(L_preimp, 2, median, na.rm = TRUE),
  prot_count   = colSums(is.finite(L_preimp) & !is.na(L_preimp))
)

## ===== 4) PCA =====
pc <- prcomp(t(Lnorm), center = TRUE, scale. = TRUE)
var_expl <- (pc$sdev^2) / sum(pc$sdev^2)

pc_df <- data.frame(pc$x[, 1:2], Sample = samples, Group = group,
                    PC1 = pc$x[,1], PC2 = pc$x[,2])

## ===== 5) Figure：PC1–PC2 + ellipses（68%）=====
fs <- 12  
p_pca <- ggplot(pc_df, aes(PC1, PC2, color = Group, shape = Group, label = Sample)) +
  stat_ellipse(type = "norm", level = 0.68, linewidth = 0.8, show.legend = FALSE) +
  geom_point(size = 3, stroke = 0.7) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 30) +
  labs(
    title = "PCA of Proteomics (log2, imputed, median-normalized)",
    x = sprintf("PC1 (%.1f%%)", 100*var_expl[1]),
    y = sprintf("PC2 (%.1f%%)", 100*var_expl[2])
  ) +
  theme_minimal(base_size = fs) +
  theme(
    text = element_text(size = fs, face = "plain"),
    axis.title = element_text(size = fs, face = "plain"),
    axis.text  = element_text(size = fs, face = "plain"),
    legend.title = element_text(size = fs, face = "plain"),
    legend.text  = element_text(size = fs, face = "plain"),
    plot.title   = element_text(size = fs, face = "plain", hjust = 0.5),
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.9)
  )
print(p_pca)
graph2tif(x = NULL, file='1_PCA analysis', font = "Arial", cairo = TRUE,   
          width = 6.5, height = 4, bg = "transparent")
graph2svg(x = NULL, file='1_PCA analysis', font = "Arial", cairo = TRUE,   
          width = 6.5, height = 4, bg = "transparent")


## ===== 6) Scree plot =====
k <- min(10, length(var_expl))
scree_df <- data.frame(PC = factor(seq_len(k)), Var = var_expl[1:k])
p_scree <- ggplot(scree_df, aes(PC, Var)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", 100*Var)), vjust = -0.3, size = 3) +
  labs(title = "Explained variance", y = "Proportion of variance", x = NULL) +
  theme_minimal(base_size = fs) +
  theme(
    text = element_text(size = fs, face = "plain"),
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.9)
  )
print(p_scree)






## Different expression analysis
# Clean → log2 transform → left-censored (MNAR) imputation for missing values →
# between-sample normalization → limma linear modeling → 
# pairwise comparisons among the three groups, with example volcano plots and heatmaps provided.
## ===== 0) Packages =====
pkgs <- c("dplyr","stringr","limma","ggplot2","ggrepel","pheatmap","tibble")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

## ===== 1) Extract the expression matrix and group labels (Protein as row names)=====
prot <- Prot1
sample_cols <- grep("^(Veh|CR|GDF15)\\.", names(prot), value = TRUE)

name_raw <- as.character(prot$Protein)
name_raw <- trimws(name_raw)

# Use Protein as row names; automatically make duplicates unique
# use PG.Genes as the first gene
mask <- is.na(name_raw) | name_raw == ""
name_raw[mask] <- trimws(sub(";.*$", "", prot$`PG.Genes`[mask]))

rn <- make.unique(name_raw)
stopifnot(!any(is.na(rn) | rn == ""))

## make new matrix and scaling 
X  <- as.matrix(prot[, sample_cols])
mode(X) <- "numeric"
X[is.nan(X)] <- NA_real_
rownames(X) <- rn

samples <- colnames(X)
group   <- factor(sub("\\..*$", "", samples), levels = c("Veh","CR","GDF15"))

# Annonation
annot <- prot %>%
  transmute(Protein = rn,
            PG.Genes = `PG.Genes`,
            PG.ProteinAccessions = `PG.ProteinAccessions`) %>%
  tibble::remove_rownames()

## ===== 2) Filtering (Robust Detection) =====
# Retain proteins quantified in at least 3 samples per group in at least 2 groups.
min_n <- 3
by_grp_ok <- sapply(levels(group), function(g){
  rowSums(!is.na(X[, group==g, drop=FALSE])) >= min_n
})
keep <- rowSums(by_grp_ok) >= 2
X <- X[keep, , drop=FALSE]
annot <- annot[keep, , drop=FALSE]

## ===== 3) Log2 Transformation + Left-Censored Imputation + Between-Sample Median Normalization =====
# log2(x + pseudo)
# Use a global low-abundance value (half of the 1st percentile) as the pseudocount
pos_vals <- X[X > 0 & is.finite(X)]
pseudo   <- as.numeric(quantile(pos_vals, 0.01, na.rm = TRUE))/2
L <- log2(X + pseudo)

# Left-censored imputation: impute missing values with the 1st percentile of each column − 1
# (conservative low-abundance imputation)
for (j in seq_len(ncol(L))) {
  colj  <- L[, j]
  q1    <- suppressWarnings(quantile(colj, 0.01, na.rm = TRUE))
  fillv <- if (is.finite(q1)) q1 - 1 else min(colj, na.rm = TRUE) - 1
  colj[!is.finite(colj) | is.na(colj)] <- fillv
  L[, j] <- colj
}

Lnorm <- sweep(L, 2, apply(L, 2, median, na.rm = TRUE), FUN = "-")

nzv <- apply(Lnorm, 1, sd) > 0
Lnorm <- Lnorm[nzv, , drop=FALSE]
annot <- annot[nzv, , drop=FALSE]
write.csv(Lnorm, "Prot_clean.csv", row.names = FALSE)


## ===== 4) limma: Three Pairwise Comparisons =====
design <- model.matrix(~ 0 + group); colnames(design) <- levels(group)
fit  <- limma::lmFit(Lnorm, design)
cont <- limma::makeContrasts(
  GDF15_vs_Veh = GDF15 - Veh,
  CR_vs_Veh    = CR    - Veh,
  GDF15_vs_CR  = GDF15 - CR,
  levels = design
)
fit2 <- limma::contrasts.fit(fit, cont)
fit2 <- limma::eBayes(fit2, robust = TRUE, trend = TRUE)

# Function to generate result tables (merge annotations and sort by P value)
get_tab <- function(coef_name){
  tt <- limma::topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")
  tt <- tt %>%
    tibble::rownames_to_column("Protein") %>%
    left_join(annot, by = "Protein") %>%
    relocate(Protein, PG.Genes, PG.ProteinAccessions, .before = 1)
  tt$contrast <- coef_name
  tt$abs_logFC <- abs(tt$logFC)
  tt
}
res_GDF15_Veh <- get_tab("GDF15_vs_Veh")
res_CR_Veh    <- get_tab("CR_vs_Veh")
res_GDF15_CR  <- get_tab("GDF15_vs_CR")
res_all <- bind_rows(res_GDF15_Veh, res_CR_Veh, res_GDF15_CR)

## ===== 5) Statistical Significance Filtering and Export =====
fdr_cut <- 0.1; lfc_cut <- log2(1.5)  # 1.5x 阈值，可按需调整
sig_all <- res_all %>% filter(adj.P.Val < fdr_cut, abs_logFC >= lfc_cut)

# export
write.csv(res_GDF15_Veh, "DE_GDF15_vs_Veh_all.csv", row.names = FALSE)
write.csv(res_CR_Veh,    "DE_CR_vs_Veh_all.csv",    row.names = FALSE)
write.csv(res_GDF15_CR,  "DE_GDF15_vs_CR_all.csv",  row.names = FALSE)
write.csv(sig_all,       "DE_all_contrasts_sig.csv", row.names = FALSE)

## ===== 6) Visualization: Volcano Plot & Heatmap =====
## Thresholds
library(ggplot2)
library(ggrepel)

volcano_plot <- function(tab, title, fdr_cut = 0.1, lfc_cut = log2(1.5), fs = 12, label_shrink = 0.85){
  
  tab <- tab %>%
    dplyr::mutate(
      cat = dplyr::case_when(
        adj.P.Val < fdr_cut & logFC >  lfc_cut  ~ "Up",
        adj.P.Val < fdr_cut & logFC < -lfc_cut  ~ "Down",
        TRUE ~ "NS"
      ),
      sig = cat != "NS"
    )
  
  up   <- sum(tab$cat == "Up",   na.rm = TRUE)
  down <- sum(tab$cat == "Down", na.rm = TRUE)
  total <- up + down
  
  label_size <- (fs/3.3) * label_shrink
  
  ## === Top20 for annotation ===
  tab_label <- tab %>%
    filter(sig) %>%
    arrange(desc(abs(logFC))) %>%              
    head(20)
  
  g <- ggplot(tab, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = cat), alpha = 0.9, size = 1.8) +
    scale_color_manual(values = c(Up = "#E64B35", Down = "#3B8EEA", NS = "grey75")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2, color = "grey50") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = 2, color = "grey50") +
    
    ggrepel::geom_text_repel(
      data = tab_label,
      aes(label = Protein, color = cat),
      size = label_size, max.overlaps = Inf,
      box.padding = 0.35, point.padding = 0.25, segment.size = 0.3,
      show.legend = FALSE
    ) +
    labs(
      title = title,
      subtitle = sprintf("Up = %d | Down = %d | Total = %d  (FDR < %.2g; |log2FC| ≥ %.2f)",
                         up, down, total, fdr_cut, lfc_cut),
      x = "Log2FoldChange", y = expression(-log[10]("FDR"))
    ) +
    theme_minimal(base_size = fs) +
    theme(
      text         = element_text(size = fs, face = "plain"),
      axis.title   = element_text(size = fs, face = "plain"),
      axis.text    = element_text(size = fs, face = "plain"),
      legend.title = element_text(size = fs, face = "plain"),
      legend.text  = element_text(size = fs, face = "plain"),
      plot.title   = element_text(size = fs, face = "plain", hjust = 0.5),
      plot.subtitle= element_text(size = fs, face = "plain"),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.9)
    ) +
    guides(color = "none")
  
  
  xmax <- max(tab$logFC, na.rm = TRUE); xmin <- min(tab$logFC, na.rm = TRUE)
  ymax <- max(-log10(tab$adj.P.Val), na.rm = TRUE)
  g + annotate("label",
               x = xmax - 0.02*(xmax - xmin), y = ymax,
               hjust = 1, vjust = 1,
               label = sprintf("Up: %d\nDown: %d\nTotal: %d", up, down, total),
               size = fs/3.3, label.size = 0.25, fill = "white")
}

p1 <- volcano_plot(res_GDF15_Veh, "GDF15 vs Veh")
p2 <- volcano_plot(res_CR_Veh,    "CR vs Veh")
p3 <- volcano_plot(res_GDF15_CR,  "GDF15 vs CR")
print(p1)
print(p2)
print(p3)

# save
graph2tif(x = p1, file='2_volcano_GDF15 vs Veh_t20', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")

graph2tif(x = p2, file='2_volcano_CR vs Veh_t20', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")

graph2tif(x = p3, file='2_volcano_GDF15 vs CR_t20', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")

#svg file
graph2svg(x = p1, file='2_volcano_GDF15 vs Veh_t20', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")

graph2svg(x = p2, file='2_volcano_CR vs Veh_t20', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")

graph2svg(x = p3, file='2_volcano_GDF15 vs CR_t20', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")


# Heatmap: take the union of significant proteins from the three comparisons
# and display row-wise Z-scores
top_prots <- res_all %>% arrange(adj.P.Val) %>% pull(Protein) %>% unique()
top_prots <- head(top_prots, 40)  
H <- Lnorm[intersect(rownames(Lnorm), top_prots), , drop=FALSE]
ann_col <- data.frame(Group = group); rownames(ann_col) <- colnames(H)

pheatmap::pheatmap(H, scale = "row", annotation_col = ann_col,
                   show_rownames = TRUE, show_colnames = FALSE,
                   main = "Differential proteins (row Z-score)")

graph2tif(x = NULL, file='2_Heatmap', font = "Arial", cairo = TRUE,   
          width = 7, height = 7, bg = "transparent")




###=========== Plot Venn ================
## ===== Parameters and Gene Sets =====
library(dplyr)
fdr_cut <- 0.1
lfc_cut <- log2(1.5)

if (!"Protein" %in% names(res_GDF15_Veh)) res_GDF15_Veh$Protein <- rownames(res_GDF15_Veh)
if (!"Protein" %in% names(res_CR_Veh))    res_CR_Veh$Protein    <- rownames(res_CR_Veh)

sig_gdf15 <- subset(res_GDF15_Veh, adj.P.Val < fdr_cut & abs(logFC) >= lfc_cut)
sig_cr    <- subset(res_CR_Veh,    adj.P.Val < fdr_cut & abs(logFC) >= lfc_cut)

# Optional: if the Protein column contains multiple names separated by semicolons,
# use the following two lines to split them into individual entries.
# split_prot <- function(s) unique(trimws(unlist(strsplit(s, ";"))))
# set_gdf15 <- unique(unlist(lapply(sig_gdf15$Protein, split_prot))); set_cr <- unique(unlist(lapply(sig_cr$Protein, split_prot)))

set_gdf15 <- unique(sig_gdf15$Protein)
set_cr    <- unique(sig_cr$Protein)

# Statistics
n_gdf15 <- length(set_gdf15)
n_cr    <- length(set_cr)
n_inter <- length(intersect(set_gdf15, set_cr))
cat("# GDF15 vs Veh:", n_gdf15, "\n# CR vs Veh:", n_cr, "\n# overlap:", n_inter, "\n")

# List object (can be exported using write.csv if needed)
common_proteins <- intersect(set_gdf15, set_cr)
unique_gdf15    <- setdiff(set_gdf15, set_cr)
unique_cr       <- setdiff(set_cr, set_gdf15)

## ===== Plot Venn =====
library(VennDiagram)
library(grid)

if (names(dev.cur()) %in% c("png","pdf","jpeg","tiff","bmp","cairo_pdf")) dev.off()

grid.newpage()
vp <- viewport(width = 0.80, height = 0.60)  
pushViewport(vp)

g <- venn.diagram(
  x = list(`GDF15 vs Veh` = set_gdf15, `CR vs Veh` = set_cr),
  filename = NULL,
  fill = c("#9FCBFF", "#FF8F8F"),
  alpha = 0.70, lty = "blank",
  cex = 1.4, cat.cex = 0
)

grid.draw(g) 

fs <- 12
grid.text(sprintf("GDF15 vs Vehicle (n = %d)", length(unique(set_gdf15))),
          x = 0.35, y = 0.90, gp = gpar(fontsize = fs, col = "#2B6E9E"))
grid.text(sprintf("CR vs Vehicle (n = %d)", length(unique(set_cr))),
          x = 0.78, y = 0.80, gp = gpar(fontsize = fs, col = "#FF6F6F"))
grid.text(sprintf("Cutoffs: FDR < %.2g, |log2FC| ≥ %.2f", fdr_cut, lfc_cut),
          x = 0.50, y = 0.06, gp = gpar(fontsize = fs, col = "grey25"))

popViewport()

graph2png(x = NULL, file='4_VennPlot_commonProtein', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 4.5, bg = "transparent")

graph2svg(x = NULL, file='4_VennPlot_commonProtein', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 4.5, bg = "transparent")






###=========== Enrichment analysis ============###
## ==== 0) Packages ====
if (!requireNamespace("clusterProfiler", quietly=TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Mm.eg.db", quietly=TRUE))   BiocManager::install("org.Mm.eg.db")
if (!requireNamespace("enrichplot", quietly=TRUE))      BiocManager::install("enrichplot")
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)
library(stringr)
library(ggplot2)

## ===== 1) Common Parameters and Utility Functions =====
fdr_cut <- 0.1
lfc_cut <- log2(1.2)

# Extract the representative gene name from each row
# (prefer the first entry in PG.Genes; otherwise use Protein)
get_symbol <- function(df){
  has_pg <- "PG.Genes" %in% names(df)
  sym <- if (has_pg) sub(";.*$", "", df$PG.Genes) else df$Protein
  sym <- trimws(sym)
  sym[sym=="" | is.na(sym)] <- NA
  sym
}

# Standardize: extract significant proteins and split them into up- and downregulated sets
# (returns SYMBOL vectors)
pick_sig_sets <- function(res_tab){
  stopifnot(all(c("adj.P.Val","logFC") %in% names(res_tab)))
  res_tab$SYMBOL <- get_symbol(res_tab)
  res_tab <- res_tab %>% filter(!is.na(SYMBOL))
  universe <- unique(res_tab$SYMBOL)    
  up   <- res_tab %>% filter(adj.P.Val < fdr_cut, logFC >  lfc_cut) %>% pull(SYMBOL) %>% unique()
  down <- res_tab %>% filter(adj.P.Val < fdr_cut, logFC < -lfc_cut) %>% pull(SYMBOL) %>% unique()
  list(up = up, down = down, universe = universe)
}

# Run enrichGO (BP/MF/CC optional) and simplify the results
run_go <- function(genes, universe, ont = "BP", title = "") {
  if (length(genes) < 5) return(NULL)  
  ego <- enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
}

plot_dot <- function(ego, top = 20, title = "") {
  if (is.null(ego)) return(NULL)
  p <- dotplot(ego, showCategory = top, font.size = 12) +
    labs(title = title) +
    theme_minimal(base_size = 12) +
    theme(panel.border = element_rect(colour="grey40", fill=NA, linewidth=0.9))
  print(p)
}

## ==== 2) GDF15 vs Veh：up-/down-regulation ====
sets_g <- pick_sig_sets(res_GDF15_Veh)
ego_g_up_BP   <- run_go(sets_g$up,   sets_g$universe, ont="BP", title="GDF15 vs Veh (Up, BP)")
ego_g_down_BP <- run_go(sets_g$down, sets_g$universe, ont="BP", title="GDF15 vs Veh (Down, BP)")


## ==== 3) CR vs Veh：up-/down-regulation ====
sets_c <- pick_sig_sets(res_CR_Veh)
ego_c_up_BP   <- run_go(sets_c$up,   sets_c$universe, ont="BP", title="CR vs Veh (Up, BP)")
ego_c_down_BP <- run_go(sets_c$down, sets_c$universe, ont="BP", title="CR vs Veh (Down, BP)")

## ==== 4) run CC、MF ====
# ego_g_up_CC   <- run_go(sets_g$up,   sets_g$universe, ont="CC", title="GDF15 vs Veh (Up, CC)")
# ego_g_up_MF   <- run_go(sets_g$up,   sets_g$universe, ont="MF", title="GDF15 vs Veh (Up, MF)")


## ==== 5) visulization（Dotplot）====
plot_dot(ego_g_up_BP,   top=12, title="GDF15 vs Veh — Upregulated (GO:BP)")
graph2tif(x = NULL, file='5_GDF15 vs Veh — Upregulated', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 6, bg = "transparent")

plot_dot(ego_g_down_BP, top=12, title="GDF15 vs Veh — Downregulated (GO:BP)")
graph2tif(x = NULL, file='5_GDF15 vs Veh — Downregulated', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 6, bg = "transparent")

plot_dot(ego_c_up_BP,   top=12, title="CR vs Veh — Upregulated (GO:BP)")
plot_dot(ego_c_down_BP, top=12, title="CR vs Veh — Downregulated (GO:BP)")


## ==== 6) export results ====
save_ego <- function(ego, file){
  if (!is.null(ego)) write.csv(as.data.frame(ego), file, row.names = FALSE)
}
save_ego(ego_g_up_BP,   "GO_BP_GDF15_vsVeh_Up.csv")
save_ego(ego_g_down_BP, "GO_BP_GDF15_vsVeh_Down.csv")
save_ego(ego_c_up_BP,   "GO_BP_CR_vsVeh_Up.csv")
save_ego(ego_c_down_BP, "GO_BP_CR_vsVeh_Down.csv")




### only show immune related pathways in GDF15 vs Veh — Downregulated (GO:BP)
library(dplyr)
library(enrichplot)
library(ggplot2)

## 1) keep following pathways
keep_terms <- c(
  "L-amino acid metabolic process",
  "regulation of cell killing",
  "carboxylic acid catabolic process",
  "T cell mediated immunity",
  "leukocyte mediated cytotoxicity",
  "positive regulation of leukocyte mediated immunity",
  "regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains"
)

## ===== 2) Filter enrichResult: Retain Only Selected Terms =====
filter_ego_by_terms <- function(ego, terms) {
  if (is.null(ego)) return(NULL)
  df <- as.data.frame(ego)
  df2 <- df %>% filter(Description %in% terms)
  
  missing_terms <- setdiff(terms, df2$Description)
  if (length(missing_terms)) message("Not found: ", paste(missing_terms, collapse = "; "))
  
  if (nrow(df2) == 0) return(NULL)
  
  # Preserve the custom order
  df2$Description <- factor(df2$Description, levels = terms[terms %in% df2$Description])
  
  ego2 <- ego
  ego2@result <- df2        
  ego2
}

## 3) Object: ego_g_down_BP → Filter and Plot
ego_sel <- filter_ego_by_terms(ego_g_down_BP, keep_terms)

fs <- 12
wrap_width <- 80   # adjust if needed (bigger = fewer lines)

p_down_selected <-
  enrichplot::dotplot(
    ego_sel,
    showCategory = Inf,
    orderBy = "x",
    font.size = fs - 1,
    label_format = wrap_width   # <-- only wrap here
  ) +
  labs(title = "GDF15 vs Veh_Down (GO)", x = "Ratio") +
  theme_minimal(base_size = fs) +
  theme(
    plot.title   = element_text(size = fs, hjust = 0.5),
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.9),
    legend.position = "bottom",           # more horizontal space for labels
    axis.text.y  = element_text(size = fs - 1, lineheight = 1.1),
    plot.margin  = margin(t = 8, r = 12, b = 8, l = 16)
  ) +
  coord_cartesian(clip = "off")            # let long labels extend if needed

p_down_selected

graph2tif(x = NULL, file='5_GDF15 vs Veh — Down_Immune', font = "Arial", cairo = TRUE,   
          width = 8.5, height = 3.5, bg = "transparent")

graph2svg(x = NULL, file='5_GDF15 vs Veh — Down_Immune', font = "Arial", cairo = TRUE,   
          width = 8.5, height = 3.5, bg = "transparent")
p=p_down_selected
ggsave("5_GDF15 vs Veh — Down_Immune_1.svg", p, width = 8.5, height = 3.5, device = svglite::svglite)


### Run enrichment analysis using both up and down proteins
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(enrichplot)
library(ggplot2)
library(patchwork)

fs   <- 12  
topN <- 20  


mk_dot <- function(ego, title, top = 20, fs = 12){
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    return(ggplot() + theme_void() +
             labs(title = paste0(title, " (no significant terms)")) +
             theme(plot.title = element_text(size = fs, hjust = 0.5)))
  }
  enrichplot::dotplot(ego, showCategory = top, font.size = fs - 1) +
    labs(title = title) +
    theme_minimal(base_size = fs) +
    theme(
      plot.title   = element_text(size = fs, hjust = 0.5),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.9)
    )
}

## —— GDF15 vs Veh：merge Up/Down —— 
sets_g <- pick_sig_sets(res_GDF15_Veh)
sig_g_all <- unique(c(sets_g$up, sets_g$down))
ego_g_all_BP <- run_go(sig_g_all, sets_g$universe, ont = "BP")

## —— CR vs Veh：merge Up/Down —— 
sets_c <- pick_sig_sets(res_CR_Veh)
sig_c_all <- unique(c(sets_c$up, sets_c$down))
ego_c_all_BP <- run_go(sig_c_all, sets_c$universe, ont = "BP")

## ===== Two Dot Plots Stacked Vertically (Top = GDF15, Bottom = CR) =====
p_g_all <- mk_dot(ego_g_all_BP,
                  sprintf("GDF15 vs Veh — All DE (n = %d)  (GO:BP)", length(sig_g_all)),
                  topN, fs)
p_c_all <- mk_dot(ego_c_all_BP,
                  sprintf("CR vs Veh — All DE (n = %d)  (GO:BP)", length(sig_c_all)),
                  topN, fs)

p_all <- (p_g_all / p_c_all) + plot_layout(guides = "collect") & theme(legend.position = "right")

p_all

# save
ggsave("GO_BP_AllDE_GDF15_on_top_CR_bottom.png", p_all,
       width = 9, height = 12, dpi = 300, bg = "white")








####===== Dotplot showing immune proteins =======##
######### Dotplot for protein related to inflammation
# Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# target protein
target_prot <- c(
  "PGLYRP2","H2-Q4","SERPINF2","SAA1","CTSE","H2-D1", "H2-Q10",
  "GSN","HSP90AB1","GPX3","CES1B","CES1C"
)

prep_df <- function(df, label, targets){
  stopifnot(all(c("Protein","logFC","adj.P.Val") %in% names(df)))
  df %>%
    dplyr::mutate(
      Protein = toupper(.data$Protein),
      Log2FC  = .data$logFC,
      NegFDR  = -log10(pmax(.data$adj.P.Val, .Machine$double.xmin))
    ) %>%
    
    dplyr::group_by(Protein) %>%
    dplyr::slice_min(order_by = .data$adj.P.Val, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$Protein %in% toupper(targets)) %>%
    dplyr::select(Protein, Log2FC, NegFDR) %>%
    dplyr::mutate(Dataset = label)
}

d1 <- prep_df(res_GDF15_Veh, "GDF15 vs Vehicle", target_prot)
d2 <- prep_df(res_CR_Veh,    "CR vs Vehicle",    target_prot)
d3 <- prep_df(res_GDF15_CR,  "GDF15 vs CR",      target_prot)

dotdat <- dplyr::bind_rows(d1, d2, d3) %>%
  dplyr::mutate(
    Protein = factor(Protein, levels = toupper(target_prot)),
    Dataset = factor(Dataset, levels = c("GDF15 vs Vehicle","CR vs Vehicle","GDF15 vs CR"))
  ) %>%
  tidyr::complete(Dataset, Protein)

missing <- dotdat %>% dplyr::filter(is.na(Log2FC)) %>% dplyr::distinct(Dataset, Protein)
if (nrow(missing) > 0) {
  message("Missing entries:\n", paste0(missing$Dataset, " : ", missing$Protein, collapse = "\n"))
}

# Plot
ggplot(dotdat, aes(x = Protein, y = Dataset, color = Log2FC, size = NegFDR)) +
  geom_point(na.rm = TRUE) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "log2FC") +
  scale_size_continuous(range = c(2.5, 9), name = expression(-log[10]("FDR"))) +
  labs(title = "Protein expression across contrasts", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title         = element_text(face = "bold", hjust = 0),
    panel.border       = element_rect(colour = "grey", fill = NA, linewidth = 1.2)
  )

graph2tif(x = NULL, file='2_Protein_inflam_2', font = "Arial", cairo = TRUE,   
          width = 6, height = 3.9, bg = "transparent")

graph2svg(x = NULL, file='2_Protein_inflam_2', font = "Arial", cairo = TRUE,   
          width = 6, height = 3.9, bg = "transparent")

p=ggplot(dotdat, aes(x = Protein, y = Dataset, color = Log2FC, size = NegFDR)) +
  geom_point(na.rm = TRUE) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "log2FC") +
  scale_size_continuous(range = c(2.5, 9), name = expression(-log[10]("FDR"))) +
  labs(title = "Protein expression across contrasts", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title         = element_text(face = "bold", hjust = 0),
    panel.border       = element_rect(colour = "grey", fill = NA, linewidth = 1.2)
  )
ggsave("2_Protein_inflam_2.svg", p, width = 6, height = 3.9, device = svglite::svglite)






####============= liver-secreted proteome ===============####
## Prepare liver-secreted proteome database
ProtDB_H <- read.csv("HumanLiverSecretedProtDB.csv", header=TRUE)
head(ProtDB_H)

ProtDB_M <- read.csv("MouseLiverSecretedProtDB.csv", header=TRUE)
head(ProtDB_M)


## Extract gene or protein names column from 'ProtDB_H' and convert to mouse gene protein
## =========== 0) Packages ===========
# BiocManager::install("biomaRt")
## ========= Human PG.Genes -> Mouse MGI =========
library(dplyr)
library(stringr)
library(biomaRt)

## ===== 0) Clean Human Gene Symbols =====
human <- ProtDB_H %>%
  transmute(sym_raw = as.character(PG.Genes)) %>%
  filter(!is.na(sym_raw), sym_raw != "") %>%
  mutate(sym = str_split(sym_raw, ";", simplify = TRUE)[,1] |> str_trim()) %>%
  distinct(sym, .keep_all = TRUE) %>%
  filter(!str_detect(sym, "^IG[HKL]"))  

## ===== 1) Connect to the Human BioMart (Use a Mirror if the Main Server Fails) =====
bm_connect <- function(mirror = c("www","uswest","asia","useast")){
  mirror <- match.arg(mirror)
  if (mirror == "www") {
    useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  } else {
    useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = mirror)
  }
}
humMart <- tryCatch(bm_connect("www"),
                    error = function(e) tryCatch(bm_connect("uswest"),
                                                 error = function(e) tryCatch(bm_connect("asia"),
                                                                              error = function(e) bm_connect("useast"))))

## ===== 2) Step 1: Map HGNC Symbols to Ensembl Gene IDs
##          (Using Attributes from the "Genes" Dataset Only) =====
id_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters    = "external_gene_name",
  values     = human$sym,
  mart       = humMart
) %>%
  as_tibble() %>%
  rename(ensembl = ensembl_gene_id, human = external_gene_name) %>%
  distinct(human, .keep_all = TRUE)

## ===== 3) Step 2: Retrieve Mouse Orthologs Using Ensembl Gene IDs
##          (Using Attributes from the "Homology" Dataset Only) =====
hom <- getBM(
  attributes = c("ensembl_gene_id",
                 "mmusculus_homolog_associated_gene_name",
                 "mmusculus_homolog_orthology_type",
                 "mmusculus_homolog_orthology_confidence"),
  filters    = "ensembl_gene_id",
  values     = id_map$ensembl,
  mart       = humMart
) %>%
  as_tibble() %>%
  rename(ensembl = ensembl_gene_id,
         mouse   = mmusculus_homolog_associated_gene_name,
         type    = mmusculus_homolog_orthology_type,
         conf    = mmusculus_homolog_orthology_confidence) %>%
  filter(mouse != "")

## ===== 4) Merge the Two Mapping Steps and Prioritize High-Confidence One-to-One Orthologs =====
map_h2m <- id_map %>%
  inner_join(hom, by = "ensembl") %>%
  arrange(desc(conf), desc(type == "ortholog_one2one")) %>%
  distinct(human, .keep_all = TRUE) %>%
  dplyr::select(human, mouse, type, conf, ensembl)

## ===== 5) Add Mouse Orthologs Back to the Original Table
##          (Add a mouse column to ProtDB_H) =====
ProtDB_H_map <- ProtDB_H %>%
  left_join(map_h2m, by = c("PG.Genes" = "human"))

## view
ProtDB_H_map %>%
  dplyr::select(PG.Genes, mouse, type, conf) %>%
  head()

write.csv(ProtDB_H_map, "ProtDB_H_map2M.csv", row.names = FALSE)



####### load our proteomics data
Prot_A <- read.csv("PlasmaProt_Mouse_GDF15CR_clean.csv", header=TRUE)
head(Prot_A)

## according to ProtDB_H_map$mouse, select all data from (Prot_All)
library(dplyr)
library(tidyr)
library(stringr)

## ===== 1) Extract the Mouse Gene List (Remove NAs, Empty Values, and Duplicates) =====
mouse_genes <- ProtDB_H_map %>%
  dplyr::filter(!is.na(mouse), mouse != "") %>%
  dplyr::distinct(mouse) %>%
  dplyr::pull(mouse)
head(mouse_genes)

## ===== 2) Find Rows in Prot_A Where PG.Genes Contains the Mouse Genes =====
# First split PG.Genes into a long-format table and assign the original row ID
Prot_A_long <- Prot_A %>%
  dplyr::mutate(.row_id = dplyr::row_number()) %>%
  tidyr::separate_rows(PG.Genes, sep = ";") %>%
  dplyr::mutate(PG.Genes = stringr::str_trim(PG.Genes)) %>%
  dplyr::filter(PG.Genes != "")

hit_rows <- Prot_A_long %>%
  dplyr::filter(PG.Genes %in% mouse_genes) %>%
  dplyr::distinct(.row_id)

# 3) Extract the corresponding rows from the original table
# (i.e., rows where PG.Genes contains at least one mouse gene)
Prot_A_matched <- Prot_A %>%
  dplyr::mutate(.row_id = dplyr::row_number()) %>%
  dplyr::semi_join(hit_rows, by = ".row_id") %>%
  dplyr::select(-.row_id)
head(Prot_A_matched)
dim(Prot_A_matched)
write.csv(Prot_A_matched, "Prot_liverSecretedProt.csv", row.names = FALSE)


### PCA analysis for liver-secreted protein
# apply a log2 transformation,
# perform left-censored (MNAR) imputation for missing values,
# and carry out median normalization across samples;
# then plot a PC1–PC2 scatter with a 68% confidence ellipse and
# a variance-explained plot for the top 10 principal components.

## ===== 0) Dependence =====
pkgs <- c("dplyr","ggplot2","ggrepel")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

## ===== 1) Extract the Expression Matrix and Group Labels
##          (Using Protein as Row Names) =====
prot <- Prot_A_matched
stopifnot("Protein" %in% names(prot))

# Sample column：Veh./CR./GDF15.
sample_cols <- grep("^(Veh|CR|GDF15)\\.", names(prot), value = TRUE)

# Use Protein as row names; automatically make duplicate row names unique
rn <- make.unique(prot$Protein)
X  <- as.matrix(prot[, sample_cols])
mode(X) <- "numeric"
X[is.nan(X)] <- NA_real_
rownames(X) <- rn

samples <- colnames(X)
group   <- factor(sub("\\..*$", "", samples), levels = c("Veh","CR","GDF15"))

## ===== 2) Filter by Detection Rate
##          (Quantified in ≥3 Samples in at Least Two Groups) =====
min_n <- 3
by_grp_ok <- sapply(levels(group), function(g){
  rowSums(!is.na(X[, group==g, drop=FALSE])) >= min_n
})
keep <- rowSums(by_grp_ok) >= 2
X <- X[keep, , drop=FALSE]

## ===== 3) Log2 Transformation + Left-Censored Imputation + Between-Sample Median Normalization =====
pos_vals <- X[X > 0 & is.finite(X)]
pseudo   <- as.numeric(quantile(pos_vals, 0.01, na.rm = TRUE))/2
L <- log2(X + pseudo)

L_preimp <- L

for (j in seq_len(ncol(L))) {
  colj  <- L[, j]
  q1    <- suppressWarnings(quantile(colj, 0.01, na.rm = TRUE))
  fillv <- if (is.finite(q1)) q1 - 1 else min(colj, na.rm = TRUE) - 1
  colj[!is.finite(colj) | is.na(colj)] <- fillv
  L[, j] <- colj
}

Lnorm <- sweep(L, 2, apply(L, 2, median, na.rm = TRUE), FUN = "-")

nzv <- apply(Lnorm, 1, sd) > 0
Lnorm <- Lnorm[nzv, , drop=FALSE]

## ===== 3.5) Quality Control (QC) Metrics
##          (Using the Pre-Imputation Matrix: L_preimp) =====
qc_tab <- data.frame(
  Sample       = samples,
  Group        = group,
  missing_frac = colMeans(!is.finite(L_preimp) | is.na(L_preimp)),
  median_int   = apply(L_preimp, 2, median, na.rm = TRUE),
  prot_count   = colSums(is.finite(L_preimp) & !is.na(L_preimp))
)

## ===== 4) PCA =====
pc <- prcomp(t(Lnorm), center = TRUE, scale. = TRUE)
var_expl <- (pc$sdev^2) / sum(pc$sdev^2)

pc_df <- data.frame(pc$x[, 1:2], Sample = samples, Group = group,
                    PC1 = pc$x[,1], PC2 = pc$x[,2])

## ===== 5) Plot：PC1–PC2 =====
fs <- 12  
p_pca <- ggplot(pc_df, aes(PC1, PC2, color = Group, shape = Group, label = Sample)) +
  stat_ellipse(type = "norm", level = 0.68, linewidth = 0.8, show.legend = FALSE) +
  geom_point(size = 3, stroke = 0.7) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 30) +
  labs(
    title = "Liver-secreted proteins (log2)",
    x = sprintf("PC1 (%.1f%%)", 100*var_expl[1]),
    y = sprintf("PC2 (%.1f%%)", 100*var_expl[2])
  ) +
  theme_minimal(base_size = fs) +
  theme(
    text = element_text(size = fs, face = "plain"),
    axis.title = element_text(size = fs, face = "plain"),
    axis.text  = element_text(size = fs, face = "plain"),
    legend.title = element_text(size = fs, face = "plain"),
    legend.text  = element_text(size = fs, face = "plain"),
    plot.title   = element_text(size = fs, face = "plain", hjust = 0.5),
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.9)
  )
print(p_pca)
graph2tif(x = NULL, file='6_PCA analysis_LiverSecPro', font = "Arial", cairo = TRUE,   
          width = 6.5, height = 4, bg = "transparent")

graph2svg(x = NULL, file='6_PCA analysis_LiverSecPro', font = "Arial", cairo = TRUE,   
          width = 6.5, height = 4, bg = "transparent")

## ===== 6) Scree plot =====
k <- min(10, length(var_expl))
scree_df <- data.frame(PC = factor(seq_len(k)), Var = var_expl[1:k])
p_scree <- ggplot(scree_df, aes(PC, Var)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.1f%%", 100*Var)), vjust = -0.3, size = 3) +
  labs(title = "Explained variance", y = "Proportion of variance", x = NULL) +
  theme_minimal(base_size = fs) +
  theme(
    text = element_text(size = fs, face = "plain"),
    panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.9)
  )
print(p_scree)






## Different expression analysis_liver secreted protein
# Clean → log2 transform → left-censored (MNAR) imputation for missing values →
# between-sample normalization → limma linear modeling → 
# pairwise comparisons among the three groups, with example volcano plots and heatmaps provided.
## ===== 0) Dependence =====
pkgs <- c("dplyr","stringr","limma","ggplot2","ggrepel","pheatmap","tibble")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

## ===== 1) Extract the Expression Matrix and Group Labels
##          (Using Protein as the Feature Name) =====
prot <- Prot_A_matched
sample_cols <- grep("^(Veh|CR|GDF15)\\.", names(prot), value = TRUE)

name_raw <- as.character(prot$Protein)
name_raw <- trimws(name_raw)

# Use Protein as row names (automatically make duplicates unique if necessary)
# Use the first gene listed in PG.Genes
mask <- is.na(name_raw) | name_raw == ""
name_raw[mask] <- trimws(sub(";.*$", "", prot$`PG.Genes`[mask]))

rn <- make.unique(name_raw)
stopifnot(!any(is.na(rn) | rn == ""))

## Reconstruct the expression matrix and normalize
X  <- as.matrix(prot[, sample_cols])
mode(X) <- "numeric"
X[is.nan(X)] <- NA_real_
rownames(X) <- rn

samples <- colnames(X)
group   <- factor(sub("\\..*$", "", samples), levels = c("Veh","CR","GDF15"))

annot <- prot %>%
  transmute(Protein = rn,
            PG.Genes = `PG.Genes`,
            PG.ProteinAccessions = `PG.ProteinAccessions`) %>%
  tibble::remove_rownames()

## ===== 2) Filtering (Robust Detection) =====
# Retain proteins quantified in ≥3 samples per group in at least 2 groups
min_n <- 3
by_grp_ok <- sapply(levels(group), function(g){
  rowSums(!is.na(X[, group==g, drop=FALSE])) >= min_n
})
keep <- rowSums(by_grp_ok) >= 2
X <- X[keep, , drop=FALSE]
annot <- annot[keep, , drop=FALSE]

## ===== 3) Log2 Transformation + Left-Censored Imputation + Between-Sample Median Normalization =====
# Apply log2(x + pseudocount), where the pseudocount is set to a global low-abundance value
# (half of the 1st percentile across all quantified values)
pos_vals <- X[X > 0 & is.finite(X)]
pseudo   <- as.numeric(quantile(pos_vals, 0.01, na.rm = TRUE))/2
L <- log2(X + pseudo)

# Left-censored imputation: replace missing values with the 1st percentile of each column minus 1
# (a more conservative imputation toward low-abundance values)
for (j in seq_len(ncol(L))) {
  colj  <- L[, j]
  q1    <- suppressWarnings(quantile(colj, 0.01, na.rm = TRUE))
  fillv <- if (is.finite(q1)) q1 - 1 else min(colj, na.rm = TRUE) - 1
  colj[!is.finite(colj) | is.na(colj)] <- fillv
  L[, j] <- colj
}

Lnorm <- sweep(L, 2, apply(L, 2, median, na.rm = TRUE), FUN = "-")

nzv <- apply(Lnorm, 1, sd) > 0
Lnorm <- Lnorm[nzv, , drop=FALSE]
annot <- annot[nzv, , drop=FALSE]

write.csv(Lnorm, "Prot_liverSecreted_clean.csv", row.names = FALSE)


## ===== 4) limma: Three Pairwise Comparisons =====
design <- model.matrix(~ 0 + group); colnames(design) <- levels(group)
fit  <- limma::lmFit(Lnorm, design)
cont <- limma::makeContrasts(
  GDF15_vs_Veh = GDF15 - Veh,
  CR_vs_Veh    = CR    - Veh,
  GDF15_vs_CR  = GDF15 - CR,
  levels = design
)
fit2 <- limma::contrasts.fit(fit, cont)
fit2 <- limma::eBayes(fit2, robust = TRUE, trend = TRUE)

# Function to generate result tables (merge annotations and sort by P value)
get_tab <- function(coef_name){
  tt <- limma::topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")
  tt <- tt %>%
    tibble::rownames_to_column("Protein") %>%
    left_join(annot, by = "Protein") %>%
    relocate(Protein, PG.Genes, PG.ProteinAccessions, .before = 1)
  tt$contrast <- coef_name
  tt$abs_logFC <- abs(tt$logFC)
  tt
}
res_GDF15_Veh <- get_tab("GDF15_vs_Veh")
res_CR_Veh    <- get_tab("CR_vs_Veh")
res_GDF15_CR  <- get_tab("GDF15_vs_CR")
res_all <- bind_rows(res_GDF15_Veh, res_CR_Veh, res_GDF15_CR)

## ===== 5) Statistical Significance Filtering and Export =====
fdr_cut <- 0.1; lfc_cut <- log2(1.5)  
sig_all <- res_all %>% filter(adj.P.Val < fdr_cut, abs_logFC >= lfc_cut)

# Export
write.csv(res_GDF15_Veh, "DE_GDF15_vs_Veh_Secr.csv", row.names = FALSE)
write.csv(res_CR_Veh,    "DE_CR_vs_Veh_all_Secr.csv",    row.names = FALSE)
write.csv(res_GDF15_CR,  "DE_GDF15_vs_CR_all_Secr.csv",  row.names = FALSE)
write.csv(sig_all,       "DE_all_contrasts_sig_Secr.csv", row.names = FALSE)

## ===== 6) Visualization: Volcano Plot & Heatmap =====
## Significance thresholds
library(ggplot2)
library(ggrepel)

volcano_plot <- function(tab, title, fdr_cut = 0.1, lfc_cut = log2(1.5), fs = 12, label_shrink = 0.85){
  
  tab <- tab %>%
    dplyr::mutate(
      cat = dplyr::case_when(
        adj.P.Val < fdr_cut & logFC >  lfc_cut  ~ "Up",
        adj.P.Val < fdr_cut & logFC < -lfc_cut  ~ "Down",
        TRUE ~ "NS"
      ),
      sig = cat != "NS"
    )
  
  up   <- sum(tab$cat == "Up",   na.rm = TRUE)
  down <- sum(tab$cat == "Down", na.rm = TRUE)
  total <- up + down
  
  label_size <- (fs/3.3) * label_shrink
  
  tab_label <- tab %>%
    filter(sig) %>%
    arrange(desc(abs(logFC))) %>%             
    head(20)
  
  g <- ggplot(tab, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = cat), alpha = 0.9, size = 1.8) +
    scale_color_manual(values = c(Up = "#E64B35", Down = "#3B8EEA", NS = "grey75")) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2, color = "grey50") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = 2, color = "grey50") +
    
    ggrepel::geom_text_repel(
      data = tab_label,
      aes(label = Protein, color = cat),
      size = label_size, max.overlaps = Inf,
      box.padding = 0.35, point.padding = 0.25, segment.size = 0.3,
      show.legend = FALSE
    ) +
    labs(
      title = title,
      subtitle = sprintf("Up = %d | Down = %d | Total = %d  (FDR < %.2g; |log2FC| ≥ %.2f)",
                         up, down, total, fdr_cut, lfc_cut),
      x = "Log2FoldChange", y = expression(-log[10]("FDR"))
    ) +
    theme_minimal(base_size = fs) +
    theme(
      text         = element_text(size = fs, face = "plain"),
      axis.title   = element_text(size = fs, face = "plain"),
      axis.text    = element_text(size = fs, face = "plain"),
      legend.title = element_text(size = fs, face = "plain"),
      legend.text  = element_text(size = fs, face = "plain"),
      plot.title   = element_text(size = fs, face = "plain", hjust = 0.5),
      plot.subtitle= element_text(size = fs, face = "plain"),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.9)
    ) +
    guides(color = "none")
  
  
  xmax <- max(tab$logFC, na.rm = TRUE); xmin <- min(tab$logFC, na.rm = TRUE)
  ymax <- max(-log10(tab$adj.P.Val), na.rm = TRUE)
  g + annotate("label",
               x = xmax - 0.02*(xmax - xmin), y = ymax,
               hjust = 1, vjust = 1,
               label = sprintf("Up: %d\nDown: %d\nTotal: %d", up, down, total),
               size = fs/3.3, label.size = 0.25, fill = "white")
}

p1 <- volcano_plot(res_GDF15_Veh, "GDF15 vs Veh_liver secreted")
p2 <- volcano_plot(res_CR_Veh,    "CR vs Veh_liver secreted")
p3 <- volcano_plot(res_GDF15_CR,  "GDF15 vs CR_liver secreted")
print(p1)
print(p2)
print(p3)

# save
graph2tif(x = p1, file='7_volcano_GDF15 vs Veh_LivSecr', font = "Arial", cairo = TRUE,   
          width = 5.5, height = 4, bg = "transparent")

graph2tif(x = p2, file='7_volcano_CR vs Veh_LivSecr', font = "Arial", cairo = TRUE,   
          width = 5.5, height = 4, bg = "transparent")

graph2tif(x = p3, file='7_volcano_GDF15 vs CR_LivSecr', font = "Arial", cairo = TRUE,   
          width = 5.5, height = 4, bg = "transparent")

## svg
# save
graph2svg(x = p1, file='7_volcano_GDF15 vs Veh_LivSecr', font = "Arial", cairo = TRUE,   
          width = 5.5, height = 4, bg = "transparent")

graph2svg(x = p2, file='7_volcano_CR vs Veh_LivSecr', font = "Arial", cairo = TRUE,   
          width = 5.5, height = 4, bg = "transparent")

graph2svg(x = p3, file='7_volcano_GDF15 vs CR_LivSecr', font = "Arial", cairo = TRUE,   
          width = 5.5, height = 4, bg = "transparent")




# Heatmap: take the union of significant proteins from the three comparisons
# and display row-wise Z-scores
top_prots <- res_all %>% arrange(adj.P.Val) %>% pull(Protein) %>% unique()
top_prots <- head(top_prots, 40)  
H <- Lnorm[intersect(rownames(Lnorm), top_prots), , drop=FALSE]
ann_col <- data.frame(Group = group); rownames(ann_col) <- colnames(H)

pheatmap::pheatmap(H, scale = "row", annotation_col = ann_col,
                   show_rownames = TRUE, show_colnames = FALSE,
                   main = "Liver secreted proteins (row Z-score)")

graph2tif(x = NULL, file='8_Heatmap_livSecr', font = "Arial", cairo = TRUE,   
          width = 7, height = 7, bg = "transparent")





###=========== Plot Venn ================
library(dplyr)
fdr_cut <- 0.1
lfc_cut <- log2(1.5)

if (!"Protein" %in% names(res_GDF15_Veh)) res_GDF15_Veh$Protein <- rownames(res_GDF15_Veh)
if (!"Protein" %in% names(res_CR_Veh))    res_CR_Veh$Protein    <- rownames(res_CR_Veh)

sig_gdf15 <- subset(res_GDF15_Veh, adj.P.Val < fdr_cut & abs(logFC) >= lfc_cut)
sig_cr    <- subset(res_CR_Veh,    adj.P.Val < fdr_cut & abs(logFC) >= lfc_cut)

set_gdf15 <- unique(sig_gdf15$Protein)
set_cr    <- unique(sig_cr$Protein)

# Statistic
n_gdf15 <- length(set_gdf15)
n_cr    <- length(set_cr)
n_inter <- length(intersect(set_gdf15, set_cr))
cat("# GDF15 vs Veh:", n_gdf15, "\n# CR vs Veh:", n_cr, "\n# overlap:", n_inter, "\n")

common_proteins <- intersect(set_gdf15, set_cr)
unique_gdf15    <- setdiff(set_gdf15, set_cr)
unique_cr       <- setdiff(set_cr, set_gdf15)

## ===== Plot Venn =====
library(VennDiagram)
library(grid)

if (names(dev.cur()) %in% c("png","pdf","jpeg","tiff","bmp","cairo_pdf")) dev.off()

grid.newpage()
vp <- viewport(width = 0.90, height = 0.60) 
pushViewport(vp)

g <- venn.diagram(
  x = list(`GDF15 vs Veh` = set_gdf15, `CR vs Veh` = set_cr),
  filename = NULL,
  fill = c("#9FCBFF", "#FF8F8F"),
  alpha = 0.70, lty = "blank",
  cex = 1.4, cat.cex = 0
)

grid.draw(g)

fs <- 12
grid.text(sprintf("GDF15 vs Vehicle (n = %d)", length(unique(set_gdf15))),
          x = 0.40, y = 0.95, gp = gpar(fontsize = fs, col = "#2B6E9E"))
grid.text(sprintf("CR vs Vehicle (n = %d)", length(unique(set_cr))),
          x = 0.78, y = 0.78, gp = gpar(fontsize = fs, col = "#FF6F6F"))
grid.text(sprintf("Cutoffs: FDR < %.2g, |log2FC| ≥ %.2f", fdr_cut, lfc_cut),
          x = 0.50, y = 0.04, gp = gpar(fontsize = fs, col = "grey25"))

popViewport()

graph2png(x = NULL, file='9_VennPlot_commonProtein_livSecr', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 4.5, bg = "transparent")

graph2svg(x = NULL, file='9_VennPlot_commonProtein_livSecr', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 4.5, bg = "transparent")





###=========== Enrichment analysis ============###
## ==== 0) Dependence ====
if (!requireNamespace("clusterProfiler", quietly=TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Mm.eg.db", quietly=TRUE))   BiocManager::install("org.Mm.eg.db")
if (!requireNamespace("enrichplot", quietly=TRUE))      BiocManager::install("enrichplot")
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(dplyr)
library(stringr)
library(ggplot2)

## ==== 1) Common Parameters and Utility Functions ====
fdr_cut <- 0.1
lfc_cut <- log2(1.1)

get_symbol <- function(df){
  has_pg <- "PG.Genes" %in% names(df)
  sym <- if (has_pg) sub(";.*$", "", df$PG.Genes) else df$Protein
  sym <- trimws(sym)
  sym[sym=="" | is.na(sym)] <- NA
  sym
}

# Extract significant genes and split into up- and downregulated sets
# (returns SYMBOL vectors)
pick_sig_sets <- function(res_tab){
  stopifnot(all(c("adj.P.Val","logFC") %in% names(res_tab)))
  res_tab$SYMBOL <- get_symbol(res_tab)
  res_tab <- res_tab %>% filter(!is.na(SYMBOL))
  universe <- unique(res_tab$SYMBOL)    
  up   <- res_tab %>% filter(adj.P.Val < fdr_cut, logFC >  lfc_cut) %>% pull(SYMBOL) %>% unique()
  down <- res_tab %>% filter(adj.P.Val < fdr_cut, logFC < -lfc_cut) %>% pull(SYMBOL) %>% unique()
  list(up = up, down = down, universe = universe)
}

# Run enrichGO
run_go <- function(genes, universe, ont = "BP", title = "") {
  if (length(genes) < 5) return(NULL)  
  ego <- enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(NULL)
  simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
}

plot_dot <- function(ego, top = 20, title = "") {
  if (is.null(ego)) return(NULL)
  p <- dotplot(ego, showCategory = top, font.size = 12) +
    labs(title = title) +
    theme_minimal(base_size = 12) +
    theme(panel.border = element_rect(colour="grey40", fill=NA, linewidth=0.9))
  print(p)
}

## ==== 2) GDF15 vs Veh：up-/down-regulation ====
sets_g <- pick_sig_sets(res_GDF15_Veh)
ego_g_up_BP   <- run_go(sets_g$up,   sets_g$universe, ont="BP", title="GDF15 vs Veh (Up, BP)")
ego_g_down_BP <- run_go(sets_g$down, sets_g$universe, ont="BP", title="GDF15 vs Veh (Down, BP)")


## ==== 3) CR vs Veh：up-/down-regulation ====
sets_c <- pick_sig_sets(res_CR_Veh)
ego_c_up_BP   <- run_go(sets_c$up,   sets_c$universe, ont="BP", title="CR vs Veh (Up, BP)")
ego_c_down_BP <- run_go(sets_c$down, sets_c$universe, ont="BP", title="CR vs Veh (Down, BP)")

## ==== 5) Dotplot====
plot_dot(ego_g_up_BP,   top=12, title="GDF15 vs Veh — Upregulated (GO:BP)")
graph2tif(x = NULL, file='5_GDF15 vs Veh — Upregulated_LivSecr', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 6, bg = "transparent")

plot_dot(ego_g_down_BP, top=15, title="GDF15 vs Veh — Downregulated (GO:BP)")
graph2tif(x = NULL, file='5_GDF15 vs Veh — Downregulated_LivSecr', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 6, bg = "transparent")

plot_dot(ego_c_up_BP,   top=12, title="CR vs Veh — Upregulated (GO:BP)")
plot_dot(ego_c_down_BP, top=12, title="CR vs Veh — Downregulated (GO:BP)")


## ==== 6) Export results ====
save_ego <- function(ego, file){
  if (!is.null(ego)) write.csv(as.data.frame(ego), file, row.names = FALSE)
}
save_ego(ego_g_up_BP,   "GO_BP_GDF15_vsVeh_Up_LivSecr.csv")
save_ego(ego_g_down_BP, "GO_BP_GDF15_vsVeh_Down_LivSecr.csv")
save_ego(ego_c_up_BP,   "GO_BP_CR_vsVeh_Up_LivSecr.csv")
save_ego(ego_c_down_BP, "GO_BP_CR_vsVeh_Down_LivSecr.csv")






####===== Dotplot showing immune proteins_liver secreted =======##
######### Dotplot for protein related to inflammation
# Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# Target protein
target_prot <- c(
  "ECM1", "ECM1.1", "COL1A1",
  "FGL1", "EPHX2", "GSN",
  "C3", "C8A", "C8B", "C8G",
  "AHSG", "AMBP", "ITIH4", "ITIH2", "ITIH3", "CTSD"
)

prep_df <- function(df, label, targets){
  stopifnot(all(c("Protein","logFC","adj.P.Val") %in% names(df)))
  df %>%
    dplyr::mutate(
      Protein = toupper(.data$Protein),
      Log2FC  = .data$logFC,
      NegFDR  = -log10(pmax(.data$adj.P.Val, .Machine$double.xmin))
    ) %>%
    dplyr::group_by(Protein) %>%
    dplyr::slice_min(order_by = .data$adj.P.Val, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$Protein %in% toupper(targets)) %>%
    dplyr::select(Protein, Log2FC, NegFDR) %>%
    dplyr::mutate(Dataset = label)
}

d1 <- prep_df(res_GDF15_Veh, "GDF15 vs Vehicle", target_prot)
d2 <- prep_df(res_CR_Veh,    "CR vs Vehicle",    target_prot)
d3 <- prep_df(res_GDF15_CR,  "GDF15 vs CR",      target_prot)

dotdat <- dplyr::bind_rows(d1, d2, d3) %>%
  dplyr::mutate(
    Protein = factor(Protein, levels = toupper(target_prot)),
    Dataset = factor(Dataset, levels = c("GDF15 vs Vehicle","CR vs Vehicle","GDF15 vs CR"))
  ) %>%
  tidyr::complete(Dataset, Protein)   

missing <- dotdat %>% dplyr::filter(is.na(Log2FC)) %>% dplyr::distinct(Dataset, Protein)
if (nrow(missing) > 0) {
  message("Missing entries:\n", paste0(missing$Dataset, " : ", missing$Protein, collapse = "\n"))
}

# Plot
ggplot(dotdat, aes(x = Protein, y = Dataset, color = Log2FC, size = NegFDR)) +
  geom_point(na.rm = TRUE) +
  scale_x_discrete(limits = toupper(target_prot), drop = FALSE) +  
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "log2FC") +
  scale_size_continuous(range = c(2.5, 9), name = expression(-log[10]("FDR"))) +
  labs(title = "Protein (liver_secreted) expression", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title         = element_text(face = "bold", hjust = 0),
    panel.border       = element_rect(colour = "grey", fill = NA, linewidth = 1.2)
  )

graph2tif(x = NULL, file='10_Protein_inflam_LivSecr', font = "Arial", cairo = TRUE,   
          width = 7, height = 4, bg = "transparent")
graph2svg(x = NULL, file='10_Protein_inflam_LivSecr', font = "Arial", cairo = TRUE,   
          width = 7, height = 4, bg = "transparent")

p=ggplot(dotdat, aes(x = Protein, y = Dataset, color = Log2FC, size = NegFDR)) +
  geom_point(na.rm = TRUE) +
  scale_x_discrete(limits = toupper(target_prot), drop = FALSE) +  
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "log2FC") +
  scale_size_continuous(range = c(2.5, 9), name = expression(-log[10]("FDR"))) +
  labs(title = "Protein (liver_secreted) expression", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title         = element_text(face = "bold", hjust = 0),
    panel.border       = element_rect(colour = "grey", fill = NA, linewidth = 1.2)
  )
ggsave("10_Protein_inflam_LivSecr.svg", p, width = 7, height = 4, device = svglite::svglite)

