########################### Spatial Transcriptomics data analysis for Liver from GDF15/CR. Author: Dongdong Wang (McMaster University) ##############
### the reference website: https://satijalab.org/seurat/articles/spatial_vignette

rm(list = ls())
options(stringsAsFactor = F)


# Creat a result folder
# set as the folder saving all results generated from the following coding

setwd("C:/xxx") # Generate your own folder


## update R
install.packages("installr")
library(installr)
updateR()
sessionInfo()

## This can take approximately 25 min. Install and load necessary packages for Seurat:
install.packages("remotes")
if (!require("hdf5r")) install.packages("hdf5r")
install.packages("Rtools")
install.packages('Seurat')
install.packages("SeuratData")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
install.packages("Signac")
BiocManager::install("SingleR")
BiocManager::install("scRNAseq")
BiocManager::install("celldex")
install.packages('SeuratObject')
install.packages("spatstat.utils")
devtools::install_github("thomasp85/patchwork")
devtools::install_github("RubD/GiottoClass")

remove.packages(grep("spatstat", installed.packages(), value = T))
.rs.restartR()
devtools::install_version("spatstat", version = "3.0.5")


update.packages(oldPkgs = c("withr", "rlang"))
remotes::install_github('satijalab/azimuth', ref = 'master')
BiocManager::install("tximport")
install.packages("tidyverse")
BiocManager::install("tximportData")
install.packages("devtools")
devtools::install_github("arleyc/PCAtest")

## Install LoupeR:
remotes::install_github("10XGenomics/loupeR")
loupeR::setup()

## install GPTCelltype package for Automatic cell type annotation
install.packages("openai")
remotes::install_github("Winnie09/GPTCelltype")

help(spatstat.utils)

BiocManager::install(c("edgeR","limma"))                   # edgeR + dependency


## install scMetabolism for quantifying metabolism activity
install.packages("data.table")
install.packages("wesanderson")
install.packages("AUCell")
install.packages("GSEABase")
install.packages("GSVA")
install.packages("VISION")
remove.packages("promises")
install.packages("promises")
install.packages("devtools")
library(promises)
library(GiottoClass)

remove.packages("VISION")
remove.packages("AUCell")
remove.packages("GSVA")
library(VISION)
library(AUCell)
library(GSVA)

devtools::install_github("YosefLab/VISION@v2.1.0") #Please note that the version would be v2.1.0
devtools::install_github("wu-yc/scMetabolism")


remotes::install_github("10XGenomics/loupeR")
loupeR::setup()
devtools::install_github('satijalab/seurat-data')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("pathview")
BiocManager::install("tximport")
BiocManager::install("tximportData")
devtools::install_github("arleyc/PCAtest")
BiocManager::install("scRNAseq")
BiocManager::install("celldex")
remotes::install_github("Winnie09/GPTCelltype")
BiocManager::install("GSVA")
BiocManager::install("GSEABase")
BiocManager::install("AUCell")
install.packages("promises")
BiocManager::install('glmGamPoi')

#devtools::install_github("YosefLab/VISION")
devtools::install_github("cran/loe")
devtools::install_github("YosefLab/VISION@v2.1.0") #Please note that the version would be v2.1.0
library(GSVA)
library(GSEABase)
library(AUCell)
library(VISION)
library(glmGamPoi)
devtools::install_github("wu-yc/scMetabolism")

BiocManager::install(c(
  "org.Hs.eg.db",   
  "SingleR",        
  "scRNAseq",       
  "celldex"         
))



#library(Azimuth)
library('SeuratObject')
library("spatstat.utils")
library(spatstat.utils)
library(Seurat)
library("hdf5r")

library(patchwork)

library(GEOquery)
library(dplyr)
library(pathview)
library("tximport")
library("readr")
library("tximportData")
library("export")
library(patchwork)
library(tidyverse)

library(SingleR)
library(scRNAseq)
library(celldex)

library(GPTCelltype)
library(openai)

library(ggplot2)
library(rsvd)

library(scMetabolism)
library(loupeR)
library(PCAtest)
library(SeuratData)

sessionInfo()



########### Overview
# this tutorial will cover the following tasks,
# which we believe will be common for many spatial analyses.

#A. Analyze two samples separately
#1.Load two data sets (control and treat) from 10X Space Ranger
#2.Normalization: SCTransform
#3.Compare the difference between two samples

#B. Analyze two samples together (After Integratation)
#1.Normalization: SCTransform
#2.Merge two different data/Working with multiple slices
#3.Normalization: SCTransform
#4.Integration two data sets 
#5.Dimensional reduction and clustering
#6.Chat-GPT4 annotation, Detecting spatially-variable features
#7.Performing any differential expression analysis
#8.Pathway analysis.


############# A. Analyze two samples separately
############# Load two data sets (control and treat) from 10X Space Ranger

############# Working with multiple slices in Seurat

### read it in and perform the same initial normalization
# liver1
liver1=Load10X_Spatial(
  data_dir <- 'C:/BaiduSyncdisk/R analysis/Spatial transcriptomics/GDF15CR/Sam19',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "CR1",
  filter.matrix = TRUE,
  to.upper = FALSE
)


# liver2
liver2=Load10X_Spatial(
  data_dir <- 'C:/BaiduSyncdisk/R analysis/Spatial transcriptomics/GDF15CR/Sam10',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "GDF15a",
  filter.matrix = TRUE,
  to.upper = FALSE
)


# liver3
liver3=Load10X_Spatial(
  data_dir <- 'C:/BaiduSyncdisk/R analysis/Spatial transcriptomics/GDF15CR/Sam18',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "CR2",
  filter.matrix = TRUE,
  to.upper = FALSE
)


# liver4
liver4=Load10X_Spatial(
  data_dir <- 'C:/BaiduSyncdisk/R analysis/Spatial transcriptomics/GDF15CR/Sam11',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "GDF15b",
  filter.matrix = TRUE,
  to.upper = FALSE
)

liver1
liver2
liver3
liver4


# Change orig.ident name and 
head(liver1[[]])
head(liver2[[]])
head(liver3[[]])
head(liver4[[]])

liver1[[]]$orig.ident='CR1'
liver2[[]]$orig.ident='GDF15a'

liver3[[]]$orig.ident='CR2'
liver4[[]]$orig.ident='GDF15b'

#Rename cell identity classes
levels(liver1)
levels(liver2)
levels(liver3)
levels(liver4)
liver1 <- RenameIdents(liver1, 'SeuratProject' = 'CR1')
liver2 <- RenameIdents(liver2, 'SeuratProject' = 'GDF15a')

liver3 <- RenameIdents(liver3, 'SeuratProject' = 'CR2')
liver4 <- RenameIdents(liver4, 'SeuratProject' = 'GDF15b')

#Name project.name
liver1@project.name = "CR1_liver"
liver2@project.name = "GDF15a_liver"
liver3@project.name = "CR2_liver"
liver4@project.name = "GDF15b_liver"



######## Data pre-processing (QC and selecting cells for further analysis)
# for liver1
plot1 <- VlnPlot(liver1, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver1, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv1_gene', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

plot1 <- VlnPlot(liver1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv1_counts', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")



# for liver2
plot1 <- VlnPlot(liver2, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver2, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv2_gene', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

plot1 <- VlnPlot(liver2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver2, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv2_counts', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")



# for liver3
plot1 <- VlnPlot(liver3, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver3, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv3_gene', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

plot1 <- VlnPlot(liver3, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver3, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv3_counts', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")



# for liver4
plot1 <- VlnPlot(liver4, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver4, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv4_gene', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

plot1 <- VlnPlot(liver4, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver4, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv4_counts', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")




# We filter cells that have unique feature counts over 7,500 or less than 1000
liver1 <- subset(liver1, subset = nFeature_Spatial > 1000 & nFeature_Spatial < 7500)
liver2 <- subset(liver2, subset = nFeature_Spatial > 1000 & nFeature_Spatial < 7500)
liver3 <- subset(liver3, subset = nFeature_Spatial > 1000 & nFeature_Spatial < 7500)
liver4 <- subset(liver4, subset = nFeature_Spatial > 1000 & nFeature_Spatial < 7500)


# Apply sctransform normalization
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform.
# Replaces NormalizeData(), ScaleData(), and FindVariableFeatures()
options(future.globals.maxSize = 4 * 1024^3)
liver1SCT <- SCTransform(liver1, method = "glmGamPoi", assay = "Spatial", verbose = TRUE)
liver2SCT <- SCTransform(liver2, method = "glmGamPoi", assay = "Spatial", verbose = TRUE)
liver3SCT <- SCTransform(liver3, method = "glmGamPoi", assay = "Spatial", verbose = TRUE)
liver4SCT <- SCTransform(liver4, method = "glmGamPoi", assay = "Spatial", verbose = TRUE)

# for liver2: show counts and features after SCTransform
plot1 <- VlnPlot(liver2SCT, features = "nFeature_SCT", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver2SCT, features = "nFeature_SCT") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2tif(x = NULL, file='1-DataCheckLiv2_gene_SCT', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")

plot1 <- VlnPlot(liver2SCT, features = "nCount_SCT", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(liver2SCT, features = "nCount_SCT") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
graph2svg(x = NULL, file='1-DataCheckLiv2_counts_SCT', font = "Arial", cairo = TRUE,   
          width = 6, height = 4.5, bg = "transparent")



## save the object at this point
saveRDS(liver1SCT, file = "liver1SCT.rds")
saveRDS(liver2SCT, file = "liver2SCT.rds")
saveRDS(liver3SCT, file = "liver3SCT.rds")
saveRDS(liver4SCT, file = "liver4SCT.rds")


## Load the data
liver1SCT <- readRDS("liver1SCT.rds")
liver2SCT <- readRDS("liver2SCT.rds")
liver3SCT <- readRDS("liver3SCT.rds")
liver4SCT <- readRDS("liver4SCT.rds")


## Dimension reduction, clustering, and visualization
# run dimension reduction and clustering on the RNA expression data

# liver1
obj=liver1SCT
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE, resolution = 1) #resolution between 0.4-1.2 typically returns good results
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
# Look at cluster IDs of the first 5 cells
head(Idents(obj), 5)

liver1SCT1= obj

# visualize the results of the clustering
p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 3)
p1 + p2

graph2tif(x = NULL, file='2-cluster_liv1', font = "Arial", cairo = TRUE,   
          width = 10, height = 6, bg = "transparent")


# liver2
obj=liver2SCT
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)

# Clustering and UMAP
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE, resolution = 1) #resolution between 0.4-1.2 typically returns good results
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
head(Idents(obj), 5)

liver2SCT1= obj

# visualize the results of the clustering
p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 3)
p1 + p2

graph2tif(x = NULL, file='2-cluster_liv2', font = "Arial", cairo = TRUE,   
          width = 10, height = 6, bg = "transparent")



# liver3
obj=liver3SCT
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE, resolution = 1) #resolution between 0.4-1.2 typically returns good results
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
head(Idents(obj), 5)

liver3SCT1= obj

# visualize the results of the clustering
p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 3)
p1 + p2

graph2tif(x = NULL, file='2-cluster_liv3', font = "Arial", cairo = TRUE,   
          width = 10, height = 6, bg = "transparent")



# liver4
obj=liver4SCT
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE, resolution = 1) #resolution between 0.4-1.2 typically returns good results
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
head(Idents(obj), 5)

liver4SCT1= obj

# visualize the results of the clustering
p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 3)
p1 + p2

graph2tif(x = NULL, file='2-cluster_liv4', font = "Arial", cairo = TRUE,   
          width = 10, height = 6, bg = "transparent")




## DimHeatmap(): easy exploration of the primary sources of heterogeneity in a dataset,
# useful when trying to decide which PCs to include for further downstream analyses
# liver1
DimHeatmap(liver1SCT1, dims = 1:6, cells = 500, balanced = TRUE)
graph2svg(x = NULL, file='3-DimHeatmap_liv1', font = "Arial", cairo = TRUE,   
          width = 8, height = 7, bg = "transparent")

# liver2
DimHeatmap(liver2SCT1, dims = 1:6, cells = 500, balanced = TRUE)
graph2svg(x = NULL, file='3-DimHeatmap_liv2', font = "Arial", cairo = TRUE,   
          width = 8, height = 7, bg = "transparent")

# liver3
DimHeatmap(liver3SCT1, dims = 1:6, cells = 500, balanced = TRUE)
graph2svg(x = NULL, file='3-DimHeatmap_liv3', font = "Arial", cairo = TRUE,   
          width = 8, height = 7, bg = "transparent")

# liver4
DimHeatmap(liver4SCT1, dims = 1:6, cells = 500, balanced = TRUE)
graph2svg(x = NULL, file='3-DimHeatmap_liv4', font = "Arial", cairo = TRUE,   
          width = 8, height = 7, bg = "transparent")



## save the object at this point                     
saveRDS(liver1SCT1, file = "Liver1_AftCluster.rds")
saveRDS(liver2SCT1, file = "Liver2_AftCluster.rds")
saveRDS(liver3SCT1, file = "Liver3_AftCluster.rds")
saveRDS(liver4SCT1, file = "Liver4_AftCluster.rds")


## Load the data
liver1SCT1 <- readRDS("Liver1_AftCluster.rds")
liver2SCT1 <- readRDS("Liver2_AftCluster.rds")
liver3SCT1 <- readRDS("Liver3_AftCluster.rds")
liver4SCT1 <- readRDS("Liver4_AftCluster.rds")





########### Integration: In order to work with multiple slices (with Merge and integration), and also compare difference btw diff conditions
## Load the data
liver1SCT <- readRDS("liver1SCT.rds")
liver2SCT <- readRDS("liver2SCT.rds")
liver3SCT <- readRDS("liver3SCT.rds")
liver4SCT <- readRDS("liver4SCT.rds")

### we provide the merge function.
liver.merge <- merge(liver1SCT, y=c(liver2SCT, liver3SCT, liver4SCT), add.cell.ids = c("CR1", "GDF15a", "CR2", "GDF15b"), project = "MASHpaper")
levels(liver.merge)

table(liver.merge$orig.ident)
liver.merge@active.ident
head(liver.merge[[]])
tail(liver.merge[[]])


# this function uses minimum of the median UMI (calculated using the raw UMI counts) of individual objects
# to reverse the individual SCT regression model using minimum of median UMI as the sequencing depth covariate. 
liver.merge=PrepSCTFindMarkers(liver.merge, assay = "SCT", verbose = TRUE)


########## integrate data from the two conditions (control and treatment)
# When aligning two genome sequences together,
# identification of shared/homologous regions can help to interpret differences
# between the sequences as well.
liver.mergeIntegr = liver.merge
liver.mergeIntegr
liver.mergeIntegr[["SCT"]]

# run standard analysis workflow
DefaultAssay(liver.mergeIntegr) <- "SCT"
VariableFeatures(liver.mergeIntegr) <- c(VariableFeatures(liver1SCT), VariableFeatures(liver2SCT),
                                         VariableFeatures(liver3SCT), VariableFeatures(liver4SCT))
#VariableFeatures(liver.mergeIntegr) <- c(VariableFeatures(liver1SCT1), VariableFeatures(liver2SCT1))

liver.mergeIntegr <- RunPCA(liver.mergeIntegr, npcs = 30, verbose = FALSE)

# integration (HarmonyIntegration)
liver.mergeIntegr <- IntegrateLayers(object = liver.mergeIntegr, method = HarmonyIntegration, orig.reduction = "pca",
                                     normalization.method = "SCT", new.reduction = "harmony", verbose = T)

# we can now visualize and cluster the datasets.
liver.mergeIntegr <- FindNeighbors(liver.mergeIntegr, reduction = "harmony", dims = 1:30)
liver.mergeIntegr <- FindClusters(liver.mergeIntegr, verbose = FALSE, resolution = 1.2, cluster.name = "harmony_clusters") # 17 cluster
liver.mergeIntegr <- RunUMAP(liver.mergeIntegr, reduction = "harmony",
                             dims = 1:30, reduction.name = "umap.harmony")

### save file
saveRDS(liver.mergeIntegr, file = "liver_merge_integration.rds")

## Load the data
liver.mergeIntegr <- readRDS("liver_merge_integration.rds")


# Visualization
DimPlot(liver.mergeIntegr, reduction = "umap.harmony", label = TRUE, group.by = c("ident", "group"))
graph2tif(x = NULL, file='17-clusterMergeIntegr_DimPlot', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")
#without labeling
DimPlot(liver.mergeIntegr, reduction = "umap.harmony", label = F, group.by = c("ident", "orig.ident"))
graph2tif(x = NULL, file='17_1-clusterMergeIntegr_DimPlotwoLabel', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")


SpatialDimPlot(liver.mergeIntegr, label = TRUE, label.size = 3)
graph2tif(x = NULL, file='17_2-clusterMergeIntegr_SpatialDimPlot', font = "Arial", cairo = TRUE,   
          width = 15, height = 6, bg = "transparent")
#without labeling
SpatialDimPlot(liver.mergeIntegr, label = F, label.size = 3)
graph2tif(x = NULL, file='17_2-clusterMergeIntegr_SpatialDimPlotwoLabel', font = "Arial", cairo = TRUE,   
          width = 15, height = 6, bg = "transparent")





#### Finding differentially expressed features (cluster biomarkers)--top15 changed genes (p_val_adj<0.001) in each cluster

# this function uses minimum of the median UMI (calculated using the raw UMI counts) of individual objects
# to reverse the individual SCT regression model using minimum of median UMI as the sequencing depth covariate. 
obj = liver.mergeIntegr
obj=PrepSCTFindMarkers(obj, assay = "SCT", verbose = TRUE)


# find markers for every cluster compared to all remaining cells,
# report only the positive ones
obj.markers <- FindAllMarkers(obj, only.pos = TRUE)

### Find top15 changed genes (p_val_adj<0.001) in each cluster
top15 <- obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.001) %>%
  arrange(desc(abs(avg_log2FC))) %>%  # arrange abs(log2FC)
  slice_head(n = 15) %>%
  ungroup()

write.csv(top15, 'Marker15forClusters.csv')




######## Annotation ########
liver.mergeIntegrAnno=obj

celltype <- c("0:Periportal hepatocytes (inflammatory phenotype)", "1: Periportal hepatocytes (steady-state)", "2: Pericentral hepatocytes (lipid accumulation)", "3: Pericentral hepatocytes (amino acid metabolism)", "4: Pericentral hepatocytes (stress-responsive traits)",
              "5: Pericentral hepatocytes (Wnt/β-catenin signaling)","6: Kupffer cells (activated, lipid-handling)", "7: Cholangiocytes", 
              "8: Hepatocytes (high metabolic, mixed-zone)", "9: Inflammatory macrophages (lipid or fibrosis-associated)", "10: Cholangiocytes/periportal epithelial cells",
              "11: Erythrocytes", "12: Cholangiocytes (progenitor-like or reactive phenotype)","13: Pro-inflammatory macrophages",
              "14: Periportal hepatocytes (TGF-beta signaling)", "15: Plasma B cells", "16: Pericentral hepatocytes (lipid metabolism)",
              "17: Cholangiocytes, HSCs, or ductal progenitors (ECM remodeling)", "18: IFN-activated macrophages")


names(celltype) <- levels(liver.mergeIntegrAnno)
liver.mergeIntegrAnno <- RenameIdents(liver.mergeIntegrAnno, celltype)

Idents(liver.mergeIntegrAnno)

liver.mergeIntegrAnno$celltype <- Idents(liver.mergeIntegrAnno)
head(liver.mergeIntegrAnno@meta.data)

# Visualize cell type annotation on UMAP
DimPlot(liver.mergeIntegrAnno)
graph2tif(x = NULL, file='18-Integrat_DimPlot_labeling', font = "Arial", cairo = TRUE,   
          width = 8, height = 5, bg = "transparent")

head(liver.mergeIntegrAnno[[]])

### save file
saveRDS(liver.mergeIntegrAnno, file = "liver_merge_integr_Ann.rds")

## Load the data
liver.mergeIntegrAnno <- readRDS("liver_merge_integr_Ann.rds")








### show cell type location in SpatialDimPlot
library(Seurat)
library(patchwork)
library(stringr)  
library(export)   

# export folder
outdir <- "spatialDimPlot_celltype"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Sample, groups
samples <- levels(factor(liver.mergeIntegrAnno$orig.ident))

# celltype list
celltypes <- sort(unique(liver.mergeIntegrAnno$celltype))

# Generate and save in a loop
for (ct in celltypes) {
    plots <- lapply(samples, function(sp) {
    obj_sp <- subset(liver.mergeIntegrAnno, orig.ident == sp)
    cells_use <- colnames(obj_sp)[obj_sp$celltype == ct]
    
    p <- if (length(cells_use) > 0) {
      SpatialDimPlot(
        obj_sp,
        cells.highlight = list(cells_use),
        cols = c("grey90", "red"),
        pt.size.factor = 1.6
      )
    } else {
      SpatialDimPlot(obj_sp, pt.size.factor = 1.6)
    }
    
    p + ggtitle(sp) + theme(legend.position = "none")
  })
  
  
  panel <- wrap_plots(plots, ncol = length(samples)) +
    plot_annotation(title = ct)
  
  
  file_safe <- str_replace_all(ct, "[^A-Za-z0-9]+", "_")
  outfile   <- file.path(outdir, paste0(file_safe, ".tif"))
  
  
  graph2tif(
    x = panel,
    file = outfile,
    width = 15, height = 6, dpi = 300,
    font = "Arial", cairo = TRUE, bg = "transparent"
  )
}






### Delete Cluster 16, because it is located around the edge of slides
to_drop <- "16: Pericentral hepatocytes (lipid metabolism)"

obj = liver.mergeIntegrAnno
obj2 <- subset(obj, subset = celltype != to_drop)

obj2$celltype <- droplevels(factor(obj2$celltype))
liver.mergeIntegrAnno1 = obj2

### save file
saveRDS(liver.mergeIntegrAnno1, file = "liver_merge_integr_Ann1.rds")




####### SpatialDimPlot by groups
## extract meta data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")

# SpatialDimPlot
head(liver.mergeIntegrAnno1@meta.data)
obj <- liver.mergeIntegrAnno1

# 1) make sure identities are what you want to label (celltype)
Idents(obj) <- "seurat_clusters"   # or "seurat_clusters" if you prefer

DimPlot(obj, reduction = "umap.harmony", label = F, group.by = c("ident", "group"))
graph2tif(x = NULL, file='17-clusterMergeIntegr_DimPlot_Group', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")
graph2svg(x = NULL, file='17-clusterMergeIntegr_DimPlot_Group', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")

DimPlot(obj, reduction = "umap.harmony", label = T, group.by = c("ident", "group"))
graph2tif(x = NULL, file='17-clusterMergeIntegr_DimPlot_Label', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")
graph2svg(x = NULL, file='17-clusterMergeIntegr_DimPlot_Label', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")


## only show cluster 9, 13, 15
obj <- liver.mergeIntegrAnno1
Idents(obj) <- "seurat_clusters"
clusters_keep <- "15"
DimPlot(
  obj,
  reduction = "umap.harmony",
  group.by = c("ident", "group"),
  label = FALSE,
  cells = WhichCells(obj, idents = clusters_keep)
)





Idents(obj) <- "celltype"   # or "seurat_clusters" if you prefer
Idents(obj) <- "seurat_clusters"
Idents(obj) <- factor(Idents(obj), levels = as.character(0:18))
DimPlot(obj, reduction = "umap.harmony", label = T, group.by = "celltype")
graph2tif(x = NULL, file='17-clusterMergeIntegr_DimPlot_Label', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")

Idents(obj) <- "seurat_clusters"   # or "seurat_clusters" if you prefer
SpatialDimPlot(obj, label = TRUE, label.size = 3)
graph2tif(x = NULL, file='17_2-clusterMergeIntegr_SpatialDimPlot1', font = "Arial", cairo = TRUE,   
          width = 15, height = 6, bg = "transparent")
graph2svg(x = NULL, file='17_2-clusterMergeIntegr_SpatialDimPlot1', font = "Arial", cairo = TRUE,   
          width = 15, height = 6, bg = "transparent")
graph2pdf(x = NULL, file='17_2-clusterMergeIntegr_SpatialDimPlot1', font = "Arial", cairo = TRUE,   
          width = 15, height = 6, bg = "transparent")

p=SpatialDimPlot(obj, label = TRUE, label.size = 3)
ggsave("17_2-clusterMergeIntegr_SpatialDimPlot2.svg", p, width = 15, height = 6, device = svglite::svglite)


## get order labeling
library(Seurat)
library(stringr)

obj <- liver.mergeIntegrAnno1   # or your object name

## 1) normalize "N:" formatting: ensure exactly one space after colon
obj$celltype <- as.character(obj$celltype)
obj$celltype <- str_replace(obj$celltype, "^\\s*(\\d+):\\s*", "\\1: ")

## 2) sort levels by numeric prefix (robust to spaces)
levs <- sort(
  unique(obj$celltype),
  index.return = FALSE,
  method = "auto"
)

# order by number before the colon
ord_idx <- order(as.integer(str_match(levs, "^(\\d+):")[,2]))
levs_sorted <- levs[ord_idx]

# keep only those actually present (defensive) and in 0..18 order
levs_sorted <- levs_sorted[order(as.integer(str_match(levs_sorted, "^(\\d+):")[,2]))]

# apply levels
obj$celltype <- factor(obj$celltype, levels = levs_sorted)

## 3) make identities = celltype so labels use this order
Idents(obj) <- obj$celltype

## 4) plot (labels + legend now 0 → 18)
p <- DimPlot(
  obj,
  reduction  = "umap.harmony",
  group.by   = "celltype",
  label      = TRUE,
  repel      = TRUE,
  label.size = 3   # size of labels on the UMAP itself
) +
  theme_bw() +
  theme(
    legend.text  = element_text(size = 12),  # increase legend text
    legend.title = element_text(size = 14, face = "bold")  # optional, bigger/bold title
  )

p
graph2tif(x = NULL, file='17-DimPlot_Labeling', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")
graph2svg(x = NULL, file='17-DimPlot_Labeling', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")





########### count how many cells in different cell types #################
## extract meta data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")

# SpatialDimPlot
head(liver.mergeIntegrAnno1@meta.data)
DimPlot(liver.mergeIntegrAnno1, label = TRUE, label.size = 3)
graph2tif(x = NULL, file='17_2-clusterMergeIntegr_SpatialDimPlot', font = "Arial", cairo = TRUE,   
          width = 15, height = 6, bg = "transparent")


liver.merge1 = liver.mergeIntegrAnno1
head(liver.merge1@meta.data)
library(data.table)
md = liver.merge1@meta.data %>% as.data.table()

## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
head(md)
md_count = md[, .N, by = c("group", "celltype")]
md_count = as.data.frame(md_count)
head(md_count)
table(md_count$group)
md_count = md_count[order(md_count$celltype),]
write.csv(md_count, '2_celltype_num2.csv')

# load packages
library(data.table)
library(dplyr)
library(ggplot2)

# counts per group × celltype
md <- as.data.table(liver.merge1@meta.data)
md_count <- md[, .N, by = .(group, celltype)] %>% as.data.frame()
head(md_count)

# df_pct has: celltype, group (CR/GDF15), N, Percent
# make sure the stacking order is what you want:
#   left = GDF15, right = CR  --> levels = c("GDF15","CR")
#   left = CR,    right = GDF15 --> levels = c("CR","GDF15")
df_pct <- md_count %>%
  mutate(group = factor(group, levels = c("CR","GDF15")))

head(df_pct)



# First calculate the percentage (Percent) within each group
library(dplyr)
library(tidyr)

df_pct <- md_count %>%
  group_by(celltype) %>%                     
  mutate(Percent = N / sum(N) * 100) %>%     
  ungroup() %>%
  mutate(group = factor(group, levels = c("CR","GDF15")))

head(df_pct)
write.csv(df_pct, '2_celltype_perctage.csv')


# centers for each segment
centers <- df_pct %>%
  select(celltype, group, Percent) %>%
  tidyr::pivot_wider(names_from = group, values_from = Percent) %>%
  mutate(x_GDF15 = `GDF15`/2,          # middle of teal (left) segment
         x_CR    = 100 - `CR`/2) %>%   # middle of red (right) segment
  select(celltype, x_GDF15, x_CR)

df_plot <- df_pct %>%
  left_join(centers, by = "celltype") %>%
  mutate(
    # use the center of THIS segment (no swap)
    pos = ifelse(group == "GDF15", x_GDF15, x_CR)
  )


library(dplyr)
library(tidyr)
library(forcats)

## Order cell types by their percentages in the GDF15 group
order_ct <- df_plot %>% 
  filter(group == "GDF15") %>%       
  arrange(Percent) %>%                
  pull(celltype)

df_plot <- df_plot %>% 
  mutate(celltype = factor(celltype, levels = order_ct))
head(df_plot)
table(df_plot$celltype)



## 2. Plot
ggplot(df_plot, aes(x = Percent, y = celltype, fill = group)) +
  geom_col(width = 0.9) +
  geom_text(aes(x = pos, label = sprintf("%.1f%%", Percent)),
            color = "black", size = 3.5) +
  geom_vline(xintercept = 50, linetype = "dotted", color = "black") +
  scale_x_continuous(breaks = seq(0, 100, 10),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "Percent", y = "Cell Type") +
  theme_bw(base_size = 14) +
  theme(
    legend.title = element_blank(),
    axis.text.x  = element_text(size = 12, hjust = 0.5),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.margin  = margin(10, 10, 10, 10)
  )

graph2tif(x = NULL, file='11_CellTypePerc2', font = "Arial", cairo = TRUE,   
          width = 12, height = 5.2, bg = "transparent")
graph2svg(x = NULL, file='11_CellTypePerc2', font = "Arial", cairo = TRUE,   
          width = 12, height = 5.2, bg = "transparent")


## make figures separately using the parenchymal cells and non-parenchymal cells
#non-parenchymal cells
library(dplyr)
library(ggplot2)
cells_keep <- c(
  "9: Activated HSCs/fibroblasts",
  "13: Pro-inflammatory macrophages",
  "15: Plasma B cells",
  "17: Cholangiocytes, HSCs, or ductal progenitors (ECM remodeling)",
  "18: IFN-activated macrophages",
  "6: Kupffer cells (activated, lipid-handling)"
)

df_plot_sub <- df_plot %>%
  filter(celltype %in% cells_keep) %>%

  mutate(celltype = factor(celltype, levels = order_ct))

## 2. Plot fig
ggplot(df_plot_sub, aes(x = Percent, y = celltype, fill = group)) +
  geom_col(width = 0.9) +
  geom_text(aes(x = pos, label = sprintf("%.1f%%", Percent)),
            color = "black", size = 3.5) +
  geom_vline(xintercept = 50, linetype = "dotted", color = "black") +
  scale_x_continuous(breaks = seq(0, 100, 10),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "Percent", y = "Cell Type") +
  theme_bw(base_size = 14) +
  theme(
    legend.title = element_blank(),
    axis.text.x  = element_text(size = 12, hjust = 0.5),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.margin  = margin(10, 10, 10, 10)
  )

graph2tif(x = NULL, file='11_CellTypePerc_immune', font = "Arial", cairo = TRUE,   
          width = 12, height = 3, bg = "transparent")
graph2svg(x = NULL, file='11_CellTypePerc_immune', font = "Arial", cairo = TRUE,   
          width = 12, height = 3, bg = "transparent")



#parenchymal cells
library(dplyr)
library(ggplot2)

cells_exclude <- c(
  "9: Activated HSCs/fibroblasts",
  "13: Pro-inflammatory macrophages",
  "15: Plasma B cells",
  "17: Cholangiocytes, HSCs, or ductal progenitors (ECM remodeling)",
  "18: IFN-activated macrophages",
  "6: Kupffer cells (activated, lipid-handling)"
)

df_plot_sub <- df_plot %>%
  filter(!celltype %in% cells_exclude) %>%
  
  mutate(celltype = factor(celltype, levels = order_ct))

## 2. Plot fig
ggplot(df_plot_sub, aes(x = Percent, y = celltype, fill = group)) +
  geom_col(width = 0.9) +
  geom_text(aes(x = pos, label = sprintf("%.1f%%", Percent)),
            color = "black", size = 3.5) +
  geom_vline(xintercept = 50, linetype = "dotted", color = "black") +
  scale_x_continuous(breaks = seq(0, 100, 10),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(x = "Percent", y = "Cell Type") +
  theme_bw(base_size = 14) +
  theme(
    legend.title = element_blank(),
    axis.text.x  = element_text(size = 12, hjust = 0.5),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.margin  = margin(10, 10, 10, 10)
  )

graph2tif(x = NULL, file='11_CellTypePerc_parenchymal', font = "Arial", cairo = TRUE,   
          width = 12, height = 4, bg = "transparent")
graph2svg(x = NULL, file='11_CellTypePerc_parenchymal', font = "Arial", cairo = TRUE,   
          width = 12, height = 4, bg = "transparent")




#### Add one column "group" in metadata
obj = liver.mergeIntegrAnno1

obj$group <- case_when(
  obj$orig.ident %in% c("CR1","CR2") ~ "CR",
  obj$orig.ident %in% c("GDF15a","GDF15b") ~ "GDF15",
  TRUE ~ obj$orig.ident   # keep original if not in list
)
head(obj@meta.data[, c("orig.ident","group")])
liver.mergeIntegrAnno1 = obj
### save file
saveRDS(liver.mergeIntegrAnno1, file = "liver_merge_integr_Ann1.rds")






## Show ident "13: Pro-inflammatory macrophages", "9: Activated HSCs/fibroblasts", "15: Plasma B cells" in spatial
## load data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
obj = liver.mergeIntegrAnno1

keep <- c("13: Pro-inflammatory macrophages",
          "9: Activated HSCs/fibroblasts",
          "15: Plasma B cells")

# subset and re-label
obj_sub <- subset(obj, subset = celltype %in% keep)

Idents(obj_sub) <- "celltype"
# Check the levels
levels(Idents(obj_sub))

# 1) Make sure identities are these three and in the order you want in the legend
lvl <- c("9: Activated HSCs/fibroblasts",
         "13: Pro-inflammatory macrophages",
         "15: Plasma B cells")
obj_sub$celltype <- factor(obj_sub$celltype, levels = lvl)
Idents(obj_sub) <- "celltype"

# 2) Plot with explicitly named colors (avoids the "Insufficient values" error)
cols3 <- c("9: Activated HSCs/fibroblasts"    = "#009E73",
           "13: Pro-inflammatory macrophages" = "#FF4D4D",
           "15: Plasma B cells"               = "#0072B2")

SpatialDimPlot(
  obj_sub,
  cols = cols3,              # named mapping
  pt.size.factor = 3,
  label = F,
  label.size = 5
) + theme(legend.position = "bottom")



graph2tif(x = NULL, file='17_SpaDimPlot_HSC_Infla_B', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")

graph2svg(x = NULL, file='17_SpaDimPlot_HSC_Infla_B', font = "Arial", cairo = TRUE,   
          width = 15, height = 5, bg = "transparent")


p=SpatialDimPlot(
  obj_sub,
  cols = cols3,              # named mapping
  pt.size.factor = 3,
  label = F,
  label.size = 5
) + theme(legend.position = "bottom")
ggsave("17_SpaDimPlot_HSC_Infla_B.svg", p, width = 15, height = 6, device = svglite::svglite)




#### Finding differentially expressed features (cluster biomarkers)
## load data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
obj = liver.mergeIntegrAnno1

# find markers for every cluster compared to all remaining cells,
# report only the positive ones
obj.markers <- FindAllMarkers(obj, only.pos = TRUE)
# Filter by fold change > 1 if desired
obj.markers <- obj.markers %>% dplyr::filter(avg_log2FC > 1)

# Select top 10 markers per cluster
top10 <- obj.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

DefaultAssay(obj) <- "Spatial"
features_use <- intersect(unique(top10$gene), rownames(obj[["Spatial"]]))

obj <- NormalizeData(obj, assay = "Spatial")
obj <- ScaleData(obj, assay = "Spatial", features = features_use)

DoHeatmap(obj, features = top10$gene) + NoLegend()

# save it in png file









##### Plot Markers for "13: Pro-inflammatory macrophages", "9: Activated HSCs/fibroblasts", "15: Plasma B cells"
##### Add Kupffer cells

## load data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
obj = liver.mergeIntegrAnno1
Idents(obj) <- "celltype"
table(Idents(obj))

## Find markers for "9: Activated HSCs/fibroblasts"
# 1) Switch to RNA (log-normalized) assay
DefaultAssay(obj) <- "Spatial"

# If needed (do once)
obj <- NormalizeData(obj)            # LogNormalize
obj <- FindVariableFeatures(obj)     # optional but good practice
# obj <- ScaleData(obj)              # not required for FindMarkers

# 2) Set the identities to your celltype column
Idents(obj) <- "celltype"

# 3) Run DE for the cluster of interest
demarkerHSC <- FindMarkers(
  obj, ident.1 = "9: Activated HSCs/fibroblasts",
  test.use = "wilcox",
  min.pct = 0.10,          # seen in ≥10% of cells in ident.1
  min.diff.pct = 0.05,     # +5% prevalence over the rest
  logfc.threshold = 0.25,  # avoid tiny effect sizes
  verbose = FALSE
)

# 4) Filter & sort as you wanted
library(dplyr)
top_genes <- demarkerHSC %>%
  filter(p_val < 0.01) %>%
  arrange(desc(avg_log2FC))
head(top_genes, 20)

write.csv(demarkerHSC, 'DE_Top15_9HSC.csv')



## Find markers for "13: Pro-inflammatory macrophages"
# 1) Switch to RNA (log-normalized) assay
DefaultAssay(obj) <- "Spatial"

# If needed (do once)
obj <- NormalizeData(obj)            # LogNormalize
obj <- FindVariableFeatures(obj)     # optional but good practice
# obj <- ScaleData(obj)              # not required for FindMarkers

# 2) Set the identities to your celltype column
Idents(obj) <- "celltype"

# 3) Run DE for the cluster of interest
demarkerHSC <- FindMarkers(
  obj, ident.1 = "13: Pro-inflammatory macrophages",
  test.use = "wilcox",
  min.pct = 0.10,          # seen in ≥10% of cells in ident.1
  min.diff.pct = 0.05,     # +5% prevalence over the rest
  logfc.threshold = 0.25,  # avoid tiny effect sizes
  verbose = FALSE
)


# 4) Filter & sort as you wanted
library(dplyr)
top_genes <- demarkerHSC %>%
  filter(p_val < 0.01) %>%
  arrange(desc(avg_log2FC))
head(top_genes, 20)

write.csv(demarkerHSC, 'DE_Top15_13ProInfl.csv')



## Find markers for "15: Plasma B cells"
# 1) Switch to RNA (log-normalized) assay
DefaultAssay(obj) <- "Spatial"

# If needed (do once)
obj <- NormalizeData(obj)            # LogNormalize
obj <- FindVariableFeatures(obj)     # optional but good practice
# obj <- ScaleData(obj)              # not required for FindMarkers

# 2) Set the identities to your celltype column
Idents(obj) <- "celltype"

# 3) Run DE for the cluster of interest

demarkerHSC <- FindMarkers(
  obj, ident.1 = "15: Plasma B cells",
  test.use = "wilcox",
  min.pct = 0.10,          # seen in ≥10% of cells in ident.1
  min.diff.pct = 0.05,     # +5% prevalence over the rest
  logfc.threshold = 0.25,  # avoid tiny effect sizes
  verbose = FALSE
)



# 4) Filter & sort as you wanted
library(dplyr)
top_genes <- demarkerHSC %>%
  filter(p_val < 0.01) %>%
  arrange(desc(avg_log2FC))
head(top_genes, 15)

write.csv(demarkerHSC, 'DE_Top15_15Bcell.csv')




## Find markers for "6: Kupffer cells (activated, lipid-handling)"
# 1) Switch to RNA (log-normalized) assay
DefaultAssay(obj) <- "Spatial"

# If needed (do once)
obj <- NormalizeData(obj)            # LogNormalize
obj <- FindVariableFeatures(obj)     # optional but good practice
# obj <- ScaleData(obj)              # not required for FindMarkers

# 2) Set the identities to your celltype column
Idents(obj) <- "celltype"

# 3) Run DE for the cluster of interest

demarkerHSC <- FindMarkers(
  obj, ident.1 = "6: Kupffer cells (activated, lipid-handling)",
  test.use = "wilcox",
  min.pct = 0.10,          # seen in ≥10% of cells in ident.1
  min.diff.pct = 0.05,     # +5% prevalence over the rest
  logfc.threshold = 0.25,  # avoid tiny effect sizes
  verbose = FALSE
)



# 4) Filter & sort as you wanted
library(dplyr)
top_genes <- demarkerHSC %>%
  filter(p_val < 0.01) %>%
  arrange(desc(avg_log2FC))
head(top_genes, 15)

write.csv(top_genes, 'DE_Top15_6Kupffer.csv')


HSC <- subset(obj, idents = c("9: Activated HSCs/fibroblasts", "13: Pro-inflammatory macrophages", "15: Plasma B cells"))
# HSC <- subset(obj, idents = c("6: Kupffer cells (activated, lipid-handling)",
#                               "9: Activated HSCs/fibroblasts", "13: Pro-inflammatory macrophages", "15: Plasma B cells"))

VlnPlot(
  HSC,
  features = c(#"Cd5l", "C1qc","Clec4f", "Timd4", "Marco", "Vsig4" ,
                "Ly6d","Cln6","Cidec","Col1a1","Gbp2",
               "Cxcl9","Cxcl10","Il1rn",
               "Ighj4","Ighj1","Jchain","Igkc","Ighm"),
  stack = TRUE,
  flip = TRUE,
  split.by = "group"
) + 
  theme(axis.text.x = element_text(angle = 70,    # change 45 → 90 for vertical
                                   hjust = 1, 
                                   vjust = 1))

graph2tif(x = NULL, file='17_HSC_B_InflExp', font = "Arial", cairo = TRUE,   
          width = 4, height = 7, bg = "transparent")
graph2svg(x = NULL, file='17_HSC_B_InflExp', font = "Arial", cairo = TRUE,   
          width = 4, height = 7, bg = "transparent")







### only showing the "13: Pro-inflammatory macrophages" in violin plot and color the M1 gene marker

## load data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
obj = liver.mergeIntegrAnno1
Idents(obj) <- "celltype"
table(Idents(obj))

## Find markers for "13: Pro-inflammatory macrophages"
# 1) Switch to RNA (log-normalized) assay
DefaultAssay(obj) <- "Spatial"

# If needed (do once)
obj <- NormalizeData(obj)            # LogNormalize
obj <- FindVariableFeatures(obj)     # optional but good practice
# obj <- ScaleData(obj)              # not required for FindMarkers

# 2) Set the identities to your celltype column
Idents(obj) <- "celltype"

# 3) Run DE for the cluster of interest
demarkerHSC <- FindMarkers(
  obj, ident.1 = "13: Pro-inflammatory macrophages",
  test.use = "wilcox",
  min.pct = 0.10,          # seen in ≥10% of cells in ident.1
  min.diff.pct = 0.05,     # >5% prevalence over the rest
  logfc.threshold = 0.25,  # avoid tiny effect sizes
  verbose = FALSE
)


# 4) Filter & sort as you wanted
library(dplyr)
top_genes <- demarkerHSC %>%
  filter(p_val < 0.01) %>%
  arrange(desc(avg_log2FC))
head(top_genes, 20)

write.csv(top_genes, 'DE_Top15_13ProInfl.csv')


## only show one cluster in UMAP
# ## only show cluster 13
# obj <- liver.mergeIntegrAnno1
# Idents(obj) <- "seurat_clusters"
# clusters_keep <- "13"
# DimPlot(
#   obj,
#   reduction = "umap.harmony",
#   group.by = c("ident", "group"),
#   label = FALSE,
#   cells = WhichCells(obj, idents = clusters_keep)
# )

#TBC
#FeaturePlot(obj, features = c("Ccl4", "Cxcl9", "Cxcl10", "Nos2", "Il1b", "Tnf", "Il6", "Ccr2"))

HSC <- subset(obj, idents = "13: Pro-inflammatory macrophages")
# HSC <- subset(obj, idents = c("6: Kupffer cells (activated, lipid-handling)",
#                               "9: Activated HSCs/fibroblasts", "13: Pro-inflammatory macrophages", "15: Plasma B cells"))

VlnPlot(
  HSC,
  features = c("Cxcl9","Cxcl10","Mmp12", "Gbp2",#for M1
               "Cd163",  "Retnla",
               "Il10", "Clec10a"),
  stack = TRUE,
  flip = TRUE,
  split.by = "group"
) + 
  theme(axis.text.x = element_text(angle = 70,    # change 45 → 90 for vertical
                                   hjust = 1, 
                                   vjust = 1))

graph2tif(x = NULL, file='17_InflExp_M1_M2', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 6, bg = "transparent")
graph2svg(x = NULL, file='17_InflExp_M1_M2', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 6, bg = "transparent")










############# Adrenergic receptor expression in the liver
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1

DotPlot(obj, features = c("Adrb1", "Adrb2", "Adrb3"),split.by = "group",assay = "SCT",scale = T)


# 1) Get DotPlot data (aggregated by celltype × gene × group)
p <- DotPlot(
  obj,
  features = c("Adrb1", "Adrb2", "Adrb3"),
  assay = "SCT",
  group.by = "celltype",
  split.by = "group",
  scale = FALSE
)

df <- p$data
# For older Seurat, split label is inside `id` -> pull it out:
if (!"split" %in% names(df)) {
  df$split    <- sub(".*_(CR|GDF15)$", "\\1", df$id)
  df$celltype <- sub("_(CR|GDF15)$",  "", df$id)
} else {
  df$celltype <- df$id
}
df$split <- factor(df$split, levels = c("CR","GDF15"))

library(ggplot2)
library(grid)  # for unit()

# set celltype levels explicitly
df$celltype <- factor(df$celltype, levels = c(
  "0:Periportal hepatocytes (inflammatory phenotype)",
  "1: Periportal hepatocytes (steady-state)",
  "2: Pericentral hepatocytes (lipid accumulation)",
  "3: Pericentral hepatocytes (amino acid metabolism)",
  "4: Pericentral hepatocytes (stress-responsive traits)",
  "5: Pericentral hepatocytes (Wnt/β-catenin signaling)",
  "6: Kupffer cells (activated, lipid-handling)",
  "7: Cholangiocytes",
  "8: Hepatocytes (high metabolic, mixed-zone)",
  "9: Inflammatory macrophages (lipid or fibrosis-associated)",
  "10: Cholangiocytes/periportal epithelial cells",
  "11: Erythrocytes",
  "12: Cholangiocytes (progenitor-like or reactive phenotype)",
  "13: Pro-inflammatory macrophages",
  "14: Periportal hepatocytes (TGF-beta signaling)",
  "15: Plasma B cells",
  "17: Cholangiocytes, HSCs, or ductal progenitors (ECM remodeling)",
  "18: IFN-activated macrophages"
))

# redraw
ggplot(df, aes(x = features.plot, y = celltype)) +
  geom_point(aes(size = avg.exp, color = split), alpha = 0.9) +
  scale_size_continuous(name = "Avg exp (log1p)", range = c(0, 8)) +
  scale_color_manual(values = c(CR = "red", GDF15 = "blue"), name = NULL) +
  facet_grid(. ~ split) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9)
  )

graph2tif(x = NULL, file='Adrb1-3expr', font = "Arial", cairo = TRUE,   
          width = 8, height = 5, bg = "transparent")






############# GR-associated gene expression in different liver cell types
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1

DotPlot(obj, features = c("Adrb1", "Dusp1", "Zfp36",
                          "Pck1", "G6pc", "Fbp1", "Tat",
                          "Angptl4", "Sgk1", "Per1", "Klf15","Pdk4","Hes1"),split.by = "group",assay = "SCT",scale = T)


# 1) Get DotPlot data (aggregated by celltype × gene × group)
p <- DotPlot(
  obj,
  features = c("Adrb1", "Dusp1", "Zfp36",
               "Pck1", "G6pc", "Fbp1", "Tat",
               "Angptl4", "Sgk1", "Per1", "Klf15","Pdk4","Hes1"),
  assay = "SCT",
  group.by = "celltype",
  split.by = "group",
  scale = FALSE
)


df <- p$data
# For older Seurat, split label is inside `id` -> pull it out:
if (!"split" %in% names(df)) {
  df$split    <- sub(".*_(CR|GDF15)$", "\\1", df$id)
  df$celltype <- sub("_(CR|GDF15)$",  "", df$id)
} else {
  df$celltype <- df$id
}
df$split <- factor(df$split, levels = c("CR","GDF15"))

library(ggplot2)
library(grid)  # for unit()

# set celltype levels explicitly
df$celltype <- factor(df$celltype, levels = c(
  "0:Periportal hepatocytes (inflammatory phenotype)",
  "1: Periportal hepatocytes (steady-state)",
  "2: Pericentral hepatocytes (lipid accumulation)",
  "3: Pericentral hepatocytes (amino acid metabolism)",
  "4: Pericentral hepatocytes (stress-responsive traits)",
  "5: Pericentral hepatocytes (Wnt/β-catenin signaling)",
  "6: Kupffer cells (activated, lipid-handling)",
  "7: Cholangiocytes",
  "8: Hepatocytes (high metabolic, mixed-zone)",
  "9: Inflammatory macrophages (lipid or fibrosis-associated)",
  "10: Cholangiocytes/periportal epithelial cells",
  "11: Erythrocytes",
  "12: Cholangiocytes (progenitor-like or reactive phenotype)",
  "13: Pro-inflammatory macrophages",
  "14: Periportal hepatocytes (TGF-beta signaling)",
  "15: Plasma B cells",
  "17: Cholangiocytes, HSCs, or ductal progenitors (ECM remodeling)",
  "18: IFN-activated macrophages"
))



# redraw (gene panels; within each gene show CR vs GDF15)
ggplot(df, aes(x = split, y = celltype)) +
  geom_point(aes(size = avg.exp, color = split), alpha = 0.9) +
  scale_size_continuous(name = "Avg exp (log1p)", range = c(0, 8)) +
  scale_color_manual(values = c(CR = "red", GDF15 = "blue"), name = NULL) +
  facet_wrap(~ features.plot, nrow = 1) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9)
  )

graph2tif(x = NULL, file='GRsignaling_geneExpr_1', font = "Arial", cairo = TRUE,   
          width = 10, height = 5, bg = "transparent")


## Normalized the corresponding CR 
library(dplyr)

# Calculate the expression value for each gene × cell type in the CR group
df <- df %>%
  group_by(features.plot, celltype) %>%
  mutate(
    CR_exp = avg.exp[split == "CR"][1],
    avg.exp.norm = case_when(
      split == "CR" ~ 1,
      split == "GDF15" ~ avg.exp / CR_exp,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

ggplot(df, aes(x = split, y = celltype)) +
  geom_point(aes(size = avg.exp.norm, color = split), alpha = 0.9) +
  scale_size_continuous(name = "Relative to CR", range = c(0, 8)) +
  scale_color_manual(values = c(CR = "red", GDF15 = "blue"), name = NULL) +
  facet_wrap(~ features.plot, nrow = 1) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9)
  )

graph2tif(x = NULL, file='GRsignaling_geneExpr_allGene_2', font = "Arial", cairo = TRUE,   
          width = 12, height = 5, bg = "transparent")




## Plot figures seperately "infla" and "nonInfla"
# "9: Inflammatory macrophages (lipid or fibrosis-associated)",
# "13: Pro-inflammatory macrophages",
# "15: Plasma B cells",

#another gene list
p <- DotPlot(
  obj,
  features = c("Dusp1", "Zfp36",
               "Pck1", "G6pc", "Fbp1", "Tat",
               "Sgk1", "Per1", "Klf15","Hes1"),
  assay = "SCT",
  group.by = "celltype",
  split.by = "group",
  scale = FALSE
)

df <- p$data
# For older Seurat, split label is inside `id` -> pull it out:
if (!"split" %in% names(df)) {
  df$split    <- sub(".*_(CR|GDF15)$", "\\1", df$id)
  df$celltype <- sub("_(CR|GDF15)$",  "", df$id)
} else {
  df$celltype <- df$id
}
df$split <- factor(df$split, levels = c("CR","GDF15"))

library(ggplot2)
library(grid)  # for unit()

# set celltype levels explicitly
df$celltype <- factor(df$celltype, levels = c(
  "0:Periportal hepatocytes (inflammatory phenotype)",
  "1: Periportal hepatocytes (steady-state)",
  "2: Pericentral hepatocytes (lipid accumulation)",
  "3: Pericentral hepatocytes (amino acid metabolism)",
  "4: Pericentral hepatocytes (stress-responsive traits)",
  "5: Pericentral hepatocytes (Wnt/β-catenin signaling)",
  "6: Kupffer cells (activated, lipid-handling)",
  "7: Cholangiocytes",
  "8: Hepatocytes (high metabolic, mixed-zone)",
  "9: Inflammatory macrophages (lipid or fibrosis-associated)",
  "10: Cholangiocytes/periportal epithelial cells",
  "11: Erythrocytes",
  "12: Cholangiocytes (progenitor-like or reactive phenotype)",
  "13: Pro-inflammatory macrophages",
  "14: Periportal hepatocytes (TGF-beta signaling)",
  "15: Plasma B cells",
  "17: Cholangiocytes, HSCs, or ductal progenitors (ECM remodeling)",
  "18: IFN-activated macrophages"
))


# ---- split cell types into two datasets ----
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1

DotPlot(obj, features = c("Dusp1", "Zfp36",
                          "Pck1", "G6pc", "Fbp1", "Tat",
                          "Sgk1", "Per1", "Klf15","Hes1"),split.by = "group",assay = "SCT",scale = T)

keep_ct <- c(
  "15: Plasma B cells",
  "13: Pro-inflammatory macrophages",
  "6: Kupffer cells (activated, lipid-handling)"
)


df_keep  <- df %>% dplyr::filter(as.character(celltype) %in% keep_ct)
df_other <- df %>%
  dplyr::filter(
    !celltype %in% keep_ct,
    !is.na(celltype)
  )

# Calculate the expression value for each gene × cell type in the CR group
df_keep <- df_keep %>%
  group_by(features.plot, celltype) %>%
  mutate(
    CR_exp = avg.exp[split == "CR"][1],
    avg.exp.norm = case_when(
      split == "CR" ~ 1,
      split == "GDF15" ~ avg.exp / CR_exp,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

df_other <- df_other %>%
  group_by(features.plot, celltype) %>%
  mutate(
    CR_exp = avg.exp[split == "CR"][1],
    avg.exp.norm = case_when(
      split == "CR" ~ 1,
      split == "GDF15" ~ avg.exp / CR_exp,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup()

# ---- helper: same plot settings, unchanged ----
make_plot <- function(dat) {
  ggplot(dat, aes(x = split, y = celltype)) +
    geom_point(aes(size = avg.exp.norm, color = split), alpha = 0.9) +
    scale_size_continuous(name = "Relative to CR", range = c(0, 8)) +
    scale_color_manual(values = c(CR = "red", GDF15 = "blue"), name = NULL) +
    facet_wrap(~ features.plot, nrow = 1) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 9)
    )
}

p_keep  <- make_plot(df_keep)
p_other <- make_plot(df_other)

print(p_keep)
graph2tif(x = NULL, file='GRsignaling_geneExpr_Infla_3', font = "Arial", cairo = TRUE,   
          width = 8, height = 3, bg = "transparent")

graph2svg(x = NULL, file='GRsignaling_geneExpr_Infla_3', font = "Arial", cairo = TRUE,   
          width = 8, height = 3, bg = "transparent")

print(p_other)
graph2tif(x = NULL, file='GRsignaling_geneExpr_nonInfla_4', font = "Arial", cairo = TRUE,   
          width = 9.5, height = 5, bg = "transparent")

graph2svg(x = NULL, file='GRsignaling_geneExpr_nonInfla_4', font = "Arial", cairo = TRUE,   
          width = 9.5, height = 5, bg = "transparent")



############# Compare DEGs in liver cell types between CR vs GDF15 ###########
## Load the data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1
table(liver.mergeIntegrAnno1$celltype)

##### Find the DE btw CR and GDF15 in cluster 0
head(obj@meta.data)
tail(obj@meta.data)
obj$celltype.stim <- paste(Idents(obj), obj$group, sep = "_")
Idents(obj) <- "celltype.stim"
obj$treat = obj$orig.ident
### save file
liver.mergeIntegrAnno1 = obj 
saveRDS(liver.mergeIntegrAnno1, file = "liver_merge_integr_Ann1.rds")



#### for loop for DEG analysis and volcano plot:

## ---------- CONFIG ----------
assay_counts <- "Spatial"  # or "RNA"
map <- c(CR1="CR", CR2="CR", GDF15a="GDF15", GDF15b="GDF15")  
lfc_cut <- 0.5
fdr_cut <- 0.05
min_count <- 10
outdir_main <- "DE_edgeR_by_celltype"

safelabel <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

suppressPackageStartupMessages({
  library(edgeR); library(ggplot2); library(ggpubr); library(ggrepel); library(ggthemes)
})

## ---------- 0) Preprocessing: Generate Clean Cell-Type Names and Aggregate Counts ----------
obj$celltype_clean <- trimws(sub("^\\s*\\d+:\\s*", "", obj$celltype))

pb_list <- Seurat::AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident", "celltype_clean")
)
pb <- pb_list[[assay_counts]]   # genes x (sample_celltype) Matrix

# Parse column names into sample and cell type
cn <- colnames(pb)
us <- regexpr("_", cn)  # 第一处下划线
samples   <- substring(cn, 1, us - 1)
celltypes <- trimws(substring(cn, us + 1))
col_annot <- data.frame(col = cn, sample = samples, celltype = celltypes, stringsAsFactors = FALSE)

col_annot$group <- unname(map[col_annot$sample])
keep_cols <- !is.na(col_annot$group)
pb <- pb[, keep_cols, drop = FALSE]
col_annot <- col_annot[keep_cols, , drop = FALSE]

celltype_levels <- c(
  "Periportal hepatocytes (inflammatory phenotype)",
  "Periportal hepatocytes (steady-state)",
  "Pericentral hepatocytes (lipid accumulation)",
  "Pericentral hepatocytes (amino acid metabolism)",
  "Pericentral hepatocytes (stress-responsive traits)",
  "Pericentral hepatocytes (Wnt/β-catenin signaling)",
  "Kupffer cells (activated, lipid-handling)",
  "Cholangiocytes",
  "Hepatocytes (high metabolic, mixed-zone)",
  "Inflammatory macrophages (lipid or fibrosis-associated)",
  "Cholangiocytes/periportal epithelial cells",
  "Erythrocytes",
  "Cholangiocytes (progenitor-like or reactive phenotype)",
  "Pro-inflammatory macrophages",
  "Periportal hepatocytes (TGF-beta signaling)",
  "Plasma B cells",
  "Cholangiocytes, HSCs, or ductal progenitors (ECM remodeling)",
  "IFN-activated macrophages"
)

celltypes_to_run <- intersect(celltype_levels, unique(col_annot$celltype))

dir.create(outdir_main, showWarnings = FALSE, recursive = TRUE)

combined_results <- list()

## ---------- 1) Loop Through Each Cell Type ----------
for (ct in celltypes_to_run) {
  message(">>> Running: ", ct)
  
  cols_ct <- col_annot$col[col_annot$celltype == ct]
  if (length(cols_ct) < 2) { message("  skip: <2 columns"); next }
  
  counts_ct <- pb[, cols_ct, drop = FALSE]
  
  smp_ct <- col_annot$sample[col_annot$celltype == ct]
  grp_ct <- factor(col_annot$group [col_annot$celltype == ct])
  
  if (length(unique(grp_ct)) < 2) { message("  skip: single group"); next }
  
  
  grp_ct <- stats::relevel(grp_ct, ref = "CR")
  design  <- model.matrix(~ grp_ct)
  coef_ix <- match(paste0("grp_ct", "GDF15"), colnames(design))
  if (is.na(coef_ix)) { message("  skip: coef not found"); next }
  
  # edgeR 
  dge <- DGEList(counts_ct, samples = data.frame(sample = smp_ct, group = grp_ct))
  keep <- filterByExpr(dge, group = grp_ct, min.count = min_count)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  if (nrow(dge) == 0) { message("  skip: no genes after filter"); next }
  
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = coef_ix)
  
  # DE list
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene <- rownames(tab)
  
  # logCPM
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  grp_means <- sapply(c("CR","GDF15"), function(g)
    rowMeans(logcpm[, grp_ct == g, drop = FALSE])
  )
  tab$mean_logCPM_CR    <- grp_means[tab$gene, "CR"]
  tab$mean_logCPM_GDF15 <- grp_means[tab$gene, "GDF15"]
  
  # save CSV
  subdir <- file.path(outdir_main, safelabel(ct))
  dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(tab, file.path(subdir, paste0("DEG_", safelabel(ct), "_", assay_counts, ".csv")),
                   row.names = FALSE)
  
  # ---------- Volcano Plot (with Top 10 Upregulated/Downregulated Labels and Counts) ----------
  deg.data <- tab
  if (!"log2FC" %in% names(deg.data)) deg.data$log2FC <- deg.data$logFC
  deg.data$FDR[is.na(deg.data$FDR)] <- 1
  deg.data$logQ <- -log10(pmax(deg.data$FDR, 1e-300))
  
  deg.data$Group <- "normal"
  deg.data$Group[deg.data$FDR < fdr_cut & deg.data$log2FC >  lfc_cut] <- "up"
  deg.data$Group[deg.data$FDR < fdr_cut & deg.data$log2FC < -lfc_cut] <- "down"
  deg.data$Group <- factor(deg.data$Group, levels = c("down","normal","up"))
  if (!"gene" %in% names(deg.data)) deg.data$gene <- rownames(deg.data)
  
  
  N <- 10
  up_idx   <- which(deg.data$Group == "up")
  down_idx <- which(deg.data$Group == "down")
  up_ord   <- up_idx[order(deg.data$log2FC[up_idx],   decreasing = TRUE)]
  down_ord <- down_idx[order(deg.data$log2FC[down_idx], decreasing = FALSE)]
  top_ids  <- c(head(up_ord, N), head(down_ord, N))
  
  deg.data$Label <- ""
  if (length(top_ids) > 0) deg.data$Label[top_ids] <- deg.data$gene[top_ids]
  
  
  cnt <- table(deg.data$Group)
  n_up   <- ifelse("up"    %in% names(cnt), as.integer(cnt[["up"]]),    0L)
  n_down <- ifelse("down"  %in% names(cnt), as.integer(cnt[["down"]]),  0L)
  n_ns   <- ifelse("normal"%in% names(cnt), as.integer(cnt[["normal"]]),0L)
  
  p <- ggscatter(
    deg.data, x = "log2FC", y = "logQ",
    color = "Group",
    palette = c(down = "#2C7BB6", normal = "#BBBBBB", up = "#D7191C"),
    size = 1.6,
    label = "Label", font.label = 8,
    repel = TRUE,
    xlab = "log2FC (GDF15 vs CR)", ylab = "-log10(FDR)"
  ) +
    theme_base() +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed") +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
    labs(title = ct,
         subtitle = sprintf("Up: %d   Down: %d   Normal: %d", n_up, n_down, n_ns))+
    theme(
      plot.title = element_text(size = 10, face = "bold"), # 改这里
      plot.subtitle = element_text(size = 8),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10)
    )
  
  tif_path <- file.path(outdir_main, paste0("Volcano_", safelabel(ct)))
  
  graph2tif(x = p, file=tif_path, font = "Arial", cairo = TRUE,   
            width = 6.5, height = 4, bg = "transparent")
  graph2svg(x = p, file=tif_path, font = "Arial", cairo = TRUE,   
            width = 6.5, height = 4, bg = "transparent")
  
  
  sig <- subset(deg.data, FDR < fdr_cut & abs(log2FC) >= lfc_cut)
  utils::write.csv(sig, file.path(subdir, paste0("SIG_", safelabel(ct), ".csv")),
                   row.names = FALSE)
  
  
  tab$celltype <- ct
  combined_results[[ct]] <- tab

  
  combined <- if (length(combined_results)) do.call(rbind, combined_results) else data.frame()
  utils::write.csv(combined, file.path(outdir_main, "ALL_celltypes_DEG_combined.csv"),
                 row.names = FALSE)

  message("Done. Results in: ", normalizePath(outdir_main, mustWork = FALSE))
}










#################### show different types of KC ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

## Input data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1
head(obj)

## Extract "6: Kupffer cells"
Idents(obj) <- "celltype"
kup <- subset(obj, idents = "6: Kupffer cells (activated, lipid-handling)")
head(kup)

## Expression Macro, Timd4, Cd5l, Vsig4 and MoKC-specific gene Clec4f between treatments

# Kupffer / MoKC markers
#markers <- c("Marco", "Timd4", "Cd5l", "Vsig4", "Clec4f", "Spp1", "Trem2")
markers <- c("Marco", "Timd4", "Spp1", "Trem2")

# Keep only genes present in the object (avoids errors if some are missing)
features_use <- intersect(markers, rownames(kup))
if (length(features_use) == 0) stop("None of the requested genes are in this object.")

# Stacked, flipped violins split by group
p <- VlnPlot(
  kup,
  features = features_use,
  stack = TRUE, flip = TRUE,
  split.by = "group",
  pt.size = 0
)
print(p)


p <- DotPlot(
  kup,
  features = c("Marco", "Timd4", "Spp1", "Trem2", "Cd5l", "Vsig4", "Clec4f","Adgre1"),#"Cd5l", "Vsig4", "Clec4f"
  assay = "SCT",
  group.by = "group",
  scale = FALSE
)
print(p)


## LAM related KC
p <- DotPlot(
  kup,
  features = c("Trem2","Gpnmb"),
  assay = "SCT",
  group.by = "group",
  scale = FALSE
)
print(p)

# Show and save
print(p)
ggplot2::ggsave("Kupffer_markers_violin.png", p, width = 8, height = 6, dpi = 300)

graph2svg(x = NULL, file='54_1_AclyAcss2_VehEvt', font = "Arial", cairo = TRUE,   
          width = 4, height = 4, bg = "transparent")




## compare Trem2+ cells between CR and GDf15 in KC
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

# -------- 1) Subset the target Kupffer cluster --------
kup <- subset(
  obj,
  subset = celltype == "6: Kupffer cells (activated, lipid-handling)"
)

# Choose an assay that holds raw counts for detection
DefaultAssay(kup) <- "Spatial"

# Find Trem2 gene name in Spatial (handles case differences)
genes_sp <- rownames(kup[["Spatial"]])
idx <- which(toupper(genes_sp) == "TREM2")
if (length(idx) == 0) stop("Trem2 not found in the Spatial assay.")
gene <- genes_sp[idx[1]]

# --- Call Trem2+ cells using raw counts from Spatial ---
expr_counts <- GetAssayData(kup, assay = "Spatial", slot = "counts")
kup$Trem2_pos <- Matrix::colSums(expr_counts[gene, , drop = FALSE]) > 0

# --- Summarize by group (e.g., CR vs GDF15) ---
by_group <- kup@meta.data %>%
  group_by(group) %>%
  summarise(
    Trem2_pos = sum(Trem2_pos),
    Total     = n(),
    Trem2_pct = 100 * Trem2_pos / Total,
    .groups = "drop"
  ) %>% arrange(group)

print(by_group)

# --- Plot counts ---
p_count <- ggplot(by_group, aes(x = group, y = Trem2_pos)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = Trem2_pos), vjust = -0.2, size = 3) +
  labs(x = NULL, y = "Trem2+ cells (Kupffer cluster 6)",
       title = "Trem2+ Kupffer cells by group (Spatial counts)") +
  theme_classic(base_size = 12)
p_count
ggsave("Trem2_positive_counts_by_group_Spatial.png", p_count, width = 6, height = 4, dpi = 300)

# --- Plot percentages ---
p_pct <- ggplot(by_group, aes(x = group, y = Trem2_pct)) +
  geom_col(width = 0.7,fill = "#ADD8E6") +
  geom_text(aes(label = sprintf("%.1f%%", Trem2_pct)), vjust = -0.2, size = 6) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  labs(x = NULL, y = "% Trem2+ within Kupffer cells",
       title = "Fraction of Trem2+ Kupffer cells by group (Spatial counts)") +
  theme_classic(base_size = 12)
p_pct

p_pct <- p_pct +
  scale_y_continuous(breaks = seq(0, 30, 10),
                     labels = scales::percent_format(scale = 1)) +
  coord_cartesian(ylim = c(0, 30)) +
  theme(
    axis.text.y  = element_text(size = 14),  
    axis.text.x  = element_text(size = 14),  
    axis.title.y = element_text(size = 14), 
    axis.title.x = element_text(size = 14), 
    legend.text  = element_text(size = 14), 
    legend.title = element_text(size = 14),
    plot.title   = element_text(size = 14) 
  )

p_pct


ggsave("17_Trem2_Perc_KC.png", p_pct, width = 3, height = 4, dpi = 300)

graph2svg(x = NULL, file='17_Trem2_Perc_KC', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 5, bg = "transparent")



## compare Trem2+ cells between CR and GDf15 in 13: Pro-inflammatory macrophages
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

# -------- 1) Subset the target Kupffer cluster --------
kup <- subset(
  obj,
  subset = celltype == "13: Pro-inflammatory macrophages"
)

# Choose an assay that holds raw counts for detection
DefaultAssay(kup) <- "Spatial"

# Find Trem2 gene name in Spatial (handles case differences)
genes_sp <- rownames(kup[["Spatial"]])
idx <- which(toupper(genes_sp) == "TREM2")
if (length(idx) == 0) stop("Trem2 not found in the Spatial assay.")
gene <- genes_sp[idx[1]]

# --- Call Trem2+ cells using raw counts from Spatial ---
expr_counts <- GetAssayData(kup, assay = "Spatial", slot = "counts")
kup$Trem2_pos <- Matrix::colSums(expr_counts[gene, , drop = FALSE]) > 0

# --- Summarize by group (e.g., CR vs GDF15) ---
by_group <- kup@meta.data %>%
  group_by(group) %>%
  summarise(
    Trem2_pos = sum(Trem2_pos),
    Total     = n(),
    Trem2_pct = 100 * Trem2_pos / Total,
    .groups = "drop"
  ) %>% arrange(group)

print(by_group)

# --- Plot counts ---
p_count <- ggplot(by_group, aes(x = group, y = Trem2_pos)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = Trem2_pos), vjust = -0.2, size = 3) +
  labs(x = NULL, y = "Trem2+ cells (Kupffer cluster 6)",
       title = "Trem2+ Kupffer cells by group (Spatial counts)") +
  theme_classic(base_size = 12)
p_count
ggsave("Trem2_positive_counts_by_group_Spatial.png", p_count, width = 6, height = 4, dpi = 300)

# --- Plot percentages ---
p_pct <- ggplot(by_group, aes(x = group, y = Trem2_pct)) +
  geom_col(width = 0.7,fill = "#ADD8E6") +
  geom_text(aes(label = sprintf("%.1f%%", Trem2_pct)), vjust = -0.2, size = 6) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  labs(x = NULL, y = "% Trem2+ within Pro-inflammatory macrophages",
       title = "Fraction of Trem2+ Kupffer cells by group (Spatial counts)") +
  theme_classic(base_size = 12)
p_pct

p_pct <- p_pct +
  scale_y_continuous(breaks = seq(0, 30, 10),
                     labels = scales::percent_format(scale = 1)) +
  coord_cartesian(ylim = c(0, 30)) +
  theme(
    axis.text.y  = element_text(size = 14),  
    axis.text.x  = element_text(size = 14),  
    axis.title.y = element_text(size = 14), 
    axis.title.x = element_text(size = 14), 
    legend.text  = element_text(size = 14),  
    legend.title = element_text(size = 14),
    plot.title   = element_text(size = 14)  
  )

p_pct

ggsave("17_Trem2_Perc_InflaMacrophage.png", p_pct, width = 3, height = 4, dpi = 300)

graph2svg(x = NULL, file='17_Trem2_Perc_InflaMacrophage', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 5, bg = "transparent")





## compare LAM cells (Trem2/Cd9/Gpnmb) between CR and GDf15 in 13: Pro-inflammatory macrophages
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

# -------- 1) Subset the target cell type --------
mac <- subset(
  obj,
  subset = celltype == "13: Pro-inflammatory macrophages"
  # If your label varies slightly, use:
  # subset = grepl("^\\s*13:\\s*Pro[- ]inflammatory macrophages", celltype)
)

# -------- 2) Use Spatial assay --------
DefaultAssay(mac) <- "Spatial"

# Robust, case-insensitive symbol matching
genes_sp <- rownames(mac[["Spatial"]])
get_gene <- function(sym){
  idx <- which(toupper(genes_sp) == toupper(sym))
  if (!length(idx)) stop(sym, " not found in the Spatial assay.")
  genes_sp[idx[1]]
}

markers  <- c("Trem2","Gpnmb")
gene_ids <- vapply(markers, get_gene, character(1))

# -------- 3) Call triple-positives using raw counts (>0) --------
expr_counts <- GetAssayData(mac, assay = "Spatial", slot = "counts")
pos_mat <- expr_counts[gene_ids, , drop = FALSE] > 0
mac$Trem2_Cd9_Gpnmb_pos <- Matrix::colSums(pos_mat) == length(gene_ids)
# (For ≥2-of-3, use: Matrix::colSums(pos_mat) >= 2)

# -------- 4) Summarize by group (CR vs GDF15) --------
by_group <- mac@meta.data %>%
  group_by(group) %>%
  summarise(
    Triple_pos = sum(Trem2_Cd9_Gpnmb_pos),
    Total      = n(),
    Triple_pct = 100 * Triple_pos / Total,
    .groups = "drop"
  ) %>% arrange(group)

print(by_group)

# -------- 5) Plot counts --------
p_count <- ggplot(by_group, aes(x = group, y = Triple_pos)) +
  geom_col(width = 0.7, fill = "#ADD8E6", color = "black") +
  geom_text(aes(label = Triple_pos), vjust = -0.2, size = 3) +
  labs(x = NULL, y = "Trem2+Cd9+Gpnmb+ cells",
       title = "Triple-positive pro-inflammatory macrophages by group (Spatial)") +
  theme_classic(base_size = 12)
p_count
ggsave("ProInflamMac_Trem2_Cd9_Gpnmb_counts_by_group_Spatial.png",
       p_count, width = 6, height = 4, dpi = 300)

# -------- 6) Plot percentages --------
p_pct <- ggplot(by_group, aes(x = group, y = Triple_pct)) +
  geom_col(width = 0.7, fill = "#ADD8E6", color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", Triple_pct)), vjust = -0.2, size = 3) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  labs(x = NULL, y = "% LAM within pro-inflammatory macrophages",
       title = "Fraction of triple-positive pro-inflammatory macrophages (Spatial)") +
  theme_classic(base_size = 12)
p_pct

# Optional: cap Y at 30%
p_pct <- p_pct +
  scale_y_continuous(breaks = seq(0, 25, 5),
                     labels = scales::percent_format(scale = 1)) +
  coord_cartesian(ylim = c(0, 25))
p_pct

ggsave("17_LAM_Trem2_Gpnmb_Perc.png",
       p_pct, width = 3, height = 4, dpi = 300)










## compare moKCs (TIM4- VSIG4+) and ResKCs (TIM4+ VSIG4+) between CR and GDf15
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

# --- 1) Subset Kupffer cluster 6 ---
kup <- subset(
  obj,
  subset = celltype == "6: Kupffer cells (activated, lipid-handling)"
)

# --- 2) Use Spatial assay + find Timd4 / Vsig4 gene IDs robustly ---
DefaultAssay(kup) <- "Spatial"
genes_sp <- rownames(kup[["Spatial"]])

get_gene <- function(sym){
  idx <- which(toupper(genes_sp) == toupper(sym))
  if (length(idx) == 0) stop(paste0(sym, " not found in Spatial assay."))
  genes_sp[idx[1]]
}
g_Timd4 <- get_gene("Timd4")
g_Vsig4 <- get_gene("Vsig4")

# --- 3) Call positivity on raw counts (>0) ---
cnt <- GetAssayData(kup, assay = "Spatial", slot = "counts")
timd4_pos <- Matrix::colSums(cnt[g_Timd4, , drop = FALSE]) > 0
vsig4_pos <- Matrix::colSums(cnt[g_Vsig4, , drop = FALSE]) > 0

# MoKC = Timd4− & Vsig4+ ; EmKC = Timd4+ & Vsig4+
kup$KC_type <- ifelse(!timd4_pos &  vsig4_pos, "MoKC (Timd4− Vsig4+)",
                      ifelse( timd4_pos &  vsig4_pos, "EmKC (Timd4+ Vsig4+)",
                              "Other"))

# --- 4) Summarize by group (e.g., CR vs GDF15) ---
sum_tbl <- kup@meta.data %>%
  group_by(group, KC_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

# keep only MoKC & EmKC for plotting
plot_tbl <- sum_tbl %>%
  filter(KC_type %in% c("MoKC (Timd4− Vsig4+)", "EmKC (Timd4+ Vsig4+)")) %>%
  arrange(group, KC_type) %>%
  mutate(
    KC_type = factor(KC_type,
                     levels = c("MoKC (Timd4− Vsig4+)", "EmKC (Timd4+ Vsig4+)"))
  )


print(plot_tbl)

# --- 5a) Plot counts (side-by-side bars) ---
p_counts <- ggplot(plot_tbl, aes(x = group, y = n, fill = KC_type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.75), vjust = -0.2, size = 3) +
  scale_fill_manual(name = NULL, drop = FALSE,
                    values = c("MoKC (Timd4− Vsig4+)" = "#00BFC4",
                               "EmKC (Timd4+ Vsig4+)" = "#F8766D")) +
  labs(x = NULL, y = "Cell count (Kupffer cluster 6)",
       title = "MoKCs vs EmKCs by group") +
  theme_classic(base_size = 12)
p_counts
ggsave("KC_counts_by_group.png", p_counts, width = 6.5, height = 4.5, dpi = 300)

# --- 5b) Plot percentages within group ---

# Normalize labels (turn Unicode minus into ASCII hyphen)
plot_tbl <- plot_tbl %>%
  mutate(KC_type = gsub("\u2212", "-", KC_type)) %>%   # "−" -> "-"
  mutate(KC_type = factor(
    KC_type,
    levels = c("MoKC (Timd4- Vsig4+)", "EmKC (Timd4+ Vsig4+)")
  ))

p_pct <- ggplot(plot_tbl, aes(x = group, y = pct, fill = KC_type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct), group = KC_type),
    position = position_dodge(width = 0.75), vjust = -0.2, size = 5
  ) +
  # y-axis to 65%
  scale_y_continuous(breaks = seq(0, 65, 10),
                     labels = scales::percent_format(scale = 1)) +
  coord_cartesian(ylim = c(0, 65), clip = "off") +
  # force both legend keys + labels
  scale_fill_manual(
    name   = NULL,
    breaks = c("MoKC (Timd4- Vsig4+)", "EmKC (Timd4+ Vsig4+)"),
    labels = c("MoKC (Timd4- Vsig4+)", "EmKC (Timd4+ Vsig4+)"),
    values = c("MoKC (Timd4- Vsig4+)" = "#00BFC4",
               "EmKC (Timd4+ Vsig4+)" = "#F8766D"),
    drop   = FALSE
  ) +
  labs(x = NULL, y = "% within Kupffer cells",
       title = "MoKCs vs EmKCs (fraction within group)") +
  theme_classic(base_size = 12) +
  theme(legend.text = element_text(size = 10))
p_pct

p_pct <- p_pct +
  scale_y_continuous(breaks = seq(0, 60, 10),
                     labels = scales::percent_format(scale = 1)) +
  coord_cartesian(ylim = c(0, 65)) +
  theme(text = element_text(size = 14))

p_pct


ggsave("17_KC_Perc_by_group.png", p_pct, width = 3.5, height = 5, dpi = 300)

graph2svg(x = NULL, file='17_KC_Perc_by_group', font = "Arial", cairo = TRUE,   
          width = 5, height = 5, bg = "transparent")







#################### show different types of 13: Pro-inflammatory macrophages ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

## Input data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1
head(obj)
table(obj$celltype)

## Extract "13: Pro-inflammatory macrophages"
Idents(obj) <- "celltype"
Inflam <- subset(obj, idents = "13: Pro-inflammatory macrophages")
head(Inflam)

## Expression Macro, Timd4, Cd5l, Vsig4 and MoKC-specific gene Clec4f between treatments

# Kupffer / MoKC markers
#markers <- c("Marco", "Timd4", "Cd5l", "Vsig4", "Clec4f", "Spp1", "Trem2")
markers <- c("Marco", "Timd4", "Spp1", "Trem2")

# Keep only genes present in the object (avoids errors if some are missing)
features_use <- intersect(markers, rownames(Inflam))
if (length(features_use) == 0) stop("None of the requested genes are in this object.")

# Stacked, flipped violins split by group
p <- VlnPlot(
  Inflam,
  features = features_use,
  stack = TRUE, flip = TRUE,
  split.by = "group",
  pt.size = 0
)
print(p)


p <- DotPlot(
  Inflam,
  features = c("Marco", "Timd4", "Spp1", "Trem2","Cd5l", "Vsig4", "Clec4f", "Adgre1"),#"Cd5l", "Vsig4", "Clec4f"
  assay = "SCT",
  group.by = "group",
  scale = FALSE
)
print(p)

## LAM markers
p <- DotPlot(
  Inflam,
  features = c("Trem2","Gpnmb"),
  assay = "SCT",
  group.by = "group",
  scale = FALSE
)
print(p)


# Show and save
print(p)
ggplot2::ggsave("Kupffer_markers_violin.png", p, width = 8, height = 6, dpi = 300)









#################### show different types of 9: Activated HSCs/fibroblasts ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

## Input data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1
head(obj)
table(obj$celltype)

## Extract "9: Activated HSCs/fibroblasts"
Idents(obj) <- "celltype"
HSC <- subset(obj, idents = "9: Activated HSCs/fibroblasts")
head(HSC)

## Expression Macro, Timd4, Cd5l, Vsig4 and MoKC-specific gene Clec4f between treatments

# Kupffer / MoKC markers
#markers <- c("Marco", "Timd4", "Cd5l", "Vsig4", "Clec4f", "Spp1", "Trem2")
markers <- c("Marco", "Timd4", "Spp1", "Trem2")

# Keep only genes present in the object (avoids errors if some are missing)
features_use <- intersect(markers, rownames(HSC))
if (length(features_use) == 0) stop("None of the requested genes are in this object.")

# Stacked, flipped violins split by group
p <- VlnPlot(
  HSC,
  features = features_use,
  stack = TRUE, flip = TRUE,
  split.by = "group",
  pt.size = 0
)
print(p)


p <- DotPlot(
  HSC,
  features = c("Marco", "Timd4", "Spp1", "Trem2","Cd5l", "Vsig4", "Clec4f", "Col1a1","Col3a1", "Adgre1"),#"Cd5l", "Vsig4", "Clec4f"
  assay = "SCT",
  group.by = "group",
  scale = FALSE
)
print(p)

## LAM markers
p <- DotPlot(
  HSC,
  features = c("Trem2","Gpnmb"),
  assay = "SCT",
  group.by = "group",
  scale = FALSE
)
print(p)


# Show and save
print(p)
ggplot2::ggsave("Kupffer_markers_violin.png", p, width = 8, height = 6, dpi = 300)











## compare LAM-like cells (Trem2/Gpnmb) between CR and GDF15 in 9: Activated HSCs/fibroblasts
library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)

# -------- 1) Subset the target cell type --------
hscfib <- subset(
  obj,
  subset = celltype == "9: Activated HSCs/fibroblasts"
  # If the label varies slightly, use:
  # subset = grepl("^\\s*9:\\s*Activated HSCs/fibroblasts", celltype)
)

# -------- 2) Use Spatial assay --------
DefaultAssay(hscfib) <- "Spatial"

# Robust, case-insensitive symbol matching
genes_sp <- rownames(hscfib[["Spatial"]])
get_gene <- function(sym){
  idx <- which(toupper(genes_sp) == toupper(sym))
  if (!length(idx)) stop(sym, " not found in the Spatial assay.")
  genes_sp[idx[1]]
}

# LAM-like markers used here: Trem2 + Gpnmb (add "Cd9" if desired)
markers  <- c("Trem2","Gpnmb")
gene_ids <- vapply(markers, get_gene, character(1))

# -------- 3) Call double-positives using raw counts (>0) --------
expr_counts <- GetAssayData(hscfib, assay = "Spatial", slot = "counts")
pos_mat <- expr_counts[gene_ids, , drop = FALSE] > 0
hscfib$LAM_pos <- Matrix::colSums(pos_mat) == length(gene_ids)
# For ≥2-of-3 when adding Cd9: Matrix::colSums(pos_mat) >= 2

# -------- 4) Summarize by group (CR vs GDF15) --------
by_group <- hscfib@meta.data %>%
  group_by(group) %>%
  summarise(
    LAM_pos   = sum(LAM_pos),
    Total     = n(),
    LAM_pct   = 100 * LAM_pos / Total,
    .groups   = "drop"
  ) %>% arrange(group)

print(by_group)

# -------- 5) Plot counts --------
p_count <- ggplot(by_group, aes(x = group, y = LAM_pos)) +
  geom_col(width = 0.7, fill = "#ADD8E6", color = "black") +
  geom_text(aes(label = LAM_pos), vjust = -0.2, size = 3) +
  labs(x = NULL, y = "Trem2+Gpnmb+ cells",
       title = "LAM-like (Trem2+Gpnmb+) cells in Activated HSCs/fibroblasts") +
  theme_classic(base_size = 12)
p_count
ggsave("HSCfib_LAM_counts_by_group_Spatial.png", p_count, width = 6, height = 4, dpi = 300)

# -------- 6) Plot percentages --------
p_pct <- ggplot(by_group, aes(x = group, y = LAM_pct)) +
  geom_col(width = 0.7, fill = "#ADD8E6", color = "black") +
  geom_text(aes(label = sprintf("%.1f%%", LAM_pct)), vjust = -0.2, size = 3) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), limits = c(0, 100)) +
  labs(x = NULL, y = "% LAM within Activated HSCs/Myofibroblasts",
       title = "Fraction of LAM-like (Trem2+Gpnmb+) cells by group (Spatial)") +
  theme_classic(base_size = 12)
p_pct

# Optional: cap Y at 25%
p_pct <- p_pct +
  scale_y_continuous(breaks = seq(0, 15, 5),
                     labels = scales::percent_format(scale = 1)) +
  coord_cartesian(ylim = c(0, 15))
p_pct

ggsave("17_HSCfib_LAM_Trem2_Gpnmb_Perc.png", p_pct, width = 3, height = 4, dpi = 300)










#################### show inflammatory gene/lipidMeta expression in KC ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Tsc22d3","Bcl6","Irf2bp2","Zbtb16","Ifngr2","Vnn1","Bhlhe41",
           "Fasn", "Acly", "Lpin1", "Mlxipl", "Hmgcr", "Lss", "Plin5", "Chka")
ct    <- "6: Kupffer cells (activated, lipid-handling)"  # your target cell type
assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
ct_kc <- "6: Kupffer cells (activated, lipid-handling)"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
genes_interest <- c("Tsc22d3","Bcl6","Irf2bp2","Zbtb16","Ifngr2","Vnn1","Bhlhe41",
                    "Fasn", "Acly", "Lpin1", "Mlxipl", "Hmgcr", "Lss", "Plin5", "Chka")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_InflammLipMetGene_KC', font = "Arial", cairo = TRUE,   
          width = 9, height = 4, bg = "transparent")
graph2svg(x = NULL, file='Panel_InflammLipMetGene_KC', font = "Arial", cairo = TRUE,   
         width = 9, height = 4, bg = "transparent")
















#################### show inflammatory gene expression in KC ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Tsc22d3","Bcl6","Irf2bp2","Zbtb16","Ppargc1a","Tat","Ifngr2","Vnn1","Bhlhe41")
ct    <- "6: Kupffer cells (activated, lipid-handling)"  # your target cell type
assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
ct_kc <- "6: Kupffer cells (activated, lipid-handling)"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
genes_interest <- c("Tsc22d3","Bcl6","Irf2bp2","Zbtb16","Ppargc1a","Tat","Ifngr2","Vnn1","Bhlhe41")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_InflammatoryGene_KC', font = "Arial", cairo = TRUE,   
          width = 5, height = 4, bg = "transparent")








#################### show lipid-handling gene expression in KC ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Fasn", "Acly", "Lpin1", "Mlxipl", "Hmgcr", "Lss", "Plin5", "Chka")
ct    <- "6: Kupffer cells (activated, lipid-handling)"  # your target cell type
assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
ct_kc <- "6: Kupffer cells (activated, lipid-handling)"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
genes_interest <- c("Fasn", "Acly", "Lpin1", "Mlxipl", "Hmgcr", "Lss", "Plin5", "Chka")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_LipidHandGene_KC', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 4, bg = "transparent")












#################### show Circadian re-entrainment gene expression in KC ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Arntl", "Dbp", "Per2", "Per3", "Usp2", "Bhlhe41")
ct    <- "6: Kupffer cells (activated, lipid-handling)"  # your target cell type
assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
ct_kc <- "6: Kupffer cells (activated, lipid-handling)"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
genes_interest <- c("Arntl", "Dbp", "Per2", "Per3", "Usp2", "Bhlhe41")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_CircadianGene_KC', font = "Arial", cairo = TRUE,   
          width = 4.2, height = 4, bg = "transparent")









#################### show Circadian re-entrainment-antiinflammatory gene expression in in cholangiocytes ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Arntl", "Usp2", "Dbp", "Cry1", "Nfil3", "Tsc22d3", "Bcl6", "Zbtb16",  "Fkbp5", 
           "Irf2bp2", "Nrg4", "Casp6")

genes_interest <- c("Arntl", "Usp2", "Dbp", "Cry1", "Nfil3", "Tsc22d3", "Bcl6", "Zbtb16",  "Fkbp5", 
                    "Irf2bp2", "Nrg4", "Casp6")

ct    <- "7: Cholangiocytes"  # your target cell type
ct_kc <- "7: Cholangiocytes"

assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
#ct_kc <- "7: Cholangiocytes"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
#genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", "Chka", "Insig1", "Plin5", "Scarb1", "Gramd1c")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
# readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_Circadian-CREB-GC-Anti-inflammatory_Cholangio', font = "Arial", cairo = TRUE,   
          width = 6.5, height = 4, bg = "transparent")








#################### show Lipid handling and membrane/sterol remodeling gene expression in cholangiocytes ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", 
           "Chka", "Plin5", "Insig1", "Gramd1c", "Scarb1")

genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", 
                    "Chka", "Plin5", "Insig1", "Gramd1c", "Scarb1")

ct    <- "7: Cholangiocytes"  # your target cell type
ct_kc <- "7: Cholangiocytes"

assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
#ct_kc <- "7: Cholangiocytes"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
#genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", "Chka", "Insig1", "Plin5", "Scarb1", "Gramd1c")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
# readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_LipidTuningGene_Cholangio', font = "Arial", cairo = TRUE,   
          width = 5.7, height = 4, bg = "transparent")









#################### show Lipid handling and membrane/sterol remodeling gene expression in cholangiocytes ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Cldn1", "Cgn", "Celsr1", "Mmp7", 
           "Herpud1", "Stip1", "Sigmar1", "Ulk1", "Coq8a", "Clpx", "Depp1")

genes_interest <- c("Cldn1", "Cgn", "Celsr1", "Mmp7", 
                    "Herpud1", "Stip1", "Sigmar1", "Ulk1", "Coq8a", "Clpx", "Depp1")

ct    <- "7: Cholangiocytes"  # your target cell type
ct_kc <- "7: Cholangiocytes"

assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
#ct_kc <- "7: Cholangiocytes"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
#genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", "Chka", "Insig1", "Plin5", "Scarb1", "Gramd1c")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
# readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_BarrierPolarityStressResilient_Cholangio', font = "Arial", cairo = TRUE,   
          width = 6.2, height = 4, bg = "transparent")










#################### show Lipid handling gene expression in Pericentral hepatocytes (stress-responsive) ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
table(liver.mergeIntegrAnno1$celltype)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Fasn", "Acly", "Elovl3", "Lpin1", "Chka", "Gpam", "Plin5", "Pde9a", 
           "Cidec", "Fitm2", "Acot2", "Acot3")

genes_interest <- c("Fasn", "Acly", "Elovl3", "Lpin1", "Chka", "Gpam", "Plin5", "Pde9a", 
                    "Cidec", "Fitm2", "Acot2", "Acot3")

ct    <- "4: Pericentral hepatocytes (stress-responsive traits)"  # your target cell type
ct_kc <- "4: Pericentral hepatocytes (stress-responsive traits)"

assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
#ct_kc <- "7: Cholangiocytes"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
#genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", "Chka", "Insig1", "Plin5", "Scarb1", "Gramd1c")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
# readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_lipidHandling_PericentralhepaStressResp', font = "Arial", cairo = TRUE,   
          width = 6.5, height = 4, bg = "transparent")












#################### show mitochondrial QC gene expression in Pericentral hepatocytes (stress-responsive) ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
table(liver.mergeIntegrAnno1$celltype)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Hyou1", "Manf", "Ero1lb", "Herpud1", "Coq8a", "Clpx", "Slc25a44", "Shld2", "Paxx")

genes_interest <- c("Hyou1", "Manf", "Ero1lb", "Herpud1", "Coq8a", "Clpx", "Slc25a44", "Shld2", "Paxx")

ct    <- "4: Pericentral hepatocytes (stress-responsive traits)"  # your target cell type
ct_kc <- "4: Pericentral hepatocytes (stress-responsive traits)"

assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
#ct_kc <- "7: Cholangiocytes"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
#genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", "Chka", "Insig1", "Plin5", "Scarb1", "Gramd1c")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
# readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_MitoQC_PericentralhepaStressResp', font = "Arial", cairo = TRUE,   
          width = 5, height = 4, bg = "transparent")











#################### show re-time-infla-cell death gene expression in Pericentral hepatocytes (stress-responsive) ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
table(liver.mergeIntegrAnno1$celltype)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Arntl", "Bcl6", "Zbtb16", "Tsc22d3", "Fkbp5", "Irf2bp2", "Dbp", "Bhlhe41", "Ifngr2", "Cish", "Casp6", "Nrg4", "Por")

genes_interest <- c("Arntl", "Bcl6", "Zbtb16", "Tsc22d3", "Fkbp5", "Irf2bp2", "Dbp", "Bhlhe41", "Ifngr2", "Cish", "Casp6", "Nrg4", "Por")

ct    <- "4: Pericentral hepatocytes (stress-responsive traits)"  # your target cell type
ct_kc <- "4: Pericentral hepatocytes (stress-responsive traits)"

assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
#ct_kc <- "7: Cholangiocytes"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
#genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", "Chka", "Insig1", "Plin5", "Scarb1", "Gramd1c")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
# readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_ClockInfl_PericentralhepaStressResp', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")












#################### show gene expression profile in other Pericentral hepatocytes (lipid accumulation state) ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
table(liver.mergeIntegrAnno1$celltype)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Fasn", "Acly", "Lpin1", "Chka", "Mlxipl", "Gpam", "Plin5",
           "Arntl", "Bcl6", "Zbtb16", "Tsc22d3", "Irf2bp2", "Fkbp4", "Dbp", "Per3", "Usp2",
           "Ifngr2", "Casp6", "Raet1d", "Nrg4")

genes_interest <- c("Fasn", "Acly", "Lpin1", "Chka", "Mlxipl", "Gpam", "Plin5",
                    "Arntl", "Bcl6", "Zbtb16", "Tsc22d3", "Irf2bp2", "Fkbp4", "Dbp", "Per3", "Usp2",
                    "Ifngr2", "Casp6", "Raet1d", "Nrg4")

ct    <- "2: Pericentral hepatocytes (lipid accumulation)"  # your target cell type
ct_kc <- "2: Pericentral hepatocytes (lipid accumulation)"

assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
#ct_kc <- "7: Cholangiocytes"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
#genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", "Chka", "Insig1", "Plin5", "Scarb1", "Gramd1c")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
# readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_genProf_PericentralhepaLipAccu', font = "Arial", cairo = TRUE,   
          width = 9, height = 4, bg = "transparent")













#################### show gene expression profile in other Pericentral hepatocytes (Wnt/β-catenin ) ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
table(liver.mergeIntegrAnno1$celltype)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Zbtb16", "Bcl6", "Tsc22d3", "Nfil3", "Arntl", "Dbp", "Per2", "Per3",
           "Fasn", "Acly", "Lpin1", "Chka", "Insig1", 
           "Hspa5", "Hyou1", "Herpud1", "Clpx", "Etfbkmt", "Depp1",
           "Wee1", "Alas1", "Pla2g6")

genes_interest <- c("Zbtb16", "Bcl6", "Tsc22d3", "Nfil3", "Arntl", "Dbp", "Per2", "Per3",
                    "Fasn", "Acly", "Lpin1", "Chka", "Insig1", 
                    "Hspa5", "Hyou1", "Herpud1", "Clpx", "Etfbkmt", "Depp1",
                    "Wee1", "Alas1", "Pla2g6")

ct    <- "5: Pericentral hepatocytes (Wnt/β-catenin signaling)"  # your target cell type
ct_kc <- "5: Pericentral hepatocytes (Wnt/β-catenin signaling)"

assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
#ct_kc <- "7: Cholangiocytes"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
#genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", "Chka", "Insig1", "Plin5", "Scarb1", "Gramd1c")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
# readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_genProf_PericentralhepaWntBeta-catenin ', font = "Arial", cairo = TRUE,   
          width = 9.5, height = 4, bg = "transparent")











#################### show gene expression profile in Periportal hepatocytes (inflammatory phenotype) ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
table(liver.mergeIntegrAnno1$celltype)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Fasn", "Acly", "Lpin1", "Chka", "Gpam", "Mlxipl", "Plin5", "Cpt1a", "Acot2", "Nceh1", 
           "Arntl", "Bcl6", "Zbtb16", "Tsc22d3", "Dbp", "Per2", "Per3", "Usp2",
           "Hyou1", "Ero1lb", "Clpx", "Coq8a", "Slc25a44", "Zfand4", "Shld2", "Paxx", "Depp1")

genes_interest <- c("Fasn", "Acly", "Lpin1", "Chka", "Gpam", "Mlxipl", "Plin5", "Cpt1a", "Acot2", "Nceh1", 
                    "Arntl", "Bcl6", "Zbtb16", "Tsc22d3", "Dbp", "Per2", "Per3", "Usp2",
                    "Hyou1", "Ero1lb", "Clpx", "Coq8a", "Slc25a44", "Zfand4", "Shld2", "Paxx", "Depp1")

ct    <- "0:Periportal hepatocytes (inflammatory phenotype)"  # your target cell type
ct_kc <- "0:Periportal hepatocytes (inflammatory phenotype)"

assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
#ct_kc <- "7: Cholangiocytes"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
#genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", "Chka", "Insig1", "Plin5", "Scarb1", "Gramd1c")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
# readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_genProf_PeriportalhepaInfla', font = "Arial", cairo = TRUE,   
          width = 13, height = 4, bg = "transparent")











#################### show gene expression profile in Periportal hepatocytes (TGF-β signaling state) ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
table(liver.mergeIntegrAnno1$celltype)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Cdkn1a", "Bcl6", "Zbtb16", "Tsc22d3", "Irf2bp2",
           "Fasn", "Acly", "Lpin1", "Gpam", "Cpt1a", "Acot3", "Cidec",
           "Hyou1", "Manf", "Clpx", "Coq8a", "Shld2",
           "Dbp", "Usp2", "Nrg4", "Paqr7", "Casp6")

genes_interest <- c("Cdkn1a", "Bcl6", "Zbtb16", "Tsc22d3", "Irf2bp2",
                    "Fasn", "Acly", "Lpin1", "Gpam", "Cpt1a", "Acot3", "Cidec",
                    "Hyou1", "Manf", "Clpx", "Coq8a", "Shld2",
                    "Dbp", "Usp2", "Nrg4", "Paqr7", "Casp6")

ct    <- "14: Periportal hepatocytes (TGF-beta signaling)"  # your target cell type
ct_kc <- "14: Periportal hepatocytes (TGF-beta signaling)"

assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
#ct_kc <- "7: Cholangiocytes"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
#genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", "Chka", "Insig1", "Plin5", "Scarb1", "Gramd1c")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
# readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_genProf_PeriportalhepaTGFb', font = "Arial", cairo = TRUE,   
          width = 11.5, height = 4, bg = "transparent")











#################### show gene expression profile in activated HSCs and myofibroblast ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
table(liver.mergeIntegrAnno1$celltype)
obj = liver.mergeIntegrAnno1
head(obj)

## here we change cell type cluster 9 name:
liver.mergeIntegrAnno1$celltype <- gsub(
  pattern = "9: Inflammatory macrophages \\(lipid or fibrosis-associated\\)",
  replacement = "9: Activated HSCs/fibroblasts",
  x = liver.mergeIntegrAnno1$celltype
)
#celltype.stim
liver.mergeIntegrAnno1$celltype.stim <- paste(obj$celltype, obj$group, sep = "_")
table(liver.mergeIntegrAnno1$celltype.stim)

###### save file
saveRDS(liver.mergeIntegrAnno1, file = "liver_merge_integr_Ann1.rds")


# ---- inputs ----
genes <- c("Bcl6", "Zbtb16", "Tsc22d3", "Fkbp5", "Arntl", "Dbp", "Usp2",
           "Fasn", "Acly", "Lpin1","Chka", "Thrsp", "Insig1", "Me1", "Gpam",
           "Hspa5", "Hyou1", "Manf", "Clpx", "Coq8a", "Slc25a44", "Zfand4", "Sec23b", "Depp1",
           "Col1a1", "Cpxm1")

genes_interest <- c("Bcl6", "Zbtb16", "Tsc22d3", "Fkbp5", "Arntl", "Dbp", "Usp2",
                    "Fasn", "Acly", "Lpin1","Chka", "Thrsp", "Insig1", "Me1", "Gpam",
                    "Hspa5", "Hyou1", "Manf", "Clpx", "Coq8a", "Slc25a44", "Zfand4", "Sec23b", "Depp1",
                    "Col1a1", "Cpxm1")

ct    <- "9: Activated HSCs/fibroblasts"  # your target cell type
ct_kc <- "9: Activated HSCs/fibroblasts"

assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
#ct_kc <- "7: Cholangiocytes"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
#genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", "Chka", "Insig1", "Plin5", "Scarb1", "Gramd1c")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
# readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_genProf_HSCs', font = "Arial", cairo = TRUE,   
          width = 11.5, height = 4, bg = "transparent")











#################### show gene expression/only inflam gene profile in Plasma B cells ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
table(liver.mergeIntegrAnno1$celltype)
obj = liver.mergeIntegrAnno1
head(obj)



# ---- inputs ----
genes <- c("Herpud1", "Coq8a", "Clpx", "Igkv19-93" )
genes_interest <- c("Herpud1", "Coq8a", "Clpx", "Igkv19-93")

#All gene panel
#genes_interest <- c("Herpud1", "Coq8a", "Clpx", "Igkv19-93", "Fasn", "Acly", "Hmgcr", "Lpin1", "Chka", "Mlxipl")


ct    <- "15: Plasma B cells"  # your target cell type
ct_kc <- "15: Plasma B cells"

assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
#ct_kc <- "7: Cholangiocytes"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
#genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", "Chka", "Insig1", "Plin5", "Scarb1", "Gramd1c")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
# readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_genProf2_PlasmaB', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 4, bg = "transparent")

graph2svg(x = NULL, file='Panel_genProf2_PlasmaB', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 4, bg = "transparent")






#################### show gene expression profile in pro-inflammatory cells ######################
## Load packages and data
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
table(liver.mergeIntegrAnno1$celltype)
obj = liver.mergeIntegrAnno1
head(obj)


# ---- inputs ----
genes <- c("Bcl6", "Tsc22d3", "Irf2bp2",  "Cdkn1a", "Acly", "Fasn", "Lpin1", "Chka", "Scarb1", "Ppp1r3b", 
           "Hyou1", "Manf", "Sigmar1", "Clpx", "Coq8a", "Atp2b2", "Tars")

genes_interest <- c("Bcl6", "Tsc22d3", "Irf2bp2",  "Cdkn1a", "Acly", "Fasn", "Lpin1", "Chka", "Scarb1", "Ppp1r3b", 
                    "Hyou1", "Manf", "Sigmar1", "Clpx", "Coq8a", "Atp2b2", "Tars")

ct    <- "13: Pro-inflammatory macrophages"  # your target cell type
ct_kc <- "13: Pro-inflammatory macrophages"

assay_use <- "SCT"   # or "RNA"
slot_use  <- "data"  # for SCT/RNA normalized values

# ---- 1) subset to the cell type; set assay ----
obj_ct <- subset(obj, subset = celltype == ct)
DefaultAssay(obj_ct) <- assay_use

# ---- 2) fetch per-cell expression + sample/group metadata ----
vars <- c(genes, "orig.ident", "group")
dat  <- FetchData(obj_ct, vars = vars, slot = slot_use)
dat$cell <- rownames(dat)
head(dat)

# ---- 3) per-sample means (average cells of same sample within this cell type) ----
per_sample <- dat |>
  pivot_longer(all_of(genes), names_to = "gene", values_to = "expr") |>
  group_by(orig.ident, group, gene) |>
  summarise(expr = mean(expr, na.rm = TRUE), .groups = "drop") |>
  mutate(
    gene  = factor(gene, levels = genes),
    group = factor(group, levels = c("CR","GDF15"))
  )

# ---- 4) group means ± s.e.m. across samples ----
sum_df <- per_sample |>
  group_by(gene, group) |>
  summarise(
    n    = dplyr::n(),
    mean = mean(expr, na.rm = TRUE),
    se   = sd(expr,  na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ---- 5) plot: bars = mean±sem; points = each sample ----
p <- ggplot(sum_df, aes(x = gene, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.18) +
  geom_point(data = per_sample,
             aes(x = gene, y = expr, color = group),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.08),
             size = 2, alpha = 0, inherit.aes = FALSE) +
  scale_fill_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6")) +
  scale_color_manual(values = c(CR = "#FF9999", GDF15 = "#2C7BB6"), guide = "none") +
  labs(
    title = ct,
    subtitle = sprintf("%s assay, %s slot, FDR from pseudobulk", assay_use, slot_use),
    x = NULL, y = "Normalized expression", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 10, face = "bold"), # 改这里
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45,hjust = 1),
    axis.text.y = element_text(size = 10)
  )
p



## get p-value using pseudobulking analysis:
# ========= packages =========
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
for (pkg in c("edgeR","Seurat","dplyr","tidyr","readr","stringr")) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE)
  library(pkg, character.only=TRUE)
}



# ========= choose raw-count assay (Visium usually "Spatial"; scRNA "RNA") =========
assay_counts <- if ("Spatial" %in% Seurat::Assays(obj)) "Spatial" else "RNA"

# ========= aggregate raw counts per (sample × celltype) =========
pb <- AggregateExpression(
  obj, assays = assay_counts, slot = "counts",
  group.by = c("orig.ident","celltype")
)[[assay_counts]]            # genes x columns(sample_celltype)

## ---- build clean annotation ----
cn <- colnames(pb)
us <- regexpr("_", cn)

annot <- data.frame(
  col           = cn,
  sample        = substring(cn, 1, us - 1),
  celltype_raw  = trimws(substring(cn, us + 1)),
  stringsAsFactors = FALSE
)

# remove trailing _<digits> appended by AggregateExpression (e.g. "_6")
annot$celltype_clean <- trimws(sub("_[0-9]+$", "", annot$celltype_raw))

# (check)
head(annot[, c("col","sample","celltype_raw","celltype_clean")], 3)


# Sample -> group map pulled directly from your metadata
sample_group <- obj@meta.data |>
  dplyr::select(orig.ident, group) |>
  dplyr::distinct() |>
  tibble::deframe()           # named vector: names=orig.ident, values=group (CR/GDF15)
head (sample_group)

# ========= edgeR helper for ONE cell type =========
pseudobulk_edger_one <- function(ct_name, pb_mat, annot, sample_group, min.count = 10L) {
  # remove leading "N:" index in your celltype string
  target <- trimws(sub("^\\d+\\s*:\\s*", "", ct_name))
  
  # pick columns for that cell type (cleaned)
  cols <- annot$col[annot$celltype_clean == target]
  message("Matched ", length(cols), " columns for: ", target)
  if (length(cols) < 2) return(NULL)
  
  counts <- pb_mat[, cols, drop = FALSE]
  sample <- annot$sample[match(cols, annot$col)]
  group  <- factor(unname(sample_group[sample]), levels = c("CR","GDF15"))
  meta   <- data.frame(sample = sample, group = group, row.names = cols)
  
  dge <- edgeR::DGEList(counts, samples = meta)
  keep <- edgeR::filterByExpr(dge, group = group, min.count = min.count)
  if (sum(keep) < 10) return(NULL)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  
  design <- model.matrix(~ group, data = dge$samples)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupGDF15")
  
  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene     <- rownames(tab)
  tab$FDR      <- p.adjust(tab$PValue, "BH")
  tab$celltype <- target
  tab
}


# ========= EXAMPLE: run for your KC cluster =========
#ct_kc <- "7: Cholangiocytes"
res_kc <- pseudobulk_edger_one(ct_kc, pb, annot, sample_group)
head(res_kc)

unique(annot$celltype_clean)

# (optional) only your 9 genes
#genes_interest <- c("Fasn", "Acly", "Lpin1", "Acacb", "Agpat2", "Gpam", "Chka", "Insig1", "Plin5", "Scarb1", "Gramd1c")
res_kc_subset <- res_kc |>
  dplyr::filter(gene %in% genes_interest) |>
  dplyr::select(celltype, gene, logFC, logCPM, F, PValue, FDR) |>
  arrange(PValue)
print(res_kc_subset)

# ========= Loop over ALL cell types =========
all_cts <- sort(unique(annot$celltype))
res_list <- lapply(all_cts, pseudobulk_edger_one, pb_mat = pb, annot = annot, sample_group = sample_group)
res_all  <- dplyr::bind_rows(res_list)
# readr::write_csv(res_all, "DE_edgeR_pseudobulk_all_celltypes.csv")


# put the label a bit above the tallest bar for each gene
ypos <- sum_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop")

# keep exactly 3 decimals (very small FDRs shown as "<0.001")
lab <- res_kc_subset %>%
  dplyr::mutate(label = ifelse(is.na(FDR), "",
                               ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)))) %>%
  dplyr::left_join(ypos, by = "gene")


# add labels to your plot 'p'
p + geom_text(data = lab,
              aes(x = gene, y = y, label = label),
              inherit.aes = FALSE, size = 3.8, vjust = 0) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

graph2tif(x = NULL, file='Panel_genProf_ProInflamma', font = "Arial", cairo = TRUE,   
          width = 9, height = 4, bg = "transparent")

graph2svg(x = NULL, file='Panel_genProf_ProInflamma', font = "Arial", cairo = TRUE,   
          width = 9, height = 4, bg = "transparent")























########################### Enrichment analysis ######################
## Load the data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1


### Metabolic analysis 
## This part includes following:
# 1. Metabolic pathway activity
# 2. Metabolic interactions
# 3. Metabolites analysis_Flux balance analysis (FBA)
#### Compare metabolic activities between GDF15 vs CR
## We test activity using 85 metabolism pathways gene sets
BiocManager::install("GSVA")
BiocManager::install("clusterProfiler")
devtools::install_github("ncborcherding/escape")
BiocManager::install("dittoSeq")
BiocManager::install("org.Mm.eg.db")


## load library
library(GSVA)
library(msigdbr)
library(escape)
library(GSEABase)
library(dittoSeq)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GSEABase)
library(org.Dm.eg.db)
library(presto)
library(msigdbr)
library(doMC)
library(dplyr)
library(Seurat)
library(irGSEA)

# Load geneset pathways
MetaPW = read.gmt("MetabolicPathway3.gmt")
str(MetaPW)
table(MetaPW$term)
head(MetaPW)

#load DEG from Pseudobulking analysis
Deg = read.csv("ALL_celltypes_DEG_combined.csv", header=TRUE)
head(Deg)


## =======================
## 0) Packages
## =======================
pkgs <- c("dplyr","tidyr","tibble","purrr","fgsea","pheatmap","stringr")
inst <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(inst)) install.packages(inst, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))


## =======================
## 1) Build pathway list from MetaPW
##    (MetaPW already in memory: term, gene)
## =======================
# clean a bit
MetaPW <- MetaPW |>
  dplyr::mutate(term = stringr::str_trim(term),
                gene = stringr::str_trim(gene)) |>
  dplyr::filter(nzchar(term), nzchar(gene))

# list: names = pathway terms, values = character vector of genes
PATHS_raw <- split(MetaPW$gene, MetaPW$term)
# (optional) drop duplicated genes within a term
PATHS <- lapply(PATHS_raw, function(v) unique(v[!is.na(v)]))


## =======================
## 2) Helper to run fgsea for ONE cell type
## =======================
rm(list = c("run_fgsea_one"))  

run_fgsea_one <- function(df_ct, pathways = PATHS,
                          minSize = 5, maxSize = 5000) {
  
  pcol <- if ("PValue" %in% names(df_ct)) "PValue" else "FDR"
  
  
  df_ct <- df_ct |>
    dplyr::filter(!is.na(logFC), is.finite(logFC),
                  !is.na(.data[[pcol]]), is.finite(.data[[pcol]])) |>
    dplyr::mutate(rank_stat = sign(logFC) * -log10(pmax(.data[[pcol]], 1e-300))) |>
    dplyr::arrange(dplyr::desc(rank_stat)) |>
    dplyr::group_by(gene) |>
    dplyr::slice_max(order_by = abs(rank_stat), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  
  stats <- df_ct$rank_stat
  names(stats) <- df_ct$gene
  
  
  PATHS_use <- lapply(pathways, function(g) intersect(g, names(stats)))
  lens <- vapply(PATHS_use, length, 1L)
  PATHS_use <- PATHS_use[lens >= minSize & lens <= maxSize]
  if (!length(PATHS_use)) return(tibble::tibble())
  
  
  if ("fgseaMultilevel" %in% getNamespaceExports("fgsea")) {
    fg <- fgsea::fgseaMultilevel(pathways = PATHS_use, stats = stats,
                                 minSize = minSize, maxSize = maxSize)
  } else {
   
    fg <- fgsea::fgsea(pathways = PATHS_use, stats = stats,
                       minSize = minSize, maxSize = maxSize)
  }
  
  tibble::as_tibble(fg) |>
    dplyr::select(pathway, size, NES, padj, leadingEdge)
}


## =======================
## 3) Run fgsea per cell type
##    (Deg already in memory)
## =======================
stopifnot(all(c("celltype","gene","logFC") %in% names(Deg)))
ct_list   <- split(Deg, Deg$celltype)
head(ct_list)

fgsea_list <- purrr::imap(ct_list, ~ run_fgsea_one(.x, PATHS) |>
                            dplyr::mutate(celltype = .y))
fgsea_tbl  <- dplyr::bind_rows(fgsea_list)

# ############%%%%%%%%% for only show 3 clusters %%%%%%%%%#############
# fgsea_tbl <- fgsea_tbl %>%
#   dplyr::filter(celltype %in% c(
#     "15: Plasma B cells",
#     "13: Pro-inflammatory macrophages",
#     "6: Kupffer cells"
#   ))


dir.create("GSEA_MetaPW_by_celltype_3Clusters", showWarnings = FALSE)
readr::write_csv(fgsea_tbl, "fgsea_all_celltypes_MetaPW_3Clusters.csv")
head(fgsea_tbl)

## =======================
## 4.1) Plot data (cell type separately)
## =======================
## =====  dotplot (color=NES, size=-log10(padj)) =====
pkgs <- c("dplyr","tidyr","stringr","forcats","ggplot2","purrr")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# 参数
alpha_padj <- 0.05   
top_k_up   <- 8      
top_k_down <- 8      
wrap_width <- 40     
outdir     <- "GSEA_MetaPW_by_celltype/dotplots_ct"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

stopifnot(all(c("pathway","NES","padj","celltype") %in% names(fgsea_tbl)))

# Plot ------------------------------------------------------
plot_fgsea_dot_ct <- function(fgsea_tbl, ct,
                              alpha_padj = 0.05, top_k_up = 8, top_k_down = 8,
                              wrap_width = 40, outdir = ".") {
  
  df <- fgsea_tbl %>%
    dplyr::filter(celltype == ct, !is.na(padj))
  
  
  sig <- df %>% dplyr::filter(padj < alpha_padj)
  
  
  up   <- sig %>% dplyr::filter(NES > 0)  %>%
    dplyr::slice_max(order_by = NES,      n = top_k_up,   with_ties = FALSE)
  down <- sig %>% dplyr::filter(NES < 0)  %>%
    dplyr::slice_min(order_by = NES,      n = top_k_down, with_ties = FALSE)
  
  keep <- dplyr::bind_rows(up, down) %>% dplyr::distinct(pathway, .keep_all = TRUE)
  
  
  if (nrow(keep) == 0) {
    keep <- df %>% dplyr::slice_max(order_by = abs(NES), n = top_k_up + top_k_down, with_ties = FALSE)
  }
  
  keep <- keep %>%
    dplyr::mutate(
      negLog10Padj   = -log10(pmax(padj, 1e-300)),
      pathway_wrapped = stringr::str_wrap(pathway, width = wrap_width)
    ) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(pathway_wrapped = forcats::fct_reorder(pathway_wrapped, NES))
  
  
  p <- ggplot(keep, aes(x = 1, y = pathway_wrapped)) +
    geom_point(aes(size = negLog10Padj, color = NES), alpha = 0.9) +
    scale_size_continuous(name = expression(-log[10](padj))) +
    scale_color_gradient2(low = "#2C7BB6", mid = "grey80", high = "#D7191C",
                          midpoint = 0, name = "NES") +
    scale_x_continuous(labels = NULL, breaks = NULL) +
    labs(
      title = ct,
      subtitle = sprintf("padj < %.2f | Top %d up / %d down", alpha_padj, top_k_up, top_k_down),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.y = element_text(size = 9, lineheight = 0.7),
      panel.grid.minor = element_blank()
    )
  
  
  h <- max(4, 0.28 * nrow(keep) + 1.5)
  fn_safe <- gsub("[^A-Za-z0-9]+", "_", ct)
  
  ggsave(file.path(outdir, paste0("fgsea_dot_", fn_safe, ".png")), p, width = 5.5, height = h, dpi = 300)
  
  p
}

# Generate plots for each cell type in batch -----------------------------------------------------
plots_by_ct <- purrr::map(unique(fgsea_tbl$celltype),
                          ~ plot_fgsea_dot_ct(fgsea_tbl, ct = .x,
                                              alpha_padj = alpha_padj,
                                              top_k_up = top_k_up,
                                              top_k_down = top_k_down,
                                              wrap_width = wrap_width,
                                              outdir = outdir))




## =======================
## 4.2) Plot data (all cell type together)
## =======================

# ---- packages ----
pkgs <- c("dplyr","tidyr","stringr","forcats","ggplot2","purrr")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- parameters ----
alpha_padj <- 0.05   
top_k_up   <- 8    
top_k_down <- 8   
wrap_width <- 40    

# fgsea_tbl: including pathway, NES, padj, celltype
stopifnot(all(c("pathway","NES","padj","celltype") %in% names(fgsea_tbl)))

# ---- Select significant genes and identify the top up-/downregulated genes within each cell type ----
fg_sig <- fgsea_tbl %>%
  filter(!is.na(padj), padj < alpha_padj)

top_by_ct <- fg_sig %>%
  group_by(celltype) %>%
  
  arrange(desc(NES), .by_group = TRUE) %>%
  mutate(rank_up = ifelse(NES > 0, row_number(), NA_integer_)) %>%
  arrange(abs(NES), .by_group = TRUE) %>% 
  ungroup() %>%
  group_by(celltype) %>%
 
  slice_max(order_by = if_else(NES > 0, abs(NES), -Inf),
            n = top_k_up, with_ties = FALSE) %>%
  ungroup() %>%
  
  bind_rows(
    fg_sig %>%
      group_by(celltype) %>%
      slice_max(order_by = if_else(NES < 0, abs(NES), -Inf),
                n = top_k_down, with_ties = FALSE) %>%
      ungroup()
  ) %>%
  distinct(celltype, pathway, .keep_all = TRUE) %>%   # 去重
  mutate(
    negLog10Padj = -log10(padj),
    dir = if_else(NES >= 0, "Up (NES>0)", "Down (NES<0)"),
    pathway_wrapped = stringr::str_wrap(pathway, width = wrap_width)
  )


top_by_ct <- top_by_ct %>%
  arrange(desc(NES)) %>%
  mutate(pathway_wrapped = forcats::fct_reorder(pathway_wrapped, NES))

# ----  dotplot ----
p <- ggplot(top_by_ct,
            aes(x = celltype, y = pathway_wrapped)) +
  geom_point(aes(size = negLog10Padj, color = NES), alpha = 0.9) +
  scale_size_continuous(name = expression(-log[10](padj))) +
  scale_color_gradient2(low = "#2C7BB6", mid = "grey80", high = "#D7191C",
                        midpoint = 0, name = "NES") +
  labs(x = "Cell type", y = NULL,
       title = "Top regulated pathways per cell type",
       subtitle = paste0("Filtered by padj < ", alpha_padj,
                         "; size = -log10(padj), color = NES")) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 9, lineheight = 0.7),
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p)

## save figure
graph2tif(x = NULL, file='Enrichment_MetabolicActivity_myGeneset', font = "Arial", cairo = TRUE,   
          width = 9, height = 9.5, bg = "transparent")










################ Gene enrichment for 3 celltypes '15: Plasma B cells', "13: Pro-inflammatory macrophages" and "6: Kupffer cells"
## Load the data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1


### Metabolic analysis 
## This part includes following:
# 1. Metabolic pathway activity
# 2. Metabolic interactions
# 3. Metabolites analysis_Flux balance analysis (FBA)
#### Compare metabolic activities between GDF15 vs CR
## We test activity using 85 metabolism pathways gene sets
BiocManager::install("GSVA")
BiocManager::install("clusterProfiler")
devtools::install_github("ncborcherding/escape")
BiocManager::install("dittoSeq")
BiocManager::install("org.Mm.eg.db")
# install irGSEA packages from Github
if (!requireNamespace("irGSEA", quietly = TRUE)) { 
  devtools::install_github("chuiqin/irGSEA", force =T)
}


## load library
library(GSVA)
library(msigdbr)
library(escape)
library(GSEABase)
library(dittoSeq)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GSEABase)

library(presto)
library(msigdbr)
library(dplyr)
library(Seurat)
library(irGSEA)

library(org.Dm.eg.db)
library(doMC)

# Load geneset pathways
MetaPW = read.gmt("MetabolicPathway3.gmt")
str(MetaPW)
table(MetaPW$term)
head(MetaPW)

#load DEG from Pseudobulking analysis
Deg = read.csv("ALL_celltypes_DEG_combined.csv", header=TRUE)
head(Deg)


## =======================
## 0) Packages
## =======================
pkgs <- c("dplyr","tidyr","tibble","purrr","fgsea","pheatmap","stringr")
inst <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(inst)) install.packages(inst, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))


## =======================
## 1) Build pathway list from MetaPW
##    (MetaPW already in memory: term, gene)
## =======================
# clean a bit
MetaPW <- MetaPW |>
  dplyr::mutate(term = stringr::str_trim(term),
                gene = stringr::str_trim(gene)) |>
  dplyr::filter(nzchar(term), nzchar(gene))

# list: names = pathway terms, values = character vector of genes
PATHS_raw <- split(MetaPW$gene, MetaPW$term)
# (optional) drop duplicated genes within a term
PATHS <- lapply(PATHS_raw, function(v) unique(v[!is.na(v)]))


## =======================
## 2) Helper to run fgsea for ONE cell type
## =======================
rm(list = c("run_fgsea_one")) 

run_fgsea_one <- function(df_ct, pathways = PATHS,
                          minSize = 5, maxSize = 5000) {
  
  pcol <- if ("PValue" %in% names(df_ct)) "PValue" else "FDR"
  
  
  df_ct <- df_ct |>
    dplyr::filter(!is.na(logFC), is.finite(logFC),
                  !is.na(.data[[pcol]]), is.finite(.data[[pcol]])) |>
    dplyr::mutate(rank_stat = sign(logFC) * -log10(pmax(.data[[pcol]], 1e-300))) |>
    dplyr::arrange(dplyr::desc(rank_stat)) |>
    dplyr::group_by(gene) |>
    dplyr::slice_max(order_by = abs(rank_stat), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  
  stats <- df_ct$rank_stat
  names(stats) <- df_ct$gene
  
  
  PATHS_use <- lapply(pathways, function(g) intersect(g, names(stats)))
  lens <- vapply(PATHS_use, length, 1L)
  PATHS_use <- PATHS_use[lens >= minSize & lens <= maxSize]
  if (!length(PATHS_use)) return(tibble::tibble())
  
  
  if ("fgseaMultilevel" %in% getNamespaceExports("fgsea")) {
    fg <- fgsea::fgseaMultilevel(pathways = PATHS_use, stats = stats,
                                 minSize = minSize, maxSize = maxSize)
  } else {
    
    fg <- fgsea::fgsea(pathways = PATHS_use, stats = stats,
                       minSize = minSize, maxSize = maxSize)
  }
  
  tibble::as_tibble(fg) |>
    dplyr::select(pathway, size, NES, padj, leadingEdge)
}


## =======================
## 3) Run fgsea per cell type
##    (Deg already in memory)
## =======================
stopifnot(all(c("celltype","gene","logFC") %in% names(Deg)))
ct_list   <- split(Deg, Deg$celltype)
head(ct_list)

fgsea_list <- purrr::imap(ct_list, ~ run_fgsea_one(.x, PATHS) |>
                            dplyr::mutate(celltype = .y))
fgsea_tbl  <- dplyr::bind_rows(fgsea_list)
table (fgsea_tbl$celltype)

############%%%%%%%%% for only show 3 clusters %%%%%%%%%#############
fgsea_tbl_1 <- fgsea_tbl %>%
  dplyr::filter(celltype %in% c(
    "Kupffer cells (activated, lipid-handling)",
    "Pro-inflammatory macrophages",
    "Plasma B cells"
  ))
fgsea_tbl <- fgsea_tbl_1

dir.create("GSEA_MetaPW_by_celltype_3Clusters", showWarnings = FALSE)
readr::write_csv(fgsea_tbl, "fgsea_all_celltypes_MetaPW_3Clusters.csv")
head(fgsea_tbl)



## =======================
## 4.2) Plot data (all cell type together)
## =======================

# ---- packages ----
pkgs <- c("dplyr","tidyr","stringr","forcats","ggplot2","purrr")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Parameters ----
alpha_padj <- 0.05   
top_k_up   <- 8     
top_k_down <- 8    
wrap_width <- 40    

# fgsea_tbl: including pathway, NES, padj, celltype
stopifnot(all(c("pathway","NES","padj","celltype") %in% names(fgsea_tbl)))

fg_sig <- fgsea_tbl %>%
  filter(!is.na(padj), padj < alpha_padj)

top_by_ct <- fg_sig %>%
  group_by(celltype) %>%
  
  arrange(desc(NES), .by_group = TRUE) %>%
  mutate(rank_up = ifelse(NES > 0, row_number(), NA_integer_)) %>%
  arrange(abs(NES), .by_group = TRUE) %>% # 备用排序
  ungroup() %>%
  group_by(celltype) %>%
  
  slice_max(order_by = if_else(NES > 0, abs(NES), -Inf),
            n = top_k_up, with_ties = FALSE) %>%
  ungroup() %>%
  
  bind_rows(
    fg_sig %>%
      group_by(celltype) %>%
      slice_max(order_by = if_else(NES < 0, abs(NES), -Inf),
                n = top_k_down, with_ties = FALSE) %>%
      ungroup()
  ) %>%
  distinct(celltype, pathway, .keep_all = TRUE) %>%   # 去重
  mutate(
    negLog10Padj = -log10(padj),
    dir = if_else(NES >= 0, "Up (NES>0)", "Down (NES<0)"),
    pathway_wrapped = stringr::str_wrap(pathway, width = wrap_width)
  )


top_by_ct <- top_by_ct %>%
  arrange(desc(NES)) %>%
  mutate(pathway_wrapped = forcats::fct_reorder(pathway_wrapped, NES))

# ----  dotplot ----
p <- ggplot(top_by_ct,
            aes(x = celltype, y = pathway_wrapped)) +
  geom_point(aes(size = negLog10Padj, color = NES), alpha = 0.9) +
  scale_size_continuous(name = expression(-log[10](padj))) +
  scale_color_gradient2(low = "#2C7BB6", mid = "grey80", high = "#D7191C",
                        midpoint = 0, name = "NES") +
  labs(x = "Cell type", y = NULL,
       title = "Top regulated pathways",
       subtitle = paste0("Filtered by padj < ", alpha_padj,
                         "; size = -log10(padj), color = NES")) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 9, lineheight = 0.7),
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p)

## save figure
graph2tif(x = NULL, file='Enrichment_MetabolicActivity_myGeneset_3Clusters', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 6, bg = "transparent")

graph2svg(x = NULL, file='Enrichment_MetabolicActivity_myGeneset_3Clusters', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 6, bg = "transparent")


ggsave("Enrichment_MetabolicActivity_myGeneset_3Clusters.svg", p,
       width = 4.5, height = 6, device = svglite::svglite)








##### Another option: GO term
## =======================
## 0) Packages
## =======================
pkgs <- c("dplyr","tidyr","tibble","purrr","fgsea","pheatmap","stringr","msigdbr")
inst <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(inst)) install.packages(inst, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

## =======================
## A) Load DEG table
## =======================
Deg <- read.csv("ALL_celltypes_DEG_combined.csv", header = TRUE)
stopifnot(all(c("celltype","gene","logFC") %in% names(Deg)))

## =======================
## B) Build GO pathways (replace the MetaPW .gmt step)
##    Choose one of the following:
## =======================


## ---- Option 2: union of BP/CC/MF (uncomment if you want all three) ----
go_bp <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") |>
  dplyr::select(term = gs_name, gene = gene_symbol)
go_cc <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:CC") |>
  dplyr::select(term = gs_name, gene = gene_symbol)
go_mf <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:MF") |>
  dplyr::select(term = gs_name, gene = gene_symbol)

go_all <- dplyr::bind_rows(go_bp, go_cc, go_mf)
use_sets <- go_all



# (optional) tidy names a bit for plotting
tidy_go_name <- function(x) {
  x <- as.character(x)
  x <- gsub('[“”"]', '', x, fixed = FALSE)  
  x <- gsub('\\s+', ' ', x, perl = TRUE)    
  x <- trimws(x)
  x
}

use_sets$term <- tidy_go_name(use_sets$term)

# Build PATHS list the same way as before
PATHS_raw <- split(use_sets$gene, use_sets$term)
PATHS     <- lapply(PATHS_raw, function(v) unique(v[!is.na(v)]))

## =======================
## C) fgsea helper (unchanged)
## =======================
run_fgsea_one <- function(df_ct, pathways = PATHS,
                          minSize = 5, maxSize = 5000) {
  pcol <- if ("PValue" %in% names(df_ct)) "PValue" else "FDR"
  df_ct <- df_ct |>
    dplyr::filter(!is.na(logFC), is.finite(logFC),
                  !is.na(.data[[pcol]]), is.finite(.data[[pcol]])) |>
    dplyr::mutate(rank_stat = sign(logFC) * -log10(pmax(.data[[pcol]], 1e-300))) |>
    dplyr::arrange(dplyr::desc(rank_stat)) |>
    dplyr::group_by(gene) |>
    dplyr::slice_max(order_by = abs(rank_stat), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  
  stats <- df_ct$rank_stat; names(stats) <- df_ct$gene
  
  pw_use <- lapply(pathways, function(g) intersect(g, names(stats)))
  lens   <- vapply(pw_use, length, 1L)
  pw_use <- pw_use[lens >= minSize & lens <= maxSize]
  if (!length(pw_use)) return(tibble::tibble())
  
  fg <- if ("fgseaMultilevel" %in% getNamespaceExports("fgsea")) {
    fgsea::fgseaMultilevel(pathways = pw_use, stats = stats,
                           minSize = minSize, maxSize = maxSize)
  } else {
    fgsea::fgsea(pathways = pw_use, stats = stats,
                 minSize = minSize, maxSize = maxSize)
  }
  
  tibble::as_tibble(fg) |>
    dplyr::select(pathway, size, NES, padj, leadingEdge)
}

## =======================
## D) Run fgsea per cell type (unchanged)
## =======================
ct_list   <- split(Deg, Deg$celltype)
fgsea_list <- purrr::imap(ct_list, ~ run_fgsea_one(.x, PATHS) |>
                            dplyr::mutate(celltype = .y))
fgsea_tbl  <- dplyr::bind_rows(fgsea_list)

dir.create("GSEA_GO_by_celltype", showWarnings = FALSE)
readr::write_csv(fgsea_tbl, "fgsea_all_celltypes_GO.csv")
getwd()
head(fgsea_tbl)



## =======================
## 4.1) Plot data (cell type separately)
## =======================
pkgs <- c("dplyr","tidyr","stringr","forcats","ggplot2","purrr")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# Parameters
alpha_padj <- 0.05  
top_k_up   <- 8     
top_k_down <- 8     
wrap_width <- 40     
outdir     <- "GSEA_GO_by_celltype/dotplots_ct"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

stopifnot(all(c("pathway","NES","padj","celltype") %in% names(fgsea_tbl)))

plot_fgsea_dot_ct <- function(fgsea_tbl, ct,
                              alpha_padj = 0.05, top_k_up = 8, top_k_down = 8,
                              wrap_width = 40, outdir = ".") {
  
  df <- fgsea_tbl %>%
    dplyr::filter(celltype == ct, !is.na(padj))
  
  
  sig <- df %>% dplyr::filter(padj < alpha_padj)
  
  
  up   <- sig %>% dplyr::filter(NES > 0)  %>%
    dplyr::slice_max(order_by = NES,      n = top_k_up,   with_ties = FALSE)
  down <- sig %>% dplyr::filter(NES < 0)  %>%
    dplyr::slice_min(order_by = NES,      n = top_k_down, with_ties = FALSE)
  
  keep <- dplyr::bind_rows(up, down) %>% dplyr::distinct(pathway, .keep_all = TRUE)
  
  if (nrow(keep) == 0) {
    keep <- df %>% dplyr::slice_max(order_by = abs(NES), n = top_k_up + top_k_down, with_ties = FALSE)
  }
  
  keep <- keep %>%
    dplyr::mutate(
      negLog10Padj   = -log10(pmax(padj, 1e-300)),
      pathway_wrapped = stringr::str_wrap(pathway, width = wrap_width)
    ) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(pathway_wrapped = forcats::fct_reorder(pathway_wrapped, NES))
  
  
  p <- ggplot(keep, aes(x = 1, y = pathway_wrapped)) +
    geom_point(aes(size = negLog10Padj, color = NES), alpha = 0.9) +
    scale_size_continuous(name = expression(-log[10](padj))) +
    scale_color_gradient2(low = "#2C7BB6", mid = "grey80", high = "#D7191C",
                          midpoint = 0, name = "NES") +
    scale_x_continuous(labels = NULL, breaks = NULL) +
    labs(
      title = ct,
      subtitle = sprintf("padj < %.2f | Top %d up / %d down", alpha_padj, top_k_up, top_k_down),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9, lineheight = 0.65),
      panel.grid.minor = element_blank()
    )
  
 
  h <- max(4, 0.28 * nrow(keep) + 1.5)
  fn_safe <- gsub("[^A-Za-z0-9]+", "_", ct)
  
  ggsave(file.path(outdir, paste0("fgsea_dot_", fn_safe, ".png")), p, width = 10, height = h, dpi = 300)
  
  p
}

plots_by_ct <- purrr::map(unique(fgsea_tbl$celltype),
                          ~ plot_fgsea_dot_ct(fgsea_tbl, ct = .x,
                                              alpha_padj = alpha_padj,
                                              top_k_up = top_k_up,
                                              top_k_down = top_k_down,
                                              wrap_width = wrap_width,
                                              outdir = outdir))


## =======================
## 4.2) Plot data (all cell type together)
## =======================

# ---- packages ----
pkgs <- c("dplyr","tidyr","stringr","forcats","ggplot2","purrr")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Parameters ----
alpha_padj <- 0.01   
top_k_up   <- 5    
top_k_down <- 5     
wrap_width <- 40     

stopifnot(all(c("pathway","NES","padj","celltype") %in% names(fgsea_tbl)))

fg_sig <- fgsea_tbl %>%
  filter(!is.na(padj), padj < alpha_padj)

top_by_ct <- fg_sig %>%
  group_by(celltype) %>%
  
  arrange(desc(NES), .by_group = TRUE) %>%
  mutate(rank_up = ifelse(NES > 0, row_number(), NA_integer_)) %>%
  arrange(abs(NES), .by_group = TRUE) %>% # 备用排序
  ungroup() %>%
  group_by(celltype) %>%
  
  slice_max(order_by = if_else(NES > 0, abs(NES), -Inf),
            n = top_k_up, with_ties = FALSE) %>%
  ungroup() %>%
  
  bind_rows(
    fg_sig %>%
      group_by(celltype) %>%
      slice_max(order_by = if_else(NES < 0, abs(NES), -Inf),
                n = top_k_down, with_ties = FALSE) %>%
      ungroup()
  ) %>%
  distinct(celltype, pathway, .keep_all = TRUE) %>%   # 去重
  mutate(
    negLog10Padj = -log10(padj),
    dir = if_else(NES >= 0, "Up (NES>0)", "Down (NES<0)"),
    pathway_wrapped = stringr::str_wrap(pathway, width = wrap_width)
  )

top_by_ct <- top_by_ct %>%
  arrange(desc(NES)) %>%
  mutate(pathway_wrapped = forcats::fct_reorder(pathway_wrapped, NES))

# ---- dotplot ----
p <- ggplot(top_by_ct,
            aes(x = celltype, y = pathway_wrapped)) +
  geom_point(aes(size = negLog10Padj, color = NES), alpha = 0.9) +
  scale_size_continuous(name = expression(-log[10](padj))) +
  scale_color_gradient2(low = "#2C7BB6", mid = "grey80", high = "#D7191C",
                        midpoint = 0, name = "NES") +
  labs(x = "Cell type", y = NULL,
       title = "Top regulated pathways per cell type",
       subtitle = paste0("Filtered by padj < ", alpha_padj,
                         "; size = -log10(padj), color = NES")) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9, lineheight = 0.65),
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p)

## save figure
graph2tif(x = NULL, file='Enrichment_MetabolicActivityxxxxx', font = "Arial", cairo = TRUE,   
          width = 9, height = 9.5, bg = "transparent")










########################### Common genes across all cell types ######################
## Load data
DegCutoff <- read.csv("DEG_summary_by_celltype_FDRabsLogFC.csv", header = TRUE)
head(DegCutoff)

Deg <- read.csv("ALL_celltypes_DEG_combined.csv", header = TRUE)

Deg <- Deg %>%
  mutate(celltype = ifelse(
    celltype == "Inflammatory macrophages (lipid or fibrosis-associated)",
    "Activated HSCs/fibroblasts",
    celltype
  ))
write.csv(Deg, "ALL_celltypes_DEG_combined.csv")

# Check
table(Deg$celltype)

head(Deg)
table(Deg$celltype)
obj <- readRDS("liver_merge_integr_Ann1.rds")
table(obj$celltype)

## --- packages ---
pkgs <- c("dplyr","tibble","tidyr","ggplot2","Seurat")
inst <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(inst)) install.packages(inst, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))


## --- inputs ---
# DEG columns: logFC, FDR, gene, celltype (as in your screenshot)
Deg <- read.csv("ALL_celltypes_DEG_combined.csv", header = TRUE)
table(Deg$celltype)
obj <- readRDS("liver_merge_integr_Ann1.rds")
head(obj@meta.data)
DefaultAssay(obj)
#DefaultAssay(obj) <- if ("SCT" %in% Assays(obj)) "SCT" else DefaultAssay(obj)
Idents(obj) <- "celltype"   # make sure this meta column exists

alpha   <- 0.05
lfc_min <- 0.5
min_cts <- 10   # ≥ 6 cell types

## 1) Significant pairs and common genes -----------------------------------
Deg_sig <- Deg %>%
  mutate(
    gene     = trimws(as.character(gene)),
    celltype = trimws(as.character(celltype))
  ) %>%
  filter(!is.na(FDR), FDR < alpha, !is.na(logFC), abs(logFC) > lfc_min)

# gene frequency across distinct cell types
# gene_freq <- Deg_sig %>% distinct(celltype, gene) %>% count(gene, name = "n_ct")

common_genes_tbl <- gene_freq %>%
  filter(n_ct >= min_cts) %>%
  left_join(
    Deg_sig %>% group_by(gene) %>% summarise(med_abs_lfc = median(abs(logFC))), 
    by = "gene"
  ) %>%
  arrange(desc(n_ct), desc(med_abs_lfc))

gene_order <- common_genes_tbl$gene
if (length(gene_order) == 0) stop("No common genes (≥6 cell types) at FDR<0.05 & |logFC|>0.5.")


# --- robust key normalizer ---
norm_key <- function(x) {
  x %>%
    as.character() %>%
    str_replace("_\\d+$", "") %>%               # drop trailing "_number"
    str_replace("^\\s*\\d+\\s*[:\\-\\.]?\\s*", "") %>%  # drop leading "12:" / "12 -" / "12 " etc.
    str_replace("^\\s*[A-Za-z]\\s+", "") %>%    # drop leading single-letter prefixes like "g "
    str_replace_all("[\u2013\u2014]", "-") %>%  # normalize dashes
    str_replace_all("[βΒ]", "beta") %>%         # normalize Greek beta
    str_replace_all("/", " ") %>%               # slashes -> space
    str_replace_all("[[:punct:]]", " ") %>%     # punctuation -> space
    str_squish() %>% tolower()
}

# --- rebuild avg_long ---
avg_mat_list <- AverageExpression(
  obj, features = gene_order, group.by = "celltype",
  assays = DefaultAssay(obj), slot = "data"
)
avg_mat <- avg_mat_list[[DefaultAssay(obj)]]

avg_long <- avg_mat %>%
  as.data.frame() %>% rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "id_pretty", values_to = "avg.exp") %>%
  mutate(
    gene      = trimws(as.character(gene)),
    id_pretty = trimws(as.character(id_pretty)),
    id_key    = norm_key(id_pretty)
  )
avg_long <- avg_long %>%
  mutate(id_pretty = str_remove(id_pretty, "^g\\s+"))

## ---- 1) set the exact order you want (six blocks) ----
gene_groups <- list(
  `circadian clock/TF`      = c("Dbp","Arntl","Tsc22d3", 'Usp2',"Bcl6", "Zbtb16"),
  `ER/mito function`       = c("Hyou1","Manf","Sigmar1","Clpx","Coq8a","Rsn7"),
  `Lipogenesis/droplet`    = c("Acly","Fasn","Lpin1","Pcyt1a", "Chka","Obsbp3", "Ppp1r3c"),
  `Sterol/BA`       = c('Abcg5',"Abcg8","Slco1a1", "Slc10a2","Cyp7a1"),
  `Glu/PL`   = c("Gck","Slc39a10","Aldh1a1","Depp1"),
  `AA/xeno`          = c("Tat","Slc7a2","Cyp2c70")
)
gene_levels <- unlist(gene_groups, use.names = FALSE)
lookup <- stack(gene_groups)                  # values=gene, ind=group name
colnames(lookup) <- c("gene","module")


# --- rebuild sig_pairs (significant pairs from your DEG) ---
sig_pairs <- Deg_sig %>%                       # from your script
  filter(gene %in% gene_levels) %>%
  transmute(
    id_pretty_deg = celltype,
    id_key        = norm_key(celltype),
    gene          = trimws(gene),
    direction     = ifelse(logFC > 0, "Up", "Down"),
    logFC         = logFC
  ) %>% distinct()

# --- diagnostics ---
ct_seurat_keys <- avg_long %>% distinct(id_key, id_pretty) %>% arrange(id_pretty)
ct_deg_keys    <- sig_pairs  %>% distinct(id_key, id_pretty_deg) %>% arrange(id_pretty_deg)
head(ct_seurat_keys)
head(ct_deg_keys)

message("# unique Seurat celltypes: ", nrow(ct_seurat_keys))
message("# unique DEG celltypes:    ", nrow(ct_deg_keys))
message("# exact key overlap:       ", length(intersect(ct_seurat_keys$id_key, ct_deg_keys$id_key)))




# --- join & plot ---
plot_df <- sig_pairs %>%
  inner_join(avg_long, by = c("id_key","gene")) %>%   # joins Seurat celltype labels
  left_join(lookup, by = "gene") %>%                  # add the module for faceting
  mutate(
    gene      = factor(gene, levels = gene_levels),   # exact x-order
    direction = factor(direction, levels = c("Down","Up")),
    size_val  = pmax(0, abs(logFC))                   # dot size = |log2FC|
  )

message("# genes to plot: ", length(levels(plot_df$gene)))
message("# celltypes to plot: ", length(unique(plot_df$id_pretty)))
message("# significant pairs: ", nrow(plot_df))



## ---- 3) draw: size = |log2FC|; optional faceting by module ----
lfc_cap <- 2.5   # (optional) cap sizes so very large effects don’t dominate
plot_df$size_val <- pmin(plot_df$size_val, lfc_cap)


if (nrow(plot_df) == 0) {
  unmapped <- anti_join(ct_deg_keys, ct_seurat_keys, by = "id_key")
  as_tibble(unmapped) %>% print(n = 99)
  stop("Still no matches — check 'unmapped' above (after removing the 'g ' prefix).")
}
table(plot_df$gene)

g <- ggplot(plot_df, aes(x = gene, y = id_pretty)) +
  geom_point(aes(size = size_val, color = direction), alpha = 0.9) +
  scale_size_continuous(name = expression("|log"[2]*"FC|"),
                        range = c(1.5, 6), limits = c(0, lfc_cap)) +
  scale_color_manual(values = c(Down = "#2C7BB6", Up = "#D7191C"),
                     name = "Direction", drop = FALSE) +
  facet_grid(~ module, scales = "free_x", space = "free_x") +   # remove this line if you don’t want blocks
  labs(
    title = sprintf("Common genes (≥%d cell types) at FDR<%.2f & |logFC|>%.1f",
                    min_cts, alpha, lfc_min),
    subtitle = "Up: higher expr in GDF15; dot size = |log2FC|",
    x = "Gene", y = "Cell type"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 9),
    panel.grid.minor = element_blank()
  )

print(g)

## save figure
graph2tif(x = NULL, file='CommonGene_clock', font = "Arial", cairo = TRUE,   
          width = 12, height = 6, bg = "transparent")










######## Enrichment analysis for ADRA, ADRB, CREB, and GR signaling ##########

## Load the data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1


### Metabolic analysis 
## This part includes following:
# 1. Metabolic pathway activity
# 2. Metabolic interactions
# 3. Metabolites analysis_Flux balance analysis (FBA)
#### Compare metabolic activities between GDF15 vs CR
## We test activity using 85 metabolism pathways gene sets
BiocManager::install("GSVA")
BiocManager::install("clusterProfiler")
devtools::install_github("ncborcherding/escape")
BiocManager::install("dittoSeq")
BiocManager::install("org.Mm.eg.db")


## load library
library(GSVA)
library(msigdbr)
library(escape)
library(GSEABase)
library(dittoSeq)
library(ggplot2)
library(clusterProfiler)
library(GSEABase)
library(org.Dm.eg.db)
library(presto)
library(msigdbr)
library(doMC)
library(dplyr)
library(Seurat)
library(irGSEA)

library(org.Mm.eg.db)

# Load geneset pathways
library(clusterProfiler)
MetaPW = read.gmt("ADRB_A_CREB_GR_Mouse_Liver.gmt")
str(MetaPW)
table(MetaPW$term)
head(MetaPW)

#load DEG from Pseudobulking analysis
Deg = read.csv("ALL_celltypes_DEG_combined.csv", header=TRUE)
head(Deg)


## =======================
## 0) Packages
## =======================
pkgs <- c("dplyr","tidyr","tibble","purrr","fgsea","pheatmap","stringr")
inst <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(inst)) install.packages(inst, repos="https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))


## =======================
## 1) Build pathway list from MetaPW
##    (MetaPW already in memory: term, gene)
## =======================
# clean a bit
MetaPW <- MetaPW |>
  dplyr::mutate(term = stringr::str_trim(term),
                gene = stringr::str_trim(gene)) |>
  dplyr::filter(nzchar(term), nzchar(gene))

# list: names = pathway terms, values = character vector of genes
PATHS_raw <- split(MetaPW$gene, MetaPW$term)
# (optional) drop duplicated genes within a term
PATHS <- lapply(PATHS_raw, function(v) unique(v[!is.na(v)]))


## =======================
## 2) Helper to run fgsea for ONE cell type
## =======================
rm(list = c("run_fgsea_one"))  

run_fgsea_one <- function(df_ct, pathways = PATHS,
                          minSize = 5, maxSize = 5000) {
  
  pcol <- if ("PValue" %in% names(df_ct)) "PValue" else "FDR"
  
  
  df_ct <- df_ct |>
    dplyr::filter(!is.na(logFC), is.finite(logFC),
                  !is.na(.data[[pcol]]), is.finite(.data[[pcol]])) |>
    dplyr::mutate(rank_stat = sign(logFC) * -log10(pmax(.data[[pcol]], 1e-300))) |>
    dplyr::arrange(dplyr::desc(rank_stat)) |>
    dplyr::group_by(gene) |>
    dplyr::slice_max(order_by = abs(rank_stat), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  
  stats <- df_ct$rank_stat
  names(stats) <- df_ct$gene
  
 
  PATHS_use <- lapply(pathways, function(g) intersect(g, names(stats)))
  lens <- vapply(PATHS_use, length, 1L)
  PATHS_use <- PATHS_use[lens >= minSize & lens <= maxSize]
  if (!length(PATHS_use)) return(tibble::tibble())
  
 
  if ("fgseaMultilevel" %in% getNamespaceExports("fgsea")) {
    fg <- fgsea::fgseaMultilevel(pathways = PATHS_use, stats = stats,
                                 minSize = minSize, maxSize = maxSize)
  } else {
  
    fg <- fgsea::fgsea(pathways = PATHS_use, stats = stats,
                       minSize = minSize, maxSize = maxSize)
  }
  
  tibble::as_tibble(fg) |>
    dplyr::select(pathway, size, NES, padj, leadingEdge)
}


## =======================
## 3) Run fgsea per cell type
##    (Deg already in memory)
## =======================
stopifnot(all(c("celltype","gene","logFC") %in% names(Deg)))
ct_list   <- split(Deg, Deg$celltype)
head(ct_list)

fgsea_list <- purrr::imap(ct_list, ~ run_fgsea_one(.x, PATHS) |>
                            dplyr::mutate(celltype = .y))
fgsea_tbl  <- dplyr::bind_rows(fgsea_list)

dir.create("GSEA_ADR_GR_by_celltype", showWarnings = FALSE)
readr::write_csv(fgsea_tbl, "fgsea_all_celltypes_ADR_GR.csv")
head(fgsea_tbl)

## =======================
## 4) Plot data
## =======================
pkgs <- c("dplyr","tidyr","stringr","forcats","ggplot2","purrr")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- Parameter ----
alpha_padj <- 0.1
top_k_up   <- 8
top_k_down <- 8
wrap_width <- 40

stopifnot(all(c("pathway","NES","padj","celltype") %in% names(fgsea_tbl)))

# 1) Significant results used only for selecting candidate genes
fg_sig <- fgsea_tbl %>% dplyr::filter(!is.na(padj), padj < 0.9)

# 2) Select the top up- and downregulated pathways for each cell type
#    (pathways to be displayed)
sel_paths_by_ct <- fg_sig %>%
  dplyr::group_by(celltype) %>%
  dplyr::slice_max(order_by = if_else(NES > 0,  abs(NES), -Inf),
                   n = top_k_up, with_ties = FALSE) %>%
  dplyr::bind_rows(
    fg_sig %>%
      dplyr::group_by(celltype) %>%
      dplyr::slice_max(order_by = if_else(NES < 0,  abs(NES), -Inf),
                       n = top_k_down, with_ties = FALSE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(pathway)                

plot_df <- fgsea_tbl %>%
  dplyr::filter(pathway %in% sel_paths_by_ct$pathway) %>%  
  dplyr::mutate(
    negLog10Padj   = -log10(pmax(padj, 1e-300)),
    pathway_wrapped = stringr::str_wrap(pathway, width = wrap_width),
    sig_group      = ifelse(!is.na(padj) & padj < alpha_padj, "sig", "nonsig"),
    NES_plot       = dplyr::if_else(sig_group == "sig", NES, NA_real_)
  ) %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(pathway_wrapped = forcats::fct_reorder(pathway_wrapped, NES))
table(plot_df$celltype)


# =======================
# NEW: split into two plots by celltype (all other settings unchanged)
# =======================
keep_ct <- c("Plasma B cells",
             "Pro-inflammatory macrophages",
             "Kupffer cells (activated, lipid-handling)")

plot_df_keep  <- plot_df %>% dplyr::filter(celltype %in% keep_ct)
plot_df_other <- plot_df %>% dplyr::filter(!celltype %in% keep_ct)

# helper to build the same plot (parameters unchanged)
make_plot <- function(dat) {
  ggplot(dat, aes(x = celltype, y = pathway_wrapped)) +
    geom_point(
      data = subset(dat, sig_group == "nonsig"),
      aes(size = negLog10Padj),
      color = "grey70", alpha = 0.7
    ) +
    geom_point(
      data = subset(dat, sig_group == "sig"),
      aes(size = negLog10Padj, color = NES_plot),
      alpha = 0.9
    ) +
    scale_size_continuous(name = expression(-log[10](padj))) +
    scale_color_gradient2(
      low = "#2C7BB6", mid = "grey80", high = "#D7191C",
      midpoint = 0, name = "NES"
    ) +
    labs(
      x = "Cell type", y = NULL,
      title = "ADR_GR Pathway activity per cell type",
      subtitle = paste0("padj < ", alpha_padj, " colored by NES; padj ≥ ", alpha_padj, " in grey")
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 10, lineheight = 0.7),
      axis.text.x = element_text(size = 10, angle = 55, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

p_keep  <- make_plot(plot_df_keep)
p_other <- make_plot(plot_df_other)

print(p_keep)
print(p_other)


## save figure
# install.packages("gdtools", type = "binary")
library("export")
graph2tif(x = p_keep, file='Enrichment_MetabolicActivity_ADRA_B_CREB_GR_Infla', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 5, bg = "transparent")

graph2svg(x = p_keep, file='Enrichment_MetabolicActivity_ADRA_B_CREB_GR_Infla', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 5, bg = "transparent")


graph2tif(x = p_other, file='Enrichment_MetabolicActivity_ADRA_B_CREB_GR_nonInfla', font = "Arial", cairo = TRUE,   
          width = 6, height = 5, bg = "transparent")

graph2svg(x = p_other, file='Enrichment_MetabolicActivity_ADRA_B_CREB_GR_nonInfla', font = "Arial", cairo = TRUE,   
          width = 6, height = 5, bg = "transparent")


## Export svg figure
# install.packages("svglite")
library(svglite)

svglite("Enrichment_MetabolicActivity_ADRA_B_CREB_GR_Infla_2.svg", width = 4.5, height = 5.5)
print(p_keep)
dev.off()

svglite("Enrichment_MetabolicActivity_ADRA_B_CREB_GR_nonInfla_2.svg", width = 8, height = 5.5)
print(p_other)
dev.off()








################### TBC: Correlation, and Compass in immune cells################

#### Do correlation analysis in 3 celltypes '15: Plasma B cells', "13: Pro-inflammatory macrophages" and "6: Kupffer cells" ####
# ## Load the data
# liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
# head(liver.mergeIntegrAnno1)
# obj = liver.mergeIntegrAnno1


### Metabolic analysis 
## This part includes following:
# 1. Metabolic pathway activity
# 2. Metabolic interactions
# 3. Metabolites analysis_Flux balance analysis (FBA)
#### Compare metabolic activities between GDF15 vs CR
## We test activity using 85 metabolism pathways gene sets
BiocManager::install("GSVA")
BiocManager::install("clusterProfiler")
devtools::install_github("ncborcherding/escape")
BiocManager::install("dittoSeq")
BiocManager::install("org.Mm.eg.db")
# install irGSEA packages from Github
if (!requireNamespace("irGSEA", quietly = TRUE)) { 
  devtools::install_github("chuiqin/irGSEA", force =T)
}


## load library
library(GSVA)
library(msigdbr)
library(escape)
library(GSEABase)
library(dittoSeq)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GSEABase)
library(presto)
library(msigdbr)
library(dplyr)
library(Seurat)

library(irGSEA)
library(org.Dm.eg.db)
library(doMC)

# Load geneset pathways
MetaPW = read.gmt("ImmuneMetabolicPathwaysCor.gmt")
head(MetaPW)
str(MetaPW)
table(MetaPW$term)

MetaPW2 <- list()
for (i in rownames(table(MetaPW$term))){
  print(i)
  MetaPW2[[i]] <- MetaPW[MetaPW$term==i,]$gene
  
}
head(MetaPW2)
gene_sets = MetaPW2

#load DEG from Pseudobulking analysis
Deg = read.csv("ALL_celltypes_DEG_combined.csv", header=TRUE)
head(Deg)

# Load data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
seurat_obj = liver.mergeIntegrAnno1
head(seurat_obj)


## ------------------------------------------------------------
## ===== 0) Setup: Gene Sets, Objects, and Utility Functions =====
## ------------------------------------------------------------
library(Seurat)
library(GSVA)
library(clusterProfiler)
library(ggplot2)
library(corrplot)
library(svglite)

# read GMT -> list
MetaPW = read.gmt("ImmuneMetabolicPathwaysCor.gmt")
gene_sets <- split(MetaPW$gene, MetaPW$term)

# read Seurat object & set Assay
seurat_obj <- readRDS("liver_merge_integr_Ann1.rds")
#if ("SCT" %in% Assays(seurat_obj)) DefaultAssay(seurat_obj) <- "SCT"

# Function to calculate the p-value matrix (Spearman correlation)
cor_mtest <- function(mat, method = "spearman") {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA_real_, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      tmp <- suppressWarnings(cor.test(mat[, i], mat[, j], method = method))
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

ct_levels_raw <- unique(as.character(seurat_obj$celltype))
ct_nums <- suppressWarnings(as.integer(sub("^\\s*(\\d+):.*$", "\\1", ct_levels_raw)))
ct_levels_ordered <- ct_levels_raw[order(ct_nums, na.last = TRUE)]



# ---- Place at the beginning of the script:
#      Compatibility wrapper for ssGSEA across different GSVA versions ----
run_ssgsea <- function(expr, gene_sets, min.sz = 2L) {
  expr <- as.matrix(expr)
  
  
  gs_in <- lapply(gene_sets, function(gs) intersect(rownames(expr), unique(gs)))
  keep  <- vapply(gs_in, length, 1L) >= min.sz
  gs_in <- gs_in[keep]
  if (length(gs_in) == 0) stop("No usable gene sets after intersecting with expression rows.")
  
  exps <- try(getNamespaceExports("GSVA"), silent = TRUE)
  if (inherits(exps, "try-error")) stop("GSVA not loaded or not installed.")
  
  
  if ("ssgsea" %in% exps) {
    return(GSVA::ssgsea(expr, gs_in, ssgsea.norm = TRUE))
  }
  
  
  if ("ssgseaParam" %in% exps) {
    
    param_try <- try(GSVA::ssgseaParam(expr, gs_in), silent = TRUE)
    if (!inherits(param_try, "try-error")) {
      res_try <- try(GSVA::gsva(param_try), silent = TRUE)
      if (!inherits(res_try, "try-error")) return(res_try)
    }
  }
  
  
  if ("gsva" %in% exps) {
    res_try2 <- try(GSVA::gsva(expr, gs_in, method = "ssgsea", ssgsea.norm = TRUE), silent = TRUE)
    if (!inherits(res_try2, "try-error")) return(res_try2)
  }
  
  stop("Could not find a compatible GSVA/ssGSEA entry point for your GSVA version.")
}


## ------------------------------------------------------------
## ===== 1) Analyze Each Cell Type in a Loop =====
## ------------------------------------------------------------
out_dir <- "per_celltype_ssGSEA_corr_ImmuMetab"
if (!dir.exists(out_dir)) dir.create(out_dir)

for (ct in ct_levels_ordered) {
  message(">>> Processing celltype: ", ct)
  
  
  ct_obj <- subset(seurat_obj, subset = celltype == ct)
  if (ncol(ct_obj) < 5) {  
    message("   (skip: too few cells: ", ncol(ct_obj), ")")
    next
  }
  
 
  expr <- as.matrix(GetAssayData(ct_obj, slot = "data", assay = DefaultAssay(ct_obj)))
  
  
  ssgsea_scores <- run_ssgsea(expr, gene_sets, min.sz = 2L)
  
  
  scr_ssgsea <- t(ssgsea_scores)
  
  
  cor_mat <- cor(scr_ssgsea, method = "spearman")
  p_mat   <- cor_mtest(scr_ssgsea, method = "spearman")
  
  
  safe_name <- gsub("[^A-Za-z0-9_\\-]+", "_", ct)   
  write.csv(cor_mat, file = file.path(out_dir, paste0("corr_", safe_name, ".csv")))
  write.csv(p_mat,   file = file.path(out_dir, paste0("pval_", safe_name, ".csv")))
  
  
  svglite(file.path(out_dir, paste0("corrplot_", safe_name, ".svg")),
      width = 8, height = 8)
  
  par(mfrow = c(1,1))
  #        
  par(mar = c(10, 0, 3, 6))
  
  # --- 1st layer: correlation heatmap---
  corrplot(
    cor_mat,
    method   = "color",
    type     = "upper",
    diag     = FALSE,
    tl.pos   = "lt",
    tl.col   = "black",
    tl.cex   = 1.5, 
    tl.srt   = 90,
    tl.offset = 0.8,
    col      = colorRampPalette(c("blue", "white", "red"))(200),
    col.lim  = c(-max(abs(cor_mat), na.rm = TRUE),
                 max(abs(cor_mat), na.rm = TRUE)),
    cl.pos   = "b",
    cl.cex   = 1.5,
    cl.ratio = 0.18,
    cl.length = 5,
    
    addCoef.col = NULL,
    number.cex  = 0,  
    is.corr     = TRUE,
    addgrid.col = NA,
    mar         = c(0,0,0,0)
  )
  
  # ---------- 2) * matrix ----------
  stars_mat <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat),
                      dimnames = dimnames(p_mat))
  
  stars_mat[p_mat >= 0.05]                      <- "ns"
  stars_mat[p_mat < 0.05  & p_mat >= 0.01]      <- "*"
  stars_mat[p_mat < 0.01  & p_mat >= 0.001]     <- "**"
  stars_mat[p_mat < 0.001]                      <- "***"
  
  
  # ---------- 3) text() add * ----------
  n <- nrow(stars_mat)
  for (i in 1:n) {
    for (j in 1:n) {
      if (j > i && stars_mat[i, j] != "") {
        text(x = j,
             y = n - i + 1,
             labels = stars_mat[i, j],
             cex = 1.5,
             font = 2)
      }
    }
  }
  
  
  mtext(ct, side = 1, line = 10, cex = 1.5, font = 2) 
  
  dev.off()
  
  message("Finish, export to:", out_dir)
}







#------ python compass analysis- FBA - input data preparation -----------------------------
# For these cell types: '15: Plasma B cells', "13: Pro-inflammatory macrophages" and "6: Kupffer cells" 

# python compass analysis- FBA
library(Seurat)
library(Matrix)
update.packages("Matrix")

# loading data
seurat_obj <- readRDS("liver_merge_integr_Ann1.rds")
DimPlot(seurat_obj, reduction = "umap.harmony", group.by = "celltype.stim", label = TRUE)


# Downsample 300 cells per identity group (cell type) (don't do this step )
head(seurat_obj@meta.data)
# Idents(seurat_obj) <- "celltype.stim"
table(seurat_obj@meta.data$celltype.stim)
# seurat_downsampled <- subset(seurat_obj, downsample = 100)
# table(seurat_obj@meta.data$celltype.stim)


# Only choose cell type: '15: Plasma B cells', "13: Pro-inflammatory macrophages" and "6: Kupffer cells"
# the target cell types
target_cells <- c(
  "15: Plasma B cells",
  "13: Pro-inflammatory macrophages",
  "6: Kupffer cells (activated, lipid-handling)",
  "9: Activated HSCs/fibroblasts"
)

# subset Seurat object
seurat_subset <- subset(
  seurat_obj,
  subset = celltype %in% target_cells
)

# check results
table(seurat_subset$celltype)


# export matrix: change 'slot' as 'data'
compassFBA1 = seurat_subset
counts_matrix <- GetAssayData(compassFBA1, assay = "SCT", layer = "data")
class(counts_matrix)
counts_matrix <- as(counts_matrix, "dgCMatrix")

# export .mtx format
writeMM(counts_matrix, file = "expression.mtx")

# export row（gene）and column（cells/spot）as tsv files
write.table(rownames(counts_matrix), file = "genes.tsv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

write.table(colnames(counts_matrix), file = "sample_names.tsv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# extract and save metadata
metadata_df <- compassFBA1@meta.data
head(metadata_df)

write.csv(metadata_df, file = "cell_metadata.csv", quote = TRUE)











#################### cell-cell communication using ‘CellChat’ #############
install.packages('NMF')
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
install.packages("parallelly")
devtools::install_github("jinworks/CellChat")

update.packages(ask = FALSE, checkBuilt = TRUE)

install.packages("pak")
pak::pak(c("ggplot2","CellChat"))

library(CellChat)
library(patchwork)
library(Seurat)



## Load the data
liver.mergeIntegrAnno1 <- readRDS("liver_merge_integr_Ann1.rds")
head(liver.mergeIntegrAnno1)
obj = liver.mergeIntegrAnno1
Idents(obj)
levels(obj)

# split samples
Liv1 <- subset(obj, orig.ident == 'CR1')
table(Liv1$orig.ident)

Liv2 <- subset(obj, orig.ident == 'GDF15a')
table(Liv2$orig.ident)

### add data
seu1 = Liv1
seu2 = Liv2

# Prepare input data for CellChat analysis
data.input1 = Seurat::GetAssayData(seu1, layer = "data", assay = "SCT") # normalized data matrix
data.input2 = Seurat::GetAssayData(seu2, layer = "data", assay = "SCT") 

# define the meta data
# a column named `samples` should be provided for spatial transcriptomics analysis,
# which is useful for analyzing cell-cell communication by aggregating multiple samples/replicates.
# Of note, for comparison analysis across different conditions, 
# users still need to create a CellChat object seperately for each condition.
meta1 = data.frame(labels = Idents(seu1), samples = "CR") # manually create a dataframe consisting of the cell labels
meta2 = data.frame(labels = Idents(seu2), samples = "GDF15") 

# a factor level should be defined for the `meta$labels` and `meta$samples`
meta1$labels <- factor(meta1$labels, levels = levels(Idents(seu1)))
meta1$samples <- factor(meta1$samples, levels = "CR")

meta2$labels <- factor(meta2$labels, levels = levels(Idents(seu2)))
meta2$samples <- factor(meta2$samples, levels = "GDF15")
unique(meta1$labels) # check the cell labels
unique(meta2$labels) # check the cell labels

# load spatial transcriptomics information
# Spatial locations of spots from full (NOT high/low) resolution images are required. 
# For 10X Visium, this information is in `tissue_positions.csv`. 
spatial.locs1 = Seurat::GetTissueCoordinates(seu1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs2 = Seurat::GetTissueCoordinates(seu2, scale = NULL, cols = c("imagerow", "imagecol")) 

nrow(spatial.locs1)
length(colnames(data.input1))
dim(spatial.locs1)
dim(data.input1)


# Scale factors of spatial coordinates
# For 10X Visium, the conversion factor of converting spatial coordinates
# from Pixels to Micrometers can be computed as the ratio of the theoretical spot size (i.e., 65um)
# over the number of pixels that span the diameter of a theoretical spot size
# in the full-resolution image (i.e., 'spot_diameter_fullres' in pixels in the 'scalefactors_json.json' file).
scalefactors1 = jsonlite::fromJSON(txt = file.path("C:/BaiduSyncdisk/R analysis/Spatial transcriptomics/GDF15CR/Sam19/spatial", 'scalefactors_json.json'))
spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor1 = spot.size/scalefactors1$spot_diameter_fullres
spatial.factors1 = data.frame(ratio = conversion.factor1, tol = spot.size/2)


scalefactors2 = jsonlite::fromJSON(txt = file.path("C:/BaiduSyncdisk/R analysis/Spatial transcriptomics/GDF15CR/Sam10/spatial", 'scalefactors_json.json'))
conversion.factor2 = spot.size/scalefactors2$spot_diameter_fullres
spatial.factors2 = data.frame(ratio = conversion.factor2, tol = spot.size/2)


#### Create a CellChat object
head(spatial.locs1)
head(spatial.locs2)
spatial.locs1=spatial.locs1[,-3]
spatial.locs2=spatial.locs2[,-3]

cellchat1 <- createCellChat(object = data.input1, meta = meta1, group.by = "labels",
                            datatype = "spatial", coordinates = spatial.locs1, spatial.factors = spatial.factors1)
cellchat2 <- createCellChat(object = data.input2, meta = meta2, group.by = "labels",
                            datatype = "spatial", coordinates = spatial.locs2, spatial.factors = spatial.factors2)

cellchat1
cellchat2


#### Set the ligand-receptor interaction database
# Use the CellChatDB with metabolic signaling
##$$$$ for all pathways
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
cellchat1@DB <- CellChatDB.use #important
cellchat2@DB <- CellChatDB.use #important

unique(CellChatDB.use$interaction$annotation)
unique(CellChatDB.use$interaction$pathway_name)


#### Preprocessing the expression data for cell-cell communication analysis
# To infer the cell state-specific communications, we identify over-expressed ligands
# or receptors in one cell group and then identify over-expressed
# ligand-receptor interactions if either ligand or receptor is over-expressed.

# subset the expression data of signaling genes for saving computation cost
#devtools::install_github('immunogenomics/presto')

cellchat1 <- subsetData(cellchat1) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat1 <- identifyOverExpressedGenes(cellchat1)
cellchat1 <- identifyOverExpressedInteractions(cellchat1)

cellchat2 <- subsetData(cellchat2) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)


#### Part II: Inference of cell-cell communication network
# Compute the communication probability and infer cellular communication network
#ptm = Sys.time()


#For certain functions, each worker needs access to certain global variables.
#If these are larger than the default limit, you will see this error. 
#To get around this, you can set options(future.globals.maxSize = X),where X is the maximum allowed size in bytes.
#So to set it to 1GB, you would run options(future.globals.maxSize = 1000 * 1024^2). 
#Note that this will increase your RAM usage so set this number mindfully.
library(future)
workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = workers)   # Windows
options(future.globals.maxSize = 24 * 1024^3, future.rng.onMisuse = "ignore")

# Remove lowly expressed genes/interactions per group to increase speed
cellchat1 <- subsetData(cellchat1)   # uses defaults inside CellChat
cellchat2 <- subsetData(cellchat2)

cellchat1@images$coordinates <- as.matrix(cellchat1@images$coordinates)
cellchat2@images$coordinates <- as.matrix(cellchat2@images$coordinates)

cellchat1 <- computeCommunProb(cellchat1, type = "truncatedMean", trim = 0.1, 
                               distance.use = FALSE, interaction.range = 120, scale.distance = NULL,
                               contact.dependent = TRUE, contact.range = 50)

cellchat2 <- computeCommunProb(cellchat2, type = "truncatedMean", trim = 0.1, 
                               distance.use = FALSE, interaction.range = 120, scale.distance = NULL,
                               contact.dependent = TRUE, contact.range = 50)

#Users can filter out the cell-cell communication if there are only few cells in certain cell groups. 
#By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)

## Infer the cell-cell communication at a signaling pathway level
#CellChat computes the communication probability on signaling pathway level by summarizing the communication
#probabilities of all ligands-receptors interactions associated with each signaling pathway.

#NB: The inferred intercellular communication network of each ligand-receptor pair
#and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat2 <- computeCommunProbPathway(cellchat2)

##### Calculate the aggregated cell-cell communication network
#We can calculate the aggregated cell-cell communication network by
#counting the number of links or summarizing the communication probability.
cellchat1 <- aggregateNet(cellchat1)
cellchat2 <- aggregateNet(cellchat2)

## netAnalysis_computeCentrality
cellchat1 <- netAnalysis_computeCentrality(cellchat1, slot.name = "netP")
cellchat2 <- netAnalysis_computeCentrality(cellchat2, slot.name = "netP")

#We can also visualize the aggregated cell-cell communication network. 
#For example, showing the number of interactions or the total interaction strength (weights)
#between any two cell groups using circle plot or heatmap plot.
#ptm = Sys.time()


### merge data
object.list <- list(ctl = cellchat1, treat = cellchat2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


# Users can now export the merged CellChat object and the list of the two separate objects for later use
save(object.list, file = "cellchat_list1.RData")
save(cellchat, file = "cellchat_obj1.RData")

# load the data
load("cellchat_list.RData")
load("cellchat_obj.RData") 

# change celltype name (from "9: Inflammatory macrophages (lipid or fibrosis-associated) to "9: Activated HSCs/Myofibroblasts")
cellchat
object.list
class(cellchat@idents)
head(cellchat@idents)
cellchat1=cellchat

class(object.list$ctl@idents)
head(object.list$ctl@idents)


### change cellchat celltype names
# 1) Collapse to one factor vector
idv <- factor(unlist(cellchat@idents, use.names = TRUE))
# 2) Rename the level (handle both CR and GDF15 suffixes if present)
old  <- "9: Inflammatory macrophages \\(lipid or fibrosis-associated\\)"
lvl  <- levels(idv)
lvl  <- sub(paste0("^", old), "9: Activated HSCs/Myofibroblasts", lvl)
levels(idv) <- lvl
# 3) (optional) order levels by the numeric prefix "0:", "1:", ..., "18:"
levs <- levels(idv)
ord  <- order(as.integer(sub(":.*", "", levs)))
idv  <- factor(idv, levels = levs[ord])
# 4) Write back to the CellChat object
cellchat1@idents <- idv
# Check
levels(cellchat1@idents)[grepl("^9:", levels(cellchat1@idents))]
cellchat=cellchat1

### change object.list celltype names
object.list1=object.list
# Helper function to rename one CellChat object
rename_cellchat_idents <- function(obj, old, new){
  idv <- factor(unlist(obj@idents, use.names = TRUE))
  
  # update levels
  lv <- levels(idv)
  lv <- sub(paste0("^", old), new, lv)
  levels(idv) <- lv
  
  # write back
  obj@idents <- idv
  return(obj)
}
# Define old and new name
old_name <- "9: Inflammatory macrophages \\(lipid or fibrosis-associated\\)"
new_name <- "9: Activated HSCs/Myofibroblasts"
# Apply to both list elements
object.list1$ctl   <- rename_cellchat_idents(object.list1$ctl,   old_name, new_name)
object.list1$treat <- rename_cellchat_idents(object.list1$treat, old_name, new_name)
# Check
levels(object.list1$ctl@idents)[grepl("^9:", levels(object.list$ctl@idents))]
levels(object.list1$treat@idents)[grepl("^9:", levels(object.list$treat@idents))]
object.list=object.list1





## Part I: Identify altered interactions and cell populations
# Whether the cell-cell communication is enhanced or not?
# Compare the total number of interactions and interaction strength

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
graph2tif(x = NULL, file='CC-interactionNumStr', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 3.5, bg = "transparent")
graph2svg(x = NULL, file='CC-interactionNumStr', font = "Arial", cairo = TRUE,   
          width = 3.5, height = 3.5, bg = "transparent")


# Compare the number of interactions and interaction strength among different cell populations
# (A) Circle plot showing differential number of interactions or interaction strength
#    among different cell populations across two datasets
# red(or blue) colored edges represent increased (or decreased) signaling in the second dataset compared to the first one.
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
graph2tif(x = NULL, file='CC-compInterat', font = "Arial", cairo = TRUE,   
          width = 20, height = 5, bg = "transparent")

# (B) Heatmap showing differential number of interactions or interaction strength
#    among different cell populations across two datasets
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
graph2tif(x = NULL, file='CC-3-HeatmapIntera', font = "Arial", cairo = TRUE,   
          width = 16, height = 13, bg = "transparent")


# (C) Circle plot showing the number of interactions or interaction strength
#     among different cell populations across multiple datasets
# The above differential network analysis only works for pairwise datasets.
# If there are more datasets for comparison, CellChat can directly show results
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
graph2tif(x = NULL, file='CC-4-CirclePlot_ctlTrt', font = "Arial", cairo = TRUE,   
          width = 20, height = 10, bg = "transparent")



# (D) Circle plot To simplify the complicated network and gain insights by showing any two/3 cell types
group.cellType <- c(rep("6: Kupffer cells (activated, lipid-handling)", 4), 
                    rep("13: Pro-inflammatory macrophages", 4), rep("15: Plasma B cells", 4), rep("9: Activated HSCs/Myofibroblasts", 4))
group.cellType <- factor(group.cellType, levels = c("6: Kupffer cells (activated, lipid-handling)", 
                                                    "13: Pro-inflammatory macrophages", "15: Plasma B cells",
                                                    "9: Activated HSCs/Myofibroblasts"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
#table(obj$celltype)
graph2tif(x = NULL, file='CC-5-3clusterInter_ctlTrt1', font = "Arial", cairo = TRUE,   
          width = 10, height = 6, bg = "transparent")

graph2svg(x = NULL, file='CC-5-3clusterInter_ctlTrt1', font = "Arial", cairo = TRUE,   
          width = 10, height = 6, bg = "transparent")



# between any two cell types using circle plot.
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)


###### Compare the major sources and targets in a 2D space
## Identify cell populations with significant changes in
## sending or receiving signals between different datasets
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
graph2tif(x = NULL, file='CC-5-SigChangeCellType1', font = "Arial", cairo = TRUE,   
          width = 12, height = 8, bg = "transparent")


## keep both scatter plots on the same axes and label only those 3 cell types
library(ggrepel)
library(patchwork)

# 1) Compute link numbers and shared size range
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link))

# 2) Build the two base plots
gg <- lapply(names(object.list), function(nm) {
  netAnalysis_signalingRole_scatter(
    object.list[[nm]],
    title        = nm,
    weight.MinMax = weight.MinMax
  )
})


# helper: which column holds the cell-type labels in the plot data
find_label_col <- function(df) {
  cand <- c("labels","group","cellgroup","cell.type","celltype","idents")
  cand[cand %in% names(df)][1]
}
targets <- c(
  "6: Kupffer cells (activated, lipid-handling)",
  "13: Pro-inflammatory macrophages",
  "15: Plasma B cells",
  "9: Activated HSCs/Myofibroblasts"
)
gg_fixed <- lapply(gg, function(p) {
  # remove any existing text/label layers
  p$layers <- Filter(function(l) !inherits(l$geom, c("GeomText","GeomLabel","GeomTextRepel")), p$layers)
  
  df <- p$data
  lab_col <- find_label_col(df)
  
  # keep only target rows for labelling (use base subsetting, not .data pronoun)
  df_lab <- df[df[[lab_col]] %in% targets, , drop = FALSE]
  
  # build plot with fixed axes and only target labels
  p +
    geom_text_repel(
      data = df_lab,
      aes_string(label = lab_col),
      size = 4.8,
      max.overlaps = Inf
    ) +
    coord_cartesian(xlim = c(-1, 6), ylim = c(-0.5, 8)) +
    theme(legend.position = "none")
})
patchwork::wrap_plots(plots = gg_fixed)
graph2tif(x = NULL, file='CC-5-SigChangeCellType2_noLabl', font = "Arial", cairo = TRUE,   
          width = 6, height = 5, bg = "transparent")
graph2svg(x = NULL, file='CC-5-SigChangeCellType2_noLabl', font = "Arial", cairo = TRUE,   
          width = 6, height = 5, bg = "transparent")


### (B) Identify the signaling changes of specific cell populations
# "6: Kupffer cells (activated, lipid-handling)", 
# "13: Pro-inflammatory macrophages", 
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "6: Kupffer cells (activated, lipid-handling)_GDF15")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "13: Pro-inflammatory macrophages_GDF15")
patchwork::wrap_plots(plots = list(gg1,gg2))
graph2tif(x = NULL, file='CC-6-signalingChanges6_13clu', font = "Arial", cairo = TRUE,   
          width = 12, height = 7, bg = "transparent")

# "15: Plasma B cells"
# "9: Activated HSCs/Myofibroblasts"
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "15: Plasma B cells_GDF15")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "9: Activated HSCs/Myofibroblasts_GDF15")
patchwork::wrap_plots(plots = list(gg1,gg2))
graph2tif(x = NULL, file='CC-6-signalingChanges9_15clu', font = "Arial", cairo = TRUE,   
          width = 12, height = 7, bg = "transparent")



### Identify altered signaling with distinct interaction strength
# (A) Compare the overall information flow of each signaling pathway or ligand-receptor pair
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2

graph2tif(x = NULL, file='CC-7-SignalInter', font = "Arial", cairo = TRUE,   
          width = 12, height = 5, bg = "transparent")

### (B) Compare outgoing (or incoming) signaling patterns associated with each cell population
library(ComplexHeatmap)
# combining all the identified signaling pathways from different datasets 
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 9)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 9)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

graph2tif(x = NULL, file='CC-8-OutInc', font = "Arial", cairo = TRUE,   
          width = 16, height = 14, bg = "transparent")


#### Part III: Identify the up-gulated and down-regulated signaling ligand-receptor pairs
# Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat, sources.use = 4,  targets.use = c(1:2), comparison = c(1, 2), angle.x = 45)
graph2tif(x = NULL, file='CC-9-ligand-receptor', font = "Arial", cairo = TRUE,
          width = 14, height = 10, bg = "transparent") # I did not use this figure


#### the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs 
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in treatment", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in treatment", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
graph2tif(x = NULL, file='4-10-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,
          width = 8, height = 5, bg = "transparent")





# Chord diagram for pathways
pathways.show <- c("MIF", "SPP1") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
graph2tif(x = NULL, file='4-11-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,
          width = 14, height = 7, bg = "transparent")  






## Chord figures
levels(object.list[[1]]@idents) 
# [1] "Hepatocytes1"    "Hepatocytes2"    "Hepatocytes3"    "Hepatocytes4"    "Hepatocytes5"   
# [6] "Hepatocytes6"    "HSCs"            "Macrophages"     "Hepatocytes7"    "Hepatocytes8"   
# [11] "B Cells"         "Erythroid cells" "Immune Cells"
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 8, targets.use = 7, lab.cex = 0.5,
                       title.name = paste0("Signaling from macrophage - ", names(object.list)[i]))
}

graph2tif(x = NULL, file='4-12-LipidcomparismMeta_r1', font = "Arial", cairo = TRUE,
          width = 14, height = 7, bg = "transparent")  











#----------------------------------------------------
###### Integration: In order to work with multiple slices (with integration)
### we provide the merge function.
liver1SCT1 <- readRDS("liver1_final.rds")
liver2SCT1 <- readRDS("liver2_final.rds")
liver.merge <- merge(liver1SCT1, y= liver2SCT1, add.cell.ids = c("control1", "treat1"), project = "CCCmeta")
levels(liver.merge)

# liver.merge <- merge(liver1SCT1, y=liver2SCT1, add.cell.ids = c("PairFed1", "GDF15_1"), project = "PairFedGDF15")
# levels(liver.merge)

table(liver.merge$orig.ident)
liver.merge@active.ident
head(liver.merge[[]])
tail(liver.merge[[]])


# this function uses minimum of the median UMI (calculated using the raw UMI counts) of individual objects
# to reverse the individual SCT regression model using minimum of median UMI as the sequencing depth covariate.
liver.merge=PrepSCTFindMarkers(liver.merge, assay = "SCT", verbose = TRUE)


########## integrate data from the two conditions (control and treatment)
# When aligning two genome sequences together,
# identification of shared/homologous regions can help to interpret differences
# between the sequences as well.
liver.mergeIntegr = liver.merge
liver.mergeIntegr
liver.mergeIntegr[["SCT"]]

# run standard analysis workflow
DefaultAssay(liver.mergeIntegr) <- "SCT"
VariableFeatures(liver.mergeIntegr) <- c(VariableFeatures(liver1SCT1),
                                         VariableFeatures(liver2SCT1))
#VariableFeatures(liver.mergeIntegr) <- c(VariableFeatures(liver1SCT1), VariableFeatures(liver2SCT1))

liver.mergeIntegr <- RunPCA(liver.mergeIntegr, npcs = 30, verbose = FALSE)

# integration (HarmonyIntegration)
liver.mergeIntegr <- IntegrateLayers(object = liver.mergeIntegr, method = HarmonyIntegration, orig.reduction = "pca",
                                     normalization.method = "SCT", new.reduction = "harmony", verbose = T)

# we can now visualize and cluster the datasets.
liver.mergeIntegr <- FindNeighbors(liver.mergeIntegr, reduction = "harmony", dims = 1:30)
liver.mergeIntegr <- FindClusters(liver.mergeIntegr, verbose = FALSE, resolution = 0.8, cluster.name = "harmony_clusters") # 17 cluster
liver.mergeIntegr <- RunUMAP(liver.mergeIntegr, reduction = "harmony",
                             dims = 1:30, reduction.name = "umap.harmony")


# Visualization
DimPlot(liver.mergeIntegr, reduction = "umap.harmony", label = TRUE, group.by = c("ident", "orig.ident"))
graph2svg(x = NULL, file='17-clusterMergeIntegr_DimPlotr1', font = "Arial", cairo = TRUE,
          width = 15, height = 5, bg = "transparent")
#without labeling
DimPlot(liver.mergeIntegr, reduction = "umap.harmony", label = F, group.by = c("ident", "orig.ident"))
graph2svg(x = NULL, file='17_1-clusterMergeIntegr_DimPlotwoLabel_2sampr1', font = "Arial", cairo = TRUE,
          width = 15, height = 5, bg = "transparent")
graph2svg(x = NULL, file='9_1-clusterMergeIntegr_DimPlotwoLabel_2sampr1', font = "Arial", cairo = TRUE,
          width = 7, height = 4, bg = "transparent")

SpatialDimPlot(liver.mergeIntegr, label = TRUE, label.size = 3)
graph2svg(x = NULL, file='17_2-clusterMergeIntegr_SpatialDimPlotr1', font = "Arial", cairo = TRUE,
          width = 15, height = 6, bg = "transparent")




### Annotation by GPT-4: the name is same as the merge before integration
# IMPORTANT! Assign your OpenAI API key. See Vignette for details
Sys.setenv(OPENAI_API_KEY = 'xxxxx') #get your own OPENAI_API_KEY

# Load packages
library(GPTCelltype)
library(openai)

# Find Markers
#liver.merge=PrepSCTFindMarkers(liver.merge, assay = "SCT", verbose = TRUE)
all_markers = FindAllMarkers(liver.mergeIntegr, assay = "SCT")
markers = all_markers

# GPT-4 annotation
res <- gptcelltype(markers, tissuename = 'liver', model = 'gpt-4')
#res <- gptcelltype(markers, tissuename = 'liver', model = 'gpt-4')
res
write.csv(res, 'Celltype_IntegCCC_r1_v6.csv')
#write.csv(res, 'AnnoInteg1.csv')

# res=["Hepatocytes1", "Hepatocytes2", "Hepatocytes3", "Hepocyte/Stellate cells", "Hepatocytes4"]
#
# res=c("Hepatocytes1", "Hepatocytes2", "Hepatocytes3", "Hepatic Stellate Cells1", "Kupffer Cells","Hepatocellular carcinoma cells1",
#       "Hepatocytes4", "B Cells", "Hepatic Stellate Cells2", "Erythrocytes", "Myofibroblasts",
#       "Hepatocellular carcinoma cells3", "Hepatocellular carcinoma cells4", "Erythrocytes", "Hepatocellular carcinoma cells5", "Hepatocellular carcinoma cells6",
#       "Cholangiocytes", "Hepatocellular carcinoma cells2", "Hepatocytes6", "Hepatocytes7")

celltype <- c("Hepatocytes1", "Hepatocytes2", "Hepatocytes3", "Hepatocytes4", "Hepatocytes5","Hepatocytes6",
              "HSCs", "Macrophages", "Hepatocytes7", "Hepatocytes8", "B Cells", "Erythroid cells", "Immune Cells")

liver.mergeIntegr1=liver.mergeIntegr
names(celltype) <- levels(liver.mergeIntegr1)
liver.mergeIntegr1 <- RenameIdents(liver.mergeIntegr1, celltype)

# Assign cell type annotation back to Seurat object
# liver.merge1=liver.merge
# liver.merge1@meta.data$celltype <- as.factor(res[as.character(Idents(liver.merge1))])

# Visualize cell type annotation on UMAP
DimPlot(liver.mergeIntegr1)
graph2svg(x = NULL, file='18-Integrat_DimPlot_labeling', font = "Arial", cairo = TRUE,
          width = 8, height = 5, bg = "transparent")

head(liver.merge1[[]])



### save file
saveRDS(liver.mergeIntegr1, file = "CCCmeta.rds")


## Load the data
liver.mergeIntegr1 <- readRDS("CCCmeta.rds")


# Chnage name and Split two datasets from integrative file
head(liver.mergeIntegr1@meta.data)
Idents(liver.mergeIntegr1)
levels(liver.mergeIntegr1)

object= liver.mergeIntegr1
table(object$orig.ident)
Idents(object)
object$celltype <- Idents(object)
head(object@meta.data)

object$orig.ident[object$orig.ident == "PairFed1"] <- "control1"
object$orig.ident[object$orig.ident == "GDF15_1"] <- "treat1"

object$stim = object$orig.ident
object$stim[object$stim == "control1"] <- "control"
object$stim[object$stim == "treat1"] <- "treat"

head(object@meta.data)


# split samples
Liv1 <- subset(object, orig.ident == 'control1')
table(Liv1$orig.ident)

Liv2 <- subset(object, orig.ident == 'treat1')
table(Liv2$orig.ident)


### add data
seu1 = Liv1
seu2 = Liv2


# show the image and annotated spots
color.use <- scPalette(nlevels(seu1)); names(color.use) <- levels(seu1)
p1 <- Seurat::SpatialDimPlot(seu1, label = F, label.size = 3, cols = color.use)
color.use <- scPalette(nlevels(seu2)); names(color.use) <- levels(seu2)
p2 <- Seurat::SpatialDimPlot(seu2, label = F, label.size = 3, cols = color.use) + NoLegend()
p1 + p2

graph2svg(x = NULL, file='10-9-ImageAnnoPlotsr1', font = "Arial", cairo = TRUE,   
          width = 8, height = 5, bg = "transparent")


# Prepare input data for CellChat analysis
data.input1 = Seurat::GetAssayData(seu1, layer = "data", assay = "SCT") # normalized data matrix
data.input2 = Seurat::GetAssayData(seu2, layer = "data", assay = "SCT") 

genes.common <- intersect(rownames(data.input1), rownames(data.input2))
colnames(data.input1) <- paste0("A1_", colnames(data.input1))
colnames(data.input2) <- paste0("A2_", colnames(data.input2))
data.input <- cbind(data.input1[genes.common, ], data.input2[genes.common, ])


# define the meta data
# a column named `samples` should be provided for spatial transcriptomics analysis,
# which is useful for analyzing cell-cell communication by aggregating multiple samples/replicates.
# Of note, for comparison analysis across different conditions, 
# users still need to create a CellChat object seperately for each condition.
meta1 = data.frame(labels = Idents(seu1), samples = "A1") # manually create a dataframe consisting of the cell labels
meta2 = data.frame(labels = Idents(seu2), samples = "A2") 

meta <- rbind(meta1, meta2)
rownames(meta) <- colnames(data.input)

# a factor level should be defined for the `meta$labels` and `meta$samples`
meta$labels <- factor(meta$labels, levels = levels(Idents(seu1)))
meta$samples <- factor(meta$samples, levels = c("A1", "A2"))
unique(meta$labels) # check the cell labels
#meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))

unique(meta$samples) # check the sample labels


# load spatial transcriptomics information
# Spatial locations of spots from full (NOT high/low) resolution images are required. 
# For 10X Visium, this information is in `tissue_positions.csv`. 
spatial.locs1 = Seurat::GetTissueCoordinates(seu1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs2 = Seurat::GetTissueCoordinates(seu2, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs <- rbind(spatial.locs1, spatial.locs2)
rownames(spatial.locs) <- colnames(data.input)

nrow(spatial.locs)
length(colnames(data.input))
dim(spatial.locs)
dim(data.input)

# Scale factors of spatial coordinates
# For 10X Visium, the conversion factor of converting spatial coordinates
# from Pixels to Micrometers can be computed as the ratio of the theoretical spot size (i.e., 65um)
# over the number of pixels that span the diameter of a theoretical spot size
# in the full-resolution image (i.e., 'spot_diameter_fullres' in pixels in the 'scalefactors_json.json' file).
scalefactors1 = jsonlite::fromJSON(txt = file.path("C:/Users/wddon/OneDrive/文档/BaiduSyncdisk/R analysis/Spatial transcriptomics/GDF15CR/Sam19/spatial", 'scalefactors_json.json'))
spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor1 = spot.size/scalefactors1$spot_diameter_fullres
spatial.factors1 = data.frame(ratio = conversion.factor1, tol = spot.size/2)


scalefactors2 = jsonlite::fromJSON(txt = file.path("C:/Users/wddon/OneDrive/文档/BaiduSyncdisk/R analysis/Spatial transcriptomics/GDF15CR/Sam10/spatial", 'scalefactors_json.json'))
conversion.factor2 = spot.size/scalefactors2$spot_diameter_fullres
spatial.factors2 = data.frame(ratio = conversion.factor2, tol = spot.size/2)

spatial.factors <- rbind(spatial.factors1, spatial.factors2)
rownames(spatial.factors) <- c("A1", "A2")


#### Create a CellChat object
head(spatial.locs)
spatial.locs=spatial.locs[,-3]
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)

cellchat


#### Set the ligand-receptor interaction database
# Use the CellChatDB with metabolic signaling

##$$$$ for all pathways
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
cellchat@DB <- CellChatDB.use
unique(CellChatDB.use$interaction$annotation)
unique(CellChatDB.use$interaction$pathway_name)
unique(CellChatDB.use$interaction$pathway_name)


# $$$$ only for metabolic signaling
# get pathways related glucose metabolism
CellChatDB <- CellChatDB.mouse
metabolic_pathways <- c("INSULIN", "IGF", "GIPR", "GCG", 
                        "LEP", "ADIPONECTIN", "RESISTIN", "ApoE")
CellChatDB.use <- subsetDB(CellChatDB, search = metabolic_pathways, key = "pathway_name")
cellchat@DB <- CellChatDB.use
unique(CellChatDB.use$interaction$annotation)
unique(CellChatDB.use$interaction$pathway_name)

###$$$ pathways for lipid metabolism
CellChatDB <- CellChatDB.mouse
lipid_pathways <- c("LEP", "ADIPONECTIN", "RESISTIN", "ANGPTL", "ApoE", "ApoA", "ApoB",
                    "LXA4", "27HC", "Cholesterol", "Calcitriol", "Desmosterol",
                    "DHEA", "DHT", "Estradiol", "Progesterone", "Testosterone")
CellChatDB.use <- subsetDB(CellChatDB, search = lipid_pathways, key = "pathway_name")
cellchat@DB <- CellChatDB.use
unique(CellChatDB.use$interaction$annotation)
unique(CellChatDB.use$interaction$pathway_name)


###$$$ pathways for amino acid metabolism
CellChatDB <- CellChatDB.mouse
aa_pathways <- c("GABA-A", "GABA-B", "Glutamate", "Glycine", "SerotoninDopamin", 
                 "Histamine", "IGFBP", "NMU", "NPY", "NTS", "VIP", 
                 "PACAP", "SOMATOSTATIN", "TAFA", "PROK")

CellChatDB.use <- subsetDB(CellChatDB, search = aa_pathways, key = "pathway_name")
cellchat@DB <- CellChatDB.use
unique(CellChatDB.use$interaction$annotation)
unique(CellChatDB.use$interaction$pathway_name)


#### Preprocessing the expression data for cell-cell communication analysis
# To infer the cell state-specific communications, we identify over-expressed ligands
# or receptors in one cell group and then identify over-expressed
# ligand-receptor interactions if either ligand or receptor is over-expressed.

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#execution.time = Sys.time() - ptm
#print(as.numeric(execution.time, units = "secs"))


#### Part II: Inference of cell-cell communication network
# Compute the communication probability and infer cellular communication network
#ptm = Sys.time()


#For certain functions, each worker needs access to certain global variables.
#If these are larger than the default limit, you will see this error. 
#To get around this, you can set options(future.globals.maxSize = X),where X is the maximum allowed size in bytes.
#So to set it to 1GB, you would run options(future.globals.maxSize = 1000 * 1024^2). 
#Note that this will increase your RAM usage so set this number mindfully.
options(future.globals.maxSize = 1050 * 1024^2)
cellchat@images$coordinates <- as.matrix(cellchat@images$coordinates)

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)

#Users can filter out the cell-cell communication if there are only few cells in certain cell groups. 
#By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
cellchat <- filterCommunication(cellchat, min.cells = 10)

## Infer the cell-cell communication at a signaling pathway level
#CellChat computes the communication probability on signaling pathway level by summarizing the communication
#probabilities of all ligands-receptors interactions associated with each signaling pathway.

#NB: The inferred intercellular communication network of each ligand-receptor pair
#and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
cellchat <- computeCommunProbPathway(cellchat)


##### Calculate the aggregated cell-cell communication network
#We can calculate the aggregated cell-cell communication network by
#counting the number of links or summarizing the communication probability.
cellchat <- aggregateNet(cellchat)


#We can also visualize the aggregated cell-cell communication network. 
#For example, showing the number of interactions or the total interaction strength (weights)
#between any two cell groups using circle plot or heatmap plot.

#ptm = Sys.time()

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count),
                 weight.scale = T, label.edge= F, title.name = "Number of interactions (Glu Met)")
graph2tif(x = NULL, file='10-10-CCCNumInte_r1', font = "Arial", cairo = TRUE,   
          width = 8, height = 8, bg = "transparent")


par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight),
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength (Glu Met)")
graph2tif(x = NULL, file='10-11-CCCInteWeig_r1', font = "Arial", cairo = TRUE,   
          width = 5, height = 5, bg = "transparent")


netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")
graph2tif(x = NULL, file='10-12-CCCHeatmap_r1', font = "Arial", cairo = TRUE,   
          width = 5, height = 5, bg = "transparent")


#### Part III: Visualization of cell-cell communication network
#Upon infering the cell-cell communication network, CellChat provides various functionality
#for further data exploration, analysis, and visualization. 
#Here we only showcase the circle plot and the new spatial plot.

#All the signaling pathways showing significant communications can be accessed by following
cellchat@netP$pathways

pathways.show <- c("IGF") 

# Circle plot
par(mfrow=c(1,1), xpd=TRUE)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
graph2tif(x = NULL, file='10-13-CCCIGF_r1', font = "Arial", cairo = TRUE,   
          width = 4, height = 4, bg = "transparent")


# Spatial plot
#for Liv1
par(mfrow=c(1,1))
# Setting `vertex.label.cex = 0` to hide the labels on the spatial plot
netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A1",
                    layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 0)
graph2tif(x = NULL, file='10-14-CCCMacrophagePathwaySpat_liv1_r1', font = "Arial", cairo = TRUE,   
          width = 5, height = 5, bg = "transparent")

#for Liv2
par(mfrow=c(1,1))
# Setting `vertex.label.cex = 0` to hide the labels on the spatial plot
netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A2",
                    layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 0)
graph2svg(x = NULL, file='10-15-CCCMacrophagePathwaySpat_liv3', font = "Arial", cairo = TRUE,   
          width = 5, height = 5, bg = "transparent")



#### Compute and visualize the network centrality scores:
# Compute the network centrality scores
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Visualize the computed centrality scores using heatmap, allowing ready
# identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

graph2tif(x = NULL, file='10-16-CCCcentrality_r1', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 4, bg = "transparent")



#### Compute the contribution of each ligand-receptor pair to the overall signaling pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)
graph2tif(x = NULL, file='10-18-CCClig-rec pairContri_r1', font = "Arial", cairo = TRUE,   
          width = 4.5, height = 2.5, bg = "transparent")


## When visualizing gene expression distribution on tissue using spatialFeaturePlot,
## users also need to provide the sample.use as an input.

# Take an input of a few genes
spatialFeaturePlot(cellchat, features = c("Itgav","Itgb3"),
                   sample.use = "A1", point.size = 0.8, color.heatmap = "Reds", direction = 1)
graph2tif(x = NULL, file='10-19-CCClig-recGeneExpSpa_s1_r1', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")

spatialFeaturePlot(cellchat, features = c("Itgav","Itgb3"),
                   sample.use = "A2", point.size = 0.8, color.heatmap = "Reds", direction = 1)
graph2tif(x = NULL, file='10-19-CCClig-recGeneExpSpa_s2_r1', font = "Arial", cairo = TRUE,   
          width = 6, height = 4, bg = "transparent")


# # Take an input of a ligand-receptor pair (not working)
# spatialFeaturePlot(cellchat, pairLR.use = "Mif-(Cd74+Cd44)", sample.use = "A1",
#                    point.size = 0.5, do.binary = FALSE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)
# 
# 
# # Take an input of a ligand-receptor pair and show expression in binary (not working)
# spatialFeaturePlot(cellchat, pairLR.use = "AREG_EGFR", sample.use = "A1",
#                    point.size = 1.5, do.binary = TRUE, cutoff = 0.05, enriched.only = F, color.heatmap = "Reds", direction = 1)


#### Part V: Save the CellChat object
saveRDS(cellchat, file = "cellchat_mouse_Liv1a2_Met.rds")







