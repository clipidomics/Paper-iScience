## Analisis datos sc-RNAseq iScience Paper
## VERSION 1.0  fecha: 20221124
## AUTOR: Oscar Pastor 
## NOTAS: 
## PMID: 34755088
# 
# Analysis notes
# Groups
# chow diet (n = 6) 
# NAFL with HFHFD for 15 weeks (n = 3)
# NASH with HFHFD for 30 (n = 2 for hepatocytes and n = 4 for NPCs) 
# NASH with HFHFD for 34 weeks (n = 1 for NPCs)

# Cells with fewer than 300 detected genes, 
# more than 60,000 detected UMI, 
# and a very high (>0.8) mitochondrial genome transcript ratio were filtered. 
# The top 2,000 variable genes (HVGs) were identified FindVariableFeatures method
# Canonical Correlation Analysis (Berry et al., 2017) was used to identify common anchors between cells among batches.

# RunCCA performs this analysis in R
# RunMultiCCA

###################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

######## LOADING FUNCTIONS
source("fnx/BASE_FUNCTION.R")
source("fnx/OUTLIERS.R")
source("fnx/ADV_STATISTICS.R")


filenames <- list.files(paste0(directory,"/output"))[grepl("(*.rds)$",list.files(paste0(directory,"/output")),ignore.case=TRUE)]

filenames <-selectOption(filenames,"Select_Data_files")

filenames

data_files<-file.path(paste0(directory,"/output"),filenames)

###Load Seurat object already processed and stored
#load file 2_datseu.time_Norm_full.rds
datseu.time<-readRDS(data_files[grep("time",data_files)])

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = datseu.time) #Not clear if Degs is present

# Verificamos si alguno de nuestros genes de interes es identificado en el analysis
my_stored_features<-readRDS(data_files[grep("my_lis",data_files)])
sum(my_stored_features %in% features) ## numero de identificados en la selección que tiene que ver con Sphingolipids
my_candidates<-my_stored_features[my_stored_features %in% features]

# We will manually add some features from Ceramide pathway to follow them over time

features2<-sort(unique(c(features,my_stored_features)))
length(features2)

#probablemente no puedo añadir genes que no estan presentes en scRNAseq realizado
all.genes <- unique(c(rownames(datseu.time[[1]]),rownames(datseu.time[[2]])))
length(rownames(datseu.time[[1]]))==length(rownames(datseu.time[[2]]))
length(all.genes)
features2.sel<-features2[features2 %in% all.genes] #select only genes in set
length(features2.sel)

# Perform integration  ****  THIS TAKES A LONG TIME !!! ****
### Subsetting Data NO SE CUANDO HACER!!!

datseu.anchors <- FindIntegrationAnchors(object.list = datseu.time, anchor.features = features2.sel)

saveRDS(datseu.time, file =   paste0( directory,"/output/FULL_post_anchored.rds"))   

rm(datseu.time) #free memory

datseu.combined <- IntegrateData(anchorset = datseu.anchors)

#Perform an integrated analysis
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(datseu.combined) <- "integrated"

rm(datseu.anchors) #free memory

#Idents(datseu.combined ) <- factor(Idents(datseu.combined ), levels = c("Chow","15weeks","30weeks"))

sum(datseu.combined$CellType=="Hepatocytes")

# Run the standard workflow for visualization and clustering
datseu.combined <- ScaleData(datseu.combined, verbose = FALSE)
datseu.combined <- RunPCA(datseu.combined, npcs = 30, verbose = FALSE)
###TOOLS to selectc how many components to choose for further analysis
### Heatmap captures the contribution of each loading to PCAs
DimHeatmap(datseu.combined, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(data.seurat, dims = 1:15, cells = 200, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
# 1. ElbowPlot  (less consuming in resources)
ElbowPlot(datseu.combined)  ### choose when the elbow changes +3-4, play with values



datseu.combined <- RunUMAP(datseu.combined, reduction = "pca", dims = 1:12)
#datseu.combined <- RunTSNE(datseu.combined, dims = 1:30, verbose = FALSE)
datseu.combined.copy<-datseu.combined

datseu.combined <- FindNeighbors(datseu.combined, reduction = "pca", dims = 1:12)

datseu.combined <- FindClusters(datseu.combined, reduction.type="cca.aligned", resolution = 0.5)

unique(datseu.combined$Time)
datseu.combined$Time<-factor(datseu.combined$Time, levels=c("Chow","15weeks","30weeks","34weeks"))

saveRDS(datseu.combined, file =   paste0( directory,"/output/datseu.combined_PrePlot.rds"))   

# Visualization
p1 <- DimPlot(datseu.combined, reduction = "umap", group.by = "Time")
p2 <- DimPlot(datseu.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
#unique(datseu.combined$CellType)
b<-DimPlot(datseu.combined, reduction = "umap", split.by = "Time",group.by = "CellType")
#DimPlot(datseu.combined, reduction = "tsne", split.by = "CellTime")

DefaultAssay(datseu.combined) <- "RNA"

FeaturePlot(datseu.combined, features = c("Degs1","Cd14","Cers2","Alb"), 
            split.by = "Time", max.cutoff = 3, cols = c("grey", "red"))

FeaturePlot(datseu.combined, features = c("Degs2","Degs1" ,"Epcam","Clec4g"), 
            split.by = "Time", max.cutoff = 3, cols = c("grey", "red"))

FeaturePlot(datseu.combined, features = c("Reln","Degs1","Degs2","Pla2g15"), 
            split.by = "Time", max.cutoff = 3, cols = c("grey", "red"))

plots <- VlnPlot(datseu.combined, features = c("Dgat2"), split.by = "Time", group.by = "CellType",
                 pt.size = 0.0001, combine = FALSE)
wrap_plots(plots = plots, ncol = 1,nrows=5,widths = 10,heights = 5)



