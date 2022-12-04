 ### Plots ####

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
library(stringr)
library(RColorBrewer)
 
######## LOADING FUNCTIONS
source("fnx/BASE_FUNCTION.R")
source("fnx/OUTLIERS.R")
source("fnx/ADV_STATISTICS.R")
####### GET DATA FILES
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory <- getwd()


filenames <- list.files(paste0(directory,"/output"))[grepl("(*.rds)$",list.files(paste0(directory,"/output")),ignore.case=TRUE)]

filenames <-selectOption(filenames,"Select_Data_files")

filenames

data_files<-file.path(paste0(directory,"/output"),filenames)

###Load Seurat object already processed and stored
#load file 2_datseu.time_Norm_full.rds
my_stored_features<-readRDS(data_files[grep("my_list",data_files)])
datseu.plots<-readRDS(data_files[grep("datseu.combined_PrePlot",data_files)])

#Preguntamos numbre de la carpeta para guardar las cosas
directory=get("directory",envir=.GlobalEnv)
print("Estas son las carpetas que hay ahora mismo creadas para analisis:")
print(matrix(gsub("(.*)\\/(.*)$","\\2",list.dirs(directory))))
carpeta_para_figuras=c("FigurasSeurat")
names(carpeta_para_figuras) <-c("Introduzca un nombre para la carpeta donde se van a guardar las figuras generadas: ")
fix(carpeta_para_figuras)
#creamos directorio nuevo
dir.create(sprintf(paste0(directory,"/",carpeta_para_figuras)),showWarnings=F);
sub_folder<-paste0(directory,"/",carpeta_para_figuras)
print(paste("Directorio",sub_folder,"creado"))


# Visualization
p1 <- DimPlot(datseu.plots, reduction = "umap", group.by = "Time")
p2 <- DimPlot(datseu.plots, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2


## Create some aliases
unique(datseu.plots@meta.data$CellType)
datseu.plots@meta.data$CellAlias<-str_sub(datseu.plots@meta.data$CellType,1,3)

Dp1<-DimPlot(datseu.plots, reduction = "umap", group.by = "CellAlias",split.by = "Time")
Dp1 + ggmin::theme_min()
select_type<-"Hepatocyte_NPC"

ggsave(filename = paste0(directory,"/",carpeta_para_figuras,"/",select_type,"_",deparse(substitute(Dp1)),".png"),
         plot=Dp1,
         width = 30, height = 10, dpi = 200, units = "cm", device='png')


# ggsave(filename = paste0(directory,"/",carpeta_para_figuras,"/",select_type,"_",deparse(substitute(Dp1)),".png"),
#        annotate_figure(
#          ggarrange(plotlist = dp1,
#                    ncol = 4, nrow = 4,
#                    common.legend = TRUE),
#          top = text_grob(paste0("Quantitation of lipid classes in ", select_type), color = "black", face = "plain", size = 14,just="center")
#        ),
#        width = 40, height = 40, dpi = 600, units = "cm", device='png')


# Get cell identity classes
Idents(datseu.plots)



#DimPlot(datseu.combined, reduction = "tsne", split.by = "CellTime")

DefaultAssay(datseu.combined) <- "RNA"

FeaturePlot(datseu.combined, features = c("Degs1","Cd14","Cers2","Alb"), 
            split.by = "Time", max.cutoff = 3, cols = c("grey", "red"))

FeaturePlot(datseu.combined, features = c("Degs2","Degs1" ,"Epcam","Clec4g"), 
            split.by = "Time", max.cutoff = 3, cols = c("grey", "red"))

FeaturePlot(datseu.combined, features = c("Reln","Degs1","Degs2","Pla2g15"), 
            split.by = "Time", max.cutoff = 3, cols = c("grey", "red"))

plots <- VlnPlot(datseu.combined, features = c("Epcam"), split.by = "Time", group.by = "CellType",
                 pt.size = 0.01, combine = FALSE)
wrap_plots(plots = plots, ncol = 1,nrows=5,widths = 10,heights = 5)



rownames(datseu.combined)





