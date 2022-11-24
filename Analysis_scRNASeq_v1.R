## Analisis datos sc-RNAseq
## VERSION 0.1  fecha: 20221001 
## AUTOR: Oscar Pastor 
## NOTAS: 
###################

library(dplyr)
library(Seurat)
library(patchwork)

####### GET DATA FILES
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
directory <- getwd()

######## LOADING FUNCTIONS
source("fnx/BASE_FUNCTION.R")
source("fnx/OUTLIERS.R")
source("fnx/ADV_STATISTICS.R")


filenames <- list.files(directory)[grepl(".*(.txt|*.tsv)$",list.files(directory),ignore.case=TRUE)]

filenames <-selectOption(filenames,"Select_Data_files")

filenames

data_files<-file.path(directory,filenames)

#data.mtx<-data.table::fread(data_files[1])
#rm(data.mtx)
#paste0(directory,"/")
data.metadata<-data.table::fread(data_files[1],data.table=F)
data.raw <- data.table::fread(data_files[2],data.table = F)

rownames(data.raw)<-data.raw$V1

FileNames<-unique(data.metadata$FileName)
CellIDs<-unique(data.metadata$CellID)
CellTypes<-unique(data.metadata$CellType)
ColumnRawNames<-unique(colnames(data.raw))
RowRawNames<-unique(data.raw$V1)

ColumnRawNames[1:100]
RowRawNames[1:100]

### Understanding data_structure

## correspondencia rawdata vs metadata CellIDS
unique(gsub("(.*)\\_(.*)","\\2",ColumnRawNames))
match(CellIDs,unique(gsub("(.*)\\_(.*)","\\2",ColumnRawNames)))

## correspondencia rawdata vs metadata FileNames=Samples
unique(gsub("(.*)\\_(.*)","\\1",ColumnRawNames))
match(FileNames,unique(gsub("(.*)\\_(.*)","\\1",ColumnRawNames)))

##Sampling to avoid lack of memmory
### Filter to reduce dataset
filter_condition<-"Hepatocyte_15|NPC_15|Hepatocyte_Chow|NPC_Chow"
filter_condition<-"Hepatocyte_|NPC_"

data.metadata.reduced<-data.metadata[grepl(filter_condition,data.metadata$FileName),]
rownames(data.metadata.reduced)<-paste(data.metadata.reduced$FileName,data.metadata.reduced$CellID,sep="_")  
unique(data.metadata.reduced$FileName)

## Hacer un from
## split vector in chuncks. TODO:seguro que hay una manera mas elegante
chunks<-3000
steps<-c(seq(from=0,to=length(rownames(data.raw)),by=chunks),length(rownames(data.raw)))
steps_fin<-steps[-c(1)]
steps_ini<-steps+1
steps_ini<-steps_ini[-c(length(steps_ini))]

data.seurat<-list()
min.cells.sel<-10
min.features.sel<-400

rm(data.seurat.combo)
for (i in seq_along(steps_ini)){
 # for (i in 8:9){
    #i=9
  data.reduced<-data.raw[steps_ini[i]:steps_fin[i],grepl(filter_condition,colnames(data.raw))]
  data.seurat[i]<-CreateSeuratObject(data.reduced, min.cells = min.cells.sel, min.features = min.features.sel, 
                                   names.field = 5,names.delim = "_",meta.data =data.metadata.reduced )
  rm(data.reduced)
  if (!length(colnames(data.seurat[[i]]))==0){ ## si no hay seleccion
     if (exists("data.seurat.combo")==F){
        data.seurat.combo<-data.seurat[[1]]
        print(paste0("Chunk:",i," data.seurat.combo tiene: ",length(rownames(data.seurat.combo))," Rows (features) X "
                     ,length(colnames(data.seurat.combo))," Cols (cells)"))
     }else {
        data.seurat.combo<-merge(data.seurat.combo,data.seurat[[i]])
        print(paste0("Chunk:",i," data.seurat.combo tiene: ",length(rownames(data.seurat.combo))," Rows (features) X "
                ,length(colnames(data.seurat.combo))," Cols (cells)"))
               
    }
  }
}

#data.reduced<-data.raw[,grepl(filter_condition,colnames(data.raw))]

summary(data.seurat.combo@meta.data$nCount_RNA) #Counts
summary(data.seurat.combo@meta.data$nFeature_RNA) #QC RNAs
unique(data.seurat.combo@meta.data$CellType)
unique(data.seurat.combo@meta.data$FileName)

## list rownames and colnames of seurat object
data.seurat.combo[["RNA"]]
rownames(data.seurat.combo)
length(rownames(data.seurat.combo))
length(colnames(data.seurat.combo))


# Load the PBMC dataset
#pbmc.data <- Read10X(data.dir = paste0( directory,"/data/pbmc3k/filtered_gene_bc_matrices/hg19/"))
# Initialize the Seurat object with the raw (non-normalized data).
#pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
#pbmc
data.seurat<-data.seurat.combo
### TODO QC control getting rid of mitocondrial genes
# Add number of genes per UMI for each cell to metadata
#Number of genes detected per UMI: 
data.seurat$log10GenesPerUMI <- log10(data.seurat$nFeature_RNA) / log10(data.seurat$nCount_RNA)
summary(data.seurat$log10GenesPerUMI)
# Proportion of transcripts mapping to mitochondrial genes. 
#data.seurat$mt<-sum(rownames(data.seurat) %in% mt)
data.seurat$mitoRatio <- PercentageFeatureSet(object = data.seurat, pattern = "^mt-")
data.seurat$mitoRatio <-data.seurat@meta.data$mitoRatio / 100



# Lets examine a few genes in the first thirty cells
data.seurat[c("Csrp2", "Atxn7l3b","Gpr182"), 1:30]
dense.size <- object.size(as.matrix(data.seurat))
dense.size


#### COMIENZO EL ANALISIS



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
data.seurat[["percent.mt"]] <- PercentageFeatureSet(data.seurat, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(data.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(data.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

##Selection of some features

data.seurat <- subset(data.seurat, subset = nFeature_RNA > 300 & nFeature_RNA < 2500 & percent.mt < 5)

###  Normalizamos

data.seurat<- NormalizeData(data.seurat, normalization.method = "LogNormalize", scale.factor = 10000)

###  Seleccion variables de interes
#By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.

data.seurat <- FindVariableFeatures(data.seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data.seurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data.seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

### Scaling of genes before PCA analysis
all.genes <- rownames(data.seurat)
#pbmc <- ScaleData(pbmc, features = all.genes)
data.seurat <- ScaleData(data.seurat,features = all.genes) #FASTER: To perform scaling on the previously identified variable features (2,000 by default)

#### Remove unwanted sources of variation
# TODO: SCTransform() es mas avanzado ver en documentacion
data.seurat <- ScaleData(data.seurat, vars.to.regress = "percent.mt")


### PCA analysis ###
data.seurat <- RunPCA(data.seurat, features = VariableFeatures(object = data.seurat))

# Examine and visualize PCA results a few different ways
print(data.seurat[["pca"]], dims = 1:5, nfeatures = 5)

#Show Loadings of two first components
VizDimLoadings(data.seurat, dims = 1:2, reduction = "pca")

#Plot PC1vsPC2
DimPlot(data.seurat, reduction = "pca")

###TOOLS to selectc how many components to choose for further analysis
### Heatmap captures the contribution of each loading to PCAs
DimHeatmap(data.seurat, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(data.seurat, dims = 1:15, cells = 200, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
# 1. ElbowPlot  (less consuming in resources)
   ElbowPlot(data.seurat)  ### choose when the elbow changes +3-4, play with values
# 2. JackStraw  (much more demanding in cpu time)
   #pbmc <- JackStraw(pbmc, num.replicate = 100)
   #pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
   #JackStrawPlot(pbmc, dims = 1:15)

###Cluster the cells
   # We chose e 10 fist components  here from ElbowPlot or JackStraw
   data.seurat <- FindNeighbors(data.seurat, dims = 1:15)
   data.seurat <- FindClusters(data.seurat, resolution = 0.5)
   
   #Look at cluster IDs of the first 5 cells
   head(Idents(data.seurat), 5) 
   
### Run non-linear dimensional reduction (UMAP/tSNE)
   
   data.seurat <- RunUMAP(data.seurat, dims = 1:30)
   # note that you can set `label = TRUE` or use the LabelClusters function to help label
   # individual clusters
   DimPlot(data.seurat, reduction = "umap")

   ### Save the object at this point so that it can easily be loaded back 
   ###TODO crear la carpeta output si no existe
 
   saveRDS(data.seurat.combo, file =   paste0( directory,"/output/hepatocyte_chow.rds"))   
   
### Finding differentially expressed features (cluster biomarkers)
   
# find all markers of cluster 2
   cluster2.markers <- FindMarkers(data.seurat, ident.1 = 2, min.pct = 0.25)
   head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
   cluster5.markers <- FindMarkers(data.seurat, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
   head(cluster5.markers, n = 5)   
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
   hepatocyte15w.markers <- FindAllMarkers(data.seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
   hepatocyte15w.markers %>%
     group_by(cluster) %>%
     slice_max(n = 2, order_by = avg_log2FC)   
   
## for evluating how good a marker is for an specific cluster use FindMarkers, Explore DE vignette
     cluster0.markers <- FindMarkers(data.seurat, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)   
     cluster1.markers <- FindMarkers(data.seurat, ident.1 = 1, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)   
    ### etc
     ### TODO make function to choose automatically
     
### Violin plots for specific marker-genes associated with clusters
   VlnPlot(data.seurat, features =rownames(cluster0.markers)[1:10])   
   VlnPlot(data.seurat, features =c("Degs1","Alb"))  ## or choose specific clustering
   # you can plot raw counts as well
   VlnPlot(data.seurat, features = rownames(cluster1.markers)[1:10], slot = "counts", log = TRUE)   
  
   
### Locate marker-genes in UMAP plots
   FeaturePlot(data.seurat, features = c("Vsig4","Degs1","Fabp5","Alb",
                                         "Clec4g","Cers2","Lyz1","Reln","Rgs5"))
   FeaturePlot(data.seurat, features = c("Degs1","Cers2","Cers5","Alb","Fabp5","Fabp4","Clec4g","Lyz1","Reln","Ccr2"))
   
### Locate marker-genes in UMAP plots
   FeaturePlot(pbmc, features = rownames(cluster1.markers)[1:10])

### Biclustering plots   
    pbmc.markers %>%
     group_by(cluster) %>%
     top_n(n = 10, wt = avg_log2FC) -> top10
   DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#### Assigning cell type identity to clusters
   
   ##This is arbitrary!!! 

  new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                        "NK", "DC", "Platelet")
   names(new.cluster.ids) <- levels(pbmc)
   pbmc <- RenameIdents(pbmc, new.cluster.ids)
   DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## Guardar los resultados del analisis 
   saveRDS(pbmc, file = paste0(directory,"/output/pbmc3k_final.rds"))
   