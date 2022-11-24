library(AnnotationHub)
library(ensembldb)
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Mus Musculus", "EnsDb"), 
              ignore.case = TRUE)
# Check versions of databases available
ahDb %>% mcols()
# Acquire the latest annotation files
id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb,  return.type = "data.frame")
# Select annotations of interest
annotations <- annotations %>% dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)
# Extract IDs for mitochondrial genes
mt <- annotations %>% dplyr::filter(seq_name == "MT") %>% dplyr::pull(gene_name)
not_mt<- annotations %>% dplyr::filter(seq_name != "MT") %>% dplyr::pull(gene_name)
# Number of UMIs assigned to mitochondrial genes
metadata$mtUMI <- Matrix::colSums(data.metadata.reduced[which(rownames(data.metadata.reduced) %in% mt),], na.rm = T)


# Calculate of mitoRatio per cell
metadata$mitoRatio <- metadata$mtUMI/metadata$nUMI