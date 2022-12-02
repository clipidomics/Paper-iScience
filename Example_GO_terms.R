library(topGO)
library(GO.db)
library(org.Mm.eg.db)
#library(annotate)
library(AnnotationDbi)
library(KEGGREST)
library(ensembldb)

##Description of the db
org.Mm.eg_dbfile()
org.Mm.eg_dbInfo()
org.Mm.eg_dbconn()
org.Mm.eg_dbschema()


###Search for GOTERM containing a specific text in this case Sphingolipids and Ceramide

goterms = unlist(Term(GOTERM))
search_patern<-"[S|s]sphingolipid|[C|c]eramide"
whmf<-grepl(search_patern, goterms)
SL_inGOterm<-goterms[whmf]  #vector que contiene los terminos GO relacionados con Cer y Sphingos
SL_inGOterm<-data.frame(cbind(GO=names(SL_inGOterm),Description=SL_inGOterm))

### List all Entrez Gen IDs related to Search pattern in GO
xx_GO_EG <- as.list(org.Mm.egGO2EG) 
xx_GO_EG[SL_inGOterm$GO]
genEGid<-unique(unlist(xx_GO_EG[SL_inGOterm$GO]))
myGenes.symbol <- data.frame(SYMBOL=mapIds(org.Mm.eg.db, as.character(unlist(genEGid)), column="SYMBOL",keytype=c("ENTREZID")))
myGenes.gename <- data.frame(GENENAME=mapIds(org.Mm.eg.db, as.character(unlist(genEGid)), "GENENAME","ENTREZID"))
myGenes<-cbind(ENTREZID=rownames(myGenes.symbol),myGenes.symbol,myGenes.gename)

# ### Add those in KEGG  #This is not working
# #Get list of KEGG databases which contain or are related to specific organism
# db_kegg<-data.frame(keggList("organism"))
# dbTOsearch<-db_kegg[grep("[M|m]us musculus",db_kegg$species),"T.number"] #seach db for organism
# kegg_terms<-cbind(KEGG_NAME=gsub("(.*)\\;\\s(.*)","\\2",keggFind(dbTOsearch,c("Ceramide"))),  ## Get KEGG_TERMS
#                   KEGG_SYMBOL=gsub("(.*)\\;\\s(.*)","\\1",keggFind(dbTOsearch,c("Ceramide"))))

## Exploring other alternatives
library (biodbKegg)
mybiodb <- biodb::newInst()
#1st consult pathways related with compounds in kegg
kegg.comp.conn <- mybiodb$getFactory()$createConn('kegg.compound')
# compounds <- mybiodb$entriesToDataframe(kegg.comp.conn$getEntry(c('C00195','C01190')), compute=FALSE)
# compounds$kegg.pathway.id
# compounds$name
kegg.comp.ids <- c('C00195','C01190','C00550') #Ceramide and glucosilceramide
#pathways <- kegg.comp.conn$getPathwayIds(kegg.comp.ids, 'mmu')
pathways <- kegg.comp.conn$getPathwayIdsPerCompound(kegg.comp.ids, 'mmu', limit = 10)
pathways<-unique(unlist(pathways,use.names=F))
library (OmnipathR)
library (orthogene)
mmu_genes<-list()
for (i in 1:length(pathways)){
  ## convierte los nombres de humano a murino
  #i=1
 #kegg_pw[[i]]<-kegg_pathway_download(sub("mmu","hsa",pathways[i]))
 kegg_pw<-kegg_pathway_download(sub("mmu","hsa",pathways[i]))
 
 if (!all(is.na(kegg_pw))){
   
   kegg_pw<-unique(kegg_pw$genesymbol_source)
   mmu_genes[[i]]<-rownames(convert_orthologs(kegg_pw,input_species = "human",
                                output_species = "mouse",
                                non121_strategy = "kbs",
                                method = "gprofiler",
                                mthreshold = 4))
   
    }
 
}
mmu_genes_list<-sort(unique(unlist(mmu_genes,use.names=F)))
mmu_genes_list

###Get paralog genes

library(biomaRt)
#mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
head(searchAttributes(mouse, pattern="region"))
#filters <- listFilters(mouse)
#attributes <- listAttributes(mouse)

results <- getBM(attributes = c("ensembl_gene_id", 
                                "external_gene_name",
                                "chromosome_name",
                                "mmusculus_paralog_ensembl_gene", 
                                "mmusculus_paralog_associated_gene_name",
                               # "mmusculus_homolog_goc_score", 
                                #"mmusculus_homolog_orthology_confidence",
                                "mmusculus_paralog_orthology_type",
                                "mmusculus_paralog_perc_id"
                                ),
                 filters = "external_gene_name",
                 values = mmu_genes_list,
                 mart = mouse)


unique(results$chromosome_name)
rm(paralogs)

results_within<-results[which(results$mmusculus_paralog_orthology_type=="within_species_paralog"),]
paralogs<-unique(sort(c(unique(results$external_gene_name),unique(results_within$mmusculus_paralog_associated_gene_name))))
length(paralogs)
paralogs

##Quick list of genes related to SL methabolism
my_list_genes_SL<-sort(unique(c(myGenes$SYMBOL,paralogs)))

saveRDS(my_list_genes_SL, file =   paste0( directory,"/output/my_list_genes_SL.rds"))   


###DEFINE MITOCHONDRIAL GENES
# Get mitocondrial annotated genes from esbml database
mt_genes<-getBM(attributes = c("external_gene_name","gene_biotype",
                               "transcript_biotype",
                               "chromosome_name", "start_position"), 
                filters = "chromosome_name",
                values = "MT",
                mart = mouse
)

saveRDS(mt_genes, file =   paste0( directory,"/output/mitochondrial_genes.rds"))   


rownames(data.seurat)[rownames(data.seurat) %in% gsub("mt-","",mt_genes$external_gene_name)]

