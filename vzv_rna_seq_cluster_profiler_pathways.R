#!/usr/bin/env Rscript

library(clusterProfiler)
library(org.Hs.eg.db)

de_from_pivot <- read.csv(file=paste0("/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019",
  "/output_hg38/deseq_results_vzv_cohrs_rna-seq_reseq_fdr_0_dot_1_09022019.csv"),header=TRUE)

genes_up <- unlist(de_from_pivot[de_from_pivot$log2FoldChange>0,1])
genes_down <- unlist(de_from_pivot[de_from_pivot$log2FoldChange<0,1])

##Useful Link: https://stat.ethz.ch/pipermail/bioconductor/2009-April/027168.html

##Get the Entrez gene IDs associated with those symbols
EG_IDs_up = mget(as.character(unlist(genes_up)), revmap(org.Hs.egSYMBOL),ifnotfound=NA)
EG_IDs_down = mget(as.character(unlist(genes_down)), revmap(org.Hs.egSYMBOL),ifnotfound=NA)

entrez_up <- unlist(unname(EG_IDs_up[!is.na(EG_IDs_up)]))
entrez_down <- unlist(unname(EG_IDs_down[!is.na(EG_IDs_down)]))

##Then get the KEGG IDs associated with those entrez genes.
##KEGG_IDs_up = unlist(mget(as.character(EG_IDs_up), org.Hs.egPATH,ifnotfound=NA))
##KEGG_IDs_down = unlist(mget(as.character(EG_IDs_down), org.Hs.egPATH,ifnotfound=NA))

##kegg_up <- unname(KEGG_IDs_up[!is.na(KEGG_IDs_up)])
##kegg_down <- unname(KEGG_IDs_up[!is.na(KEGG_IDs_down)])

typeof(unlist(entrez_up))
length(entrez_up)

kegg_paths_up <- enrichKEGG(entrez_up, 
           organism = "hsa", 
           keyType = "kegg",
           pvalueCutoff = 0.05, 
           pAdjustMethod = "BH", 
           qvalueCutoff = 0.2)

kegg_paths_down <- enrichKEGG(entrez_down, 
                            organism = "hsa", 
                            keyType = "kegg",
                            pvalueCutoff = 0.05, 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.2)

kegg_paths_up@result$Description

kegg_paths_down@result$Description

barplot(kegg_paths_up, showCategory=30)

barplot(kegg_paths_down, showCategory=30)

##https://yulab-smu.github.io/clusterProfiler-book/chapter12.html

##Gene-concept network
#### convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, foldChange=geneList)