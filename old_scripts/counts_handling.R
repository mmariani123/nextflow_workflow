#!/usr/bin/env Rscript

library(magrittr)
library(dplyr)
library(tidyr)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)

#counts <- featureCounts(bams, blablabla_restofcommand)$counts
#deseqdata <- DESeqDataSetFromMatrix(countData=counts, colData=sampleInfo, design=~condition)

##setwd("/slipstream/home/mmariani/projects/hhv6_rna/output/bowtie2_counts")
setwd("C:/Users/mmari/OneDrive/Documents/GitHub/nextflow_workflow/output/counts")

files <- list.files(path="./", pattern=".counts")
frames <- lapply(files,read.table,sep="\t",header=TRUE,skip=1,stringsAsFactors=FALSE)
geneNames <- frames[[1]]$Geneid
for(z in 1:length(files))
{
  colnames(frames[[z]])[7] <- gsub(".counts","",files[z])
  frames[[z]] <- frames[[z]][,c(7)]
}
big_frame <- do.call(cbind,frames)
rownames(big_frame) <- geneNames
colnames(big_frame) <- gsub(".counts","",files)
countsMatrix <- big_frame
rm(big_frame)

write.csv(countsMatrix,
          file="C:/Users/mmari/OneDrive/Documents/GitHub/nextflow_workflow/output/combined_counts_matrix.csv",
          row.names=FALSE,
          quote=FALSE)

coldata <- read.table("C:/Users/mmari/OneDrive/Documents/GitHub/nextflow_workflow/output/samples_for_deseq2.txt")

#Do some low count pre-filtering: (we want 3 or more samples to have 5 or more genes)
dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]

dds <- DESeqDataSetFromMatrix(countData = countsMatrix,
                              colData = coldata,
                              design= ~ condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_2h_vs_24h")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="condition_2h_vs_24h", type="apeglm")

# omit NAs from results:
res <- na.omit(res)

res[,abs(res$log2FoldChange) > 1 & res$padj <= 0.05] 

res2v24Sig <- res[which((res$padj < 0.05) & abs(res$log2FoldChange) >= 1),]

res2v24Sig <- res2v24Sig[
  with(res2v24Sig , order(log2FoldChange, padj)),
]

res2v24SigDf <- as.data.frame(res2v24Sig)
res2v24SigDfTopDown <- res2v24SigDf[res2v24SigDf$log2FoldChange<=-1,] %>% arrange(padj) %>% top_n(100)
res2v24SigDfTopUp   <- res2v24SigDf[res2v24SigDf$log2FoldChange>=1,]  %>% arrange(padj) %>% top_n(100)

geneList = sort(geneList, decreasing = TRUE)

genesDown <- gseGO(geneList = rownames(res2v24SigDfTopDown),
              OrgDb         = org.Hs.eg.db,
              ont           = "CC",
              minGSSize     = 100,
              maxGSSize     = 500,
              pvalueCutoff  = 0.05,
              verbose       = FALSE)

genesUp <- gseGO(geneList = rownames(res2v24SigDfTopUp),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   minGSSize     = 100,
                   maxGSSize     = 500,
                   pvalueCutoff  = 0.05,
                   verbose       = FALSE)
