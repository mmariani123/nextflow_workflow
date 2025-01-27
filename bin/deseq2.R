#!/usr/bin/env Rscript

library(magrittr)
library(dplyr)
library(tidyr)
library(DESeq2)
#library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)
library(patchwork)
#library(ggpubr)

#counts <- featureCounts(bams, blablabla_restofcommand)$counts
#deseqdata <- DESeqDataSetFromMatrix(countData=counts, colData=sampleInfo, design=~condition)

##setwd("/slipstream/home/mmariani/projects/hhv6_rna/output/bowtie2_counts")

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
#if (length(args)==0) {
#  stop("No arguments applied to deseq2.R", call.=FALSE)
#}

inputDir  <- args[1]

countsMatrix = matrix(
  
  # Taking sequence of elements  
  c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
  
  # No of rows
  nrow = 3,   
  
  # No of columns
  ncol = 3,         
  
  # By default matrices are in column-wise order
  # So this parameter decides how to arrange the matrix
  byrow = TRUE          
)

write.csv(countsMatrix,
          file=args[3],
          row.names=FALSE,
          quote=FALSE)

##files     <- list.files(path=inputDir, pattern=".counts", full.names = TRUE)
##frames    <- lapply(files,read.table,sep="\t",header=TRUE,skip=1,stringsAsFactors=FALSE)
##geneNames <- frames[[1]]$Geneid
##for(z in 1:length(files))
##{
##  colnames(frames[[z]])[7] <- gsub(".counts","",basename(files[z]))
##  frames[[z]] <- frames[[z]][,c(7)]
##}
##big_frame <- do.call(cbind,frames)
##rownames(big_frame) <- geneNames
##colnames(big_frame) <- gsub(".counts","",basename(files))
##countsMatrix <- big_frame
##rm(big_frame)
##
###Read sample file for DESeq2:
##
###coldata <- read.table(paste0(args[2],"/samples_for_deseq2.txt"))
##coldata <- read.table(args[2])
##
###Write counts matrix for DESeq2:
##
###write.csv(countsMatrix,
###          file=paste0(args[2],"/combined_counts_matrix.csv"),
###          row.names=FALSE,
###          quote=FALSE)
##
##countsMatrix = matrix(
##  
##  # Taking sequence of elements  
##  c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
##  
##  # No of rows
##  nrow = 3,   
##  
##  # No of columns
##  ncol = 3,         
##  
##  # By default matrices are in column-wise order
##  # So this parameter decides how to arrange the matrix
##  byrow = TRUE          
##)
##
##write.csv(countsMatrix,
##          file=args[3],
##          row.names=FALSE,
##          quote=FALSE)
##
###write.csv(countsMatrix,
###          file="combined_counts_matrix.csv",
###          row.names=FALSE,
###          quote=FALSE)
##
###Do some low count pre-filtering: (we want 3 or more samples to have 5 or more genes)
##dds <- estimateSizeFactors(dds)
##idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
##dds <- dds[idx,]
##
##dds <- DESeqDataSetFromMatrix(countData = countsMatrix,
##                              colData = coldata,
##                              design= ~ condition)
##dds <- DESeq(dds)
##resultsNames(dds) # lists the coefficients
##res <- results(dds, name="condition_2h_vs_24h")
### or to shrink log fold changes association with condition:
##res <- lfcShrink(dds, coef="condition_2h_vs_24h", type="apeglm")
##
### omit NAs from results:
##res <- na.omit(res)
##
##res[,abs(res$log2FoldChange) > 1 & res$padj <= 0.05] 
##
##res2v24Sig <- res[which((res$padj < 0.05) & abs(res$log2FoldChange) >= 1),]
##
##res2v24Sig <- res2v24Sig[
##  with(res2v24Sig , order(log2FoldChange, padj)),
##]
##
##res2v24SigDf <- as.data.frame(res2v24Sig)
##
##genes.df <- bitr(toupper(rownames(res2v24SigDf)), 
##                    fromType = "SYMBOL",
##                    toType = "ENTREZID",
##                    OrgDb = org.Hs.eg.db
##)
##rownames(genes.df) <- genes.df$SYMBOL
###Remember we may not get complete mapping
##
##rownames(res2v24SigDf) <- toupper(rownames(res2v24SigDf))
##
##merged <- merge(res2v24SigDf,
##                genes.df,
##                by="row.names",
##                all.x=TRUE,
##                all.y=TRUE)
##
##merged <- na.omit(merged)
##
##mergedDown <- merged[merged$log2FoldChange<=-1,] %>% arrange(padj) #%>% top_n(100)
##mergedUp   <- merged[merged$log2FoldChange>=1,]  %>% arrange(padj) #%>% top_n(100)
##
##genesDown <- enrichGO(gene          = mergedDown$ENTREZID,
##                      OrgDb         = org.Hs.eg.db,
##                      keyType       = 'ENTREZID',
##                      ont           = "BP",
##                      pAdjustMethod = "BH",
##                      pvalueCutoff  = 0.01,
##                      qvalueCutoff  = 0.05)
##
##genesUp <- enrichGO(gene            = mergedUp$ENTREZID,
##                      OrgDb         = org.Hs.eg.db,
##                      keyType       = 'ENTREZID',
##                      ont           = "BP",
##                      pAdjustMethod = "BH",
##                      pvalueCutoff  = 0.01,
##                      qvalueCutoff  = 0.05)
##
###clusterProfiler::goplot(genesDown, environment())
###clusterProfiler::goplot(genesUp)
##
##plotDown <- dotplot(genesDown, showCategory=10)
##plotUp   <- dotplot(genesUp, showCategory=10)
##
##cowplot::plot_grid(plotDown,
##                   plotUp)
##
##plotDown + plotUp + patchwork::plot_layout(guides = "collect") +
##  scale_fill_continuous(limits = range(c(plotDown$data$p.adjust, plotUp$data$p.adjust)))
##
##pathwaysPlotOut <- ggarrange(
##  plotDown, 
##  plotUp,
##  align = "h", 
##  labels = c("Genes down-regulated pathways", "Genes up-regulated pathways"),
##  common.legend = TRUE,
##  legend="right"
##)
##
##ggsave(plot=pathwaysPlotOut,
##       #filename = paste0(args[2],"/pathways.png"),
##       filename = args[4],
##       #filename = "pathways.png",
##       width = 12,
##       height = 8,
##       device = "png")
