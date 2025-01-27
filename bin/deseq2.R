#!/usr/bin/env Rscript

library(magrittr)
library(dplyr)
library(tidyr)
library(DESeq2)
library(apeglm)
library(pheatmap)
library(org.Hs.eg.db)
library(ggplot2)
library(cowplot)
library(patchwork)

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

files     <- list.files(path=inputDir, pattern=".counts", full.names = TRUE)
frames    <- lapply(files,read.table,sep="\t",header=TRUE,skip=1,stringsAsFactors=FALSE)
geneNames <- frames[[1]]$Geneid
for(z in 1:length(files))
{
  colnames(frames[[z]])[7] <- gsub(".counts","",basename(files[z]))
  frames[[z]] <- frames[[z]][,c(7)]
}
big_frame <- do.call(cbind,frames)
rownames(big_frame) <- geneNames
colnames(big_frame) <- gsub(".counts","",basename(files))
countsMatrix <- big_frame
rm(big_frame)

#Read sample file for DESeq2:

#coldata <- read.table(paste0(args[2],"/samples_for_deseq2.txt"))
coldata <- read.table(args[2],
                      header = TRUE,
                      sep = "\t",
                      stringsAsFactors = FALSE)

#Write counts matrix for DESeq2:

#write.csv(countsMatrix,
#          file=paste0(args[2],"/combined_counts_matrix.csv"),
#          row.names=FALSE,
#          quote=FALSE)

#write.csv(countsMatrix,
#          file="combined_counts_matrix.csv",
#          row.names=FALSE,
#         quote=FALSE)

write.csv(countsMatrix,
          file=args[3],
          row.names=TRUE,
          quote=FALSE)

dds <- DESeqDataSetFromMatrix(countData = countsMatrix,
                              colData   = coldata,
                              design    = ~ condition)

#Do some low count pre-filtering: (we want 3 or more samples to have 5 or more genes)
dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]

dds <- DESeq(dds)

#resultsNames(dds) # lists the coefficients

res <- results(dds, name="condition_2h_vs_24h")

# or to shrink log fold changes association with condition:
# NB: type='apeglm' requires installing the Bioconductor package 'apeglm'
res <- lfcShrink(dds, coef="condition_2h_vs_24h", type="apeglm")

#get results as dataframe:
res.df <- res@listData
res.df <- do.call(cbind,res.df)
# omit NAs from results:

res.df <- na.omit(res.df)

write.table(res.df,
            file = args[4],
            sep="\t")

#log2(n + 1) transform
ntd <- normTransform(dds)

#Create heatmap
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

#df <- as.data.frame(colData(dds)[,c("condition","type")])
df <- as.data.frame(colData(dds)[,"condition"])
colnames(df) <- "condition"
rownames(df) <- rownames(colData(dds))

pheatmap(assay(ntd)[select,], 
         cluster_rows=FALSE, 
         show_rownames=TRUE,
         cluster_cols=FALSE, 
         annotation_col=df,
         silent = TRUE,
         filename = args[5],
         width=6,
         height=6)

#Create PCA plot
vsd <- vst(dds, blind=FALSE)

#plotPCA(vsd, intgroup=c("condition", "type"))
plotPCA(vsd, intgroup="condition")

ggsave(filename = args[6],
       device="png",
       width=6,
       height=6)

#res.df[,abs(res.df$log2FoldChange) > 1 & res.df$padj <= 0.05] 

#res2v24Sig <- res.df[which((res.df$padj < 0.05) & abs(res.df$log2FoldChange) >= 1),]

#res2v24Sig <- res2v24Sig[
#  with(res2v24Sig , order(log2FoldChange, padj)),
#]
