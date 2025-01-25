#!/usr/bin/env Rscript


library(DESeq2)

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

dds <- DESeqDataSetFromMatrix(countData = countsMatrix,
                              colData = coldata,
                              design= ~ condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_2h_vs_24h")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="condition_2h_vs_24h", type="apeglm")

