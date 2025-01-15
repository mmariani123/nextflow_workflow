#!/usr/bin/env Rscript

##Mike Mariani 2019
##cluster profiler

#!/usr/bin/env Rscript

## Mike Mariani UVM 2019

## 11/10/2019

library(clusterProfiler)
library(org.Hs.eg.db)
library(qpcR)

##From RNA re-Seq of VZV HFL infected over mock, padj <= 0.05 lgfc>=1
deseq.in <- read.table(file="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_debrowser/up+down.csv",
                       header=TRUE,
                       stringsAsFactors = FALSE,
                       sep=",")

##Seth wants me to take top 1000 genes for each up and down:
deseq.up <- deseq.in[(deseq.in$padj <= 0.05 & deseq.in$log2FoldChange >= 1),]
deseq.down <- deseq.in[(deseq.in$padj <= 0.05 & deseq.in$log2FoldChange <= -1),]

deseq.up.top <- deseq.up[order(deseq.up$padj),][1:1000,]
deseq.down.top <- deseq.down[order(deseq.down$padj),][1:1000,]

##asign genes:
deseq.up.top.EG_IDs = mget(as.character(unlist(deseq.up.top)), revmap(org.Hs.egSYMBOL),ifnotfound=NA)
deseq.down.top.EG_IDs = mget(as.character(unlist(deseq.down.top)), revmap(org.Hs.egSYMBOL),ifnotfound=NA)

##Get genes:
deseq.up.top.entrez  <- unlist(unname(deseq.up.top.EG_IDs[!is.na(deseq.up.top.EG_IDs)]))
deseq.down.top.entrez  <- unlist(unname(deseq.down.top.EG_IDs[!is.na(deseq.down.top.EG_IDs)]))

##Can run Kegg individually or skip this part and jump below:
kegg.markers.up <- enrichKEGG(deseq.up.top.entrez, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
kegg.markers.down <- enrichKEGG(deseq.down.top.entrez, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
barplot(kegg.markers.up, showCategory=nrow(kegg.markers.up))
barplot(kegg.markers.down, showCategory=nrow(kegg.markers.down))

##formula interface, combine all upregulated genes by cluster
##into a single data.frame:
mydf <- data.frame(
  Entrez = c(deseq.up.top.entrez,
             deseq.down.top.entrez
  ),
  group = c(rep("Up" ,times=length(deseq.up.top.entrez)), 
            rep("Down" ,times=length(deseq.down.top.entrez))
  )
)

##Run Kegg Pathway Analysis Between Clusters:
xx.formula.kegg.combined <- compareCluster(Entrez~group, 
                                           data=mydf,
                                           fun='enrichKEGG'##, 
                                           ##OrgDb='org.Hs.eg.db'
)

xx.formula.go.combined <- compareCluster(Entrez~group, 
                                data=mydf,
                                fun='enrichGO', 
                                OrgDb='org.Hs.eg.db'
)

xx.formula.group.go.combined <- compareCluster(Entrez~group, 
                                   data=mydf,
                                   fun='groupGO', 
                                   OrgDb='org.Hs.eg.db'
)

##Dot Plot of Significant Pathways by up and down deseq expression:

pdf.output.dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/cluster_profiler_output"

pdf(file=paste0(pdf.output.dir,"/up.down.kegg.dot.05.pdf"),
      height=8,
      width=8)
dotplot(xx.formula.kegg.combined, 
        showCategory=length(unique(xx.formula.kegg.combined@compareClusterResult[xx.formula.kegg.combined@compareClusterResult$p.adjust<=0.05,"Description"])),
        split="Cluster",
        color="p.adjust")
dev.off()

pdf(file=paste0(pdf.output.dir,"/up.down.go.dot.05.pdf"),
    height=8,
    width=8)
dotplot(xx.formula.go.combined, 
        showCategory=length(unique(xx.formula.go.combined@compareClusterResult[xx.formula.go.combined@compareClusterResult$p.adjust<=0.05,"Description"])),
        split="Cluster",
        color="p.adjust")
dev.off()

nrow(xx.formula.group.go.combined@compareClusterResult)
pdf(file=paste0(pdf.output.dir,"/up.down.group.go.pdf"),
    height=8,
    width=8)
dotplot(xx.formula.group.go.combined, 
        showCategory=nrow(xx.formula.group.go.combined@compareClusterResult),
        split="Cluster",
        color="p.adjust")
dev.off()

##Doesn't work
##xx.formula.go.kegg <- clusterProfiler::simplify(xx.formula.kegg)

xx.formula.go.combined.simplify <-  clusterProfiler::simplify(xx.formula.go.combined)

##Doesn't work:
##xx.formula.group.simplfy <- clusterProfiler::simplify(xx.formula.group)

length(unique(xx.formula.go.combined.simplify@compareClusterResult[xx.formula.go.combined.simplify@compareClusterResult$p.adjust<=0.05,"Description"]))

pdf(file=paste0(pdf.output.dir,"/up.down.go.simplified.dot.05.pdf"),
    height=8,
    width=8)
dotplot(xx.formula.go.combined.simplify, 
        showCategory=length(unique(xx.formula.go.combined.simplify@compareClusterResult[xx.formula.go.combined.simplify@compareClusterResult$p.adjust<=0.05,"Description"])),
        split="Cluster",
        color="p.adjust")
dev.off()

##I'm not sure what enrichPathway is for may be for comparing 
##multiple sets of pathways.

##xx.formula.pathway <- compareCluster(##Entrez~group,
##  data=mydf,
##  fun='enrichPathway'##, 
##  ##OrgDb='org.Hs.eg.db'
##)
##
##pathways <- c(kegg.markers.up, kegg.markers.down)
##
##barplot(pathways, nr=2, beside=T)
##
##dotplot(pathways, 
##        showCategory=20,
##        split=TRUE)

############# For particular Pathways result, get the genes and store as a data frame ############
##################################################################################################

##List all the pathways:
cat(xx.formula.go.combined.simplify@compareClusterResult$Description,sep="\n")

##convert results to data.frame is desired (might be easier to work with):
xx.go.combined.simplify.frame <- as.data.frame(xx.formula.go.combined.simplify)

##Pick your pathway by names(s) to get the genes:
desired.pathways<-c("DNA-binding transcription activator activity, RNA polymerase II-specific",
                    "cell adhesion molecule binding")

desired.genes <- xx.go.combined.simplify.frame[xx.go.combined.simplify.frame$Description %in% desired.pathways,"geneID"]

genes.list <- list()
for(i in 1:length(desired.genes))
{
  genes.now <- unlist(strsplit(desired.genes[i],split="/"))
  translate <- bitr(genes.now, fromType="ENTREZID", toType="SYMBOL", OrgDb=org.Hs.eg.db)
  genes.list[[i]] <- translate$SYMBOL
}

genes.out.frame <- data.frame(do.call(qpcR:::cbind.na,genes.list),
                              stringsAsFactors = FALSE)

colnames(genes.out.frame) <- desired.pathways

##Check output frame:
head(genes.out.frame, n=10)
