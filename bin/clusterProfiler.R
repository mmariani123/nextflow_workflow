#!/usr/bin/env Rscript

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
