#!/usr/bin/env Rscript

##Mike Mariani UVM 2019-2020

##Input raw counts data files, 
##such as STAR --> featureCounts files
##then plot counts heatmap, boxplots,
##PCA, and tracks

#!/usr/bin/env Rscript

##Mike Mariani UVM 2019-2020

##Header material for HHV6 detection projected related materials

################################ Preprocessing ############################################

library(data.table)
library(ggfortify)
library(factoextra)
library(ggplot2)
library(cowplot)
library(patchwork)
library(dplyr)
library(rtracklayer)
library(R.4Cker)
library(seqsetvis)
library(ssvRecipes)
library(pdftools)
library(magick)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(IRanges)
library(ggdendro)
library(peakC)
library(stats)
library(ggforce)
library(DESeq2)

source("/slipstream/home/mmariani/scripts/read.files.mm.R")
source("/slipstream/home/mmariani/scripts/hhv6_detection_scripts/r_scripts/mm_bdg_clustering_method.R")
source("/slipstream/home/mmariani/scripts/bedgraph_to_wiggraph_mm.R")
source("/slipstream/home/mmariani/scripts/output_circos_mm.R")
source("/slipstream/home/mmariani/scripts/createR4CkerObjectFromFiles_mm_modify.R")
source("/slipstream/home/mmariani/scripts/createR4CkerObjectFromDFs_mm_modify.R")
source("/slipstream/home/mmariani/scripts/hhv6_detection_scripts/r_scripts/peakc_mm_bdg_to_peakc.R")
source("/slipstream/home/mmariani/scripts/domainograms_mm.R")
source("/slipstream/home/mmariani/scripts/domainograms_2_mm.R")
source("/slipstream/home/mmariani/scripts/ideograms_karyotypes_mm.R")
source("/slipstream/home/mmariani/scripts/annotate_chromhmm_mm.R")
source("/slipstream/home/mmariani/scripts/domainograms_mm.R")
source("/slipstream/home/mmariani/scripts/hhv6_detection_scripts/r_scripts/chomp_chrom_ends_mm.R")
source("/slipstream/home/mmariani/scripts/hhv6_detection_scripts/r_scripts/sum_scores_mm.R")
source("/slipstream/home/mmariani/scripts/hhv6_detection_scripts/r_scripts/chrom_test_anova_tukey_mm.R")
source("/slipstream/home/mmariani/scripts/hhv6_detection_scripts/r_scripts/chrom_test_wilcoxon_mm.R")
source("/slipstream/home/mmariani/scripts/hhv6_detection_scripts/r_scripts/transform_bw_to_data_frame_mm.R")
source("/slipstream/home/mmariani/scripts/hhv6_detection_scripts/r_scripts/peakc_mm_bdg_to_peakc.R")
source("/slipstream/home/mmariani/scripts/hhv6_detection_scripts/r_scripts/peakc_mm_peakc_to_regions.R")

set.seed(seed=1)

############################### Variables and Loading ##############################################

##For chromosome ends comparison functions defined below:
hg38.hhv6a.len.file <- "/slipstream/home/mmariani/references/ucsc_hg38_canonical_with_hhv6a_newest/ucsc_hg38_canonical_with_hhv6a_newest.len"
hg38.hhv6a.len <- read.table(file=hg38.hhv6a.len.file,
                             sep="\t",
                             stringsAsFactors = FALSE,
                             header = FALSE)
colnames(hg38.hhv6a.len) <- c("chrom", "length")
hg38.hhv6a.bac.len.file <- "/slipstream/home/mmariani/references/ucsc_hg38_canonical_with_hhv6a_gfp/ucsc_hg38_canonical_with_hhv6a_gfp.len"
hg38.hhv6a.bac.len <- read.table(file=hg38.hhv6a.bac.len.file,
                                 sep="\t",
                                 stringsAsFactors = FALSE,
                                 header = FALSE)
colnames(hg38.hhv6a.bac.len) <- c("chrom", "length")
##hg38.hhv6a.bac.len <- rbind(hg38.hhv6a.len,c("HHV-6A_BAC.",160111))
##chroms.original <- hg38.hhv6a.bac.len$chrom
chroms.to.be.labelled <- hg38.hhv6a.bac.len$chrom[!hg38.hhv6a.bac.len$chrom %in% c("chrM","HHV-6A_BAC.")]
chroms.p <- paste0(chroms.to.be.labelled,"_p")
chroms.q <- paste0(chroms.to.be.labelled,"_q")
chroms.all.tels <- c(chroms.p,chroms.q)

enz_file.smc <- read.table(file=paste0("/slipstream/home/mmariani/references/r4cker_ucsc_hg38_canonical_with_hhv6a_newest_60",
                                       "/ucsc_hg38_canonical_with_hhv6a_newest_hindiii_flanking_sites_60_2.bed"),
                           stringsAsFactors = FALSE)

enz_file.293 <- read.table(file=paste0("/slipstream/home/mmariani/references/r4cker_ucsc_hg38_canonical_with_hhv6a_gfp_60",
                                       "/ucsc_hg38_canonical_with_hhv6a_gfp_hindiii_flanking_sites_60_2.bed"),
                           stringsAsFactors = FALSE)

cytoband.hg38.hhv6a.gfp    <- "/slipstream/home/mmariani/references/cytoband_files/cytoband.ucsc.hg38.and.hhv6a.gfp.txt"
cytoband.hg38.hhv6a.newest <- "/slipstream/home/mmariani/references/cytoband_files/cytoband.ucsc.hg38.and.hhv6a.newest.txt"

chromosome.index.hg38.hhv6a.gfp=c(
  "chr1",
  "chr2",
  "chr3",
  "chr4",
  "chr5",
  "chr6",
  "chr7",
  "chr8",
  "chr9",
  "chr10",
  "chr11",
  "chr12",
  "chr13",
  "chr14",
  "chr15",
  "chr16",
  "chr17",
  "chr18",
  "chr19",
  "chr20",
  "chr21",
  "chr22",
  "chrM",
  "chrX",
  "chrY",
  "HHV-6A_BAC.")

chromosome.index.hg38.hhv6a.newest=c(
  "chr1",
  "chr2",
  "chr3",
  "chr4",
  "chr5",
  "chr6",
  "chr7",
  "chr8",
  "chr9",
  "chr10",
  "chr11",
  "chr12",
  "chr13",
  "chr14",
  "chr15",
  "chr16",
  "chr17",
  "chr18",
  "chr19",
  "chr20",
  "chr21",
  "chr22",
  "chrM",
  "chrX",
  "chrY",
  "chrHHV6A")

source("/slipstream/home/mmariani/scripts/hhv6_detection_scripts/r_scripts/hhv6_detection_header.R")

######################################## FIGURE 1 #############################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

output.dir <- "/slipstream/home/mmariani/projects/hhv6_detection/figure1_final"
setwd(output.dir)

############################ Timecourse RNA-Seq data ##########################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

hhv6a.gfp.genes.rev <- c(
  "U100",
  "U95",
  "U94",
  "U91",
  "U90",
  "U86",
  "U85",
  "U84",
  "U83",
  "U82",
  "U81",
  "U79",
  "U77",
  "U76",
  "U75",
  "U74",
  "U73",
  "U72",
  "U71",
  "U70",
  "U69",
  "U68",
  "U67",
  "U60",
  "U65",
  "U64",
  "U63",
  "U62",
  "U59",
  "U58",
  "U57",
  "U56",
  "U55",
  "U54",
  "EcoGPT",
  "GFP",
  "U53.5",
  "U53",
  "U52",
  "U51",
  "U50",
  "U49",
  "U48",
  "U47A",
  "U47",
  "U46",
  "U45",
  "U44",
  "U43",
  "U42",
  "U41",
  "U40",
  "U39",
  "U38",
  "U37",
  "U36",
  "U35",
  "U34",
  "U33",
  "U32",
  "U31",
  "U30",
  "U29",
  "U28",
  "U27",
  "U26",
  "U25",
  "U24A_1",
  "U24A",
  "U23",
  "U22",
  "U21",
  "U20",
  "U19",
  "U18",
  "U17",
  "U15",
  "U14",
  "U13",
  "U12",
  "U11",
  "U10",
  "U7",
  "U4",
  "U3",
  "U2"
)

hhv6a.gfp.genes <- rev(hhv6a.gfp.genes.rev)

input.dir <- "/slipstream/home/mmariani/projects/hhv6_detection/hhv6a_rna/hhv6a_time_course/star_output_hhv6a_gfp"
rna.counts.files <- list.files(input.dir, pattern="L002_R1_001.pe.v.Aligned.sortedByCoord.out.counts$",full.names = TRUE)

rna.counts.frames <- lapply(rna.counts.files,
                            read.table,
                            header=FALSE,
                            stringsAsFactors=FALSE,
                            sep="\t",
                            skip=2)
for(i in 1:length(rna.counts.files)){
  rna.counts.frames[[i]]$sample <- gsub("_L002_R1_001.pe.v.Aligned.sortedByCoord.out.counts","",basename(rna.counts.files[i]))
  print(ncol(rna.counts.frames[[i]]))
}

rna.big.frame.in <- do.call(rbind,rna.counts.frames[1:7])
setDT(rna.big.frame.in)

sort.samples <- gtools::mixedsort(unique(rna.big.frame.in$sample))

##countToFPKM::fpkm(counts,
##                  featureLength = ,
##                  meanFragmentLength = 

##can use picard to get mean frafment length:
##java -jar picard.jar CollectInsertSizeMetrics \
##I=input.bam \
##O=insert_size_metrics.txt \
##H=insert_size_histogram.pdf \
##M=0.5

mmm <- rna.big.frame.in[,c("V1","V7","sample")]
mmm.casted <- dcast(mmm,V1~sample,value.var="V7")
mmm.rownames <- as.character(mmm.casted$V1)
mmm.casted <- mmm.casted[,-1]
rownames(mmm.casted) <- mmm.rownames

m <- as.matrix(mmm.casted,dimnames=list(rownames(mmm.casted),colnames(mmm.casted)))

################## Raw RNA-Seq counts heatmap ##########################################
########################################################################################
########################################################################################
########################################################################################
########################################################################################

rna.counts <- m

rna.big.frame.long <- as.data.frame(rna.counts)
rna.big.frame.long$gene <- mmm.rownames
rna.big.frame <- melt(rna.big.frame.long, by="gene")
colnames(rna.big.frame) <- c("gene","sample","count")
setDT(rna.big.frame)

rna.big.frame[sample %in% sort.samples[1:2],time:="early",]
rna.big.frame[sample %in% sort.samples[3:5],time:="mid",]
rna.big.frame[sample %in% sort.samples[6:7],time:="late",]

sub.rna.frame <- subset(rna.big.frame, time %in% c("early", "mid", "late"))
sub.rna.frame <- subset(sub.rna.frame, !(gene %in% c("DR1","DR1_1","DR6","DR6_1")))
sub.rna.frame$time <- factor(sub.rna.frame$time, levels = unique(sub.rna.frame$time))
sub.rna.frame$gene <- factor(sub.rna.frame$gene, levels = rev(hhv6a.gfp.genes))

##remove U53.5 gene
##sub.rna.frame <- sub.rna.frame[!(sub.rna.frame$gene %in% "U53.5"),]

casted.rna <- dcast(sub.rna.frame, gene~sample, value.var=("count"))
dat.rna <- as.matrix(casted.rna[,2:7])  # numerical columns
rownames(dat.rna) <- casted.rna$gene
row.order <- hclust(dist(dat.rna))$order
col.order <- hclust(dist(t(dat.rna)))$order
dat.rna.new <- dat.rna[row.order, col.order] ## re-order matrix accoring to clustering
##melted.rna <- melt(dat.rna.new, id="gene")
##colnames(metled.rna) <- c("gene", "sample", "count")

hc.rna <- hclust(dist(dat.rna))
hcd.rna <- as.dendrogram(hc.rna)
col.clust.rna.plot <- plot(hcd.rna, leaflab="none")

hc2.rna <- hclust(dist(t(dat.rna)))
hcd2.rna <- as.dendrogram(hc2.rna)
plot(hcd2.rna, leaflab="none")

dendr.col.rna <- dendro_data(hcd2.rna, type="rectangle")
dendr.row.rna <- dendro_data(hcd.rna, type="rectangle") 

dendro.row.rna.plot <- ggplot() +
  geom_segment(data=ggdendro::segment(dendr.row.rna), aes(x=x, y=y, xend=xend, yend=yend)) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(plot.margin=margin(0,0,0,0,"cm")) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border=element_blank()) +
  scale_x_continuous(expand = c(0.005,0.005)) +
  theme(plot.margin=margin(0,0,0.5,0,"cm")) +
  coord_flip() 

dendro.col.rna.plot <- ggplot() +
  geom_segment(data=ggdendro::segment(dendr.col.rna), aes(x=x, y=y, xend=xend, yend=yend)) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border=element_blank()) ##+
##scale_x_continuous(expand = c(0.0,0.0)) +
##theme(plot.margin=margin(0,0,0,0,"cm"))

heatmap.rna <- ggplot(sub.rna.frame, aes(x=gsub("PBS_|_.*","",sample),y=gene,fill=ifelse(count==0,1,log2(count+1))))+
  geom_tile() +
  scale_fill_gradient2(low="blue",mid="yellow",high="red",midpoint = median(log2(sub.rna.frame$count+1))) +
  theme(axis.text.x = element_text(angle=90)) +
  ##theme(legend.position="none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5)) +
  ##scale_fill_continuous(limits=c(0, 300), breaks=seq(0,300,by=50)) +
  theme(legend.position = c(1.2, 1.1),
        panel.border=element_blank(),
        panel.grid = element_blank()) +
  ylab("gene") +
  labs(fill = "log2(count+1)") +
  xlab("sample") +
  theme(plot.margin=margin(0,0,0,0,"cm")) ##+
##theme(margin

heatmap.dendro.rna <- cowplot::plot_grid(
  dendro.col.rna.plot, 
  plot_spacer()+theme_dendro(), 
  heatmap.rna, 
  dendro.row.rna.plot, 
  nrow=2, 
  ncol=2, 
  rel_heights = c(1,5),
  rel_widths = c(4,1))##,
##scale = c(0.8,1,1,1),
##align="v",
##axis="lr")

heatmap.dendro.rna <- (dendro.col.rna.plot | plot_spacer() + theme_bw()) / (heatmap.rna | dendro.row.rna.plot) +
  plot_layout(nrow=2, heights=c(1,4))

boxplot.rna <- ggplot(sub.rna.frame, aes(x=time,y=ifelse(count==0,1,log2(count+1)),fill=time)) +
  geom_boxplot() +
  ylab("count") +
  theme_bw()

pc.casted <- dcast(sub.rna.frame,gene~sample,value.var=c("count"))
casted.pr <- prcomp(pc.casted[,c(2:8)], center=TRUE, scale=TRUE)
summary(casted.pr)
pca.frame <- data.frame(casted.pr$rotation, stringsAsFactors = FALSE)
pca.frame$time <- c("early", "early", "mid", "mid", "mid", "late", "late")
pca.frame$time <- factor(pca.frame$time, levels=c("early","mid","late"))
pca_data_perc=round(100*(casted.pr$sdev)^2/sum((casted.pr$sdev)^2),1)

library(ggforce)
pca.rna <- ggplot(pca.frame, aes(PC1,PC2, color = time)) +
  geom_point(size=5) +
  ##stat_ellipse(alpha = 0.5) +
  ##geom_ellipse(aes(x0 = pca.frame, y0 = PC2, a = 2, b = 1, angle = -pi / 3, m1 = 3)) +
  ##geom_ellipse(aes(x0 = mean(pca.frame[pca.frame$time=="early",]$PC1), y0=mean(pca.frame[pca.frame$time=="early",]$PC2), a = 0.005, b = 0.25, angle = 45*pi/180),color="red") +
  ##ggforce::geom_mark_ellipse(aes(fill=time, color=time)) ##label=time)) +
  ##ggforce::geom_mark_ellipse(aes(fill=time, label=time), expand = unit(5, "mm"), radius=2) +
  labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")")) +
  theme_bw()

figure1.rna.raw.plot <- heatmap.dendro.rna | (boxplot.rna / pca.rna)

ggsave(filename = paste0(output.dir,"/","figure1.rna.counts.raw.pdf"),
       height=20,
       width=20,
       device="pdf",
       plot=figure1.rna.raw.plot)