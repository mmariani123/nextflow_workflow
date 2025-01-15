#!/usr/bin/env Rscript

##setwd("/slipstream/home/mmariani/projects/hhv6_rna/output/bowtie2_counts")
setwd("/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38")

files <- list.files(path="./", pattern=".counts.counts")
frames <- lapply(files,read.table,sep="\t",stringsAsFactors=FALSE)
for(z in 1:length(files))
{
  frames[[z]]$sample <- gsub(".counts.counts","",files[z])
}
big_frame <- do.call(rbind,frames)
colnames(big_frame) <- c("gene","count","sample")

big_cast <- data.table::dcast(big_frame, gene ~ sample, value.var="count")

write.csv(big_cast,
          file="combined.counts.for.pivot.csv",
          row.names=FALSE,
          quote=FALSE)

frames_counts <- list()
frames_summary <- list()
for(i in 1:length(frames))
{
  frames_counts[[i]] <- head(frames[[i]],n=-5)
  frames_counts[[i]]$sample <- rep(files[i],times=nrow(frames_counts[[i]]))
  frames_summary[[i]] <- tail(frames[[i]],n=5)
  frames_summary[[i]]$sample <- rep(files[i],times=nrow(frames_summary[[i]]))
}
big_counts <- do.call(rbind,frames_counts)
big_summary <- do.call(rbind,frames_summary)

colnames(big_counts) <- c("gene","count","sample")  
colnames(big_summary) <- c("feature","count","sample")  

big_counts <- big_counts[,c(3,1,2)]
big_summary <- big_summary[,c(3,1,2)]

big_counts$sample <- gsub("\\..*","",big_counts$sample)
big_summary$sample <- gsub("\\..*","",big_summary$sample)

samples_list <- c("TZ_9_D1_8_S8_L002_R1_001",
"TZ_9_D2_9_S9_L002_R1_001",
"TZ_10_D3_12_S10_L002_R1_001",
"PBS_D1_1_S1_L002_R1_001",
"PBS_D2_2_S2_L002_R1_001",
"PBS_D3_3_S3_L002_R1_001",
"PBS_D4_4_S4_L002_R1_001",
"PBS_D5_5_S5_L002_R1_001",
"PBS_D6_6_S6_L002_R1_001",
"PBS_D7_7_S7_L002_R1_001")

genes_list <- c("DR1",
"DR1_1",
"DR6",
"DR6_1",
"U2",
"U3",
"U4",
"U7",
"U10",
"U11",
"U12",
"U13",
"U14",
"U15",
"U17",
"U18",
"U19",
"U20",
"U21",
"U22",
"U23",
"U24A",
"U24A_1",
"U25",
"U26",
"U27",
"U28",
"U29",
"U30",
"U31",
"U32",
"U33",
"U34",
"U35",
"U36",
"U37",
"U38",
"U39",
"U40",
"U41",
"U42",
"U43",
"U44",
"U45",
"U46",
"U47",
"U47A",
"U48",
"U49",
"U50",
"U51",
"U52",
"U53",
"U53.5",
"GFP",
"EcoGPT",
"U54",
"U55",
"U56",
"U57",
"U58",
"U59",
"U62",
"U63",
"U64",
"U65",
"U60",
"U67",
"U68",
"U69",
"U70",
"U71",
"U72",
"U73",
"U74",
"U75",
"U76",
"U77",
"U79",
"U81",
"U82",
"U83",
"U84",
"U85",
"U86",
"U90",
"U91",
"U94",
"U95",
"U100")

##big_summary$sample <- factor(big_summary$sample,levels(factor(samples_list)))
#big_counts <- big_counts[order(factor(big_counts$sample, levels = factor(samples_list))),]
#big_summary <- big_summary[order(factor(big_summary$sample, levels = factor(samples_list))),]
#big_counts$sample <- as.factor(big_counts$sample)
#big_summary$sample <- as.factor(big_summary$sample)
#temp$size_f = factor(temp$size, levels=c('50%','100%','150%','200%'))

big_counts$sample <- factor(big_counts$sample, levels=samples_list)
big_summary$sample <- factor(big_summary$sample, levels=samples_list)
big_counts$gene <- factor(big_counts$gene, levels=genes_list)

ggplot(big_counts, aes(x=gene,y=count)) +
  geom_col() +
  #facet_grid(rows = vars(sample)) +
  facet_grid(sample ~ .) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.text.y = element_text(angle = 0))
