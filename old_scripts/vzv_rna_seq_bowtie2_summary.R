#!/usr/bin/env Rscript 

setwd("/slipstream/home/mmariani/projects/hhv6_rna/output/bowtie2_counts")

files <- list.files(path="./",pattern=".bowtie2.log")
tables <- lapply(files,read.table,sep="\n",stringsAsFactors = FALSE)
for(i in 1:length(tables))
{
  tables[[i]]$sample <- files[i] 
}
big_table <- do.call(rbind,tables)
colnames(big_table) <- c("info","sample")

big_table$info <- gsub("\\s+reads.*","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\s+were\\s+paired.*","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\saligned\\s+","",big_table$info,perl=TRUE)
big_table$info <- gsub("exactly\\s+","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\s+aligned\\s+concordantly\\s+","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\s+aligned\\s+discordantly\\s+","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\s+times; of these:","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\s+pairs aligned\\s+","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\s+times concordantly or discordantly; of these:","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\s+mates make up the pairs; of these:","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\s+overall alignment rate","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\s+times","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\s+time","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\s+pairs","",big_table$info,perl=TRUE)
big_table$info <- gsub("concordantly\\s+","",big_table$info,perl=TRUE)
big_table$info <- gsub("discordantly\\s+","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\)0","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\)1","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\)>1","",big_table$info,perl=TRUE)
big_table$info <- gsub("\\)","",big_table$info,perl=TRUE)
big_table$info <- gsub("%","",big_table$info,perl=TRUE)

big_table <- big_table[!grepl("----",big_table$info),]

percents <- numeric(0)
splits <- strsplit(big_table$info,split="\\s+\\(",perl=TRUE)
for(i in 1:length(splits))
{
  percents[i] <- as.numeric(splits[[i]][2])
}
percents_table <- data.frame(info=percents,sample=big_table$sample)
percents_table_clean <- na.omit(percents_table)

big_table$info <- gsub("\\s+\\(.*","",big_table$info,perl=TRUE)

big_frame <- as.data.frame(matrix(as.numeric(big_table$info),nrow = 13))
percents_frame <- as.data.frame(matrix(percents_table_clean$info,nrow = 8))

colnames(big_frame) <- unique(gsub("\\.bowtie2\\.log","",big_table$sample,perl=TRUE))
colnames(percents_frame) <- unique(gsub("\\.bowtie2\\.log","",big_table$sample,perl=TRUE))

total_frame <- rbind(big_frame,percents_frame)

fields <- c(
"total_reads",
"paired_reads",
"aligned_concordantly_0_times",
"aligned_concordantly_exactly_1_time",
"aligned_concordantly_>1_times",
"pairs_aligned_concordantly_0_times",
"pairs_aligned_discordantly_1_time",
"pairs_aligned_0_times_concordantly_or_discordantly",
"mates_make_up_the_pairs",
"mates_aligned_0_times",
"mates_aligned_exactly_1_time",
"mates_aligned_>1_times",
"percent_overall_alignment_rate",
"percent_paired",
"percent_aligned_concordantly_0_times",
"percent_aligned_concordantly_exactly_1_time",
"percent_aligned_concordantly_>1_times",
"percent_pairs_aligned_discordantly_1_time",
"percent_mates_aligned_0_times",
"percent_mates_aligned_exactly_1_time",
"percent_mates_aligned_>1_times"
)

field_levels <- c(
  "total_reads",
  "percent_overall_alignment_rate",
  "paired_reads",
  "percent_paired",
  "aligned_concordantly_0_times",
  "percent_aligned_concordantly_0_times",
  "aligned_concordantly_exactly_1_time",
  "percent_aligned_concordantly_exactly_1_time",
  "aligned_concordantly_>1_times",
  "percent_aligned_concordantly_>1_times",
  "pairs_aligned_concordantly_0_times",
  "pairs_aligned_discordantly_1_time",
  "percent_pairs_aligned_discordantly_1_time",
  "pairs_aligned_0_times_concordantly_or_discordantly",
  "mates_make_up_the_pairs",
  "mates_aligned_0_times",
  "percent_mates_aligned_0_times",
  "mates_aligned_exactly_1_time",
  "percent_mates_aligned_exactly_1_time",
  "mates_aligned_>1_times",
  "percent_mates_aligned_>1_times"
)

sample_levels <- c(
  "TZ_9_D1_8_S8_L002_R1_001",
  "TZ_9_D2_9_S9_L002_R1_001", 
  "TZ_10_D3_12_S10_L002_R1_001", 
  "PBS_D1_1_S1_L002_R1_001",
  "PBS_D2_2_S2_L002_R1_001",  
  "PBS_D3_3_S3_L002_R1_001",
  "PBS_D4_4_S4_L002_R1_001",
  "PBS_D5_5_S5_L002_R1_001",  
  "PBS_D6_6_S6_L002_R1_001", 
  "PBS_D7_7_S7_L002_R1_001"
)

total_frame$fields <- fields

total_frame$fields <- factor(total_frame$fields, levels=field_levels)

total_frame <- total_frame[,c(11,9,10,8,1,2,3,4,5,6,7)]

total_melted <- melt(total_frame, id = c("fields")) 

total_melted$variable <- factor(total_melted$variable, levels=sample_levels)  

cat(big_table$info,sep="\n")

pdf(file="/slipstream/home/mmariani/projects/hhv6_rna/output/bowtie2.summary.pdf",
    height=24,
    width=24)
ggplot(data=total_melted,aes(x=variable,y=value,fill=variable)) +
  geom_col() +
  facet_grid(rows=vars(fields),scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(strip.text.y = element_text(angle = 0)) +
  ylab("value") +
  xlab("sample") +
  ggtitle("HHV6-GFP RNA-Seq Time-course Bowtie2 alignment summary (--no-unal)") +
  theme(legend.position="none")
dev.off()

#################################### Example Bowtie2 Metrics ######################################

#11342370 reads; of these:

#11342370 (100.00%) were paired; of these:

#10865188 (95.79%) aligned concordantly 0 times

#461204 (4.07%) aligned concordantly exactly 1 time

#15978 (0.14%) aligned concordantly >1 times

#----

#10865188 pairs aligned concordantly 0 times; of these:

#39560 (0.36%) aligned discordantly 1 time

#----

#10825628 pairs aligned 0 times concordantly or discordantly; of these:

#21651256 mates make up the pairs; of these:

#21596100 (99.75%) aligned 0 times

#53095 (0.25%) aligned exactly 1 time

#2061 (0.01%) aligned >1 times

#4.80% overall alignment rate

######################################## EXTRA ##################################################
#big_cast <- dcast(big_table, sample ~ info, fun=print)
#
#df.new = big_table[seq(1, nrow(big_table), 15), ]
#df.new$info
#
#write.table(big_table,
#            file="/slipstream/home/mmariani/projects/hhv6_rna/output/bowtie2.summary.csv",
#            sep=",",
#            row.names=FALSE,
#            col.names=TRUE,
#            quote=FALSE)
