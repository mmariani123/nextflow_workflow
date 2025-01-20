#!/usr/bin/env Rscript

##If counts files are separate:

setwd("/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output")

files.list <- list.files(path="./",pattern=".txt")

files.frames <- lapply(files.list,read.table,sep="\t",header=FALSE)
for(i in 1:length(files.frames))
{
  files.frames[[i]]$name <- files.list[i]
}

big.frame <- do.call(rbind,files.frames)

colnames(big.frame) <- c("gene","count","sample")

big.frame$sample <- gsub(".counts.counts","",big.frame$sample)

casted <- dcast(big.frame,gene ~ sample, value.var = "count")

write.csv(casted,file="vzv_rna-seq_resequenced.csv",col.names = TRUE,row.names = FALSE,quote = FALSE)

##11/09/2019
##If counts files are combined:

setwd("/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_featureCounts")

files.list.c <- list.files(path="./",pattern=".txt$")
files.list.c

big.frame.c <- read.table(files.list.c[[1]],
                          sep="\t",
                          header=TRUE,
                          stringsAsFactors = FALSE)

debrowser.frame <- big.frame.c[,c(1,7:12)]
colnames(debrowser.frame) <- c("gene",gsub("\\..*","",colnames(big.frame.c)[7:12]))

output.dir <- "/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_debrowser"

debrowser.frame <- debrowser.frame[,c(1,5,6,7,2,3,4)]
##Write DEBrowser counts file:
write.csv(debrowser.frame,
          file=paste0(output.dir,"/hg38_star_and_featureCounts_for_debrowser_11092019.csv"),
          col.names = TRUE,
          row.names = FALSE,
          quote = FALSE)

##Write DEBrowser meta-data file:

##Example:
##samples	treatment	batch
##exper_rep1	cond1	1
##exper_rep2	cond1	2
##exper_rep3	cond1	1
##control_rep1	cond2	2
##control_rep2	cond2	1
##control_rep3	cond2	2

meta.data.frame <- data.frame(
  samples = colnames(debrowser.frame)[2:7],
  treatment = c(rep("VZV",times=3),rep("MOCK",times=3)),
  batch = rep(1,times=6),
  stringsAsFactors = FALSE
)

write.csv(meta.data.frame,
          file=paste0(output.dir,"/hg38_star_and_featureCounts_metadata_for_debrowser_11092019.csv"),
          col.names = TRUE,
          row.names = FALSE,
          quote = FALSE)
