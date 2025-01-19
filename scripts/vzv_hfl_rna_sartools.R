#!/usr/bin/env Rscript

library("SARTools")
library("openxlsx")

work.dir <- "/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/sartools_edger"
work.dir <- "/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/sartools_deseq2"

setwd(work.dir)

##counts.df <- openxlsx::read.xlsx(paste0(work.dir,"FINAL RAW COUNTS.xlsx"), 
##                  sheet = 1, 
##                  startRow = 1, 
##                  colNames = TRUE)

##target.frame.gcap.vzvp <- openxlsx::read.xlsx(paste0(work.dir,"/RAW COUNTS SEPARATED.xlsx"), 
##                                              sheet = 2, 
##                                              startRow = 1, 
##                                              colNames = TRUE)

##target.frame.gcap.vzvn <- openxlsx::read.xlsx(paste0(work.dir,"/RAW COUNTS SEPARATED.xlsx"), 
##                                              sheet = 3, 
##                                              startRow = 1, 
##                                              colNames = TRUE)

##target.frame.norm.vzvp <- openxlsx::read.xlsx(paste0(work.dir,"/RAW COUNTS SEPARATED.xlsx"), 
##                                              sheet = 4, 
##                                              startRow = 1, 
##                                              colNames = TRUE)

##target.frame.norm.vzvn <- openxlsx::read.xlsx(paste0(work.dir,"/RAW COUNTS SEPARATED.xlsx"),
##                                              sheet = 5, 
##                                              startRow = 1, 
##                                              colNames = TRUE)

##Read in featureCounts combined file.
##featureCounts was run after STAR alignment

hg38.feature.counts.combined <- read.table(file="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_featureCounts/vzv_mock_rna_re-seq_star_feature_counts_combined_11092019.txt",
                                   sep="\t",
                                   header=TRUE,
                                   stringsAsFactors = FALSE,
                                   skip = 1
                                   )

hg38.counts.combined <- hg38.feature.counts.combined[,c(1,7:12)]

hg38.counts.mock <- hg38.counts.combined[,c(1,2:4)]
hg38.counts.infected <- hg38.counts.combined[,c(1,5:7)] 

colnames(hg38.counts.mock) <- c("gene","InfectedS1","InfectedS2","InfectedS3")
colnames(hg38.counts.infected) <- c("gene","MockS4","MockS5","MockS6")

length(c(colnames(hg38.counts.mock)[-1],
         colnames(hg38.counts.infected)[-1]
))

duplicated(c(colnames(hg38.counts.mock)[-1],
             colnames(hg38.counts.infected)[-1]
          ))

any(duplicated(hg38.counts.mock$gene))
any(duplicated(hg38.counts.infected$gene))
##FALSE
##FALSE

## Ouptut the counts matrices by treatment
## as \t text files:

##Write counts files ot both the edger and deseq2 output folders:

write.table(x=hg38.counts.infected,
            file = paste0("/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/sartools_edger",
                          "/infected.hg38.counts.txt"),
            col.names=FALSE,
            row.names=FALSE,
            quote=FALSE,
            sep="\t")

write.table(x=hg38.counts.mock,
            file = paste0("/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/sartools_edger",
                          "/mock.hg38.counts.txt"),
            col.names=FALSE,
            row.names=FALSE,
            quote=FALSE,
            sep="\t")

write.table(x=hg38.counts.infected,
            file = paste0("/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/sartools_deseq2",
                          "/infected.hg38.counts.txt"),
            col.names=FALSE,
            row.names=FALSE,
            quote=FALSE,
            sep="\t")

write.table(x=hg38.counts.mock,
            file = paste0("/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/sartools_deseq2",
                          "/mock.hg38.counts.txt"),
            col.names=FALSE,
            row.names=FALSE,
            quote=FALSE,
            sep="\t")

nrow(hg38.counts.mock)
nrow(hg38.counts.infected)
##26485 for both

##SARTools::
##first column: unique names of the samples (short but informative 
##as they will be displayed on all the figures);
##second column: name of the count files;
##third column: biological conditions;
##optional columns: further information about 
##the samples (day of library preparation for example).

target.frame.hg38.mock <- data.frame(names=colnames(hg38.counts.mock)[-1],
                                     file="mock.hg38.counts.txt",
                                     condition="mock",
                                     ##optional columns,
                                     stringsAsFactors = FALSE
                                     )
          
target.frame.hg38.infected <- data.frame(names=colnames(hg38.counts.infected)[-1],
                                         file="infected.hg38.counts.txt",
                                         condition="infected",
                                         ##optional columns
                                         stringsAsFactors = FALSE
                                        )

target.frame.hg38.total <- rbind(target.frame.hg38.mock,
                                 target.frame.hg38.infected)
                                 
colnames(target.frame.hg38.total) <- c("label","files","group")

##Write target frame to both edger and deseq2 folders:

write.table(x=target.frame.hg38.total,
            file = paste0("/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/sartools_edger/target.frame.hg38.total.txt"),
            sep="\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

write.table(x=target.frame.hg38.total,
            file = paste0("/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/sartools_deseq2/target.frame.hg38.total.txt"),
            sep="\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

##Run both edger and deseq2 versions of sartools:

source("/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/sartools_edger/hfl_script_edger.r")

##source("/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/sartools_deseq2/hfl_script_deseq2.r")
