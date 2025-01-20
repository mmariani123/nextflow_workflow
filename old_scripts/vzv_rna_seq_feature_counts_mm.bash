#!/usr/bin/env bash

##Mike Mariani UVM 2019-2021

##11/09/2019 - 02/24/2021
##qsub -cwd -pe threads 20 -j yes -o /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/logs /slipstream/home/mmariani/scripts/vzv_interactions_scripts/vzv_rna_seq_feature_counts_mm.bash

##Let's run feature counts on the 
##Hg38 aligned STAR files
##featureCounts is part of the Subread package
##I installed it via bioconda

#http://bioinf.wehi.edu.au/featureCounts/

##input_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_alignment"
##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_featureCounts"
##gtf_path="/slipstream/home/mmariani/references/hg38_ucsc_prebuilt/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
##
##cd $input_dir
##
##/slipstream/home/mmariani/programs/miniconda3/bin/featureCounts \
##-T 20 \
##-t exon \
##-g gene_id \
##-a $gtf_path \
##-o $output_dir"/"vzv_mock_rna_re-seq_star_feature_counts_combined_11092019.txt \
##H4_R4_S4_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##H5_R7_S5_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##H6_R8_S6_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##HV1_R1_S1_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##HV2_R2_S2_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##HV3_R3_S3_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam

##input_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_vzv_star_and_counts_02212021"
##output_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_vzv_star_and_counts_02212021"
##gtf_path="/slipstream/home/mmariani/references/vzv/vzv.gtf"
##
##cd $input_dir
##
##/slipstream/home/mmariani/programs/miniconda3/bin/featureCounts \
##-T 20 \
##-t CDS \
##-g gene_id \
##-a $gtf_path \
##-o $output_dir"/"vzv_mock_rna_re-seq_star_feature_counts_combined_02242021.txt \
##H4_R4_S4_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##H5_R7_S5_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##H6_R8_S6_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##HV1_R1_S1_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##HV2_R2_S2_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##HV3_R3_S3_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam
##
########################################## VIRUS #######################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

input_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_and_vzv_star_cellranger_and_counts_02242021"
vzv_path_1=$input_dir"/HV1_R1_S1_L002_R1_001.vzv.Aligned.out.sorted.bam"
vzv_path_2=$input_dir"/HV2_R2_S2_L002_R1_001.vzv.Aligned.out.sorted.bam"
vzv_path_3=$input_dir"/HV3_R3_S3_L002_R1_001.vzv.Aligned.out.sorted.bam"
mock_path_1=$input_dir"/H4_R4_S4_L002_R1_001.vzv.Aligned.out.sorted.bam"
mock_path_2=$input_dir"/H5_R7_S5_L002_R1_001.vzv.Aligned.out.sorted.bam"
mock_path_3=$input_dir"/H6_R8_S6_L002_R1_001.vzv.Aligned.out.sorted.bam"

genes_file="/slipstream/home/mmariani/references/vzv_10x_depledge/FW__current_VZV_annotation/gtf_dumas_adjusted_cellranger.gtf"

cd /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_and_vzv_star_cellranger_and_counts_02242021

/slipstream/home/mmariani/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-s 2 \
-T 16 \
-a $genes_file \
-o $(basename $vzv_path_1 ".bam")".counts" \
$vzv_path_1

/slipstream/home/mmariani/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-s 2 \
-T 16 \
-a $genes_file \
-o $(basename $vzv_path_2 ".bam")".counts" \
$vzv_path_2

/slipstream/home/mmariani/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-s 2 \
-T 16 \
-a $genes_file \
-o $(basename $vzv_path_3 ".bam")".counts" \
$vzv_path_3

/slipstream/home/mmariani/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-s 2 \
-T 16 \
-a $genes_file \
-o $(basename $mock_path_1 ".bam")".counts" \
$mock_path_1

/slipstream/home/mmariani/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-s 2 \
-T 16 \
-a $genes_file \
-o $(basename $mock_path_2 ".bam")".counts" \
$mock_path_2

/slipstream/home/mmariani/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-s 2 \
-T 16 \
-a $genes_file \
-o $(basename $mock_path_3 ".bam")".counts" \
$mock_path_3

##htseq-count \
##-s reverse \
##$vzv_path_1 \
##$genes_file \
##> $(basename $vzv_path_1 ".bam")".counts"
##
##htseq-count \
##-s reverse \
##$vzv_path_2 \
##$genes_file \
##> $(basename $vzv_path_2 ".bam")".counts"
##
##htseq-count \
##-s reverse \
##$vzv_path_3 \
##$genes_file \
##> $(basename $vzv_path_3 ".bam")".counts"
##
##htseq-count \
##-s reverse \
##$mock_path_1 \
##$genes_file \
##> $(basename $mock_path_1 ".bam")".counts"
##
##htseq-count \
##-s reverse \
##$mock_path_2 \
##$genes_file \
##> $(basename $mock_path_2 ".bam")".counts"
##
##htseq-count \
##-s reverse \
##$mock_path_3 \
##$genes_file \
##> $(basename $mock_path_3 ".bam")".counts"

