#!/usr/bin/env bash

##For Galaxy use (SGE):
##qsub -cwd -pe threads 8 -j yes -o /slipstream/home/mmariani/projects/hhv6_time_course/logs ./star_script_mm.bash
##qsub -cwd -pe threads 8 -j yes -o /slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/logs ./star_script_mm.bash
##11/09/2019 (human) and 02/21/2021 (virus)
##qsub -cwd -pe threads 16 -j yes -o /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/logs /slipstream/home/mmariani/scripts/vzv_interactions_scripts/vzv_rna_seq_star_script_mm.bash

##For VACC use (PBS-Torque):
##PBS -q poolmemq
##PBS -l nodes=1:ppn=16,mem=128gb,vmem=128gb
##PBS -l walltime=30:00:00
##PBS -N hhv6_rna_seq
##PBS -j oe
##PBS -o /gpfs1/home/s/f/sfrietze/mike_m/hhv6_rna_project/logs

##PBS -M michael.mariani@uvm.edu
##PBS -m bea

##spack load STAR

##This script I have smoothed out in August 2019 and pairs well the the 
##index_star.bash script, it can be used for single-end or paired-end star data
##there is a featureCounts block added to both SE and PE routines, and 
##bamCoverage routine and/or samtools stats routines (along with multiqc) 
##to compile can also be added, or run separately

star_path="/slipstream/home/mmariani/programs/miniconda3/bin/STAR"
##mode="single_end_h"
mode="single_end_v"
num_threads=16
t_value="exon"
g_value="gene_id"

##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_with_vzv"
##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38"
##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_vzv"
##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_then_vzv"
##output_dir="/slipstream/home/mmariani/projects/hhv6_time_course/output_hhv6"
##output_dir="/slipstream/home/mmariani/projects/hhv6_time_course/output_hg38"
##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_vzv"
##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_alignment"

##ref_dir="/slipstream/home/mmariani/references/ucsc_hg38_with_vzv/star"
##ref_dir_h="/slipstream/home/mmariani/references/ucsc_hg38_star"
##ref_dir_v="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/star"
##ref_dir_v="/slipstream/home/mmariani/references/hhv3/star"

##gtf_file="/slipstream/home/mmariani/references/ucsc_hg38_with_vzv/ucsc_hg38_with_vzv.gtf"
##gtf_file_h="/slipstream/home/mmariani/references/ucsc_hg38_with_vzv/ucsc_hg38.gtf"
##gtf_file_v="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest_mm_cds_to_exon.gtf"
##gtf_file_v="/slipstream/home/mmariani/references/hhv3/hhv3.single.exon.changed.to.cds.gtf"
##gtf_file_v="/slipstream/home/mmariani/references/hhv3/hhv3.all.changed.to.exon.gtf"

############################### Single end alignment: human #################################################################
#######################################################################################################################

if [ $mode == "single_end_h" ]
then

output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_alignment"
ref_dir_h="/slipstream/home/mmariani/references/ucsc_hg38_star"
gtf_file_h="/slipstream/home/mmariani/references/ucsc_hg38_with_vzv/ucsc_hg38.gtf"

while IFS=$'\t' read -r -a my_array
do

cd $output_dir

base_name=$(basename "${my_array[0]}" ".fastq.gz")

$star_path \
--runThreadN $num_threads \
--genomeDir $ref_dir_h \
--readFilesIn "${my_array[0]}" \
--readFilesCommand zcat \
--outFileNamePrefix $output_dir"/"$base_name".hg38." \
--outSAMtype BAM SortedByCoordinate \
> $base_name".log.txt"

file_name=$output_dir"/"$base_name".hg38.Aligned.sortedByCoord.out.bam"
samtools index $file_name
samtools idxstats $file_name > $output_dir"/"$(basename $file_name ".bam")".idxstats"
samtools stats $file_name > $output_dir"/"$(basename $file_name ".bam")".stats"
samtools flagstat $file_name > $output_dir"/"$(basename $file_name ".bam")".flagstat"

##--outReadsUnmapped Fastx
##feature_counts_output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_featureCounts"

##~/programs/miniconda3/bin/featureCounts \
##-T $num_threads \
##-t $t_value \
##-g $g_value \
##-a $gtf_file_h \
##-o $feature_counts_output_dir"/"$base_name".hg38.counts.txt" \
##$output_dir"/"$base_name".hg38.Aligned.sortedByCoord.out.bam"

##Can do a second round of alignment here (e.g. subtractive alignment)
##$star_path \
##--runThreadN 8 \
##--genomeDir $ref_dir_v \
##--readFilesIn $output_dir"/"$base_name".hUnmapped.out.mate1" \
##--outFileNamePrefix $output_dir"/"$base_name."v" \
##--outSAMtype BAM SortedByCoordinate 
##
##~/programs/miniconda3/bin/featureCounts \
##-T 8 \
##-t exon \
##-g gene_id \
##-a $gtf_file \
##-o $output_dir"/"$base_name".v.counts.txt" \
##$output_dir"/"$base_name".v.Aligned.sortedByCoord.out.bam"

done < /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/files/vzv_rna_files.txt

fi

############################### Single end alignment: virus #################################################################
#######################################################################################################################

if [ $mode == "single_end_v" ]
then

##output_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_vzv_star_and_counts_02212021"
##output_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_vzv_all_exon_star_and_counts_02242021"
##output_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_and_vzv_star_and_counts_02242021"
output_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_and_vzv_star_cellranger_and_counts_02242021"
##ref_dir_v="/slipstream/home/mmariani/references/vzv/star"
##ref_dir_v="/slipstream/home/mmariani/references/vzv/star_all_exon"
ref_dir_v="/slipstream/home/mmariani/references/vzv/vzv_depledge/star"
##gtf_file_v="/slipstream/home/mmariani/references/vzv/vzv.gtf"
##gtf_file_v="/slipstream/home/mmariani/references/vzv/vzv_all_exon.gtf"
gtf_file_v="/slipstream/home/mmariani/references/vzv/vzv_depledge/vzv_depledge_adjusted_cellranger.gtf"

while IFS=$'\t' read -r -a my_array
do

cd $output_dir

base_name=$(basename "${my_array[0]}" ".fastq.gz")

$star_path \
--runThreadN $num_threads \
--genomeDir $ref_dir_v \
--genomeSAindexNbases 7 \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript transcript_id \
--sjdbGTFtagExonParentGene gene_id \
--readFilesIn "${my_array[0]}" \
--readFilesCommand zcat \
--outFileNamePrefix $output_dir"/"$base_name".vzv." \
> $base_name".log.txt"

##--outSAMtype BAM SortedByCoordinate \

file_name=$output_dir"/"$base_name".vzv.Aligned.sortedByCoord.out.bam"
samtools index $file_name
samtools idxstats $file_name > $output_dir"/"$(basename $file_name ".bam")".idxstats"
samtools stats $file_name > $output_dir"/"$(basename $file_name ".bam")".stats"
samtools flagstat $file_name > $output_dir"/"$(basename $file_name ".bam")".flagstat"

##--outReadsUnmapped Fastx
##feature_counts_output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_featureCounts"

##~/programs/miniconda3/bin/featureCounts \
##-T $num_threads \
##-t $t_value \
##-g $g_value \
##-a $gtf_file_h \
##-o $feature_counts_output_dir"/"$base_name".hg38.counts.txt" \
##$output_dir"/"$base_name".hg38.Aligned.sortedByCoord.out.bam"

##Can do a second round of alignment here (e.g. subtractive alignment)
##$star_path \
##--runThreadN 8 \
##--genomeDir $ref_dir_v \
##--readFilesIn $output_dir"/"$base_name".hUnmapped.out.mate1" \
##--outFileNamePrefix $output_dir"/"$base_name."v" \
##--outSAMtype BAM SortedByCoordinate 
##
##~/programs/miniconda3/bin/featureCounts \
##-T 8 \
##-t exon \
##-g gene_id \
##-a $gtf_file \
##-o $output_dir"/"$base_name".v.counts.txt" \
##$output_dir"/"$base_name".v.Aligned.sortedByCoord.out.bam"

done < /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/files/vzv_rna_files.txt

fi

######################################### Paired end alignment ########################################################
#######################################################################################################################

####Testing:
##while read -ra array
##do
##echo "${array[0]} ${array[1]}"
####echo $file1 $file2
##done < /slipstream/home/mmariani/projects/hhv6_time_course/files/hhv6_rna_files.txt

##if [ $mode == "paired_end" ]
##then
##
##while read -ra array
##do
##
##cd $output_dir
##
##base_name=$(basename "${array[0]}" ".fastq.gz")
##
##$star_path \
##--runThreadN $num_threads \
##--genomeDir $ref_dir_v \
##--readFilesIn "${array[0]}" "${array[1]}" \
##--readFilesCommand zcat \
##--outFileNamePrefix $output_dir"/"$base_name".pe.v." \
##--outSAMtype BAM SortedByCoordinate \
##> $base_name".log.txt"
##
####--outReadsUnmapped Fastx
##feature_counts_output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_featureCounts"
##~/programs/miniconda3/bin/featureCounts \
##-T $num_threads \
##-t $t_value \
##-g $g_value \
##-a $gtf_file_v \
##-o $feature_counts_output_dir"/"$base_name".pe.v.counts.txt" \
##$output_dir"/"$base_name".pe.v.Aligned.sortedByCoord.out.bam"
##
####Can do a second round of alignment here (e.g. subtractive alignment)
####$star_path \
####--runThreadN 8 \
####--genomeDir $ref_dir_v \
####--readFilesIn $output_dir"/"$base_name".hUnmapped.out.mate1" \
####--outFileNamePrefix $output_dir"/"$base_name."v" \
####--outSAMtype BAM SortedByCoordinate 
####
####~/programs/miniconda3/bin/featureCounts \
####-T 8 \
####-t exon \
####-g gene_id \
####-a $gtf_file \
####-o $output_dir"/"$base_name".v.counts.txt" \
####$output_dir"/"$base_name".v.Aligned.sortedByCoord.out.bam"
##
##done < /slipstream/home/mmariani/projects/hhv6_time_course/files/hhv6_rna_files.txt
##
##fi
