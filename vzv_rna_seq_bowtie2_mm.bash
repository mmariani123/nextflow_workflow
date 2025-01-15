#!/usr/bin/env bash

##Mike Mariani UVM 2019-2021

##qsub -cwd -pe threads 16 -j yes -o /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/logs /slipstream/home/mmariani/scripts/vzv_interactions_scripts/vzv_rna_seq_bowtie2_mm.bash

bowtie2_path="/slipstream/home/mmariani/programs/miniconda3/bin/bowtie2"
bowtie2_build_path="/slipstream/home/mmariani/programs/miniconda3/bin/bowtie2-build"
ref_path="/slipstream/home/mmariani/references/vzv/vzv.fasta"
bowtie2_ref="/slipstream/home/mmariani/references/vzv/vzv"
num_threads=16
output_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_vzv_bowtie2_02222021"

$bowtie2_build_path \
--threads $num_threads \
"/slipstream/home/mmariani/references/vzv/vzv.fasta" \
"/slipstream/home/mmariani/references/vzv/vzv"

while IFS=$'\t' read -r -a my_array
do

cd $output_dir

base_name=$(basename "${my_array[0]}" ".fastq.gz")

$bowtie2_path \
-p 8 \
-x $bowtie2_ref \
-U "${my_array[0]}" \
1> $base_name".sam" \
2> $base_name".bowtie2.log"
##-1 "${my_array[0]}" \
##-2 "${my_array[1]}" \

samtools sort -@ 15 -o $base_name".sorted.bam" $base_name".sam"
samtools index $base_name".sorted.bam"
samtools stats $base_name".sorted.bam" > $base_name".sorted.stats"
samtools idxstats $base_name".sorted.bam" > $base_name".sorted.idxstats"

done < /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/files/vzv_rna_files.txt

##Create bigwigs:

cd $output_dir

for i in *.sorted.bam
do

base_name=$(basename $i ".bam")

/slipstream/home/mmariani/programs/miniconda3/bin/bamCoverage \
-p 16 \
-b $i \
-o $base_name".bw" \
--effectiveGenomeSize 124884 \

done
