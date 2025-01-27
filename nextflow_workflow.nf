#!/usr/bin/env nextflow

// Michael P. Mariani PhD, 2025
// Working with Nextflow. Working on developing RNA-Seq pipeline
// Testing on Windows 10 Ubuntu App, WSL2, root directory  

nextflow.enable.dsl=2

//Run commands on linux shell:
//conda activate nextflow_testing #activate our dedicated conda environment
//$HOME/nextflow-24.10.3-dist run -w ./nextflow_workflow/work $HOME/nextflow_workflow/nextflow_workflow.nf -with-report -with-timeline -with-dag flowchart.png
//NB: setting "work" dir explicitly above

/*
* Pipeline parameters
*/

//Primary input
projectDir                 = "/root/nextflow_workflow"
params.outdir               = "${projectDir}/output"
params.fastq                = "${projectDir}/fastq"
params.reference            = "${projectDir}/reference"
params.fasta                = "${projectDir}/reference/mm39.fa"
//params.index              = "${projectDir}/reference/mm39.fai"
params.bams                 = "${projectDir}/output/alignment"
params.counts               = "${projectDir}/output/counts"
params.deseq2               = "${projectDir}/output/deseq2"
params.deseq2_samples_file  = "${projectDir}/output/deseq2/samples_for_deseq2.txt"
params.gtf                  = "${projectDir}/reference/refGene.gtf"
params.type                 = "single"
params.frag_len             = "8"
params.threads              = "8"
params.overHang             = "75"
params.index_star           = "false"
params.run_fastqc           = "false"
params.run_star_align       = "false"
params.run_samtools         = "false"
params.run_feature_counts   = "false"
params.run_deseq2           = "true"
params.run_cluster_profiler = "false"
params.run_final_analyses   = "false"

//Accessory files
//params.reference        = "${projectDir}/data/ref/ref.fasta"
//params.reference_index  = "${projectDir}/data/ref/ref.fasta.fai"
//params.reference_dict   = "${projectDir}/data/ref/ref.dict"
//params.intervals        = "${projectDir}/data/ref/intervals.bed"

process FASTQC {
  container 'quay.io/biocontainers/fastqc:0.11.9--0'

  tag "Running FastQC on ${sampleid}"

  publishDir "${params.outdir}/fastqc"
  //, mode: 'copy'

  input: 
    tuple val(sampleid), path(fastq)

  output:
    path("*.html")

  script:
    """
    fastqc \
	-t ${params.threads} \
	${fastq}
    """
	
}

process STAR_INDEX {
	container 'quay.io/biocontainers/star:2.7.10b--h6b7c446_1'
	
	publishDir params.reference, mode: "copy"
    
	input:
    path fasta
	path gtf
	
    //output:
	//path params.reference
	//The above output caused scope issues and script failure
	
    script:
	"""
	STAR 
	--runMode genomeGenerate \
    --genomeDir ${params.reference} \
    --genomeFastaFiles ${fasta} \
    --runThreadN ${params.threads} \
	--sjdbGTFfile ${gtf} \
	--sjdbOverhang ${params.overHang}
    """
}

process STAR_ALIGNMENT {
	container 'quay.io/biocontainers/star:2.7.10b--h6b7c446_1'
	publishDir "${params.outdir}/alignment", mode:"copy"
	
	input: 
    tuple val(sampleid), path(fastq)
	path fasta
	
	output:
	//tuple path(fastq), path("*"), emit: star_alignment
	path "${fastq}.bam"     , emit: bam
	
	script:
	"""
	STAR \
	--runThreadN        ${params.threads} \
	--genomeDir         ${fasta} \
	--readFilesIn       ${fastq} \
	--readFilesCommand  "zcat" \
	--outFileNamePrefix ${params.outdir}"/alignment/"${sampleid}".mm39." \
	--outSAMtype        BAM SortedByCoordinate
	#>"star.index.mm39.log.txt"
	"""
}

process SAMTOOLS_WORK {

 	//Generate samtools files:
    
	//docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
	container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    publishDir "${params.outdir}/alignment", mode: 'copy'

	input:
    tuple val(sample), path(bam)

    output:
    path "${sample}.sorted.bam"
    path "${sample}.sorted.bam.bai"
    path "${sample}.idxstats"
	path "${sample}.stats"
	path "${sample}.flagstat"
	
    shell:
    """
	dirName="\$(dirname ${bam})"
	samtools sort     ${bam} -o "!{sample}.sorted.bam"
    samtools index    "\${dirName}/${sample}.sorted.bam"
	samtools idxstats "\${dirName}/${sample}.sorted.bam" > "${sample}.idxstats"
	samtools stats    "\${dirName}/${sample}.sorted.bam" > "${sample}.stats"
	#The output can be visualized using plot-bamstats.
	samtools flagstat "\${dirName}/${sample}.sorted.bam" > "${sample}.flagstat"
	"""
}   

process FEATURE_COUNTS {

	//Can use featureCouts container or install with CONDA as usual 
	publishDir "${params.outdir}/counts", mode: 'copy'

	input:
    path gtf
	tuple val(sample), path(bam)
	
	output:
    path "${sample}.counts"
	
	script:
	"""
	featureCounts \
	-s "2" \
	-T ${params.threads} \
	-a ${gtf} \
	-o ${sample}".counts" \
	${bam}
	#NB: remember with -s parameter want to ensure reverse strandedness of original reads
	"""

}

process DESEQ2{

    publishDir "${params.deseq2}", mode: 'copy'

	//Make sure R and libraries are installed in conda env:
	//install -c conda-forge r-essentials
	//then run R and install bioconductor
	//if (!requireNamespace("BiocManager", quietly = TRUE))
    //install.packages("BiocManager")
	//or ...
	//conda install -c conda-forge 'r-magrittr'
	//etc.
	
	//BiocManager::install('magrittr')
	//BiocManager::install('dplyr')
	//BiocManager::install('tidyr')
	//BiocManager::install('DESeq2')
	//BiocManager::isntall('fgsea')
    //BiocManager::install('DOSE')
	//BiocManager::install('enrichplot')
	//BiocManager::install('clusterProfiler')
	//BiocManager::install('org.Hs.eg.db')
	//BiocManager::install('ggplot2')
	//BiocManager::install('cowplot')
	//BiocManager::install('patchwork')
	//BiocManager::install('ggpubr')
	//BiocManager::install("apeglm") for use with DESEQ2
	//BiocManager::install("pheatmap") for use with DESEQ2
	
	conda 'r-magrittr'
	conda 'r-dplyr'
	conda 'r-tidyr'
	conda 'r-DESeq2'
	conda 'r-apeglm'
	conda 'r-pheatmap'
	//conda 'r-clusterProfiler'
	conda 'r-org.Hs.eg.db'
	conda 'r-ggplot2'
	conda 'r-cowplot'
	conda 'r-patchwork'
	//conda 'r-ggpubr'

    input:
    path countsInputDir
	//path countsOutputDir
	path inputSamplesFile
	
	output:
	//path "combined_counts_matrix.csv" into csv_ch
	//path "samples_for_deseq2.txt" into txt_ch
	//path "pathways.png" into png_ch
	path "combined_counts_matrix.csv" //, emit: csv_file
	path "deseq2_df.txt"
	path "deseq2_heatmap.png"
	path "deseq2_pca.png"
	//path "pathways.png" //, emit: png_file
	//file("*.csv") into csv_ch
	//file("*.png") into png_ch
    //path "*.csv", emit: csv_files
	//path "*.png", emit: png_files
	
	script:
    """
	deseq2.R "${countsInputDir}" "${inputSamplesFile}" "combined_counts_matrix.csv" "deseq2_df.txt" "deseq2_heatmap.png" "deseq2_pca.png"
	"""
}

process CLUSTER_PROFILER{

    publishDir "${params.deseq2}", mode: 'copy'
	
	conda 'r-magrittr'
	conda 'r-dplyr'
	conda 'r-tidyr'
	//conda 'r-clusterProfiler'
	conda 'r-org.Hs.eg.db'
	conda 'r-ggplot2'
	conda 'r-cowplot'
	conda 'r-patchwork'
	//conda 'r-ggpubr'

    input:
    path countsInputDir
	//path countsOutputDir
	path inputSamplesFile
	
	output:
	//path "combined_counts_matrix.csv" into csv_ch
	//path "samples_for_deseq2.txt" into txt_ch
	//path "pathways.png" into png_ch
	path "combined_counts_matrix.csv" //, emit: csv_file
	//path "pathways.png" //, emit: png_file
	//file("*.csv") into csv_ch
	//file("*.png") into png_ch
    //path "*.csv", emit: csv_files
	//path "*.png", emit: png_files
	
	script:
    """
	clusterProfiler.R "${countsInputDir}" "${inputSamplesFile}" combined_counts_matrix.csv
	"""
}

Channel
  .fromPath("${params.fastq}/*{.fastq.gz,.fq.gz,.fastq,.fq}")
  .map { it -> tuple( it.simpleName, it ) }
  .ifEmpty { error "Cannot find any fastq files in ${params.fastq}" }
  .set { fastq_files }
  
Channel
  .fromPath("${params.bams}/*{.bam,.sam}")
  .map { it -> tuple( it.simpleName, it ) }
  .ifEmpty { error "Cannot find any bam files in ${params.bam}" }
  .set { bam_files }
  
Channel
  .fromPath("${params.bams}/*{.sorted.bam,.sorted.sam}")
  .map { it -> tuple( it.simpleName, it ) }
  .ifEmpty { error "Cannot find any sorted bam files in ${params.bam}" }
  .set { sorted_bam_files }
  
star_fasta     = file(params.fasta)
star_gtf       = file(params.gtf)
star_reference = file(params.reference)
//deseq2_input   = file(params.counts)
//deseq2_output  = file(params.deseq2)
//deseq2_samples = file(params.deseq2_samples_file)
def deseq2_input = Channel.fromPath(params.counts)
//def deseq2_output = Channel.fromPath(params.deseq2)
def deseq2_samples = Channel.fromPath(params.deseq2_samples_file)
//deseq2_file1   = file("${params.deseq2}/combined_counts_matrix.csv")
//deseq2_file2   = file("${params.deseq2}/samples_for_deseq2.txt")
//deseq2_file3   = file("${params.deseq2}/pathways.png")
def cluster_rofiler_input = file("${params.clusterProfiler}")

workflow {
	
	if ( params.index_star == 'true' ){
	//	if !(file(params.index_star, checkIfExists: true)) { // is ! breaking it here?
			
			STAR_INDEX(star_fasta, star_gtf)
		
	//	}
	//	else {
			//echo "STAR genome already indexed!"
	//	}
	}
	
	if ( params.run_fastqc == 'true' ){
	
		FASTQC(fastq_files)
	
	}
	
	if ( params.run_star_align == 'true'){
	
		STAR_ALIGNMENT(fastq_files, star_reference)
		
	}
	
	if ( params.run_samtools == 'true' ){
	
		SAMTOOLS_WORK(bam_files)
		
	}
	
	if ( params.run_feature_counts == 'true' ){
	
		FEATURE_COUNTS(star_gtf, sorted_bam_files)
		
	}
	
	if ( params.run_deseq2 == 'true' ){
	
		//DESEQ2(deseq2_input, deseq2_output, deseq2_samples)
		DESEQ2(deseq2_input, deseq2_samples)
		
	}
	
	if ( params.run_cluster_profiler == 'true' ){
	
		CLUSTER_PROFILER(cluster_rofiler_input)
		
	}
	
	if ( params.run_final_analyses == 'true' ){
	
	
		
	}
	
}
