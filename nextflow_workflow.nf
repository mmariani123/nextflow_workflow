#!/usr/bin/env nextflow

// Michael P. Mariani PhD, 2025
// Working with Nextflow. Working on developing RNA-Seq pipeline
// Testing on Windows Ubuntu App, WSL2, root directory  

//https://github.com/nextflow-io/nextflow/issues/2082
//Note that: Nextflow does not have its own syntax parser 
//but instead uses AST manipulation techniques to extend 
//the Groovy syntax into Nextflow DS

nextflow.enable.dsl=2

/*
* Pipeline parameters
*/

//Primary input
$projectDir      = "/root/nextflow_workflow"
params.outdir    = "${projectDir}/output"
params.fastq     = "${projectDir}/fastq"
params.reference = "${projectDir}/reference"
params.fasta     = "${projectDir}/reference/mm39.fa"
//params.index     = "${projectDir}/reference/mm39.fai"
params.bams      = "${projectDir}/output/*.bam"
params.gtf       = "${projectDir}/reference/refGene.gtf"
params.type      = "single"
params.frag_len  = "8"
params.threads   = "8"
params.overHang  = "75"
params.index_star = "false"
params.run_fastqc = "false"
params.run_star_align = "true"
params.run_samtools = "false"

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
	--outFileNamePrefix ${params.outdir}"/"${sampleid}".mm39." \
	--outSAMtype        BAM SortedByCoordinate
	//"star.index.mm39.logs.txt.log.txt"
	"""
}

process SAMTOOLS_WORK {

 	//Generate samtools files:
    
	container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    publishDir params.outdir, mode: 'symlink'

    input:
        path input_bam

    output:
        path "${input_bam}.bai"
		path "${input_bam}.idxstats"
		path "${input_bam}.stats"
		path "${input_bam}.flagstats"

    script:
    """
    samtools index ${input_bam}
	samtools idxstats ${input_bam}
	samtools stats ${input_bam}
	samtools flagstat ${input_bam}
    """
}

Channel
  .fromPath("${params.fastq}/*{.fastq.gz,.fq.gz,.fastq,.fq}")
  .map { it -> tuple( it.simpleName, it ) }
  .ifEmpty { error "Cannot find any fastq files in ${params.fastq}" }
  .set { fastq_files }
  
star_fasta = file(params.fasta)
star_gtf   = file(params.gtf)
star_reference = file(params.reference)
bams_ch	   = Channel.fromPath(params.bams) 

workflow {

	// fastq_ch           = Channel.fromFilePath(params.fastq, checkIfExists: true)
	// fasta_ch      	  = Channel.fromPath(params.fasta, checkIfExists: true)
	// gtf_ch             = Channel.fromPath(params.gtf, checkIfExists: true)
	// Load the file paths for the accessory files (reference and intervals)
    // I manually set some channels above
	// ref_file           = file(params.fasta)
    // ref_index_file     = file(params.index)
	// fasta_ch = Channel.fromPath(params.fasta)
	// fasta_ch.view()
	// fasta_ch = Channel.of(params.fasta)
	// reference_ch = Channel.fromPath(params.reference)
	// reference_ch.view()
	// gtf_ch   = Channel.fromPath(params.gtf)
	// gtf_ch.view()
	// reference_ch = Channel.fromPath(params.reference)
	
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
		SAMTOOLS_WORK(bams_ch)
	}
	
}
