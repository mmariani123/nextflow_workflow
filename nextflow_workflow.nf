#!/usr/bin/env nextflow

// Michael P. Mariani PhD, 2025
// Working with Nextflow. Working on developing RNA-Seq pipeline
// Testing on Windows Ubuntu App, WSL2, root directory  

nextflow.enable.dsl=2

//Run command on linux shell:
//$HOME/nextflow-24.10.3-dist run -w ./nextflow_workflow $HOME/nextflow_workflow/nextflow_workflow.nf
//NB: setting "work" dir explicitly above

/*
* Pipeline parameters
*/

//Primary input
$projectDir                 = "/root/nextflow_workflow"
params.outdir               = "${projectDir}/output"
params.fastq                = "${projectDir}/fastq"
params.reference            = "${projectDir}/reference"
params.fasta                = "${projectDir}/reference/mm39.fa"
//params.index              = "${projectDir}/reference/mm39.fai"
params.bams                 = "${projectDir}/output/alignment"
params.gtf                  = "${projectDir}/reference/refGene.gtf"
params.type                 = "single"
params.frag_len             = "8"
params.threads              = "8"
params.overHang             = "75"
params.index_star           = "false"
params.run_fastqc           = "false"
params.run_star_align       = "false"
params.run_samtools         = "false"
params.run_feature_counts   = "true"
params.run_deseq2           = "false"
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
	
	if ( params.run_cluster_profiler == 'true' ){
	
	
		
	}
	
	if ( params.run_final_analyses == 'true' ){
	
	
		
	}
	
}
