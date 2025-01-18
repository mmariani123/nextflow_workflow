#!/usr/bin/env nextflow

/*
* Pipeline parameters
*/

//Primary input
$porjectDir = $HOME
params.outdir   = "$projectDir/output"
params.reads    = "$projectDir/data/*.fastq"
params.fasta    = "$projectDir/reference/mm39.fa"
params.index    = "$projectDir/reference/mm39.fai"
params.bams     = "${projectDir}/output/*.bam"
params.gtf      = "$projectDir/reference/mm39.gtf
params.type     = "single"
params.frag_len = "8"
params.threads  = "8"

//Accessory files
//params.reference        = "${projectDir}/data/ref/ref.fasta"
//params.reference_index  = "${projectDir}/data/ref/ref.fasta.fai"
//params.reference_dict   = "${projectDir}/data/ref/ref.dict"
//params.intervals        = "${projectDir}/data/ref/intervals.bed"

process download_reference{

	publishDir "${params.outdir}/reference", pattern: "*.fa", mode:'copy'

	output:
	file("*.fa") into fasta_downloaded

        when: !params.fasta

	script:
	"""
	wget --no-check-certificate https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
	gunzip mm39.fa.gz
	"""
}

process STAR_INDEX {
    publishDir "${params.reference", mode:"copy"

    input:
    file(fasta)

    output:
    tuple path(fasta), path("*"), emit: bwa_index

    script:
    """
	$star_path \
	--runMode genomeGenerate \
	--runThreadN params.threads \
	--genomeDir params.fasta \
	--genomeFastaFiles params.fasta \
	--sjdbGTFfile params.gtf
}

process STAR_ALIGNMENT {
	publishDir "${params.reference", mode:"copy"
	
	input:
	file(fastq)
	file(fasta)
	file(gtf)
	
	output:
	tuple path(fastq), path("*"), emit: star_alignment
	
	script:
	"""
	$star_path \
	--runThreadN        ${threads} \
	--genomeDir         ${fasta} \
	--readFilesIn       ${fastq} \
	--readFilesCommand  "zcat" \
	--outFileNamePrefix $output_dir"/${reads}.mm39." \
	--outSAMtype BAM    "SortedByCoordinate"
	"""
	
}

/*
 * Generate samtools files:
 */

process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    publishDir params.outdir, mode: 'symlink'

    input:
        path input_bam
		path input_bam
		path input_bam
		path input_bam

    output:
        path "${input_bam}.bai"
		path "${input_bam}.idxstats"
		path "${input_bam}.stats"
		path "${input_bam}.flagstats"

    script:
    """
    samtools index '$input_bam'
	samtools idxstats '$input_bam'
	samtools stats $file_name '$input_bam'
	samtools flagstat '$input_bam'
    """
}

workflow {

	fastq_ch      = Channel.fromFilePairs(params.fastq, checkIfExists: true)
	reference_ch  = ch_fasta = params.fasta ? Channel.value(file(params.fasta)) : fasta_downloaded
		//Channel.fromPath(params.fasta, checkIfExists: true)
	bam_ch		  = Channel.fromPath(params.bams checkIfExists: true)
    
	// Load the file paths for the accessory files (reference and intervals)
    ref_file        = file(params.reference)
    ref_index_file  = file(params.reference_index)
    ref_dict_file   = file(params.reference_dict)
    intervals_file  = file(params.intervals)

	STAR_INDEX(
	reference.ch
    "star.index.mm39.logs.txt"
    --runMode genomeGenerate \
	--runThreadN params.threads \
	--genomeDir params.fasta \
	--genomeFastaFiles params.fasta \
	--sjdbGTFfile params.gtf
	)
	
	STAR_ALIGNMENT(
	--readFilesIn = fastq_ch
	"star.index.mm39.logs.txt.log.txt"
	--runThreadN  = params.threads 
	--genomeDir   = params.outdir
	)
    
	// Create index file for input BAM file
    SAMTOOLS_INDEX(bam_ch)

}

