#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input
params.reads_bam = "${projectDir}/data/bam/reads_mother.bam"
params.outdir    = "results_genomics"

//https://github.com/nextflow-io/nextflow/issues/2082
//Note that: Nextflow does not have its own syntax parser 
//but instead uses AST manipulation techniques to extend 
//the Groovy syntax into Nextflow DS


/*
 * Generate BAM index file
 */
 
/*
pull samtools docker container:
docker pull community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
Then spin it up
docker run -it -v ./data:/data community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464
*/

/*
 * Generate BAM index file
 */
process SAMTOOLS_INDEX {

    container 'community.wave.seqera.io/library/samtools:1.20--b5dfbd93de237464'

    publishDir params.outdir, mode: 'symlink'

    input:
        path input_bam

    output:
        path "${input_bam}.bai"

    script:
    """
    samtools index '$input_bam'
    """
}

fastqPath = "/data/*.fastq.gz"

Channel
      .fromPath(fastqPath)
      .set{fastqChannel}

channelFastq.view()

process run_bash_script {
    
    input:
    val index_star_mm

    output:
    stdout

    script:
    """
    bash ${index_star_mm.bash}
    """
}

params.fasta = ''
params.index = ''
params.input = ''
params.input_type = ''
params.fragment_length = ''
params.standard_deviation = ''
params.cpus = ''
params.outdir = ''

def helpMessage() {
    log.info"""
    ====================================================
                    Practicing NextFlow RNA-Seq
   ====================================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow -bg run main.nf --input "/data/reads/*_r{1,2}.fastq.gz" --input_type 'paired_end' --fasta '/data/index/GRCh38.cdna.fa' --cpus 2

   Arguments:
      --input                        [path] Path to input RNA-Seq data. Use a suitable wildcard glob pattern to capture paired end
                                           or single end reads and surround in double quotes

      --input_type                    [str] Input data type
                                           Available: paired-end, single-end

      --fasta                        [path] Path to reference cDNA file
                                           Automatically downloaded if left empty

      --index                        [path] Path to reference cDNA Kallisto index file
                                           Automatically generated if left empty

      --outdir                       [path] Directory to write results to

   Kallisto arguments:
      --cpus                          [int] Number of CPU cores to use for alignment

      --fragment_length               [int] Estimated average fragment length (single end data)

      --standard_deviation            [int] Estimated standard deviation of fragment length (single end data)

    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Check parameters
if(!params.input) exit 1, "error: No input data provided."

if(!params.input_type) exit 1, "error: Input data type not provided, please provide 'paired-end' or 'single-end'."

if(!params.outdir) exit 1, "error: No output directory provided. Please provide a directory to write results to."

if(!params.cpus) exit 1, "error: Please specify the number of CPUs for Kallisto pseudo-alignment."

if(params.input_type == 'single-end' && (!params.fragment_length || !params.standard_deviation)) exit 1, "error: Single end data selected, but --fragment_length and/or --standard_deviation not specified. Please provide values for both."

// Place reads in channel
if(params.input_type == 'paired-end'){
  Channel
        .fromFilePairs(params.input)
	.set{ch_reads}
}else if(params.input_type == 'single-end'){
  Channel
	.fromPath(params.input)
	.map{ it -> [it.simpleName, it]}
	.set{ch_reads}
}

process download_reference{

	publishDir "${params.outdir}/reference", pattern: "*.fa", mode:'copy'

	output:
	file("*.fa") into fasta_downloaded

        when: !params.fasta

	script:
	"""
	wget --no-check-certificate http://ftp.ensembl.org/pub/release-103/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
	gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
	mv Homo_sapiens.GRCh38.cdna.all.fa GRCh38.cdna.fa
	"""
}

ch_fasta = params.fasta ? Channel.value(file(params.fasta)) : fasta_downloaded

process create_index{

	publishDir "${params.outdir}/index", pattern: "*.idx", mode:'copy'

	input:
	file(fasta) from ch_fasta

	output:
	file("*.idx") into index_created

	when: !params.index

	script:
	"""
	kallisto index -i GRCh38.fa.idx $fasta
	"""
}

ch_index = params.index ? Channel.value(file(params.index)) : index_created

process kallisto{

	publishDir "${params.outdir}/quant", pattern: "${base}", mode:'copy'

	input:
	tuple val(base), file(reads) from ch_reads
	file(index) from ch_index

	output:
	file("${base}") into kallisto_out

	script:
	if(params.input_type == 'paired-end'){
	"""
	kallisto quant \
	-i $index \
	-t ${params.cpus} \
	-o ${base}/ \
	--bias \
	$reads
	"""
	}else if(params.input_type == 'single-end'){
	"""
	kallisto quant \
	--single \
	-l ${params.fragment_length} \
	-s ${params.standard_deviation} \
	-i $index \
	-t ${params.cpus}/ \
	-o ${base}/ \
	--bias \
	$reads
	"""
	}
}
