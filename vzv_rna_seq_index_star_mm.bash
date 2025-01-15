#!/usr/bin/env bash

##For Galaxy use (SGE)
##qsub -cwd -pe threads 8 -j yes -o /slipstream/home/mmariani/projects/hhv6_time_course/logs ./index_star.bash
##qsub -cwd -pe threads 8 -j yes -o /slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/logs ./index_star.bash
##11/09/2019 (human) and 02/21/2021 (virus)
##qsub -cwd -pe threads 20 -j yes -o /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/logs /slipstream/home/mmariani/scripts/vzv_interactions_scripts/vzv_rna_seq_index_star_mm.bash

##For VACC use (PBS-Torque):
##PBS -q poolmemq
##PBS -l nodes=1:ppn=16,mem=64gb,vmem=64gb
##PBS -l walltime=03:00:00
##PBS -N star_index_job
##PBS -j oe   

########################################### NOTES ##############################################################
################################################################################################################

##This script will do indexing for STAR
##note that:
##"I believe you don't need --sjdbGTFtagExonParentTranscript Parent 
##to index GTF files 
##https://www.biostars.org/p/300066/

##To combine GTF of hg38 with vzv:
##gffread hhv3.gff3 -T -o hhv3.gtf
##cat ucsc_hg38.gtf hhv3.gtf > ucsc_hg38_with_vzv.gtf 

#https://www.biostars.org/p/93883/
#This is quiet close to what you have asked:
#https://groups.google.com/forum/#!topic/rna-star/J6qH9JCysZw
#Basically, sjdbOverhang should be set as readlength -1. 
#So if you have 75 bp read then it should be set to 74. 
#Whereas, alignSJDBoverhangMin option ignores the alignment 
#with a small spilce overhangs. I use the default settings for this parameter.

##For larger genome, e.g. hg38, or combined hg38+virus genome:

##If using gff3 for example with virus ref:
##From Seth: "I’ve used 7 for this value with VZV in the past. 
##From the manual for small genomes, ‘this needs to be scaled down, 
##with a typical value of min(14, log2(GenomeLength)/2 - 1). 
##For example, for 1 megaBase genome, this is equal to 9, 
##for 100 kiloBase genome, this is equal to 7.’" 

##The below parameters should work for indexing small gff3 files (or gtf files, e.g. after 
##converting .gff3 to .gtf with gffread tool) such as our herpesvirus genomes, 
##but need to set CDS as exon if there are no exons present,
##or could create lines with "exon" corresponding to the "CDS": in column 3,
##or if there are a handful of exon lines already present then could
##change the few "exons" in column 3 to "CDS", or could change all "CDS"
##in column 3 to "exon". Make sure to note what you decide to do when altering
##gff/gtf/gff3 files and indexing with STAR.

##For newest HHV6a genome for example, gffread produces gtf with 
##only "CDS" so I changed --sjdbGTFfeatureExon exon to
##--sjdbGTFfeatureExon CDS when indexing with STAR

##This works and will do counts by gene for a gtf file
##With no "exon" in 3rd column - instead has "CDS":
##--sjdbGTFfeatureExon CDS \
##--sjdbGTFtagExonParentTranscript transcript_id \
##--sjdbGTFtagExonParentGene gene_id

##This works and will do counts by gene for a gff3 file
##With no "exon" in 3rd column - instead has "CDS":
##--sjdbGTFfeatureExon CDS \
##--sjdbGTFtagExonParentTranscript Parent \
##--sjdbGTFtagExonParentGene gene \

##Can also consider setting sjdbOverhang (based on read size)
##--sjdbOverhang 49
##star_path="/gpfs1/home/s/f/sfrietze/programs/STAR/bin/Linux_x86_64/STAR"
##star_ref_dir="/slipstream/home/mmariani/references/"
##star_gtf_path="/gpfs1/home/s/f/sfrietze/references/hhv1/hhv1_gffread.gtf"
##star_fasta_path="/gpfs1/home/s/f/sfrietze/references/hhv1/hhv1.fasta"

########################################### CODE ##################################################################################
###################################################################################################################################

star_path="/slipstream/home/mmariani/programs/miniconda3/bin/STAR"

######################################### Index Small virus genomes (e.g. VZV) ####################################################
###################################################################################################################################

##HHV6 Newest
##gffread hhv6_newest.gff3 -T -o hhv6_newest.gtf
##star_path="/slipstream/home/mmariani/programs/miniconda3/bin/STAR"
##star_ref_dir="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/star"
##star_fasta_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest_mm_cds_to_exon.gtf"
##star_gtf_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest.gtf"
##star_gff3_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest.gff3"

##HHV6_GFP
##star_path="/slipstream/home/mmariani/programs/miniconda3/bin/STAR"
##star_ref_dir="/slipstream/home/mmariani/references/hhv6/hhv6_gfp/star"
##star_fasta_path="/slipstream/home/mmariani/references/hhv6/hhv6_gfp/hhv6_gfp.fasta"
##star_gff3_path="/slipstream/home/mmariani/references/hhv6/hhv6_gfp/hhv6_gfp.gff3"

##VZV
##gffread vzv.gff3 -T -o vzv.gtf
##September? 2019 attempt to index small viral genomes with STAR:
##star_ref_dir="/slipstream/home/mmariani/references/hhv3/star/"
##star_fasta_path="/slipstream/home/mmariani/references/hhv3/hhv3.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.gtf"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.exon.removed.gtf"
##star_gff3_path="/slipstream/home/mmariani/references/hhv3/hhv3.gff3"
##star_gff3_path="/slipstream/home/mmariani/references/hhv3/hhv3.exon.removed.gff3"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.single.exon.changed.to.cds.gtf"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.all.changed.to.exon.gtf"
##Here (above) I changed the single "exon" (3rd column) line in the gtf file created 
##from the genbank gff3 to "CDS", thus all 3rd column entries are "CDS" in the gtf 
##file then I indexed with STAR below

##02/21/2021 need to reindec for newer version of star
##dont know which one I used above, so I will go with 
##vzv.gtf which I think was created from 
##the gffread program, first.

##star_ref_dir="/slipstream/home/mmariani/references/vzv/star/"
##star_fasta_path="/slipstream/home/mmariani/references/vzv/vzv.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/vzv/vzv.gtf"
##
##$star_path \
##--runMode genomeGenerate \
##--runThreadN 20 \
##--genomeDir $star_ref_dir \
##--genomeFastaFiles $star_fasta_path \
##--genomeSAindexNbases 7 \
##--sjdbGTFfile $star_gtf_path \
##--sjdbGTFfeatureExon exon \
##--sjdbGTFtagExonParentTranscript transcript_id \
##--sjdbGTFtagExonParentGene gene_id

##--readFilesCommand zcat \

######################################### Index Hg38 and VZV all exons gtf ##############################################################################
###################################################################################################################################

##star_ref_dir="/slipstream/home/mmariani/references/ucsc_hg38_canonical_with_vzv/star"
##star_fasta_path="/slipstream/home/mmariani/references/ucsc_hg38_canonical_with_vzv/ucsc_hg38_canonical_with_vzv.fa"
##star_gtf_path="/slipstream/home/mmariani/references/ucsc_hg38_canonical_with_vzv/ucsc_hg38_canonical_with_vzv_all_exons.gtf"
##
##$star_path \
##--runMode genomeGenerate \
##--runThreadN 20 \
##--genomeDir $star_ref_dir \
##--genomeFastaFiles $star_fasta_path \
##--genomeSAindexNbases 14 \
##--sjdbGTFfile $star_gtf_path \
##--sjdbGTFfeatureExon exon \
##--sjdbGTFtagExonParentTranscript transcript_id \
##--sjdbGTFtagExonParentGene gene_id

##--readFilesCommand zcat \

######################################### Index Small virus genomes (e.g. VZV) ####################################################
###################################################################################################################################

##HHV6 Newest
##gffread hhv6_newest.gff3 -T -o hhv6_newest.gtf
##star_path="/slipstream/home/mmariani/programs/miniconda3/bin/STAR"
##star_ref_dir="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/star"
##star_fasta_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest_mm_cds_to_exon.gtf"
##star_gtf_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest.gtf"
##star_gff3_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest.gff3"

##HHV6_GFP
##star_path="/slipstream/home/mmariani/programs/miniconda3/bin/STAR"
##star_ref_dir="/slipstream/home/mmariani/references/hhv6/hhv6_gfp/star"
##star_fasta_path="/slipstream/home/mmariani/references/hhv6/hhv6_gfp/hhv6_gfp.fasta"
##star_gff3_path="/slipstream/home/mmariani/references/hhv6/hhv6_gfp/hhv6_gfp.gff3"

##VZV
##gffread vzv.gff3 -T -o vzv.gtf
##September? 2019 attempt to index small viral genomes with STAR:
##star_ref_dir="/slipstream/home/mmariani/references/hhv3/star/"
##star_fasta_path="/slipstream/home/mmariani/references/hhv3/hhv3.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.gtf"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.exon.removed.gtf"
##star_gff3_path="/slipstream/home/mmariani/references/hhv3/hhv3.gff3"
##star_gff3_path="/slipstream/home/mmariani/references/hhv3/hhv3.exon.removed.gff3"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.single.exon.changed.to.cds.gtf"
##star_gtf_path="/slipstream/home/mmariani/references/vzv/vzv.all.changed.to.exon.gtf"
##Here (above) I changed the single "exon" (3rd column) line in the gtf file created 
##from the genbank gff3 to "CDS", thus all 3rd column entries are "CDS" in the gtf 
##file then I indexed with STAR below

##02/24/2021 Now try with all "exon" or "CDS" as the third column

##star_ref_dir="/slipstream/home/mmariani/references/vzv/star_all_exon/"
##star_fasta_path="/slipstream/home/mmariani/references/vzv/vzv.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/vzv/vzv_all_exon.gtf"
##
##$star_path \
##--runMode genomeGenerate \
##--runThreadN 20 \
##--genomeDir $star_ref_dir \
##--genomeFastaFiles $star_fasta_path \
##--genomeSAindexNbases 7 \
##--sjdbGTFfile $star_gtf_path \
##--sjdbGTFfeatureExon exon \
##--sjdbGTFtagExonParentTranscript transcript_id \
##--sjdbGTFtagExonParentGene gene_id

##--readFilesCommand zcat \

##02/24/2021 Now try with all "CDS" as the third column

##star_ref_dir="/slipstream/home/mmariani/references/vzv/star_all_cds/"
##star_fasta_path="/slipstream/home/mmariani/references/vzv/vzv.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/vzv/vzv.single.exon.changed.to.cds.gtf"
##
##$star_path \
##--runMode genomeGenerate \
##--runThreadN 20 \
##--genomeDir $star_ref_dir \
##--genomeFastaFiles $star_fasta_path \
##--genomeSAindexNbases 7 \
##--sjdbGTFfile $star_gtf_path \
##--sjdbGTFfeatureExon CDS \
##--sjdbGTFtagExonParentTranscript transcript_id \
##--sjdbGTFtagExonParentGene gene_id
##
####--readFilesCommand zcat \

##02/24/2021
##Now try with the new depledge gtf that we were working 
##with for the scRNA data

star_ref_dir="/slipstream/home/mmariani/references/vzv/vzv_depledge/star"
star_fasta_path="/slipstream/home/mmariani/references/vzv/vzv_depledge/vzv_depledge_cellranger.fa"
star_gtf_path="/slipstream/home/mmariani/references/vzv/vzv_depledge/vzv_depledge_adjusted_cellranger.gtf"

$star_path \
--runMode genomeGenerate \
--runThreadN 20 \
--genomeDir $star_ref_dir \
--genomeFastaFiles $star_fasta_path \
--genomeSAindexNbases 7 \
--sjdbGTFfile $star_gtf_path \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript transcript_id \
--sjdbGTFtagExonParentGene gene_id

##--readFilesCommand zcat \

######################################### Index Hg38 ##############################################################################
###################################################################################################################################

##11/09/2019, Hg38 indexing:
##Use the prebuilt ucsc reference files for hg38 to build 
##a new STAR genome, for newest version of STAR:

##star_fasta_path="/slipstream/home/mmariani/references/hg38_ucsc_prebuilt/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa"
##star_ref_dir="/slipstream/home/mmariani/references/ucsc_hg38_star"
##star_gtf_path="/slipstream/home/mmariani/references/hg38_ucsc_prebuilt/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
##
##cd "/slipstream/home/mmariani/references/ucsc_hg38_star"
##
##$star_path \
##--runMode genomeGenerate \
##--runThreadN 16 \
##--genomeDir $star_ref_dir \
##--genomeFastaFiles $star_fasta_path \
##--sjdbGTFfile $star_gtf_path \
##> "star.index.hg38.11092019.logs.txt"
