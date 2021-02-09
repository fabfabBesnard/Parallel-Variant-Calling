#!/usr/bin/env nextflow

/*
# Pipeline 
#What it does: 
#Script used for Variant analysis of Mutation accumulation lines
#For a given sample, this script pipes all steps from raw reads to bam files and to a gVCF. gVCF allows subsequent joint genotyping with other samples for a better variant analysis.
# 1)map with bwa 
#		mem faster and more accurate algorithm (ideal for small genome and "long" 100-pb Illumina reads
#		-a Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments. 
#		-M=shorter split reads flagged as secondary (mem allows split alignment), necessary for PICARD compatibility
# 2) (Optional) Filter using FLAG for uniquely and properly mapped as a pair ('ufilter' stands for Unmapped filtered)
#		Note: it should not be necessery to filter like this: GATK-tools' filters will do it. This deprecated step can still be activated using the config file
#		it is not necessary to use Picard/FixMateInformation or samtools/fixmate beacause GATK/indelrealigner does it
# 3) do the preliminary procedures of GATK best practices before accurate calling 
#		->File formatting/compatibility: Add Read Groups, mark duplicated Reads, Index
#		->Sequencing data & Alignment improvement: 
#				* (optional) Realignment around Indel
#				* (Optional) BQSR. Note: with Haplotype Caller, BQSR only marginally improves false calls while this step takes a very long time. 
# 4) do a variant call with GATK Haplotype caller in a gVCF mode.
#		-> This allow to genotype together a cohort of samples and to easily add further samples into the analysis using GENOTYPEgVCF


# required: 
# bwa version: 0.7.5a-r405 or later -> PSMN version (June 2016)=Bwakit/0.7.15
# samtools (PSMN version June 2016= 0.1.18)
# Picard Version: 1.110+. PSMN version (June 2016)= 2.3.0
# GATK 3.7+. PSMN version (June 2016): 3.6 (java1.8 compatible)
#OPTION: check $filter and $BQSR option in the config file.
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
      nextflow run script/GATK_to_gVCF.nf -c script/GATK_to_gVCF.config --read /home/rmarin/V300042688_L2_AE06084935-608* --index /home/fbesnard/Reference_genomes/Moss/Physcomitrella_patens.Phypa_V3.dna_rm.toplevel* -profile psmn


    Required arguments:
      --read                Full path to directory and name of reads in fastq.gz
      --index               Full path to directory of genome

    Nextflow config:
      -c                            Path to config file: src/chip_analysis.config
      -profile                      Profil used by nextflow to run the pipeline (you have choice between singularity, docker, psmn or ccin2p3)
                                    For local utilisation use singularity or docker
    Save option:
      --outdir                      Specify where to save the output from the nextflow run (default: "./results/")

    help message:
      --help                        Print help message
    """
      .stripIndent()
  }


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////
/* --         DEFAULT PARAMETER VALUES         -- */
////////////////////////////////////////////////////

params.help = false
params.index= false
params.read = false
params.outdir = 'results'


/*
 * SET UP CONFIGURATION VARIABLES
 */
//params.help="False"
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       HEADER LOG INFO                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Header log info
log.info nfcoreHeader()
def summary = [:]
summary['read']                 = params.read ?: 'Not supplied'
summary['index']                 = params.index ?: 'Not supplied'
summary['Config Profile']         = workflow.profile
summary['Output']                 = params.outdir
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          STEP 1 : MAPPING                           -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
#Step1. MAPPING
printf "${col} **STEP1** MAPPING pe-reads with bwa (mem, -aM, default) + sort bam output ${NC} \n"
#Note: bwa allowed for 4 threads. samtools sort syntax (with -T) fixed to new samtools version.
for ((a=1; a <=$NbRUNS; a++))
do
	echo "           ------1.$a Map Sequencing Run n°$a ------"
	bwa mem -t 4 -aM $REFERENCE ${READS[$((2*$a-1))]} ${READS[$((2*$a))]} | samtools view -buS -|samtools sort - -T $OUTPUT_NAME -o ${OUTPUT_NAME}.SR$a.bam
done

*/

if (params.read) {
        Channel
            .fromFilePairs(params.read, size:-1)
            .ifEmpty { error "Cannot find any file matching: ${params.read}" }
            .set{ fastqgz }
}

if (params.index) {
        Channel
            .fromPath( params.index )
            .ifEmpty { error "Cannot find index: ${params.fasta}" }
            .set { index }
}

index.view()
fastqgz.view()


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       NF-CORE HEADER                                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

def nfcoreHeader() {
    // Log colors ANSI codes
    c_reset =  "\033[0m";
    c_dim = "\033[2m";
    c_black = "\033[0;30m";
    c_green = "\033[0;32m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_purple = "\033[0;35m";
    c_cyan = "\033[0;36m";
    c_white ="\033[0;37m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  CHIP-SEQ Pipeline          ${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}
