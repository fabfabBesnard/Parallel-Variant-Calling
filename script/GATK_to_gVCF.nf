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
      nextflow run script/GATK_to_gVCF.nf -c script/GATK_to_gVCF.config --reads "/home/rmarin/V300042688_L2_AE06084935-608*" --genomeindex "/home/fbesnard/Reference_genomes/Moss/Physcomitrella_patens.Phypa_V3.dna_rm.toplevel" -profile psmn

    Don't use dot in file name, use underscore instead
    ( bad exemple : truc.v2.machin3.fasta , good exemple : truc_v2_machin3.fasta )
    
    Required arguments:
      --reads                    Full path to directory and name of reads in fastq.gz
      --genomefasta              Full path to file of reference genome 
      --genomeindex              Full path to directory of index

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
params.genomeindex= false
params.genomefasta = false
params.reads = false
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
summary['reads']                 = params.reads ?: 'Not supplied'
summary['genomeindex']           = params.genomeindex  ?: 'Not supplied'
summary['genomefasta']           = params.genomefasta  ?: 'Not supplied'
summary['Config Profile']        = workflow.profile
summary['Output']                = params.outdir
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          STEP   : INDEX                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


if (params.reads) {
        Channel
            .fromFilePairs(params.reads, size:2)
            .ifEmpty { error "Cannot find any file matching: ${params.reads}" }
            .set{ fastqgz }
}

if (params.genomefasta) {
        Channel
            .fromPath( params.genomefasta )
            .ifEmpty { error "Cannot find any file matching: ${params.genomefasta}" }
            //.map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
            .into { fasta_file ; testfa }

        testfa.view()

        process index_fasta {
                            label "bwa"
                            tag "$fasta.simpleName"
                            publishDir "results/mapping/index/", mode: 'copy'

                            input:
                              file fasta  from fasta_file

                            output:
                              file "${fasta.baseName}.*" into index_files , testindexout
                              file "*_bwa_report.txt" into index_files_report

                            script:
                              """
                              bwa index -p ${fasta.baseName} ${fasta} \
                              &> ${fasta.baseName}_bwa_report.txt
                              """

                            }
      testindexout.view()
}

if (params.genomeindex) {
        Channel
            .fromPath( params.genomeindex )
            .ifEmpty { error "Cannot find any file matching: ${params.genomeindex}" }
            //.map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
            .toList()
            .into { index_files ; testindex }
        testindex.view()

}

/* AJOUTER POSSIBILITER D'INDEX DEJA FAIT 

./nextflow run script/GATK_to_gVCF.nf -c script/GATK_to_gVCF.config --reads "/home/rmarin/V300042688_L2_AE06084935-608_{1,2}.fq.gz" --genomeindex /home/fbesnard/Reference_genomes/Moss/Physcomitrella_patens.Phypa_V3.dna_rm.toplevel.fa.* -profile psmn

test.view()
[Physcomitrella_patens.Phypa_V3.dna_rm.toplevel.fa, [/home/fbesnard/Reference_genomes/Moss/Physcomitrella_patens.Phypa_V3.dna_rm.toplevel.fa.amb]]

Probleme tous les fichier ne sont pas recuperer seulement le premier
*/


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          STEP   : MAPPING                           -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



process mapping_fastq {
  label 'bwa'
  tag "$pair_id"
  publishDir "results/mapping/${pair_id}/sam", mode: 'copy'

  input:
  set pair_id, file(reads) from fastqgz
  file index from index_files.collect()

  output:
  set pair_id, "${pair_id}.sam" into sam_files
  set pair_id, "${pair_id}_bwa_report.txt" into mapping_repport_files

  script:
  index_id = index[0].baseName
  """
  bwa mem -t ${task.cpus} \
  -aM ${index_id} ${reads[0]} ${reads[1]} \
  -o ${pair_id}.sam &> ${pair_id}_bwa_report.txt
  """

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       STEP   : SAM to BAM                           -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



process sam_to_bam {
  label 'samtools'
  tag "$pair_id"
  publishDir "results/mapping/${pair_id}/bam", mode: 'copy'

  input:
  set pair_id, "${pair_id}.sam" from sam_files

  output:
  set pair_id, "${pair_id}.bam" into bam_files
  

  script:
"""
samtools view -buS ${pair_id}.sam | samtools sort - -o ${pair_id}.bam
"""


}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       STEP   : READ GROUPS                          -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////



/*
Adding Read-Group Informations

for ((a=1; a <=$NbRUNS; a++))
do
	echo "           ------2.$a Read Groups for Sequencing Run n°$a ------"
	java -Xmx${RAM}g -jar $PICARD AddOrReplaceReadGroups \
	I=$OUTPUT_NAME.SR$a.bam \
	RGID=${RGID[$a]} RGSM=$RGSM RGLB=$RGLB RGPU=${RGPU[$a]} RGPL=$RGPL \
	O=$OUTPUT_NAME.RG.SR$a.bam
	printf "${col} Check & Clean files before next step ${NC} \n"
	if [ -e $OUTPUT_NAME.RG.SR$a.bam ]; then
	rm $OUTPUT_NAME.SR$a.bam
	else     
	printf "\e[0;31m ##ERROR## Expected bam file was not generated. Script will abort ${NC} \n"    
	exit 1
	fi
done
*/

process picard_tools {
  label 'picardtools'
  tag "$pair_id"
  publishDir "results/picard/${pair_id}/bam", mode: 'copy'

  input:
  set pair_id, bam_file from bam_files

  output:
  set pair_id, "${pair_id}_readGroup_MarkDuplicates.bam" into bam_files_RG
  set pair_id, "${pair_id}_marked_dup_metrics.txt" into picardmetric_files

  script:
  //java -Xmx${RAM}g -jar picard AddOrReplaceReadGroups \
  """
  PicardCommandLine AddOrReplaceReadGroups \
       I=${bam_file} \
       O="${pair_id}_readGroup.bam" \
       RGID=4 \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=20

  PicardCommandLine MarkDuplicates \
      I="${pair_id}_readGroup.bam" \
      O="${pair_id}_readGroup_MarkDuplicates.bam" \
      M="${pair_id}_marked_dup_metrics.txt"
  """

}


//https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       STEP   : MULTIQC                          -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

 process MultiQC {
     label "multiQC"
     publishDir "${params.outdir}/multiQC", mode: 'copy'

     input:
     //file report_fastqc from fastqc_report.collect().ifEmpty([])
     //file report_trim from trimming_report.collect().ifEmpty([])
     //file report_adptoRemoval from adapter_removal_report.collect().ifEmpty([])
     file report_mapping from mapping_repport_files.collect()
     file report_picard from picardmetric_files.collect()

     output:
     file "*multiqc_*" into multiqc_report


     when:
     !params.skipMultiqc

     script:
     """
     multiqc -f . \\
     
     """
  }

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
  ${c_blue}    @@       @@@   @@     @@@@@@@@   @@@@     @@     @@@     @@@@ @@@@@@@@@@ 
         @     @    @  @     @@    @@    @      @  @     @ @@    @       @@ 
          @   @    @@@@@@    @@@@@,      @     @@@@@@    @   @@  @       @@ 
           @ @    @     @@   @@    @     @    @     @@   @     @ @       @@ 
            @@   @@@    @@@  @@     @  @@@@  @@     @@  @@       @       @@ 
                                                                       
                                  @@@.       @      @@@     @@@     @@@@@@@   @@@@@@@
                               @       @    @ @      @       @       @@       @     @@ 
                              @@           @   @     @       @       @@@@@@   @ @  
                               @@         @@@@@@@    @       @       @@       @@ @     
                                 @@@@@   @      @@   @@@@@   @@@@@  @@@@@@@   @   @ @
                                                    ${c_purple}  ROMUALD MARIN PIPELINE${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}




/*
 GET EXTENSION FASTA FILE TO TEST IT
 TRY .BAM OF .VCF IF READ ALLREZADY MAPPED

 fasta_file.into{ fasta_file_test_zip;
                  fasta_file_zip }

 ext=fasta_file_test_zip.getVal().getExtension()


  if(ext=="gz" || ext=="bz" || ext=="zip"){

*/
