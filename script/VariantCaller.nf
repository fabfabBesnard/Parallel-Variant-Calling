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

A COMPL2TER !!! 


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
      --genomeindex              Full path to directory of index (Optionnal)

    Nextflow config:
      -c                            Path to config file: src/chip_analysis.config (Optionnal)
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
            .into { fasta_file ; fasta_file2GATK; fasta3}

        if (!params.genomeindex){
                  process index_fasta {
                                      label "bwa"
                                      tag "$fasta.simpleName"
                                      publishDir "results/mapping/index/", mode: 'copy'

                                      input:
                                        file fasta  from fasta_file

                                      output:
                                        file "${fasta.baseName}.*" into index_files 
                                        file "*_bwa_report.txt" into index_files_report

                                      script:
                                        """
                                        bwa index -p ${fasta.baseName} ${fasta} \
                                        &> ${fasta.baseName}_bwa_report.txt
                                        """

                                      }
        }
}
/*
else { exit 1,
  log.warn "=================================================================\n" +
           "  WARNING! No genome fasta file precised.\n" +
           "  Use '--genomefasta' \n" +
           "  Or '--help' for more informations \n" +
           "======================================================================="
}
*/
if (params.genomeindex) {
        Channel
            .fromPath( params.genomeindex )
            .ifEmpty { error "Cannot find any file matching: ${params.genomeindex}" }
            //.map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
            .toList()
            .into { index_files ; testindex }
        //testindex.view()

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
  publishDir "${params.outdir}/mapping/${pair_id}/sam", mode: 'copy'

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
  publishDir "${params.outdir}/mapping/${pair_id}/bam", mode: 'copy'

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


process readGroup_MarkDuplicates {
  label 'picardtools'
  tag "$pair_id"
  publishDir "${params.outdir}/picard/${pair_id}/bam", mode: 'copy'

  input:
  set pair_id, bam_file from bam_files

  output:
  set pair_id, "${pair_id}_readGroup_MarkDuplicates.bam" into bam_files_RG_MD
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
/* --                   STEP   : ufilter and Get basics statistics        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
printf "${col} **STEP4** Generating basic stats on bam file ${NC} \n"
lastbam=`ls *.bam`
lastbam_base="${lastbam%.*}"
if [ -n "$ufilter" ]
then
	printf "Filtering out unmapped and unpaired reads from bam file \n" 
	#FLAG / meaning -> 4=unmapped ; 256=not primary alignment ; 3=paired in a proper pair
	samtools view -hu -F4 $lastbam | samtools view -hu -F256 - | samtools view -hb -f3 - > $lastbam_base.ufilter.bam
	samtools flagstat $lastbam_base.ufilter.bam > ./metrics/$lastbam_base.ufilter.flagstat.txt
	if [ -e $lastbam_base.ufilter.bam ]; then    
		rm $lastbam
		else printf "\e[0;31m ##ERROR## Expected filtered bam file was not generated during the filtering step. Script will abort ${NC} \n"
		printf "Last bam generated is: $lastbam"    
		exit 1
	fi
else 
	samtools flagstat $lastbam > ./metrics/$lastbam_base.flagstat.txt
fi
*/

process filter {
  label 'samtools'
  tag "$pair_id"
  publishDir "${params.outdir}/filter/${pair_id}/bam", mode: 'copy'

  input:
  set pair_id, bam_RD_MD from bam_files_RG_MD

  output:
  set pair_id, "${pair_id}_ufilter.bam" into bam_files_RG_MD_filter , bam2GATK
  set pair_id, "${pair_id}_ufilter_flagstat.txt" into flagstat_files
  set pair_id, "${pair_id}_ufilter.bam.bai" into bam_index_samtools

  script:
  """
  samtools view -hu -F4 ${bam_RD_MD} | samtools view -hu -F256 - | samtools view -hb -f3 - > ${pair_id}_ufilter.bam
	samtools flagstat ${pair_id}_ufilter.bam > ${pair_id}_ufilter_flagstat.txt
  #creation de l'index BAI 
  samtools index ${pair_id}_ufilter.bam ${pair_id}_ufilter.bam.bai
  """
}

//testoutfilter.view()

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       STEP   : Index BAM                            -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
#Step5. Index bam
printf "${col} **STEP5** Indexing BAM${NC} \n"
lastbam=`ls *.bam`
lastbam_base="${lastbam%.*}"
java -Xmx${RAM}g -jar $PICARD BuildBamIndex \
INPUT=$lastbam
*/

//https://gatk.broadinstitute.org/hc/en-us/articles/360037057932-BuildBamIndex-Picard-

process index_bam {
  label 'picardtools'
  tag "$pair_id"
  publishDir "${params.outdir}/index_bam/${pair_id}/bam", mode: 'copy'

  input:
  set pair_id, bam_file from bam_files_RG_MD_filter

  output:
  set pair_id, "${pair_id}.bai" into bam_index

  script:
  """
  PicardCommandLine BuildBamIndex \
  INPUT="${bam_file}" \
  OUTPUT="./${pair_id}.bai"
  """
}


/*

#Step6 (optional). Realign around Indels
if [ -n "$IndelRealign" ]
then
	printf "${col} **STEP6** Realign reads around Indel ${NC} \n"
	printf "${col}            ------6.1. Creating interval table------- ${NC} \n"
	java -Xmx${RAM}g -jar $GATKHOME/GenomeAnalysisTK.jar -T RealignerTargetCreator \
	-R $REFERENCE \
	-I $lastbam \
	-o ./metrics/indelrealigner.$OUTPUT_NAME.intervals

	printf "${col}            ------6.2. Realigning Reads in bam------- ${NC} \n"
	java -Xmx${RAM}g -jar $GATKHOME/GenomeAnalysisTK.jar -T IndelRealigner \
	-R $REFERENCE \
	-I $lastbam -o $lastbam_base.realign.bam \
	-targetIntervals ./metrics/indelrealigner.$OUTPUT_NAME.intervals --filter_bases_not_stored

	printf "${col} Check & Clean files before next step ${NC} \n"
	if [ -e $lastbam_base.realign.bam ]; then    
		rm $lastbam_base.ba*
		else printf "\e[0;31m ##ERROR## Expected realigned bam file was not generated. Script will abort ${NC} \n"
		printf "Last bam generated is: $lastbam"    
		exit 1
	fi
fi

#Step7 (optional). BQSR (With HC as caller, only one step is necessary when bootstrapping a first call)
lastbam=`ls *.bam`
lastbam_base="${lastbam%.*}"
if [ -n "$BQSR" ]; then
	printf "${col} **STEP7** BQSR ${NC} \n"
	#Create a folder to contain files related to the BSQR
	mkdir BQSR_files
	printf "${col}            ------7.1. Perform a first variant call (HC)------- ${NC} \n"
	java -jar -Xmx${RAM}g $GATKHOME/GenomeAnalysisTK.jar -T HaplotypeCaller \
	-R $REFERENCE \
	-stand_call_conf 10 \
	-I $lastbam -o ./BQSR_files/$OUTPUT_NAME.HC0.vcf
	printf "${col}            ------7.2. Analyze BaseQuality Covariates by boostrapping the previous vcf------- ${NC} \n"
	java -jar -Xmx${RAM}g $GATKHOME/GenomeAnalysisTK.jar -T BaseRecalibrator \
	-R $REFERENCE \
	-I $lastbam -knownSites ./BQSR_files/$OUTPUT_NAME.HC0.vcf \
	-o ./BQSR_files/$OUTPUT_NAME.BQSR.table
	printf "${col}            ------7.3. Analyze Remaining covariates after BQSR------- ${NC} \n"
	java -jar -Xmx${RAM}g $GATKHOME/GenomeAnalysisTK.jar -T BaseRecalibrator \
	-R $REFERENCE \
	-I $lastbam -knownSites ./BQSR_files/$OUTPUT_NAME.HC0.vcf \
	-BQSR ./BQSR_files/$OUTPUT_NAME.BQSR.table \
	-o ./BQSR_files/$OUTPUT_NAME.post-BQSR.table
	printf "${col}            ------7.4. Generate plots and stats of the BQSR ------- ${NC} \n"
	java -jar -Xmx${RAM}g $GATKHOME/GenomeAnalysisTK.jar -T AnalyzeCovariates \
	-R $REFERENCE \
	-before ./BQSR_files/$OUTPUT_NAME.BQSR.table -after ./BQSR_files/$OUTPUT_NAME.post-BQSR.table \
	-plots ./BQSR_files/$OUTPUT_NAME.BQSRplots.pdf
	printf "${col}            ------7.5. Apply BQSR to the bam ------- ${NC} \n"
	java -jar -Xmx${RAM}g $GATKHOME/GenomeAnalysisTK.jar -T PrintReads \
	-R $REFERENCE \
	-I $lastbam -BQSR ./BQSR_files/$OUTPUT_NAME.BQSR.table \
	-o $lastbam_base.BQSR.bam
	printf "${col} Check & Clean files before next step ${NC} \n"
	if [ -e $lastbam_base.BQSR.bam ]; then    
		rm $lastbam_base.ba*
	else     
		printf "\e[0;31m ##ERROR## Expected BQSR-bam file was not generated. Script will abort ${NC} \n"
		printf "Last bam generated is: $lastbam"
	exit 1
	fi
fi
*/




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       STEP   : Variant calling GATK                 -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
#Step8. Call (With HC as caller, in gVCF mode)
printf "${col} **STEP8** Variant Calling --Genotype $OUTPUT_NAME -- ${NC} \n"
#Take the last bam produced
lastbam=`ls *.bam`
lastbam_base="${lastbam%.*}"
java -jar -Xmx${RAM}g  $GATKHOME/GenomeAnalysisTK.jar -T HaplotypeCaller \
-R $REFERENCE \
--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \
-I $lastbam \
-o $OUTPUT_NAME.gVCF.vcf
*/

//https://gatk.broadinstitute.org/hc/en-us/articles/360037057932-BuildBamIndex-Picard-

process index_ref_for_GATK {
  label 'samtools'
  tag "$fasta"

  input:
  file fasta from fasta_file2GATK

  output:
  file "*" into fasta_fai

  script:
  """
  samtools faidx $fasta 
  """
}

process variant_calling {
  label 'gatk'
  tag "$pair_id"
  publishDir "${params.outdir}/gvf/${pair_id}", mode: 'copy'

  input:
  set pair_id, bam_file from bam2GATK
  set pair_id, bam_file_index from bam_index_samtools
  file fasta_fai from fasta_fai
  file fasta from fasta3


  output:
  file "*.g.vcf.gz" into final_chan
  file "*.txt" into gatkmetric_files

  script:
  """
  gatk CreateSequenceDictionary -R $fasta
  gatk --java-options "-Xmx4G" \
  HaplotypeCaller \
  -R $fasta \
  -I $bam_file \
  -O "./${pair_id}.g.vcf.gz" \
  -ERC GVCF

  #rapport multiQC
  gatk --java-options "-Xmx4G" VariantEval  \
   -R $fasta \
   -O "./${pair_id}.txt"\
   --eval "./${pair_id}.g.vcf.gz" 
  """
}


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
     file report_flagstat from flagstat_files.collect()
     //file report_trim from trimming_report.collect().ifEmpty([])
     //file report_adptoRemoval from adapter_removal_report.collect().ifEmpty([])
     file report_mapping from mapping_repport_files.collect()
     file report_picard from picardmetric_files.collect()
     file report_picard from gatkmetric_files.collect()

     output:
     file "*multiqc_*" into multiqc_report

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

    return """-${c_dim}--------------------------------------------------${c_reset}-
${c_purple}    @@      @@@   @@     @@@@@@@@   @@@@      @@    @@@     @@@@  @@@@@@@@@@ 
${c_purple}      @     @@   @  @     @@    @@    @      @  @     @ @@    @       @@ 
${c_purple}       @   @@   @@@@@@    @@@@@,      @     @@@@@@    @   @@  @       @@ 
${c_green}        @ @@   @     @@   @@    @     @    @     @@   @     @ @       @@ 
${c_green}         @@   @@@    @@@  @@     @  @@@@  @@      @@ @@       @       @@ 
${c_purple}                                                           
${c_cyan}                              @@@.       @      @@@     @@@     @@@@@@@   @@@@@@@
${c_cyan}                           @       @    @ @      @       @       @@       @     @@ 
${c_cyan}                          @@           @   @     @       @       @@@@@@   @@@@@ 
${c_blue}                           @@         @@@@@@@    @       @       @@       @@   @     
${c_blue}                             @@@@@   @      @@   @@@@@   @@@@@  @@@@@@@   @     @@
${c_yellow}                                                           ROMUALD MARIN PIPELINE${c_reset}
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
