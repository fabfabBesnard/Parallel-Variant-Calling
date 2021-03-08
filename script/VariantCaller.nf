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
log.info "-\033[2m-------------------------------------------------------------------------\033[0m-"

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
            .into { fasta_file ; fasta_file2GATK; fasta3;fasta_variantmetric; fasta_variantcalling ; fasta_joingvcf; fasta_extract ; fasta_dict}

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



process Mapping_reads {
  label 'bwa'
  tag "$pair_id"
  publishDir "${params.outdir}/mapping/", mode: 'copy'

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



process Sam_to_bam {
  label 'samtools'
  tag "$pair_id"
  publishDir "${params.outdir}/mapping/", mode: 'copy'

  input:
  set pair_id, "${pair_id}.sam" from sam_files

  output:
  set pair_id, "${pair_id}.bam" into bam_files, bam_files_breakdancer

  

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

//https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-

process ReadGroup_MarkDuplicates {
  label 'picardtools'
  tag "$pair_id"
  publishDir "${params.outdir}/picard/", mode: 'copy'

  input:
  set pair_id, bam_file from bam_files

  output:
  set pair_id, "${pair_id}_readGroup_MarkDuplicates.bam" into bam_files_RG_MD
  set pair_id, "${pair_id}_marked_dup_metrics.txt" into picardmetric_files

  script:
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


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   STEP   : ufilter and Get basics statistics        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


process Filter_index {
  label 'samtools'
  tag "$pair_id"
  publishDir "${params.outdir}/filter/", mode: 'copy'

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


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   STEP   : Variant calling GATK and get metrics     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//https://gatk.broadinstitute.org/hc/en-us/articles/360037057932-BuildBamIndex-Picard-

process Create_ref_index {
  label 'samtools'
  tag "$fasta"

  input:
  file fasta from fasta_file2GATK

  output:
  file "*" into fasta_fai , fasta_fai_variantmetric , fasta_fai_gvcftovcf , fasta_fai_extract

  script:
  """
  samtools faidx $fasta 
  """
}

process Create_ref_dictionary {
  label 'gatk'
  tag "$fasta"
  
  input:
  file fasta from  fasta_dict

  output:
  file "*.dict" into fasta_dict_Variant_calling, fasta_dict_Variant_metric , fasta_dict_Gvcf_to_vcf , fasta_dict_extract

  script:
  """
  gatk CreateSequenceDictionary -R $fasta
  """
}

process Variant_calling {
  label 'gatk'
  tag "$pair_id"
  publishDir "${params.outdir}/vcf/", mode: 'copy'

  input:
  set pair_id, bam_file from bam2GATK
  set pair_id, bam_file_index from bam_index_samtools
  file fasta_fai from fasta_fai.collect()
  file fasta from fasta_variantcalling.collect()
  file fasta_dict from fasta_dict_Variant_calling.collect()

  output:
  set pair_id, "${pair_id}.g.vcf.gz" into gvcf_before_rename

  script:
  """
  gatk --java-options "-Xmx4G" \
  HaplotypeCaller \
  -R $fasta \
  -I $bam_file \
  -O "./${pair_id}.g.vcf.gz" \
  -ERC GVCF
  """
}

process Add_sample_name {
  label 'picard'
  tag "$pair_id"
  publishDir "${params.outdir}/vcf/", mode: 'copy'

  input:
  set pair_id, filename from gvcf_before_rename

  output:
  file "*.g.vcf.gz" into gvcf, gvcf2metrics

  script:
  """
  PicardCommandLine RenameSampleInVcf \
  INPUT= $filename \
  OUTPUT= ${pair_id}.g.vcf.gz \
  NEW_SAMPLE_NAME= $pair_id
  
  """
}

process Variant_metric {
  label 'gatk'
  publishDir "${params.outdir}/vcf/", mode: 'copy'

  input:
  file file_gvcf from gvcf2metrics.collect()
  file fasta_fai from fasta_fai_variantmetric
  file fasta from fasta_variantmetric
  file fasta_dict from fasta_dict_Variant_metric

  output:
  file "*.txt" into gatkmetric_files

  script:
  """
  ls *g.vcf.gz > myvcf.list

  for input in \$(cat myvcf.list)
  do
    gatk --java-options "-Xmx4G" IndexFeatureFile -I \$input
  done 
  
  gatk --java-options "-Xmx4G" VariantEval  \
   -R $fasta \
   -O "./GATK_metrics.txt"\
   --eval myvcf.list
  """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   STEP   : Join gvcf into vcf                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

process Gvcf_to_vcf {
  label 'gatk'
  publishDir "${params.outdir}/vcf/", mode: 'copy'

  input:
  file file_gvcf from gvcf.collect()
  file fasta_fai from fasta_fai_gvcftovcf
  file fasta from fasta_joingvcf
  file fasta_dict from fasta_dict_Gvcf_to_vcf

  output:
  file "final.vcf.gz" into vcf

  script:
  """
  ## make list of input variant files
  ls *g.vcf.gz > myvcf.list

  for input in \$(cat myvcf.list)
  do
    gatk --java-options "-Xmx4G" IndexFeatureFile -I \$input
  done 
  
  gatk --java-options "-Xmx4G" CombineGVCFs \
  -R $fasta \
  --variant myvcf.list \
  -O combined.g.vcf.gz

  gatk --java-options "-Xmx4G" GenotypeGVCFs \
  -R $fasta \
  -V combined.g.vcf.gz \
  -O final.vcf.gz
  """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   STEP   : Variant filtering                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


process Extract_SNPIndel_Filtration {
  label 'gatk'
  publishDir "${params.outdir}/vcf/", mode: 'copy'

  input:
  file file_vcf from vcf
  file fasta_fai from fasta_fai_extract
  file fasta from fasta_extract
  file fasta_dict from fasta_dict_extract

  output:
  file "filtered_snps.vcf" into vcf_snp
  file "filtered_indels.vcf" into vcf_indel

  script:
  """
  gatk IndexFeatureFile \
     -I $file_vcf
  
  gatk SelectVariants \
      -R $fasta \
      -V $file_vcf \
      --select-type-to-include SNP \
      -O raw_snps.vcf

  gatk VariantFiltration \
        -R $fasta \
        -V raw_snps.vcf \
        -O filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
  
  gatk SelectVariants \
      -R $fasta \
      -V $file_vcf \
      --select-type-to-include INDEL \
      -O raw_indels.vcf

  gatk VariantFiltration \
        -R $fasta \
        -V raw_indels.vcf \
        -O filtered_indels.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0" 
  """
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   STEP   : Variant filtering                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


process Extract_SNPIndel_Filtration {
  label 'gatk'
  publishDir "${params.outdir}/vcf/", mode: 'copy'

  input:
  file file_vcf from vcf
  file fasta_fai from fasta_fai_extract
  file fasta from fasta_extract
  file fasta_dict from fasta_dict_extract

  output:
  file "filtered_snps.vcf" into vcf_snp
  file "filtered_indels.vcf" into vcf_indel

  script:
  """
  gatk IndexFeatureFile \
     -I $file_vcf
  
  gatk SelectVariants \
      -R $fasta \
      -V $file_vcf \
      --select-type-to-include SNP \
      -O raw_snps.vcf

  gatk VariantFiltration \
        -R $fasta \
        -V raw_snps.vcf \
        -O filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
  
  gatk SelectVariants \
      -R $fasta \
      -V $file_vcf \
      --select-type-to-include INDEL \
      -O raw_indels.vcf

  gatk VariantFiltration \
        -R $fasta \
        -V raw_indels.vcf \
        -O filtered_indels.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0" 
  """
}


// A CONTINUER VARIANT RECALIBRATION pour snp et indel 
// https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
// voir script VQSR_dna_raw

// Appel de variant structuraux avec pindel et breakdancer
/*
process SV_calling {
  label 'breakdancer'
  publishDir "${params.outdir}/breakdancer/", mode: 'copy'

  input:
  set pair_id, bam from bam_files_breakdancer

  output:
  file "*.txt" into gatkmetric_files

  script:
  """
  bam2cfg.pl -g -h tumor.bam normal.bam > BRC6.cfg 

  breakdancer_max -t -q 10 -d BRC6.ctx BRC6.cfg > BRC6.ctx 
  """
}

 */

/*
//https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants
process Filter_lowdp {
  label 'snpsift'
  publishDir "${params.outdir}/vcf/", mode: 'copy'

  input:
  file file_vcf from vcf_filter

  output:
  file "group_sample_filter.vcf" into vcf_filter

  script:
  """
  zcat $file_vcf | SnpSift filter "( GEN[ALL].DP >= 3 ) & isHom(GEN[ALL])" > group_sample_filter.vcf

  """
}

#1. Background SNPs
printf "**STEP1** Extract bona fide background SNPs \n"
java -Xmx48g -jar $GATKHOME/GenomeAnalysisTK.jar -T SelectVariants \
-R $REFERENCE \
--variant $input_vcf \
-select 'vc.isSNP()' \
-o out.vcf
cat out.vcf | java -Xmx48g -jar $snpEff_dir/SnpSift.jar filter "countVariant() = 6" > VQSR/$input_base.bckgd_SNP.vcf


#2. Background Indels
printf "**STEP2** Extract bona fide background Indels \n"
java -Xmx48g -jar $GATKHOME/GenomeAnalysisTK.jar -T SelectVariants \
-R $REFERENCE \
--variant $input_vcf \
-select 'vc.isIndel()' \
-o out.vcf
cat out.vcf | java -Xmx48g -jar $snpEff_dir/SnpSift.jar filter "countVariant() = 6" > VQSR/$input_base.bckgd_Indels.vcf
#clean
rm out.vcf*
*/

/*
#variable declaration
snpEff_dir=/applis/PSMN/debian9/software/Generic/snpEff/4.3t/snpEff/
REFERENCE=/home/fbesnard/Reference_genomes/Moss/Physcomitrella_patens.Phypa_V3.dna.toplevel.fa
input_vcf=VXLvsMS.dna.filter.VQSR.vcf
input_base="${input_vcf%.*}"

SS=VXL
mutantfile=/Xnfs/rdpdb/moss/VariantAnalysis_E1/mutagenized_samples.txt

#internal functions
function specific_call () {
cat $1 | java -jar $snpEff_dir/SnpSift.jar filter "(FILTER = 'PASS') & isHom(GEN[${MS[1]}]) & !(GEN[${MS[1]}].GT = './.') & !(GEN[${MS[1]}].GT = GEN[${SS}].GT) & !(GEN[${MS[1]}].GT = GEN[${MS[2]}].GT) & !(GEN[${MS[1]}].GT = GEN[${MS[3]}].GT) & !(GEN[${MS[1]}].GT = GEN[${MS[4]}].GT) & !(GEN[${MS[1]}].GT = GEN[${MS[5]}].GT)"
}

#body
cd $ExecutionDIR

declare -A MS
idx=1
while read line ; do 
                MS[$idx]=$line
                idx=$((idx + 1))
        done < <(cat $mutantfile)  #Process substitution to bypass the subshell created by the while (http://mywiki.wooledge.org/BashFAQ/024)  
echo "Nb of samples in the array: ${#MS[@]}"

printf "Compute specific variants for sample ${MS[1]} \n"
specific_call $input_vcf >> ${MS[1]}.dna.specific.vcf

for ((i=2; i<=${#MS[@]}; i++)); do
                temp=${MS[1]}
                MS[1]=${MS[$i]}
                MS[$i]=$temp
                printf "Compute specific variants for sample ${MS[1]} \n"
                specific_call $input_vcf >> ${MS[1]}.dna.specific.vcf
        done
printf "end \n"
*/


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

    return """-${c_dim}-------------------------------------------------------------------------${c_reset}-
${c_purple}    @@      @@@   @@     @@@@@@@@   @@@@      @@    @@@     @@@@  @@@@@@@@@@ 
${c_purple}      @     @@   @  @     @@    @@    @      @  @     @ @@    @       @@ 
${c_purple}       @   @@   @@@@@@    @@@@@,      @     @@@@@@    @   @@  @       @@ 
${c_green}        @ @@   @     @@   @@    @     @    @     @@   @     @ @       @@ 
${c_green}         @@   @@@    @@@  @@     @  @@@@  @@      @@ @@       @       @@ 
${c_purple}                                                           
${c_cyan}                           @@@.       @      @@@     @@@     @@@@@@@   @@@@@@@
${c_cyan}                        @       @    @ @      @       @       @@       @     @@ 
${c_cyan}                       @@           @   @     @       @       @@@@@@   @@@@@ 
${c_blue}                        @@         @@@@@@@    @       @       @@       @@   @     
${c_blue}                          @@@@@   @      @@   @@@@@   @@@@@  @@@@@@@   @     @@
${c_yellow}                                                       ROMUALD MARIN PIPELINE${c_reset}
-${c_dim}-------------------------------------------------------------------------${c_reset}-
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
