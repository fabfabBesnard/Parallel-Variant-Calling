#!/usr/bin/env nextflow

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
      ./nextflow run script/VariantCaller.nf -c script/VariantCaller.config --reads "/home/rmarin/ref/*_{1,2}.fq.gz" --genomefasta ../ref/Physcomitrella_patens_Phypa_V3_dna_rm_toplevel.fa -profile psmn --genomeindex "../ref/index/Physcomitrella_patens_Phypa_V3_dna_rm_toplevel.fa*" -resume --annotation ../ref/annotation/Physcomitrella_patens.Phypa_V3.48.gff3  "SAMPLE_V300042688_L2_AE97758923-605"

    Please avoid to use dot in file name, use underscore instead
    ( bad exemple : truc.v2.machin3.fasta , good exemple : truc_v2_machin3.fasta )
    
    Required arguments:
      --reads                    Full path to directory and name of reads in fastq.gz
      --genomefasta              Full path to file of reference genome (.fa or .fasta or .gz) 
      --annotation               Full path to file of annotation genome (.gff) 
      --scriptdir                Full projectdir path !!!!
    Optionnal arguments:
      --genomeindex              Full path to directory of index (Optionnal)
      --vqsrfile                 Variant calibration step, give a refernce fin in vcf with short indel and snp (default: false).
      --sampletable              Table in .csv who contain name of different sample (sep : ',')
      --vqsrrate                 Defaut 99.0

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

// possibiliter de precier une zone de recherche et ou genome masqué/non masqué et utiliser par gatk ou non 
// A voir pour mettre tout ca dans le fichier de config !  https://github.com/nf-core/sarek/blob/master/nextflow.config

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


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
summary['scriptdir']             = params.scriptdir ?: 'Not supplied'
summary['reads']                 = params.reads ?: 'Not supplied'
summary['genomeindex']           = params.genomeindex  ?: 'Not supplied'
summary['genomefasta']           = params.genomefasta  ?: 'Not supplied'
summary['Ploidy']                = params.ploidy
summary['Config Profile']        = workflow.profile
summary['annotationgff']         = params.annotationgff  ?: 'Not supplied'
summary['annotationname']        = params.annotationname
summary['sampletable']           = params.sampletable  ?: 'Not supplied'
summary['minglobalqual']         = params.minglobalqual
summary['mindepth']              = params.mindepth
summary['VQSR']                  = params.vqsrfile
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
            .ifEmpty{ error "Cannot find any file matching: ${params.reads}" }
            .set{ fastqgz }
        Channel
            .fromPath(params.reads)
            .set{ fastqgz_fastqc}
}

if (params.genomefasta) {
  // voir pour faire ne pas creer une channel mais juste une variable reulisable metadata = file("metadata.tsv")
        Channel
            .fromPath( params.genomefasta )
            .ifEmpty { error "Cannot find any file matching: ${params.genomefasta}" }
            .into { fasta_file ; fasta_VQSR ;
            fasta_file2GATK; 
            fasta3;fasta_variantmetric; 
            fasta_variantcalling ; fasta_joingvcf; 
            fasta_extract_Extract_SNP_VQSR ; fasta_Extract_INDEL_VQSR ;
            fasta_BaseRecalibrator ; 
            fasta_dict ; fasta_snpeff ; fasta_Snpeff_variant_effect ;
            fasta_Structural_Variant_calling_GATK ; fasta_Structural_Variant_calling_GATK_prepare ;  fasta_pindel ; fasta_cnv; fasta_metasv}

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


if (params.genomeindex) {
        Channel
            .fromPath( params.genomeindex )
            .ifEmpty { error "Cannot find any file matching: ${params.genomeindex}" }
            //.map { it -> [(it.baseName =~ /([^\.]*)/)[0][1], it]}
            .toList()
            .into { index_files  }
 

}

if (params.annotationgff) {
        Channel
            .fromPath( params.annotation )
            .ifEmpty { error "Cannot find any file matching: ${params.annotation}" }
            .set{ annotation }
}

/*
 GET EXTENSION FASTA FILE TO TEST IT
 TRY .BAM OF .VCF IF READ ALLREZADY MAPPED

 fasta_file.into{ fasta_file_test_zip;
                  fasta_file_zip }

 ext=fasta_file_test_zip.getVal().getExtension()


  if(ext=="gz" || ext=="bz" || ext=="zip"){

*/

process Fastqc {
    label 'fastqc'
    tag "$read"

    input:
    file read from fastqgz_fastqc

    output:
    file "*.{zip,html}" into fastq_repport_files

    script:
    """
    fastqc --quiet --threads ${task.cpus} --format fastq --outdir ./ \
             ${read}
    """
  }

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                          STEP   : MAPPING                           -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//Change les noms des pair_id avec les noms des echantillons directement dans la premiere channel
if (params.sampletable) {
      Channel
        .fromPath( params.sampletable )
        .splitCsv( skip: 0)
        .join(fastqgz)
        .set{fastqsamplename}

  process Mapping_reads_and_add_sample_name {
    label 'bwa'
    tag "$sample_id"
    publishDir "${params.outdir}/mapping/", mode: 'copy'

    input:
    set pair_id,sample_id, file(reads) from fastqsamplename
    file index from index_files.collect()

    output:
    set sample_id, "${sample_id}.sam" into sam_files
    set sample_id, "${sample_id}_bwa_report.txt" into mapping_repport_files

    script:
    index_id = index[0].baseName
    """
    bwa mem -t ${task.cpus} \
    -aM ${index_id} ${reads[0]} ${reads[1]} \
    -o ${sample_id}.sam &> ${sample_id}_bwa_report.txt
    """
  }
}
else{
  process Mapping_reads {
    label 'bwa'
    tag "$pair_id"
    publishDir "${params.outdir}/mapping/", mode: 'copy'

    input:
    set pair_id,file(reads) from fastqgz
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

//https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-

process Add_ReadGroup_and_MarkDuplicates_bam {
  label 'picardtools'
  tag "$pair_id"

  input:
  set pair_id, bam_file from bam_files

  output:
  set pair_id, "${pair_id}_readGroup_MarkDuplicates.bam" into bam_files_RG_MD , bam_for_strartingstrain 
  set pair_id, "${pair_id}_marked_dup_metrics.txt" into picardmetric_files

  script:
  """
  PicardCommandLine AddOrReplaceReadGroups \
       I=${bam_file} \
       O="${pair_id}_readGroup.bam" \
       RGID=${pair_id} \
       RGLB=lib1 \
       RGPL=illumina \
       RGPU=unit1 \
       RGSM=${pair_id}
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

process Filtering_and_Indexing_bam {
  label 'samtools'
  tag "$pair_id"

  input:
  set pair_id, bam_RD_MD from bam_files_RG_MD

  output:
  set pair_id, "${pair_id}_ufilter.bam", "${pair_id}_ufilter.bam.bai" into bam_files_RG_MD_filter,filtered_bam_pindel, bam2GATK , bam2GATK_SV, filtered_bam_files_breakdancer
  set pair_id, "${pair_id}_ufilter_flagstat.txt" into flagstat_files
  set pair_id, "${pair_id}_ufilter.bam.bai" into bam_index_samtools
  set pair_id, "${pair_id}.bam",  "${pair_id}.bam.bai" into non_filtered_bam_pindel , bam_for_lumpy_1,bam_for_lumpy_2, non_filtered_bam_files_cnvnator, non_filtered_bam_files_breakdancer , non_filtered_bam_files_metasv

  script:
  """
  samtools view -hu -F4 ${bam_RD_MD} | samtools view -hu -F256 - | samtools view -hb -f3 - > ${pair_id}_ufilter.bam
	samtools flagstat ${pair_id}_ufilter.bam > ${pair_id}_ufilter_flagstat.txt
  samtools index ${pair_id}_ufilter.bam ${pair_id}_ufilter.bam.bai
 
  #index pour bam avant filtration
  cp ${bam_RD_MD} ${pair_id}.bam
  samtools index ${pair_id}.bam ${pair_id}.bam.bai

  """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --               STEP   : Variant calling GATK and get metrics         -- */
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
  file "*" into fasta_fai , fasta_fai_VQSR, fasta_fai_Structural_Variant_calling_GATK,fasta_fai_Structural_Variant_calling_GATK_prepare, fasta_fai_variantmetric , fasta_fai_gvcftovcf ,fasta_fai_gvcftovcf_after_bqsr, fasta_fai_extract_Extract_SNP_VQSR ,fasta_fai_Extract_INDEL_VQSR , fasta_fai_BaseRecalibrator , fasta_fai_after_bqsr , fasta_fai_pindel , fastafai_metasv

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
  file "*.dict" into fasta_dict_Variant_calling, fasta_dict_VQSR, fasta_dict_Structural_Variant_calling_GATK,fasta_dict_Structural_Variant_calling_GATK_prepare, fasta_dict_Variant_metric , fasta_dict_Gvcf_to_vcf , fasta_dict_Gvcf_to_vcf_after_bqsr , fasta_dict_extract_Extract_SNP_VQSR,fasta_dict_Extract_INDEL_VQSR,fasta_dict_extract_after_bqsr , fasta_dict_BaseRecalibrator , fasta_dict_Variant_calling_after_bqsr

  script:
  """
  gatk  --java-options "-Xmx${task.memory.giga}g" CreateSequenceDictionary -R $fasta
  """
}

//https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode

/*  --sample-ploidy / -ploidy

Ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy).
Sample ploidy - equivalent to number of chromosomes per pool. In pooled experiments this should be = # of samples in pool * individual sample ploidy

int  2  [ [ -∞  ∞ ] ]  */


process Variant_calling {
  label 'gatk'
  tag "$pair_id"
  publishDir "${params.outdir}/variant/", mode: 'copy'

  input:
  set pair_id, bam_file, bam_file_index from bam2GATK
  file fasta_fai from fasta_fai.collect()
  file fasta from fasta_variantcalling.collect()
  file fasta_dict from fasta_dict_Variant_calling.collect()
  val p from params.ploidy

  output:
  set pair_id, "${pair_id}.g.vcf.gz" into gvcf
  file "*.g.vcf.gz" into gvcf2metrics

  script:
  """
  gatk --java-options "-Xmx${task.memory.giga}g" HaplotypeCaller \
  -R $fasta \
  -I $bam_file \
  -O "./${pair_id}.g.vcf.gz" \
  -ploidy $p \
  -ERC GVCF
  """
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   STEP   : Join gvcf into vcf                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//https://gatk.broadinstitute.org/hc/en-us/articles/360056967892-VariantEval-BETA-

// Voir doc eval pour comparoverlap ???


process Join_Gvcf_to_vcf {
  label 'gatk'
  tag "and compute metrics"
  publishDir "${params.outdir}/variant/", mode: 'copy'

  input:
  file file_gvcf from gvcf.collect()
  file fasta_fai from fasta_fai_gvcftovcf
  file fasta from fasta_joingvcf
  file fasta_dict from fasta_dict_Gvcf_to_vcf

  output:
  file "final.vcf.gz" into vcf_for_indel , vcf_for_snp
  file "*.txt" into gatkmetric_files

  script:
  """
  ## make list of input variant files
  ls *g.vcf.gz > myvcf.list

  for input in \$(cat myvcf.list)
  do
    gatk --java-options "-Xmx${task.memory.giga}g" IndexFeatureFile -I \$input
  done 
  
  gatk --java-options "-Xmx${task.memory.giga}g" CombineGVCFs \
  -R $fasta \
  --variant myvcf.list \
  -O combined.g.vcf.gz

  gatk --java-options "-Xmx${task.memory.giga}g" GenotypeGVCFs \
  -R $fasta \
  -V combined.g.vcf.gz \
  -O final.vcf.gz

  gatk --java-options "-Xmx${task.memory.giga}g" VariantEval  \
   -R $fasta \
   -O "GATK_metrics.txt"\
   -EV CompOverlap -EV IndelSummary -EV CountVariants -EV MultiallelicSummary \
   --eval final.vcf.gz
  """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   STEP   : Variant filtering                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// https://gatk.broadinstitute.org/hc/en-us/articles/360035532412-Can-t-use-VQSR-on-non-model-organism-or-small-dataset
//https://gatk.broadinstitute.org/hc/en-us/articles/360036350452-VariantFiltration
process Extract_SNP_and_filtering {
  label 'gatk'
  tag "$file_vcf"
  publishDir "${params.outdir}/variant/", mode: 'copy'

  input:
  file file_vcf from vcf_for_snp
  file fasta_fai from fasta_fai_extract_Extract_SNP_VQSR
  file fasta from fasta_extract_Extract_SNP_VQSR
  file fasta_dict from fasta_dict_extract_Extract_SNP_VQSR

  output:
  file "filtered_snps.vcf" into vcf_snp , snp_files, bqsr_vcf_snp
  file "raw_snps.vcf" into rawsnp

  script:
  """
  gatk --java-options "-Xmx${task.memory.giga}g" IndexFeatureFile \
     -I $file_vcf
  
  gatk --java-options "-Xmx${task.memory.giga}g" SelectVariants \
      -R $fasta \
      -V $file_vcf \
      --select-type-to-include SNP \
      -O raw_snps.vcf

  gatk --java-options "-Xmx${task.memory.giga}g" VariantFiltration \
        -R $fasta \
        -V raw_snps.vcf \
        -O pre_filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

  gatk --java-options "-Xmx${task.memory.giga}g" SelectVariants \
      -R $fasta \
      -V pre_filtered_snps.vcf \
      --exclude-filtered \
      -O filtered_snps.vcf
  """
}

process Extract_INDEL_and_filtering {
  label 'gatk'
  tag "$file_vcf"
  publishDir "${params.outdir}/variant/", mode: 'copy'

  input:
  file file_vcf from vcf_for_indel
  file fasta_fai from fasta_fai_Extract_INDEL_VQSR
  file fasta from fasta_Extract_INDEL_VQSR
  file fasta_dict from fasta_dict_Extract_INDEL_VQSR

  output:
  file "filtered_indels.vcf" into indel , indel_files , bqsr_vcf_indel , testgatkmetasv
  file "raw_indels.vcf" into rawindel

  script:
  """
  gatk --java-options "-Xmx${task.memory.giga}g" IndexFeatureFile \
     -I $file_vcf
  
  gatk --java-options "-Xmx${task.memory.giga}g" SelectVariants \
      -R $fasta \
      -V $file_vcf \
      --select-type-to-include INDEL \
      -O raw_indels.vcf

  gatk --java-options "-Xmx${task.memory.giga}g" VariantFiltration \
        -R $fasta \
        -V raw_indels.vcf \
        -O pre_filtered_indels.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0" 

  gatk --java-options "-Xmx${task.memory.giga}g" SelectVariants \
      -R $fasta \
      -V pre_filtered_indels.vcf \
      --exclude-filtered \
      -O filtered_indels.vcf
  """
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         STEP   : VQSR                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
 //https://gatk.broadinstitute.org/hc/en-us/articles/360051306591-ApplyVQSR
 //https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator
if (params.vqsrfile) {

  Channel
            .fromPath( params.vqsrfile )
            .ifEmpty { error "Cannot find any file matching: ${params.vqsrfile}" }
            .set{ vcfforvqsr }

  process  VariantQualityScoreRecalibration {
    label 'gatk'
    tag "snp + indel"
    publishDir "${params.outdir}/variant/", mode: 'copy'

    input:
    file rawsnp from rawsnp
    file rawindel from rawindel
    file fasta_fai from fasta_fai_VQSR
    file fasta from fasta_VQSR
    file fasta_dict from fasta_dict_VQSR
    file vqsrvcf from vcfforvqsr
    val percent from params.vqsrrate

    output:
    file "filtered_indels_VQSR.vcf" into indelVQSR
    file "filtered_snps_VQSR.vcf" into snpVQSR
    
    script:
    """
    ## make list of input variant files
    ls *.vcf* > myvcf.list

    for input in \$(cat myvcf.list)
    do
       gatk --java-options "-Xmx${task.memory.giga}g" IndexFeatureFile -I \$input
    done 

  
    gatk --java-options "-Xmx${task.memory.giga}g" SelectVariants \
      -R $fasta \
      -V $vqsrvcf \
      --select-type-to-include INDEL \
      -O indels.vcf
    
    gatk --java-options "-Xmx${task.memory.giga}g" SelectVariants \
      -R $fasta \
      -V $vqsrvcf \
      --select-type-to-include SNP \
      -O snps.vcf

    gatk --java-options "-Xmx${task.memory.giga}g" VariantRecalibrator \
        -R $fasta \
        -V $rawsnp \
        -tranche $percent \
        --resource:backgroundSNP,known=true,training=true,truth=true,prior=10.0 snps.vcf \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP \
        -O output_snps.recal \
        --tranches-file output_snps.tranches \
        --rscript-file output_snps.plots.R \
        --max-gaussians 4

    gatk --java-options "-Xmx${task.memory.giga}g" ApplyVQSR \
        -R $fasta \
        -V $rawsnp \
        -O pre-filtered_snps_VQSR.vcf \
        --truth-sensitivity-filter-level $percent \
        --tranches-file output_snps.tranches \
        --recal-file output_snps.recal \
        -mode SNP
    
    gatk --java-options "-Xmx${task.memory.giga}g" SelectVariants \
      -R $fasta \
      -V pre-filtered_snps_VQSR.vcf \
      --exclude-filtered \
      -O filtered_snps_VQSR.vcf

    gatk --java-options "-Xmx${task.memory.giga}g" VariantRecalibrator \
      -R $fasta \
      -V $rawindel \
      --trust-all-polymorphic \
      -tranche $percent \
      -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
      -mode INDEL \
      --resource:backgroundindel,known=true,training=true,truth=true,prior=10.0 indels.vcf \
      -O cohort_indels.recal \
      --max-gaussians 4 \
      --tranches-file cohort_indels.tranches

    gatk --java-options "-Xmx${task.memory.giga}g" ApplyVQSR \
        -R $fasta \
        -V $rawindel \
        -O pre-filtered_indels_VQSR.vcf \
        --truth-sensitivity-filter-level $percent \
        --tranches-file cohort_indels.tranches \
        --recal-file cohort_indels.recal \
        -mode INDEL
    
    gatk --java-options "-Xmx${task.memory.giga}g" SelectVariants \
      -R $fasta \
      -V pre-filtered_indels_VQSR.vcf \
      --exclude-filtered \
      -O filtered_indels_VQSR.vcf
    """
  }

  process Extract_specific_variant_VQSR{
        // No label ! Launch with PSMN configuration
        tag "$vcf"
        publishDir "${params.outdir}/variant/", mode: 'copy'

        input:
        val scriptpath from params.scriptdir
        val qual from params.minglobalqual
        val dpmin from params.mindepth
        //val wt from params.WT 
        file vcf from indelVQSR.concat(snpVQSR)

        output:
        file "*.vcf" into good_variant

        script:
        """
        $scriptpath/extract_specific.py $vcf $qual $dpmin
        """
    }

}
else{
  process Extract_specific_variant{
        // No label ! Launch with PSMN configuration
        tag "$vcf"
        publishDir "${params.outdir}/variant/", mode: 'copy'

        input:
        val scriptpath from params.scriptdir
        val qual from params.minglobalqual
        val dpmin from params.mindepth
        //val wt from params.WT 
        file vcf from snp_files.concat(indel_files)

        output:
        file "*.vcf" into good_variant

        script:
        """
        $scriptpath/extract_specific.py $vcf $qual $dpmin
        """
    }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   STEP   : Structural variant caller                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

process Structural_Variant_calling_breakdancer {
  label 'breakdancer'
  tag "${pair_id}"
  publishDir "${params.outdir}/structural_variant/", mode: 'copy'

  input:
  set pair_id, bam , bai from non_filtered_bam_files_breakdancer
  val scriptpath from params.scriptdir

  output:
  set  pair_id, "breakdancer_${pair_id}.ctx" into breakdancer_files , breackdancer_metasv
  file "${pair_id}_config.cfg" into config_breakdancer
  set pair_id, "breakdancer_${pair_id}.vcf" into breakdancer_vcf

  script:
  //Attention la doc est fausse il ne faut pas utiliser de chevron pour bam2cfg sinon le fichier de config n'est pas bon 
  """

  bam2cfg $bam -o ${pair_id}_config.cfg
  
  breakdancer-max ${pair_id}_config.cfg > breakdancer_${pair_id}.ctx

  # Transform ctx into vcf
  python $scriptpath/breakdancer2vcf.py -i 'breakdancer_${pair_id}.ctx' -o 'breakdancer_${pair_id}.vcf'
  """
}


// https://github.com/ALLBio/allbiotc2/blob/master/breakdancer/breakdancer2vcf.py


//http://gmt.genome.wustl.edu/packages/pindel/user-manual.html
process Structural_Variant_calling_pindel {
  label 'pindel'
  tag "${pair_id}"
  publishDir "${params.outdir}/structural_variant/", mode: 'copy'

  input:
  file fasta from fasta_pindel.collect()
  file fasta_fai from fasta_fai_pindel.collect()
  set pair_id, bam, bai from non_filtered_bam_pindel
  file config from config_breakdancer.collect()
  val p from params.ploidy

  output:
  file "${pair_id}_SV_pindel*" into pindel_metasv
  file "pindel_${pair_id}.vcf" into pindel_vcf

  script:
  // recupere la taille d'insert avec les fichier config de breackdancer dans la variable leninsert
  // Creation d'un fichier de config pour pindel avec la taille moyenne des inserts 
  """
  leninsert=\$(cut -f9 ${pair_id}_config.cfg | cut -f2 -d ":" | cut -f1 -d ".")

  echo "${bam}	\$leninsert	${pair_id}" > config

  grep ">" $fasta | awk '/^>/ -F " "{  print \$1" $p" }' | sed 's/>//' > fileploidy

  pindel -f $fasta \
  -T ${task.cpus} \
  -i config \
  --Ploidy fileploidy \
  -o ${pair_id}_SV_pindel

  pindel2vcf -r $fasta -R ${fasta.baseName} -P ${pair_id}_SV_pindel -v pindel_${pair_id}.vcf -d 20210101 -G
  """
}

process Structural_Variant_calling_Lumpy {
  label 'lumpy'
  tag "${pair_id}"
  publishDir "${params.outdir}/structural_variant/", mode: 'copy'

  input:
  val scriptpath from params.scriptdir
  set pair_id , bam , bai from bam_for_lumpy_2

  output:
  file "Lumpy_${pair_id}.vcf" into lumpy_out

  script:
  """
  # Extract the discordant paired-end alignments.
  samtools view -b -F 1294 $bam > sample.discordants.unsorted.bam

  samtools view -h $bam | $scriptpath/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > sample.splitters.unsorted.bam
  
  # Sort both alignments
  samtools sort sample.discordants.unsorted.bam -o sample.discordants.bam
  samtools sort sample.splitters.unsorted.bam -o sample.splitters.bam

  #Run Lumpy in express mode
  lumpyexpress \
    -B $bam \
    -S sample.splitters.bam \
    -D sample.discordants.bam \
    -o Lumpy_${pair_id}.vcf
  """
}

//https://github.com/abyzovlab/CNVnator
process Structural_Variant_calling_CNVnator {
  label 'cnvnator'
  tag "${pair_id}"
  publishDir "${params.outdir}/structural_variant/", mode: 'copy'

  input:
  file fasta from fasta_cnv.collect()
  set pair_id, bam, bai from non_filtered_bam_files_cnvnator

  output:
  file "${pair_id}_CNV.call" into cnvnator_out

  script:
  """
  cnvnator -root file.root -tree $bam

  cnvnator -root file.root -his 100 -fasta $fasta

  cnvnator -root file.root -stat 100

  cnvnator -root file.root -partition 100

  cnvnator -root file.root -call 100  > ${pair_id}_CNV.call
  """
}

// ---------------------
// Breakseq http://bioinform.github.io/breakseq2/ 
// https://hub.docker.com/r/szarate/breakseq2
// ---------------------

//https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4528635/
//http://bioinform.github.io/metasv/
process Group_Structural_Variant_with_Metasv{
        label 'metasv'
        tag "$pair_id"
        publishDir "${params.outdir}/structural_variant/", mode: 'copy'

        input:
        file fasta from fasta_metasv.collect()
        file fastafai from fastafai_metasv.collect()
        file pindelout from pindel_metasv.collect()  
        //file gatkindel from testgatkmetasv.collect() // #--gatk_vcf $gatkindel \
        file cnv from cnvnator_out.collect()
        file lumpy from lumpy_out.collect()
        set pair_id, breakdancerout , bam , bamindex from breackdancer_metasv.join( non_filtered_bam_files_metasv ) 

        output:
        //set pair_id, "${pair_id}_SV.vcf" into metasvout
        file "${pair_id}_SV.vcf" into vcfmetasv_withnonspecific
        file "raw_${pair_id}_SV.vcf" into raw_metasv
        val pair_id into id

        script:
        """
        grep -v "IMPRECISE" Lumpy_${pair_id}.vcf  > Lumpy.vcf


        run_metasv.py \
        --num_threads ${task.cpus} \
        --reference $fasta \
        --breakdancer_native $breakdancerout \
        --pindel_native ${pair_id}_SV_pindel* \
        --cnvnator_native ${pair_id}_CNV.call \
        --lumpy_vcf Lumpy.vcf\
        --outdir out \
        --sample $pair_id \
        --filter_gaps \
        --bam $bam \
        --minsvlen 5 \
        --disable_assembly \
        #--spades spades.py \
        #--age age_align \
        --keep_standard_contigs

        gunzip out/variants.vcf.gz
        mv out/variants.vcf raw_${pair_id}_SV.vcf
        grep -v "LowQual" raw_${pair_id}_SV.vcf | grep -v "IMPRECISE" > ${pair_id}_SV.vcf

        """
}


process Find_specific_SV{
        publishDir "${params.outdir}/structural_variant/", mode: 'copy'

        input:
        val scriptpath from params.scriptdir
        file SV from vcfmetasv_withnonspecific.collect()

        output:
        file "*filtered_SV.vcf" into vcfmetasv 

        script:
        """
        $scriptpath/extract_specific_SV.py
        """
}

/*
process Prepare_Structural_Variant_calling_GATK {
    label 'gatk'
    
    input:
    file fasta from fasta_Structural_Variant_calling_GATK_prepare.collect()
    file fasta_fai from fasta_fai_Structural_Variant_calling_GATK_prepare.collect()
    file fasta_dict from fasta_dict_Structural_Variant_calling_GATK_prepare.collect()

    output:
    file "*" into file_to_GATK_SV

    script:
    """
    gatk BwaMemIndexImageCreator \
     -I $fasta \
     -O reference.img
    
    gatk FindBadGenomicKmersSpark \
    -R $fasta \
    -O kmers_to_ignore.txt
    """
  }
  
  process Structural_Variant_calling_GATK {
    label 'gatk'
    tag "$pair_id"
    publishDir "${params.outdir}/stuctural_variant/", mode: 'copy'

    input:
    set pair_id, bam_file, bai from bam2GATK_SV
    file fasta_fai from fasta_fai_Structural_Variant_calling_GATK.collect()
    file fasta from fasta_Structural_Variant_calling_GATK.collect()
    file fasta_dict from fasta_dict_Structural_Variant_calling_GATK.collect()
    file vgdjfhhsdkbh from file_to_GATK_SV.collect()

    output:
    file "*vcf" into GATK_SV

    script:
    """
    gatk BwaMemIndexImageCreator \
     -I $fasta \
     -O reference.img
    
    gatk FindBadGenomicKmersSpark \
    -R $fasta \
    -O kmers_to_ignore.txt

    gatk FindBreakpointEvidenceSpark \
     -I $bam_file \
     --aligner-index-image reference.img \
     --kmers-to-ignore kmers_to_ignore.txt \
     -O aligned_contigs.sam

    gatk StructuralVariationDiscoveryPipelineSpark \
     -I $bam_file \
     -R $fasta \
     --aligner-index-image reference.img \
     --kmers-to-ignore kmers_to_ignore.txt \
     --contig-sam-file aligned_contigs.sam \
     -O GATK_${pair_id}.vcf 
    """
  }
*/


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   STEP   : Snp effect                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Peu de docker disponible, aucun pour la derniere version de snpEff 
// voir pour en crer un ?? 
// Voir pour faire option pour database deja crée et lui passer le chemin du rep data ? 
// voir piur utiliser la database deja crée dans snpeff 
// voir pour la localisation du config file ? possibilité d'aller chercher dans un autre repertoire


if (params.annotationgff) {
    
    //A modifier car actuellement le docker contient deja la base de données pour Physcomitrella_patens

    process Snpeff_build_database {
        label 'snpeff'
        tag "$gff"

        input:
        file gff from annotation
        file fa from fasta_snpeff

        output:
        file "snpEff.config" into configsnpeff
        file "data/*" into snpfile

        script:
        """
        mkdir data
        mkdir ./data/${fa.baseName}/

        cp $gff ./data/${fa.baseName}/genes.gff
        cp $fa ./data/${fa.baseName}/sequences.fa

        echo "#Physcomitrium (Physcomitrella) patens (${fa.baseName}, 11/03/2021)" >> snpEff.config
        echo "${fa.baseName}.genome : Physcomitrium patens"  >> snpEff.config

        snpeff build -gff3 -c snpEff.config -v ${fa.baseName}
        """
        }

      //http://pcingola.github.io/SnpEff/ss_extractfields/
    process Snpeff_variant_effect_withgff {
        label 'snpeff'
        tag "$file_vcf"
        publishDir "${params.outdir}/snpeff/$file_vcf", mode: 'copy'

        input:
        file configfile from configsnpeff.collect()
        file fa from fasta_Snpeff_variant_effect.collect()
        file file_vcf from good_variant.flatten().concat(vcfmetasv)
        file directorysnpeff from snpfile.collect()

        output:
        file "snpeff_${file_vcf}" into anno
        file "snpEff_summary.html" into summary
        file "snpEff_genes.txt" into snptxt
        file "tab_snpeff_${file_vcf}" into tabsnpeff

        script:
        """
        snpeff ${fa.baseName} \
        -c $configfile \
        -v $file_vcf > snpeff_${file_vcf}

        cat snpeff_${file_vcf} | vcfEffOnePerLine.pl | snpsift extractFields -e '.' - "ANN[*].GENEID" "ANN[*].GENE" CHROM POS REF ALT DP "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].BIOTYPE" | uniq -u > tab_snpeff_${file_vcf}
        """
      }
}


if (params.annotationname) {
    //http://pcingola.github.io/SnpEff/ss_extractfields/
    process Snpeff_variant_effect {
        label 'snpeff'
        tag "$file_vcf"
        publishDir "${params.outdir}/snpeff/$file_vcf", mode: 'copy'

        input:
        file fa from fasta_Snpeff_variant_effect.collect()
        file file_vcf from good_variant.flatten().concat(vcfmetasv.flatten())

        output:
        file "snpeff_${file_vcf}" into anno
        file "snpEff_summary.html" into summary
        file "snpEff_genes.txt" into snptxt
        file "tab_snpeff_${file_vcf}" into tabsnpeff

        script:
        """
        snpeff $params.annotationname \
        -v $file_vcf > snpeff_${file_vcf}

        cat snpeff_${file_vcf} | vcfEffOnePerLine.pl | snpsift extractFields -e '.' - "ANN[*].GENE" CHROM POS REF ALT DP "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].BIOTYPE" | uniq -u > tab_snpeff_${file_vcf}
        """
      }
}

process Final_process {
      tag "$pair_id"
      publishDir "${params.outdir}/final", mode: 'copy'

      input:
      file fi from tabsnpeff.collect()
      val pair_id from id

      output:
      file "${pair_id}.tsv" into final_files

      script:
      """
      liste_fichiers=`ls *${pair_id}*`

      for fichier in \$liste_fichiers
      do
        head -1 \$fichier > ${pair_id}.tsv
      done
      
      for i in \$liste_fichiers
      do 
       sed '1d;\$d' \$i >> ${pair_id}.tsv
      done
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
     file report_fastqc from fastq_repport_files.collect()
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
${c_purple}  @@      @@@   @@     @@@@@@@@   @@@@      @@    @@@     @@@@  @@@@@@@@@@ 
${c_purple}    @     @@   @  @     @@    @@    @      @  @     @ @@    @       @@ 
${c_purple}     @   @@   @@@@@@    @@@@@,      @     @@@@@@    @   @@  @       @@ 
${c_green}      @ @@   @     @@   @@    @     @    @     @@   @     @ @       @@ 
${c_green}       @@   @@@    @@@  @@     @  @@@@  @@      @@ @@       @       @@ 
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