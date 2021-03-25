#!/usr/bin/env nextflow

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
      ./nextflow run script/VariantCaller.nf -c script/VariantCaller.config --reads "/home/rmarin/ref/*_{1,2}.fq.gz" --genomefasta ../ref/Physcomitrella_patens_Phypa_V3_dna_rm_toplevel.fa -profile psmn --genomeindex "../ref/index/Physcomitrella_patens_Phypa_V3_dna_rm_toplevel.fa*" -resume --annotation ../ref/annotation/Physcomitrella_patens.Phypa_V3.48.gff3 --startingstrain "SAMPLE_V300042688_L2_AE97758923-605"

    Please avoid to use dot in file name, use underscore instead
    ( bad exemple : truc.v2.machin3.fasta , good exemple : truc_v2_machin3.fasta )
    
    Required arguments:
      --reads                    Full path to directory and name of reads in fastq.gz
      --genomefasta              Full path to file of reference genome (.fa or .fasta or .gz) 
      --genomeindex              Full path to directory of index (Optionnal)
      --annotation               Full path to file of annotation genome (.gff) 
    Optionnal arguments:
      --skipbqsr                 Skip base calibration step (default: activated).
      --startingstrain           Name of starting strain (for exemple file : V300042688_L2_AE97758923-605_2.fq.gz , name "V300042688_L2_AE97758923-605" )

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


////////////////////////////////////////////////////
/* --         DEFAULT PARAMETER VALUES         -- */
////////////////////////////////////////////////////

params.help = false
params.genomeindex= false
params.genomefasta = false
params.reads = false
params.annotation = false
params.skipbqsr = false
params.startingstrain = false
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

// possibilité de donner un fichier tableau avec mes noms de sample en fonction des pair_id 

// Faire option pour voir si le dossier resultst n'est pas deja crée (dans le cas mettre un log)
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
summary['annotation']            = params.annotation  ?: 'Not supplied'
summary['startingstrain']        = params.startingstrain  ?: 'Not supplied'
summary['BQSR']                  = params.skipbqsr ? 'Skipped' : 'Yes'
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
            .into { fasta_file ; 
            fasta_file2GATK; 
            fasta3;fasta_variantmetric; 
            fasta_variantcalling ; fasta_joingvcf; 
            fasta_joingvcf_after_bqsr; fasta_extract  ; 
            fasta_extract_after_bqsr; fasta_BaseRecalibrator ; 
            fasta_dict ; fasta_snpeff ; fasta_Snpeff_variant_effect ;
            fasta_Structural_Variant_calling_GATK ; fasta_Structural_Variant_calling_GATK_prepare ;
            fasta_variantcalling_after_bqsr;
            fasta_pindel}

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

if (params.annotation) {
        Channel
            .fromPath( params.annotation )
            .ifEmpty { error "Cannot find any file matching: ${params.annotation}" }
            .set{ annotation }
}

if (params.annotation) {
        Channel
            .value( params.startingstrain )
            .into{ startingstrain ; startingstrain_2 ; startingstrain_SV }
}



/*
 GET EXTENSION FASTA FILE TO TEST IT
 TRY .BAM OF .VCF IF READ ALLREZADY MAPPED

 fasta_file.into{ fasta_file_test_zip;
                  fasta_file_zip }

 ext=fasta_file_test_zip.getVal().getExtension()


  if(ext=="gz" || ext=="bz" || ext=="zip"){

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
  set pair_id, "${pair_id}.sam" into sam_files, sam_files_for_startingstrain
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
  publishDir "${params.outdir}/picard/", mode: 'copy'

  input:
  set pair_id, bam_file from bam_files

  output:
  set pair_id, "${pair_id}_readGroup_MarkDuplicates.bam" into bam_files_RG_MD , bam_files_RG_MD_for_bqsr , bam_for_strartingstrain 
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

/*
if (params.startingstrain) {
  process Get_starting_strain {
    label 'samtools'
    tag "$pair_id"
    publishDir "${params.outdir}/mapping/", mode: 'copy'

    input:
    val startingstrain from startingstrain
    set pair_id, "${pair_id}.sam" from bam_for_strartingstrain

    output:
    file "${startingstrain}.bam" into startingstrainbam

    when:
    pair_id == startingstrain

    script:
    """
    mv ${pair_id}.sam ${startingstrain}.bam
    """
  }
}
*/

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   STEP   : ufilter and Get basics statistics        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

process Filtering_and_Indexing_bam_ {
  label 'samtools'
  tag "$pair_id"
  publishDir "${params.outdir}/filter/", mode: 'copy'

  input:
  set pair_id, bam_RD_MD from bam_files_RG_MD

  output:
  set pair_id, "${pair_id}_ufilter.bam", "${pair_id}_ufilter.bam.bai" into bam_files_RG_MD_filter, bam2GATK , bam2GATK_SV, filtered_bam_files_breakdancer
  set pair_id, "${pair_id}_ufilter_flagstat.txt" into flagstat_files
  set pair_id, "${pair_id}_ufilter.bam.bai" into bam_index_samtools , bam_index_samtools_after_bqsr
  set pair_id, "${pair_id}.bam", "${pair_id}.bam.bai" into nofiltered_bam_pindel

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
  file "*" into fasta_fai , fasta_fai_Structural_Variant_calling_GATK,fasta_fai_Structural_Variant_calling_GATK_prepare, fasta_fai_variantmetric , fasta_fai_gvcftovcf ,fasta_fai_gvcftovcf_after_bqsr, fasta_fai_extract , fasta_fai_extract_after_bqsr , fasta_fai_BaseRecalibrator , fasta_fai_after_bqsr , fasta_fai_pindel

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
  file "*.dict" into fasta_dict_Variant_calling, fasta_dict_Structural_Variant_calling_GATK,fasta_dict_Structural_Variant_calling_GATK_prepare, fasta_dict_Variant_metric , fasta_dict_Gvcf_to_vcf , fasta_dict_Gvcf_to_vcf_after_bqsr , fasta_dict_extract,fasta_dict_extract_after_bqsr , fasta_dict_BaseRecalibrator , fasta_dict_Variant_calling_after_bqsr

  script:
  """
  gatk CreateSequenceDictionary -R $fasta
  """
}

//https://gatk.broadinstitute.org/hc/en-us/articles/360035890411-Calling-variants-on-cohorts-of-samples-using-the-HaplotypeCaller-in-GVCF-mode

process Variant_calling {
  label 'gatk'
  tag "$pair_id"
  publishDir "${params.outdir}/variant/", mode: 'copy'

  input:
  set pair_id, bam_file, bam_file_index from bam2GATK
  file fasta_fai from fasta_fai.collect()
  file fasta from fasta_variantcalling.collect()
  file fasta_dict from fasta_dict_Variant_calling.collect()

  output:
  set pair_id, "${pair_id}.g.vcf.gz" into gvcf
  file "*.g.vcf.gz" into gvcf2metrics

  script:
  """
  gatk HaplotypeCaller \
  -R $fasta \
  -I $bam_file \
  -O "./${pair_id}.g.vcf.gz" \
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


process Gvcf_to_vcf {
  label 'gatk'
  publishDir "${params.outdir}/variant/", mode: 'copy'

  input:
  file file_gvcf from gvcf.collect()
  file fasta_fai from fasta_fai_gvcftovcf
  file fasta from fasta_joingvcf
  file fasta_dict from fasta_dict_Gvcf_to_vcf

  output:
  file "final.vcf.gz" into vcf
  file "*.txt" into gatkmetric_files

  script:
  """
  ## make list of input variant files
  ls *g.vcf.gz > myvcf.list

  for input in \$(cat myvcf.list)
  do
    gatk  IndexFeatureFile -I \$input
  done 
  
  gatk  CombineGVCFs \
  -R $fasta \
  --variant myvcf.list \
  -O combined.g.vcf.gz

  gatk  GenotypeGVCFs \
  -R $fasta \
  -V combined.g.vcf.gz \
  -O final.vcf.gz

  gatk VariantEval  \
   -R $fasta \
   -O "GATK_metrics.txt"\
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

// Ajouter lien vers GATK filtration
process Extract_SNPIndel_Filtration {
  label 'gatk'
  tag "$file_vcf"
  publishDir "${params.outdir}/variant/", mode: 'copy'

  input:
  file file_vcf from vcf
  file fasta_fai from fasta_fai_extract
  file fasta from fasta_extract
  file fasta_dict from fasta_dict_extract

  output:
  file "filtered_snps.vcf" into vcf_snp , vcf_snp_for_snpeff
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
        -O pre_filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

  gatk SelectVariants \
      -R $fasta \
      -V pre_filtered_snps.vcf \
      --exclude-filtered \
      -O filtered_snps.vcf

  gatk SelectVariants \
      -R $fasta \
      -V $file_vcf \
      --select-type-to-include INDEL \
      -O raw_indels.vcf

  gatk VariantFiltration \
        -R $fasta \
        -V raw_indels.vcf \
        -O pre_filtered_indels.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0" 

  gatk SelectVariants \
      -R $fasta \
      -V pre_filtered_indels.vcf \
      --exclude-filtered \
      -O filtered_indels.vcf
  """
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                         STEP   : BQSR                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

if (!params.skipbqsr) {
  process BaseRecalibrator {
    label 'gatk'
    tag "$pair_id"
    publishDir "${params.outdir}/variant_after_bqsr/", mode: 'copy'

    input:
    file fasta_fai from fasta_fai_BaseRecalibrator.collect()
    file fasta from fasta_BaseRecalibrator.collect()
    file fasta_dict from fasta_dict_BaseRecalibrator.collect()
    file bqsr_snps from vcf_snp.collect()
    file bqsr_indels from vcf_indel.collect()
    //file bqsr_snps from bqsr_vcf_snp.collect()
    //file bqsr_indels from bqsr_vcf_indel.collect()
    set pair_id, bam from bam_files_RG_MD_for_bqsr

    output:
    set pair_id, "${pair_id}_recal_reads.bam" into bam_after_bqsr
    
    script:
    """
    gatk  IndexFeatureFile -I $bqsr_snps

    gatk  IndexFeatureFile -I $bqsr_indels

    gatk BaseRecalibrator \
        -R $fasta \
        -I $bam \
        --known-sites $bqsr_snps \
        --known-sites $bqsr_indels \
        -O recal_data.table 
    
    gatk ApplyBQSR \
        -R $fasta \
        -I $bam \
        -bqsr recal_data.table \
        -O ${pair_id}_recal_reads.bam
    """
  }

  process Variant_calling_after_bqsr {
    label 'gatk'
    tag "$pair_id"
    publishDir "${params.outdir}/variant_after_bqsr/", mode: 'copy'

    input:
    set pair_id, bam_file from bam_after_bqsr
    set pair_id, bam_file_index from bam_index_samtools_after_bqsr
    file fasta_fai from fasta_fai_after_bqsr.collect()
    file fasta from fasta_variantcalling_after_bqsr.collect()
    file fasta_dict from fasta_dict_Variant_calling_after_bqsr.collect()

    output:
    set pair_id, "${pair_id}_after_bqsr.g.vcf.gz" into gvcf_after_bqsr

    script:
    """
    gatk  \
    HaplotypeCaller \
    -R $fasta \
    -I $bam_file \
    -O "./${pair_id}_after_bqsr.g.vcf.gz" \
    -ERC GVCF
    """
  }

  process Gvcf_to_vcf_after_bqsr {
    label 'gatk'
    publishDir "${params.outdir}/variant_after_bqsr/", mode: 'copy'

    input:
    file file_gvcf from gvcf_after_bqsr.collect()
    file fasta_fai from fasta_fai_gvcftovcf_after_bqsr
    file fasta from fasta_joingvcf_after_bqsr
    file fasta_dict from fasta_dict_Gvcf_to_vcf_after_bqsr

    output:
    file "final_after_bqsr.vcf.gz" into vcf_after_bqsr

    script:
    """
    ## make list of input variant files
    ls *g.vcf.gz > myvcf.list

    for input in \$(cat myvcf.list)
    do
      gatk  IndexFeatureFile -I \$input
    done 
    
    gatk  CombineGVCFs \
    -R $fasta \
    --variant myvcf.list \
    -O combined.g.vcf.gz

    gatk  GenotypeGVCFs \
    -R $fasta \
    -V combined.g.vcf.gz \
    -O final_after_bqsr.vcf.gz
    """
  }


  process Extract_SNPIndel_Filtration_after_bqsr {
    label 'gatk'
    tag "$file_vcf"
    publishDir "${params.outdir}/variant_after_bqsr/", mode: 'copy'

    input:
    file file_vcf from vcf_after_bqsr
    file fasta_fai from fasta_fai_extract_after_bqsr
    file fasta from fasta_extract_after_bqsr
    file fasta_dict from fasta_dict_extract_after_bqsr

    output:
    file "filtered_snps_after_bqsr.vcf" into vcf_snp_after_bqsr , vcf_snp_for_snpeff_after_bqsr
    file "filtered_indels_after_bqsr.vcf" into vcf_indel_after_bqsr

    script:
    """
    gatk IndexFeatureFile \
      -I $file_vcf
    
    gatk SelectVariants \
        -R $fasta \
        -V $file_vcf \
        --select-type-to-include SNP \
        -O bqsr_snps.vcf

    gatk VariantFiltration \
          -R $fasta \
          -V bqsr_snps.vcf \
          -O pre_filtered_snps_after_bqsr.vcf \
          -filter-name "QD_filter" -filter "QD < 2.0" \
          -filter-name "FS_filter" -filter "FS > 60.0" \
          -filter-name "MQ_filter" -filter "MQ < 40.0" \
          -filter-name "SOR_filter" -filter "SOR > 4.0" \
          -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
          -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
    
    gatk SelectVariants \
      -R $fasta \
      -V pre_filtered_snps_after_bqsr.vcf \
      --exclude-filtered \
      -O filtered_snps_after_bqsr.vcf
    
    gatk SelectVariants \
        -R $fasta \
        -V $file_vcf \
        --select-type-to-include INDEL \
        -O bqsr_indels.vcf

    gatk VariantFiltration \
          -R $fasta \
          -V bqsr_indels.vcf \
          -O pre_filtered_indels_after_bqsr.vcf \
          -filter-name "QD_filter" -filter "QD < 2.0" \
          -filter-name "FS_filter" -filter "FS > 200.0" \
          -filter-name "SOR_filter" -filter "SOR > 10.0" 
    
    gatk SelectVariants \
      -R $fasta \
      -V pre_filtered_indels_after_bqsr.vcf \
      --exclude-filtered \
      -O filtered_indels_after_bqsr.vcf 
    """
  }
}

//
//
//
// https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-
//
//




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   STEP   : Structural variant caller                -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Essayer les verions beta de GATK pour trouver les variants struct. 

/// !!!!!!!!!!!!!!!! Attention a ne pas utiliser les channel de fichier bam filtré

// A tester sur nombre de variant trouvé avec pindel et break dans les 2 cas filter / non filtrer

// precier d'ou viennent les fichier bam (filter, ... ect )

process Structural_Variant_calling_breakdancer {
  label 'breakdancer'
  tag "${pair_id}"
  publishDir "${params.outdir}/stuctural_variant/", mode: 'copy'

  input:
  set pair_id, bam , bai from filtered_bam_files_breakdancer
  //file startingstrainbam from startingstrainbam.collect()

  output:
  file "breakdancer_${pair_id}.ctx" into breakdancer_files
  file "${pair_id}_config.cfg" into confing_breakdancer

  script:
  //Attention la doc est fausse il ne faut pas utiliser de chevron pour bam2cfg sinon le fichier de config n'est pas bon 
  """

  bam2cfg $bam -o ${pair_id}_config.cfg
  
  breakdancer-max ${pair_id}_config.cfg > breakdancer_${pair_id}.ctx
  """
}

//http://gmt.genome.wustl.edu/packages/pindel/user-manual.html
process Structural_Variant_calling_pindel {
  label 'pindel'
  tag "${pair_id}"
  publishDir "${params.outdir}/stuctural_variant/", mode: 'copy'

  input:
  file fasta from fasta_pindel.collect()
  file fasta_fai from fasta_fai_pindel.collect()
  set pair_id, bam, bai from nofiltered_bam_pindel
  file config from confing_breakdancer.collect()

  output:
  file "pindel_${pair_id}.vcf" into pindel_vcf

  script:
  // recupere la taille d'insert avec les fichier config de breackdancer dans la variable leninsert
  // Creation d'un fichier de config pour pindel avec la taille moyenne des inserts 
  """
  leninsert=\$(cut -f9 ${pair_id}_config.cfg | cut -f2 -d ":" | cut -f1 -d ".")

  echo "${bam}	\$leninsert	${pair_id}" > config
  pindel -f $fasta \
  -i config \
  -o ${pair_id}_SV_pindel

  pindel2vcf -r $fasta -R ${fasta.baseName} -P ${pair_id}_SV_pindel -v pindel_${pair_id}.vcf -d 20210101 -G
  """
}



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


if (params.annotation) {
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

    process Snpeff_variant_effect {
        label 'snpeff'
        tag "$file_vcf"
        publishDir "${params.outdir}/snpeff/", mode: 'copy'

        input:
        file configfile from configsnpeff
        file fa from fasta_Snpeff_variant_effect
        file file_vcf from vcf_snp_for_snpeff
        file directorysnpeff from snpfile

        output:
        file "snpeff_${file_vcf}" into anno
        file "snpEff_summary.html" into summary
        file "snpEff_genes.txt" into snptxt

        script:
        """
        snpeff ${fa.baseName} \
        -c $configfile \
        -v $file_vcf \
        > snpeff_${file_vcf}
        """
    }
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