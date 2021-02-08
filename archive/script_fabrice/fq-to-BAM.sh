#!/bin/bash
## script: 'fq-to-BAM_PSMN.sh'
## Fabrice Besnard, 2018 v0.1, 
# Requires a config file
# last edit: 2020-03-09

#######################
#       READ ME       #
#######################

#What it does: 
#Script used for P3.p project (Variant analysis of Mutation accumulation lines)
#For a given sample, this script pipes all steps from raw reads to bam files ready te be called for Structural Variations with PINDEL. 
#For this purpose, it is important to not filter out unmapped reads to keep information from broken read pairs (ie, inactivate the option 'ufilter' in the config file)
#GATK's option BQSR and realignment around indels are more dedicated to small Indels and SNV -> they could be inactivated too (at least the time-consuming BQSR step)

#Outline of the script:
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

# required: 
# bwa version: 0.7.5a-r405 or later -> PSMN version (June 2016)=Bwakit/0.7.15
# samtools (PSMN version June 2016= 0.1.18)
# Picard Version: 1.110+. PSMN version (June 2016)= 2.3.0
# GATK 3.3+. PSMN version (June 2016): 3.6 (java1.8 compatible)
#OPTION: check $ufilter $BQSR and $IndelRealign options in the config file.

#Usage:
#sh fq-to-BAM.sh Config-file_mysample.txt
#From where it is started, it creates a subfolder named "mysample" to store output files if it does not exist already. 
#Suggestion: create an 'analysis' folder where you perform Variant analysis for all your samples.

#######################
# Source Config file  #
#######################

CONFIGFILE=$1
. $CONFIGFILE #in tcsh shell, source function is called by '.' 

####################################
#  Inform users of the parameters  #
####################################
printf "${col} **Running fq-to-BAM_PSMN.sh** ${NC} \n"
printf "${col} Reference genome selected: ${NC} $REFERENCE \n"
printf "${col} Numbers of sequencing runs used for mapping: $NbRUNS \n"
for ((a=1; a <=$NbRUNS; a++))
	do
		echo "Sequencing Run $a, Fwd pair is: ${READS[$((2*$a-1))]}"
		echo "Sequencing Run $a, Reverse pair is: ${READS[$((2*$a))]}"
		echo "Sequencing Run $a, RGID is: ${RGID[$a]}"
		echo "Sequencing Run $a, RGPU is: ${RGPU[$a]}"
	done 
printf "${col} Final file generated will be a bam file staring with: ${NC} $OUTPUT_NAME. It contains ReadGroup, where 'sample' is $RGSM \n"
if [ -n "$ufilter" ] || [ -n "$BQSR" ] || [ -n "$IndelRealign" ]; then
printf "Options activated : $ufilter $BQSR $IndelRealign \n"
fi
if [ -n "$ufilter" ]; then 
printf "Info on BAM file: option filter activated. Reads unmapped or secondary alignments will be filtered out so that the bam file will only contain mapped reads in a proper pair.\n"
fi
if [ -n "$IndelRealign" ]; then 
printf "Info on BAM file: option IndelRealigner is activated: reads around indel will be realigned (increase accuracy for small indels)\n"
fi
if [ -n "$BQSR" ]; then 
printf "Info on BAM file: BQSR was applied: quality scores will be modified \n"
fi

printf "${col} Numbers of sequencing runs used for mapping: $NbRUNS \n"

#######################
#     Script Body     #
#######################
#Create a subfolder if necessary
if [ ! -e $OUTPUT_NAME ]
then 
	mkdir $OUTPUT_NAME
else
	echo "a folder called $OUTPUT_NAME already exists. Output data will be stored in it."
fi

cd $OUTPUT_NAME

#Step1. MAPPING
printf "${col} **STEP1** MAPPING pe-reads with bwa (mem, -aM, default) + sort bam output ${NC} \n"
#Note: bwa allowed for 4 threads. samtools sort syntax (with -T) fixed to new samtools version.
for ((a=1; a <=$NbRUNS; a++))
do
	echo "           ------1.$a Map Sequencing Run n°$a ------"
	bwa mem -t 4 -aM $REFERENCE ${READS[$((2*$a-1))]} ${READS[$((2*$a))]} | samtools view -buS -|samtools sort - -T $OUTPUT_NAME -o ${OUTPUT_NAME}.SR$a.bam
done

#Step2. READ GROUPS
printf "${col} **STEP2** Adding Read-Group Informations ${NC} \n"

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

#Step3. Merge and Mark Duplicated Reads at once
printf "${col}**STEP3** Merge all Sequencing Runs and Mark duplicated reads ${NC} \n"
mkdir metrics

#This function allows to feed the MarkDuplicate tool of Picard with all bam previously generated:
function multipleDupe () { #$1 is the number of sequencing Runs=Nber of Bam inputs
	declare -A INPUTS
	for ((a=1; a <=$1; a++))
	do
	INPUTS[$a]="INPUT=$OUTPUT_NAME.RG.SR$a.bam"  #Note: this function only works if $OUTPUT_NAME is defined in the shell env
	done
	java -jar $PICARD MarkDuplicates \
	`echo "${INPUTS[@]}"`\
	OUTPUT=$OUTPUT_NAME.RG.dedup.SRmerged.bam \
	METRICS_FILE=./metrics/dedup.$OUTPUT_NAME.SRmerged.txt
}

multipleDupe $NbRUNS

printf "${col} Check & Clean files before next step ${NC} \n"
if [ -e $OUTPUT_NAME.RG.dedup.SRmerged.bam ]
then
	for ((a=1; a <=$NbRUNS; a++))
	do
		rm $OUTPUT_NAME.RG.SR$a.bam #Delete intermediate .bam files for separate SeqRuns
	done	    
else
	printf "\e[0;31m ##ERROR##Expected merged+dedup bam file was not generated. Script will abort ${NC} \n"    
	exit 1
fi

#Step4. ufilter and Get basics statistics
printf "${col} **STEP4** Generating basic stats on bam file ${NC} \n"
lastbam=`ls *.bam`
lastbam_base="${lastbam%.*}"
if [ -n "$ufilter" ]
then
	printf "Filtering out unmapped and unpaired reads from bam file \n" 
	#FLAG / meaning -> 4=unmapped ; 256=not primary alignment ; 3=paired in a proper pair
	samtools view -hu -F4 $lastbam | samtools view -hu -F256 - | samtools view -hb -f3 - > $lastbam_base.ufilter.bam
	samtools flagstat $lastbam_base.ufilter.bam > ./metrics/$lastbam_base.ufilter.flagstat.txt
else 
	samtools flagstat $lastbam > ./metrics/$lastbam_base.flagstat.txt
fi

#Step5. Index bam
printf "${col} **STEP5** Indexing BAM${NC} \n"
lastbam=`ls *.bam`
lastbam_base="${lastbam%.*}"
java -Xmx${RAM}g -jar $PICARD BuildBamIndex \
INPUT=$lastbam

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
		rm $lastbam_base.BQSR.bam.ba*
	else     
		printf "\e[0;31m ##ERROR## Expected BQSR-bam file was not generated. Script will abort ${NC} \n"
		printf "Last bam generated is: $lastbam"
	exit 1
	fi
fi

################################
# End of script: info to users #
################################
printf "${col} **END OF SCRIPT** Your output bam is: $lastbam ${NC} \n"
#Check the bam
if [ -n "$ufilter" ] && [[ "$lastbam" != *"ufilter"* ]]; then
	printf "\e[0;31m ##Warning## ${NC} Check ufilter step: expected bam file was not generated \n"
fi
if [ -n "$IndelRealigner" ] && [[ "$lastbam" != *"realign"* ]]; then
	printf "\e[0;31m ##Warning## ${NC} Check Indel-realignment step: expected bam file was not generated \n"
fi
if [ -n "$BQSR" ] && [[ "$lastbam" != *"BQSR"* ]]; then
	printf "\e[0;31m ##Warning## ${NC} Check BQSR step: expected bam file was not generated \n"
fi
#End of script
