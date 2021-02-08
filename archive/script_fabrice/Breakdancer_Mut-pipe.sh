#! bin/bash

###########################################################################
#### Copyright (C) 2020 - INRAe (Fabrice Besnard, RDP)
#### This script is a free software: you can redistribute it
#### and/or modify it under the terms of the GNU General Public
#### License (GNU GPL) as published by the Free Software Foundation, either
#### version 3 of the License, or (at your option) any later version.
#### Distributed without any warranty.
###########################################################################
# last edit: 2020-06-25
#Version v3
#Tested on PSMN

#######################
#       READ ME       #
#######################

## What it does:
#Takes as inputs the bam file of a mutagenized strain (MS) and the bam file of the ancestral starting strain (SS) before mutagenesis and perform a breakdancer call for SV, substraction of two sets of variations (to keep only MS-specific variations) and filtering. 
#Ideally, bam files should be generated with non-masked genome version, and unmapped read pairs should be filtered out.
#It outputs one file containing the MS-specific SV in the vcf format
#If provided, It can start from a vcf for the SS (instead of a bam file)

## Required: 
#Breakdancer (tested: version 1.4.5-unstable-66-4e44b43)
#Bedtools (tested: v2.25.0)
#Python scripts (6 scripts): breakdancer2vcf.2.py / CompareVCF.py / Vcf_to_DF.py / HardFilter.py / Filter_vcf_from_DF.py / Filter_coverage.py

#Usage:
#bash Breakdancer_Mut-pipe.sh config_file

#######################
# Source Config file  #
#######################
CONFIGFILE=$1
source $CONFIGFILE 

####################################
#  Inform users of the parameters  #
####################################
printf "${col} **Running Breakdancer_Mut-pipe.sh** ${NC} \n"
if [ -n "$SS_bam" ] && [ -n "$SS_vcf" ]; then
	printf "Either choose a bam or a vcf for the starting strain. Both were selected. \n"
	printf "Script will abort. \n"
	exit 1
fi

#######################
#     Script Body     #
#######################
#Create a subfolder if necessary
if [ ! -e $MS_NAME.vs.$SS_NAME.BKD ]
then 
	mkdir $MS_NAME.vs.$SS_NAME.BKD
else
	echo "a folder called '$MS_NAME.vs.$SS_NAME.BKD' already exists. Output data will be stored in it."
fi

cd $MS_NAME.vs.$SS_NAME.BKD

#Step1: Breakdancer perl
printf "${col} Step1 (Breakdancer call) ${NC} \n"
perl $Breakdancer/bam2cfg.pl -g -h $MS_bam > $MS_NAME.cfg
breakdancer-max $MS_NAME.cfg > bkd_$MS_NAME.txt
python $path_to_pythonscripts/breakdancer2vcf.2.py -i bkd_$MS_NAME.txt -o bkd_$MS_NAME.vcf
if [ -n "$SS_bam" ]; then #SS bam also provided
	perl $Breakdancer/bam2cfg.pl -g -h $SS_bam > $SS_NAME.cfg
	breakdancer-max  $SS_NAME.cfg > bkd_$SS_NAME.txt
	python $path_to_pythonscripts/breakdancer2vcf.2.py -i bkd_$SS_NAME.txt -o bkd_$SS_NAME.vcf
fi
	
#Step2: Compare MS and SS
printf "${col} Step2 (Compare MS and SS)${NC} \n"
if [ -z "$SS_vcf" ]; then
	SS_vcf=bkd_$SS_NAME.vcf
fi
python $path_to_pythonscripts/CompareVCF.py -m specific \
-f1 bkd_$MS_NAME.vcf -f2 $SS_vcf \
-o 100 -p 0.5 \
-r

#Step3: Convert to DF to easy filter each field
printf "${col} Step3 (Convert to dataframe)${NC} \n"
python $path_to_pythonscripts/Vcf_to_DF.py \
-v bkd_${MS_NAME}_compare.vcf \
-o bkd_${MS_NAME}_compare.tsv
#output: a file named '${MS_NAME}_compare.tsv'

#Step4: HardFilter
#WARNING: make sure that the hard-coded filters are ok for analysis. Edit the script if necessary. 
printf "${col} Step4 (hard filters)${NC} \n"
python $path_to_pythonscripts/HardFilter_DF.py -f bkd_${MS_NAME}_compare.tsv
#output: a file named '${MS_NAME}_compare_filter.tsv'

#Step5: Edit vcf with filter
printf "${col} Step5 (Apply filters to vcf)${NC} \n"
python $path_to_pythonscripts/Filter_vcf_from_DF.py \
--mode keep \
--df bkd_${MS_NAME}_compare_filter.tsv --vcf bkd_${MS_NAME}_compare.vcf
#output: a file named 'bkd_${MS_NAME}_compare_from-df.vcf'

#Step6: Filter out variations lying in regions with abnormal high coverage (>4 x times mean coverage excluding Mt)
printf "${col} Step6 (Coverage filters)${NC} \n"
printf "Using bedtools to compute coverage for MS bam file $MS_bam: \n"
bedtools genomecov -ibam $MS_bam -bga > ${MS_NAME}_coverage
result=`grep -v "MtDNA" ${MS_NAME}_coverage | awk '{interval += $3-$2; cumcov += ($3-$2)*$4} END {print "mean cov="cumcov/interval" ; interval=" interval}'`
echo $result

#Use the result of coverage to filter out SV lying in regions that show at least four times more coverage than the mean nuclear genome 
tres=`echo $result | awk -F'[=; ]' '{print $3}'| awk '{print $0*4}'`

python $path_to_pythonscripts/Filter_coverage.py \
--vcf bkd_${MS_NAME}_compare_from-df.vcf \
--cov ${MS_NAME}_coverage \
--tres $tres \
--window 2000
#output: a file named 'bkd_${MS_NAME}_compare_from-df_cov-filter.vcf'

################################
# End of script: info to users #
################################
if [ -n "$print_stats" ]; then
	python $path_to_pythonscripts/Vcf_to_DF.py -sn -v bkd_${MS_NAME}_compare_from-df_cov-filter.vcf
fi

printf "info on the evolution of the number of variants: \n"
for file in *.vcf; do
	nb_variants=`grep -vc "#" $file`
	echo "$file: number of variants=$nb_variants"
done

printf "${col} END of script ${NC}\n"
#end of script
