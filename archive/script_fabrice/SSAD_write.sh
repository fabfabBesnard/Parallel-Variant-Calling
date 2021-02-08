#!/bin/bash

###########################################################################
#### Copyright (C) 2020 - INRAe (Fabrice Besnard, RDP)
#### This script is a free software: you can redistribute it
#### and/or modify it under the terms of the GNU General Public
#### License (GNU GPL) as published by the Free Software Foundation, either
#### version 3 of the License, or (at your option) any later version.
#### Distributed without any warranty.
###########################################################################
# last edit: 2020-08-31
#Version v1
version=$(echo "version v1")


#######################
#       READ ME       #
#######################
#SSAD_write script (Sample Specific Allele Depth_write)
#Input= vcf with at least 2 "samples". One sample is considered as the reference, another as the sample of interest.
#Ouput= the vcf with a new INFO field 'SSAD' indicating various metrics of 

#All vcf files are tab-delimited files ordered as follows:
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1	sample2
#1	2	3	4	5	6	7	8	9	10	11	etc... (Nbering of columns)

#usage: bash SSAD_write.sh -s sample_name -r ref_name -o vcf.vcf (optional with -o option: give outputname

##############################
#  Read & check User inputs  #
##############################
printf "SSAD_write.sh \n"
#Recall usage if bad inputs:
usage() { echo "Usage: $0 [-s <string>] [-r <string>] [-o <nothing>] file.vcf string";
	echo "	-s sample name -as given in the vcf file- for which specific variants are wanted";
	echo "	-r sample name of the reference 'starting strain' -as given in the vcf file-, which was mutagenized to yield other strains";
	echo "	-o specify whether the name of the output file is given as the second positional argument"
	exit 1; }
#Manage inputs as options:
while getopts ":s:r:o" arg; do
    case $arg in
        s)
            s=${OPTARG}
            ;;
	r)
            r=${OPTARG}
            ;;
	o)
            output=user_defined
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

#Basic checks for inputs
if [ -z $s ] || [ -z $r ] ; then
	echo "ERROR: define a sample name and reference name";
	usage 
fi
if [[ $s == *".vcf" ]] || [[ $r == *".vcf" ]]; then 
	echo "ERROR : at least one sample name is missing";
	echo "note: alternatively, sample name ends with '.vcf', which cannot work with this program";
	usage 
fi
if [ -z $1 ] ; then
	echo "ERROR: provide a vcf file" ;
	usage 
fi

#Manage output file naming
vcfbase=$(basename $1 .vcf)
if [ -z $output ]; then
#no output name given by user -> impose the following name:
	outputname=$vcfbase.ssad-W.vcf
#if an output name has been given, just check whether it comes with or without the extension .vcf (and adapt to provide the right output ending with only one '.vcf')
elif [[ -z $2 ]]; then 
	echo "Warning: if you want a specific name for output (option -o), please provide the name as second argument (after the name of the input vcf)"
	echo "Output name will have the default nomenclature = appending 'ssad-W.vcf'"
	outputname=$vcfbase.ssad-W.vcf ;
elif [[ $2 == *".vcf" ]]; then 
	outputname=$2
else
	outputname=$2.vcf
fi

#echo "the vcf file is $1"
echo "the sample name is $s"
echo "the reference name is $r"

#Advanced checks for vcf files
vcf_field=(\#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT)

counter=0
grep "^#CHROM" $1 | awk '{split($0, fields)}{for (i=1; i<=9; i++) print fields[i]}' | \
while read line; do 
	#echo $line 
	if [ "$line" != "${vcf_field[$counter]}" ]; then
		echo "the input vcf file has some formating issues."
		echo "The following field was not found: ${vcf_field[$counter]} "
	fi
	counter=$((counter+1))
done

nb_samples=$(grep "^#CHROM" $1 | awk '{print NF-9}')
if [ -z $nb_samples ]; then
	echo "Error: no sample fields were found. Check format of the input vcf";
	usage
else
	echo "Input vcf contains $nb_samples samples"
fi

get_index(){ #provided that the array sample is a gobal shell variable
value=$1
for i in "${!all_samples[@]}"; do
   if [[ "${all_samples[$i]}" = "${value}" ]]; then
       echo "${i}";
   fi
done
}

declare -A all_samples
idx=1
while read line ; do 
		all_samples[$idx]=$line
		idx=$((idx + 1))
	done < <(grep "^#CHROM" $1 | awk '{split($0, fields)}{for (i=10; i<=NF; i++) print fields[i]}')
#for i in "${all_samples[@]}"; do echo "$i"; done

#Get the idx of the sample/ref given in input
mysample=$(get_index $s)
#echo $mysample
ref=$(get_index $r)
#echo $ref

if [ -z $mysample ] || [ -z $ref ] ; then
	echo "ERROR: sample and reference sample names were not found in the vcf file";
	usage 
fi

#######################
#     Script Body     #
#######################
##1. Header: check presence of anterior SSAD, save & append with new information
cmd=$0
inputfile=$1
comment1(){
echo "##SSAD_write=$version;"  >> $1
echo "##SSAD_write_Cmd=$cmd -s $s -r $r $inputfile " >> $1
echo "##run `date '+%Y-%m-%d'`" >> $1
echo "##INFO=<ID=SSAD,Number=6,Type=String,Description=\"Comparison of allele depth between a specific sample versus reference sample and other samples: 'SSAD= Name of the specific sample | Name of the Reference Sample | Depth of major allele (Smajo) in the specific sample | Ratio of the major allele (Smajo) in the specific sample | Depth of Smajo allele in the reference sample | Ratio of Smajo allele in the reference sample | nb of other samples having reads for Smajo | total depth of Smajo in other samples\">" >> $1
}
comment2(){
echo "##SSAD_write=Overwrite previous data" >> $1
echo "##SSAD_write=$version;"  >> $1
echo "##SSAD_write_Cmd=$cmd -s $s -r $r $inputfile " >> $1
echo "##run `date '+%Y-%m-%d'`" >> $1
}

#check the presence of anterior ssad ($anterior)
anterior=$(grep -c "^##SSAD_write" $1)
#start writing the header of the outputfile
grep "^##" $1 > $outputname
#append comment depending on the presence/absence of anterior ssad
if [ $anterior -eq 0 ]; then
	comment1 $outputname
else
	echo "warning: the vcf already contains ssad computations. They will be overwritten."
	comment2 $outputname
fi

#add the last line of comment (ie columns/field names of the vcf):
grep "^#CHROM" $1 >> $outputname

##2. Parsing the vcf data (Awk code-block)
grep -v "#" $1 | \
awk -v Sidx="$mysample" -v Ridx="$ref" -v Sname="$s" -v Rname="$r" \
'
#block of awk-functions
function get_SSAD_idx(array, size){
	idx=0;
	for (i=1; i<=size; i++)
		if (substr(array[i],1,4) == "SSAD") idx=i;
	return idx
}
function get_idx_max(my_array, size){
	max=0;
	for (i=1; i<=size; i++)
	if (max<my_array[i]) max=my_array[i];
	idx=0; 
	for (i=1; i<=size; i++)
		if (my_array[i] == max) idx=i;
	return idx
}
function get_idx_second(my_array, size, max){
	val=0;
	for (i=1; i<=size; i++)
	if (val<my_array[i] && my_array[i] != max) val=my_array[i];
	#return val;
	idx=0; 
	for (i=1; i<=size; i++)
		if (my_array[i] == val) idx=i;
	return idx
}
function allele_ratio(num1, num2){
	if (num2 == 0) return 1
	return (1-num2/(num1+num2))
}
function array_idx_ratio(array, size,idx){
#return the ratio of the value of array[idx] over the sum of all the array
	sum=0
	for (i=1; i<=size; i++) sum=sum+array[i];
	if (sum == 0) return "NA"
	else return array[idx]/sum;
	
}
function nb_strain_Sallele_read(string){
	n=split(string, element);
	counter=0;
	for (i=1; i<=n; i++) if (element[i] != 0) counter++;
	return counter	
}
function total_AD_other_strain(string){
	n=split(string, element);
	sum=0
	for (i=1; i<=n; i++) sum=sum+element[i];
	return sum
}
BEGIN { OFS="\t" }
## awk-block 1: start parsing vcf 
{split($9, FORMAT, ":") ; split($(Sidx+9), mysample, ":"); split($(Ridx+9), reference, ":")} 
## awk-block 2:compute stats for the different samples
{if (FORMAT[2] != "AD") next; #case where format field is empty 
S_nad=split(mysample[2], ad_mysample, ","); 
#mysample stats
Sallele_idx=get_idx_max(ad_mysample, S_nad);
Sallele_depth=ad_mysample[Sallele_idx];
Sallele2_idx=get_idx_second(ad_mysample, S_nad, ad_mysample[Sallele_idx]);
Sratio=allele_ratio(ad_mysample[Sallele_idx], ad_mysample[Sallele2_idx]);
#ref stats
R_nad=split(reference[2], ad_ref, ",");
Rallele_idx=get_idx_max(ad_ref, R_nad);	
Rallele_depth=ad_ref[Sallele_idx] ; #In the reference sample, the depth of the allele which is max in the specific sample
Rratio=array_idx_ratio(ad_ref, R_nad, Sallele_idx); #In the reference sample, computes the ratio of the depth of the allele for which specific sample is max over all other allele depth
#write
other_ad=""; #set to null
for (i=10; i<=NF; i++) {		
	if (i != (Sidx+9) && i != (Ridx+9)) {split($i, sample, ":"); split(sample[2], ad_sample, ","); other_ad=other_ad FS ad_sample[Sallele_idx]}
	}
nb_other=nb_strain_Sallele_read(other_ad)
AD_other=total_AD_other_strain(other_ad)
}
## awk-block 3: print line by adding INFO data (n°8)
{N=split($8, INFO, ";");
SSAD_idx=get_SSAD_idx(INFO, N);
if (SSAD_idx == 0){
	for (i=1; i<=N; i++){$8=$8";"INFO[i]};
	$8=$8";SSAD="Sname"|"Rname"|"Sallele_depth"|"Sratio"|"Rallele_depth"|"Rratio"|"nb_other"|"AD_other; 
	print $0}
else {
	delete INFO[SSAD_idx];
	for (i=1; i<=N; i++){$8=$8";"INFO[i]};
	$8=$8"SSAD="Sname"|"Rname"|"Sallele_depth"|"Sratio"|"Rallele_depth"|"Rratio"|"nb_other"|"AD_other; 
	print $0}
}'>> $outputname
