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
#SSAD_filter script (Sample Specific Allele Depth_filter)
#Input= vcf with at least 2 "samples". One sample is considered as the reference, another as the sample of interest. Must have been processed with SSAD_write script.
#Ouput= filter vcf with a new INFO field 'SSAD' indicating various metrics of 

#All vcf files are tab-delimited files ordered as follows:
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1	sample2
#1	2	3	4	5	6	7	8	9	10	11	etc... (Nbering of columns)

#usage: bash SSAD_filter.sh -s specific_sample_name -r reference_sample_name -o vcf.vcf (optional: output_name if option -o is given)

##############################
#  Read & check User inputs  #
##############################
printf "SSAD_filter.sh \n"
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

echo "the vcf file is $1"
echo "the sample name is $s"
echo "the reference name is $r"

#Manage output file naming
#Expected default behavior: the program will append ssad-F just before the .vcf extension.
#	if ssad-W is present, it will be changed to ssad-F
#	if several ssad-W are present, only the lastone will be changed to ssad-F
vcfbase=$(basename $1 .vcf)

src="ssad-W" 
dest="ssad-F"

if [ -z $output ]; then
#no output name given by user -> impose the following name:
	if [[ $vcfbase = *"$src"* ]]; then
		prefix=${vcfbase%"$src"*}                  # Extract content before the last instance
		suffix=${vcfbase#"$prefix"}                # Extract content *after* our prefix
  		outputname=${prefix}${suffix/"$src"/"$dest"}.vcf  # Append unmodified prefix w/ suffix w/ replacement
	else
  		outputname=$vcfbase.ssad-F.vcf
	fi
#if option -o is activated but no outputname given
elif [[ -z $2 ]]; then
	echo "Warning: if you want a specific name for output (option -o), please provide the name as second argument (after the name of the input vcf)"
	echo "Output name will have the default nomenclature: appending 'ssad-F.vcf' and replacing the last occcurence of 'ssad-W' if present"
	if [[ $vcfbase = *"$src"* ]]; then
		prefix=${vcfbase%"$src"*}                  # Extract content before the last instance
		suffix=${vcfbase#"$prefix"}                # Extract content *after* our prefix
  		outputname=${prefix}${suffix/"$src"/"$dest"}.vcf  # Append unmodified prefix w/ suffix w/ replacement
	else
  		outputname=$vcfbase.ssad-F.vcf
	fi
#if an output name has been given, just check whether it comes with or without the extension .vcf (and adapt to provide the right output ending with only one '.vcf')
elif [[ $2 == *".vcf" ]]; then 
	outputname=$2
else
	outputname=$2.vcf
fi

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
#save & append header
grep "^##" $1 > $outputname
cmd=$0
inputfile=$1
comment(){
echo "##SSAD_filter=$version"  >> $1
echo "##SSAD_filter_Cmd=$cmd -s $s -r $r $inputfile " >> $1
}

comment $outputname
#add the last line of comment (ie columns/field names of the vcf):
grep "^#CHROM" $1 >> $outputname

#Reminder: SSAD fields
#SSAD[1]: Sample name
#SSAD[2]: Ref name
#SSAD[3]: Depth of Sample majo allele (Smajo)
#SSAD[4]: Smajo ratio
#SSAD[5]: Smajo depth in Reference
#SSAD[6]: Smajo ratio in Reference
#SSAD[7]: Nb other sample with reads for Smajo
#SSAD[8]: Total depth Smajo in other samples

#Start parsing vcf
#Remove comment lines
grep -v "#" $1 | \
awk -v Sname="$s" -v Rname="$r" \
'
#block of awk-functions
function get_SSAD_idx(array, size){
	idx=0;
	for (i=1; i<=size; i++)
		if (substr(array[i],1,4) == "SSAD") idx=i;
	return idx
}
BEGIN { OFS="\t" }
## awk-block 1: get the INFO field 
{
N=split($8, INFO, ";");
SSAD_idx=get_SSAD_idx(INFO, N);
split(INFO[SSAD_idx], SSAD, "|");
SSAD[1]=substr(SSAD[1],6); # because SSAD[1] starts with this 5 characters: "SSAD="
if (SSAD[1] != Sname || SSAD[2] != Rname ){
	next 
	}
else {
	if (SSAD[3] >= 2 && SSAD[4] >=0.9 && SSAD[5] <= 0 && SSAD[6] <= 0 && SSAD[7] <= 2 && SSAD[8] <= 2)
		print $0
	}
}
'>> $outputname
