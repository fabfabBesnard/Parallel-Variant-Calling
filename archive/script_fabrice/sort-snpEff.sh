#!/bin/bash
##This is sort-snpEff.sh
##Created by F. BESNARD, june 2016.
#Run on PSMN

#What it does
#sort and order in a human readble format the output of a vcf file annotated by snpEff

#Usage
# bash sort-snpEff.sh myvcf.snpeff.vcf myouput.txt

printf "\e[1;32m START OF SCRIPT 'sort-snpEff.sh' \e[0m \n"

printf "##### SORTING annotated variants  for $1 ##### \n" > $2.txt
echo "#`date`" >> $2.txt

#Sort candidate genes based on the functional annotation
##HIGH IMPACT mutations
printf "\n#HIGH IMPACT mutations \n" >> $2.txt
cat $1 | \
/applis/PSMN/generic/SnpEff/4.2/scripts/vcfEffOnePerLine.pl | \
java -jar /applis/PSMN/generic/SnpEff/4.2/SnpSift.jar filter "ANN[0].IMPACT has 'HIGH'" | \
java -jar /applis/PSMN/generic/SnpEff/4.2/SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].GENE" "ANN[*].FEATUREID" "ANN[*].EFFECT" "LOF[*].GENE">> $2.txt

NbHIGH=$(cat $1 | \
	/applis/PSMN/generic/SnpEff/4.2/scripts/vcfEffOnePerLine.pl | \
	java -jar /applis/PSMN/generic/SnpEff/4.2/SnpSift.jar filter "ANN[0].IMPACT has 'HIGH'" | grep -vc "#")

##MODERATE IMPACT mutations
printf "\n#MODERATE IMPACT mutations \n">> $2.txt
cat $1 | \
/applis/PSMN/generic/SnpEff/4.2/scripts/vcfEffOnePerLine.pl | \
java -jar /applis/PSMN/generic/SnpEff/4.2/SnpSift.jar filter "ANN[0].IMPACT has 'MODERATE'" | \
java -jar /applis/PSMN/generic/SnpEff/4.2/SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].GENE" "ANN[*].FEATUREID" "ANN[*].EFFECT" "LOF[*].GENE">> $2.txt

NbMODERATE=$(cat $1 | \
/applis/PSMN/generic/SnpEff/4.2/scripts/vcfEffOnePerLine.pl | \
java -jar /applis/PSMN/generic/SnpEff/4.2/SnpSift.jar filter "ANN[0].IMPACT has 'MODERATE'" | grep -vc "#")

##LOW IMPACT mutations
printf "\n#LOW IMPACT mutations \n">> $2.txt
cat $1 | \
/applis/PSMN/generic/SnpEff/4.2/scripts/vcfEffOnePerLine.pl | \
java -jar /applis/PSMN/generic/SnpEff/4.2/SnpSift.jar filter "ANN[0].IMPACT has 'LOW'" | \
java -jar /applis/PSMN/generic/SnpEff/4.2/SnpSift.jar extractFields - CHROM POS REF ALT "ANN[*].GENE" "ANN[*].FEATUREID" "ANN[*].EFFECT" "LOF[*].GENE">> $2.txt

NbLOW=$(cat $1 | \
/applis/PSMN/generic/SnpEff/4.2/scripts/vcfEffOnePerLine.pl | \
java -jar /applis/PSMN/generic/SnpEff/4.2/SnpSift.jar filter "ANN[0].IMPACT has 'LOW'" | grep -vc "#" )

##Counts of variants per functional class
echo "##Counts of variants per functional class"
echo "##Counts of variants per functional class">> $2.txt

echo "Number of High-impact mutations: $NbHIGH"
echo "Number of High-impact mutations: $NbHIGH">> $2.txt
echo "Number of Moderate-impact mutations: $NbMODERATE"
echo "Number of Moderate-impact mutations: $NbMODERATE">> $2.txt
echo "Number of Low-impact mutations: $NbLOW"
echo "Number of Low-impact mutations: $NbLOW">> $2.txt

echo "##end of file" >> $2.txt

printf "\e[1;32m END OF SCRIPT 'find_gene_PSMN.sh' \e[0m \n"
