#!/bin/bash
##This is MapVarFreq.sh
##Created by F. BESNARD 2016/06/28
# Adapted for PSMN 
#last edit: 2016-06-29

#What is does:
#From a .vcf file (it can contain several??? samples), it computes for each reported variant position the frequency of reads for the ALT allele over the sum of (REF + ALT) read counts.
#Results are stored in a file named $argument2.freq.txt
#Then plots of frequency and coverage are geneated with a Rcript and stored in a new folder named "plots".

#Requirements:
#This scripts uses the R script Plot-fqcy_PSMN.R
#a file containing the scaf names and their sizes -> get it from vcf ??

#Usage:
#sh map_mutant_PSMN.sh myvcf.vcf NameoftheSample Outfile-prefix
#Takes as 1st argument .vcf file containing the genotype for JU170 snps
#2nd argument is the sample of the vcf file to study. #Important Note: samples in vcf must be without space
#3rd argument is the outfile prefix (can be the same as sample name)

printf "\e[1;32m STARTING 'MapVarFreq_PSMN.sh' \e[0m \n"
printf "\e[0;32m computing allele frequencies from file $1 \e[0m \n"
printf "\e[0;32m Output file named $3.freq.txt will be generated \e[0m \n"
printf "\e[0;32m Frequency and coverage plots will be generated in a new folder named '$3_MapVarFreq-plots' \e[0m \n"

#Get the sample and its position
#All vcf files are tab-delimited files ordered as follows:
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1  sample2
SampleLine=$(grep "#CHROM" $1 | cut -f10-)
NbSample=`echo $SampleLine | grep -c "\s"`
NbSample=$(($NbSample+1))

#Store samples into an associative array: key=SampleName ; Values=Order of sample in the vcf
declare -A SAMPLES
for i in `seq 1 $NbSample`
do 	
	SampleName=$(echo $SampleLine | cut -d " " -f$i,$i )
	SAMPLES[$SampleName]=$i
done

#Get Chr info in a temp file ending by '.ContigID.txt'
vcfbase=$(basename $1 .vcf)
grep "##contig=<ID=" $1 > $vcfbase.ContigID.txt
sed -i -e "s/##contig=<ID=//" -e "s/,length=[0-9]*>//" $vcfbase.ContigID.txt

#Prepare a subfolder to store
mkdir $3_MapVarFreq-plots

#Compute frequencies
printf "temporary file generated \n"
grep -v "#" $1 > filter.vcf #Remove the header. 
#Isolate the infos corresponding to sample given in argument
SampleRank=${SAMPLES[$2]}
echo "Collecting vcf info for $2 which is the n°$SampleRank in $1"
#The corresponding col in the tab-delimited vfc is:
SampleCol=$(($SampleRank+9))
printf "Computing frequencies in progress. Progressing through: \n"
#Trick: use awk rather than a grep here in case contigs are just numbers
# Extract following fields: CHROM POS REF ALT sampleFORMAT
cat $vcfbase.ContigID.txt | while read line 
				do 
				echo $line  
				awk -v line=$line '$1==line {print $0}' filter.vcf |   
				cut -f 1,2,4,5,$SampleCol |
				awk '{split($4,ALT,",")}{n=split($5, FORMAT, ":")}{split(FORMAT[1], allele, "/")}{allele1=allele[1];allele2=allele[2]}{split(FORMAT[2],AD,",")} #split all fields that need to be analyzed
					{if (n==1) ##Case where no GT has been computed (Format field is empty)
						next;
							if (allele[1]==0 || allele1==allele2) ##Exclude situtations with Het of two non-ref alleles
								if (allele[1]==".")   ## Case of NA position FORMAT has only two fields, GT and AD: ./.:0,0
									print $1,$2,$3,$4,FORMAT[1], "NA", "NA", "NA";
								else if (AD[1]+AD[allele2+1]>0)  ##If statement avoids division by 0. AD[1] is always the depth for the Ref allele (GT=0 in vcf). The alt allele (either HM or Het with Ref) has an index in the AD array which is GT+1
									print $1,$2,$3,$4,FORMAT[1],FORMAT[3],AD[1]+AD[allele2+1],AD[allele2+1]/(AD[1]+AD[allele2+1]);
									else print $1,$2,$3,$4,FORMAT[1],FORMAT[3],AD[1]+AD[allele2+1],0}'>> ./$3_MapVarFreq-plots/$3.freq.txt
				done
#Resulting files should contain the following fields:
#CHROM POS REF ALT GT DP AD(ref)+AD(sample) AD(sample)/(AD(ref)+AD(sample))

printf "All frequencies computed. Cleaning temp. files \n"
#Clean temp filesWS_MapVarFreq-plots/
rm $vcfbase.ContigID.txt
rm filter.vcf

printf "Generating plots \n"
cd $3_MapVarFreq-plots
Rscript /home/fbesnard/SCRIPTS_PSMN/MapVarFreq_plot.R $3.freq.txt $3


printf "\e[1;32m END OF SCRIPT 'map_mutant_PSMN.sh' \e[0m \n"
