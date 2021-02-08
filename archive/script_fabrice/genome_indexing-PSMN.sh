#!/bin/sh
## script: 'genome_indexing-PSMN.sh'
## Fabrice Besnard, 2014 v1.0
# last edit: 2016-06-01
#version for PSMN 

# required:
# bwa index
# samtools 
# Picard Version: 1.110+. PSMN version (June 2016): 2.3.0
#For picard: the env variable $PICARD is created upon loading the module

#Warning: use it from the directory where you want have the final files, because the variable $filepath is only local 

echo "Starting program genome_indexing.sh with fasta file $1"
filename=$(basename "$1")
filepath=${1%/*}/
filenamewithoutextension=${filename%.*}

#1. Create BWA index of the assembly
#Since the genome is small, we use the '-a is' option (linear-time algorithm)
printf "\033[1;32m **STEP1** Indexing $filename with bwa -is algorithm \033[0m \n"
bwa index -a is $1
  
#2.Generate fasta file index
printf "\033[1;32m **STEP2** Indexing $filename with samtools \033[0m \n"
samtools faidx $1
  
#3.Generate a dictionary with Picard
printf "\033[1;32m **STEP3** Create a sequence dictionary for $filename \033[0m \n"
java -jar $PICARD CreateSequenceDictionary \
REFERENCE=$1 \
OUTPUT=${filenamewithoutextension}.dict
