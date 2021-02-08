#!/bin/bash

mygene=$1

#Transform gene ID into the ID used in Augustus file
ID=$(echo $mygene | sed 's/nOt.2.0.1.//' | sed 's/g0/g/')

#Retrieve the sequence
SEQUENCE=$(awk '/start gene '$ID'/,/end gene '$ID'/' /home/fbesnard/Reference_genomes/Otipulae/Annotation/nOt.2.0.1.aug.gff | \
grep -A1000 "# protein sequence" | grep -B1000 "]$" | \
sed 's/# protein sequence = \[//' | sed 's/# //' | sed 's/\]//' | perl -pe 's/\n//g')

echo ">$mygene"
echo $SEQUENCE
