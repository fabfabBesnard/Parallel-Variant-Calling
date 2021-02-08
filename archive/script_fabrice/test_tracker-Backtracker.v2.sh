#!/bin/bash


CONFIGFILE=$1
source $CONFIGFILE

#Inform users of the parameters:
printf "Starting program: `date` \n"
printf "${col} **Running test -> testing** ${NC} \n"
printf "${col} T-DNA selected: ${NC} $TDNA \n"
printf "${col} Reference genome selected: ${NC} $REFERENCE \n"
printf "${col} Reads_F selected: ${NC} $READS_F \n"
printf "${col} Reads_R selected: ${NC} $READS_R  \n"

#Step 14. Backtrack each region of the plant genome to the T-DNA
echo "Step 14"
idx=1
cat ${OUTPUT_NAME}_regions.txt | while read line
	do 
		bash ${tracker_path}/TDNA_Backtracker.sh \
				-f $CONFIGFILE \
				-R $line \
				-o TDNA.$idx \
				-d $OUTPUT_NAME.vsTDNA.Sfwd.bam \
				-i $OUTPUT_NAME.vsTDNA.Srev.bam \
				$OUTPUT_NAME.vsTDNA_coverage \
				$OUTPUT_NAME.Mates.WG.bam
		((idx++))
	done

#Step 15. Store and organize files
#Move all bam files and coverage files to data subfolder
echo "Step 15"
mv *.bam ./data/
mv *.bai ./data/
mv *_coverage ./data/

printf "${col} ** END of test 'tracker-Backtracker' **${NC}\n"
