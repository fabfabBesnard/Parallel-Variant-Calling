#!/bin/bash

##Submit the script GATK_fq-to-gVCF.v4_PSMN.sh for sample FG06_13_C4 (moss project)
#For PSMN
#By FabFab, 2020-02-12
#Description: FG06_13_C4/non-masked/BQSR/gVCF

### variables SGE
### shell du job
#$ -S /bin/bash
### nom du job 
#$ -N testsampling
### file d'attente
#$ -q h48-CLG6226Rdeb192
### charger l'environnement utilisateur pour SGE
#$ -cwd
### exporter les variables d'environnement sur tous les noeuds d'execution
#$ -V
### mails en debut et fin d'execution
#$ -m be
### where to put error file (-e) and output file (-o)
#$ -e /home/rmarin/testsampling.err
#$ -o /home/rmarin/testsampling.out

# aller dans le repertoire de travail/soumission
# important, sinon, le programme est lancé depuis ~/
## cd "${SGE_O_WORKDIR}" || { echo "cannot cd to ${SGE_O_WORKDIR}"; exit 1; }

### configurer l'environnement
module use /applis/PSMN/Modules
module load Bwakit/0.7.15
#SAMtools/1.1 is loaded automatically with previous commands
module load picard/2.3.0
module load GATK/3.6 #Java version of picard/2.3.0 and of GATK/3.6 are compatible (Java1.8) 

### execution du programme
SCRIPTDIR=/home/rmarin
ProjectDIR=/home/rmarin
#Launch the sub file from ProjectDir/GATK.1
bash ${SCRIPTDIR}/GATK_fq-to-gVCF.v4_PSMN.sh ${ProjectDIR}/my.config
# end
