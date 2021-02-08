#!/bin/bash

##Submit the script GATK_fq-to-genotypeJU170.sh
#For PSMN
#By FabFab, 2016-06-01
#
### variables SGE
### shell du job
#$ -S /bin/bash
### nom du job 
#$ -N FastQC_DRB27-WS
### file(s) d'attente(s)
#$ -q x5650lin24ibB
### charger l'environnement utilisateur pour SGE
#$ -cwd
### exporter les variables d'environnement sur tous les noeuds d'execution
#$ -V
### mails en debut et fin d'execution
#$ -m be
### where to put error file (-e) and output file (-o)
#$ -e /home/jlegrand/arabidopsis_drb27/Fab_analysis/Fab_labbook/FastQC_DRB27-WS.err
#$ -o /home/jlegrand/arabidopsis_drb27/Fab_analysis/Fab_labbook/FastQC_DRB27-WS.log

# aller dans le repertoire de travail/soumission
# important, sinon, le programme est lancé depuis ~/
cd ${SGE_O_WORKDIR}

### configurer l'environnement
source /usr/share/modules/init/bash
module use /applis/PSMN/Modules
module load Base/psmn
module load FastQC/0.11.2

### execution du programme
SCRIPTDIR=/home/jlegrand/arabidopsis_drb27/Fab_analysis/Fab_scripts
sh ${SCRIPTDIR}/FastQC_queue.sh DRB27 WS

# fin
