# Journal de Bord



## Objectif : Prise en main du script initial
Faire tourner le script existant avec un jeu de données echantionnée sur le PSMN

- Les données fastq paired end de test 
`/home/ycoude01/F19FTSEUHT1414_MOSkriR/Clean/FG06_13_C4/V300042688_L2_AE06084935-608_1.fq.gz /home/ycoude01/F19FTSEUHT1414_MOSkriR/Clean/FG06_13_C4/V300042688_L2_AE06084935-608_2.fq.gz`

- Le script fastq to gvcf  `/home/fbesnard/SCRIPTS_PSMN/GATK_fq-to-gVCF.v4_PSMN.sh`
- Exemple d'un fichier de config ` /Xnfs/rdpdb/moss/VariantAnalysis_E1/labbook/GATK1.dna/FG06_13_C4.config`

Copie des fichiers en local 
```
romuald@RDP-ZBook15-YL:~$ scp -oProxyCommand="ssh rmarin@ssh.psmn.ens-lyon.fr netcat -w1 allo-psmn %p" rmarin@allo-psmn:~/V300042688_L2_AE06084935-608_1.fq.gz stage/.
rmarin@ssh.psmn.ens-lyon.fr's password: 
V300042688_L2_AE06084935-608_1.fq.gz          100% 3255MB   5.1MB/s   10:38 
```


Connexion au psmn depuis ens 
```
ssh rmarin@allo-psmn
ssh -X cl6242comp1
```

Installation de htseq pour echantilloner les read en pair end avec un script trouvé (https://bioinfo-fr.net/rna-seq-plus-de-profondeur-ou-plus-dechantillons)

```
sudo apt install python3-pip
pip3 install HTSeq

#sur le serveur 
pip install HTSeq --user
```

lancement script d'echantillonage 
```
python3 fastq_sample.py 0.1 V300042688_L2_AE06084935-608_1.fq.gz V300042688_L2_AE06084935-608_2.fq.gz fastq_sampling1.fq fastq_sampling2.fq
```
echec : 
```
rmarin@cl6242comp1:~$ qsub sub_test.sh
Unable to run job: Job was rejected because job requests unknown queue "r815lin128ib"
Job was rejected because job requests unknown queue "monointeldeb24"
Job was rejected because job requests unknown queue "monointeldeb48"
Exiting.

```

Apres modification, le script fonctionne, le repertoire contient les fichiers 

```
rmarin@cl6242comp1:~$ tree TEST_sampling/
TEST_sampling/
├── metrics
│   ├── dedup.TEST_sampling.SRmerged.txt
│   └── TEST_sampling.RG.dedup.SRmerged.ufilter.flagstat.txt
├── TEST_sampling.RG.dedup.SRmerged.bam
├── TEST_sampling.RG.dedup.SRmerged.ufilter.bai
└── TEST_sampling.RG.dedup.SRmerged.ufilter.bam
```

## Objectif : Creation du pipeline en nextflow

Installation de nextflow 
https://www.nextflow.io/docs/latest/getstarted.html

Récuperation d'un pipeline existant avec fichier de config provenant de Laurent modolo https://gitbio.ens-lyon.fr/lmodolo

Le bug est de transformer les differentes scripts dans un pipeline nextflow : 
Modification du pipeline afin de faire comme le script GATK_fq-to-gVCF.v4_PSMN.sh


test avec read /home/ycoude01/F19FTSEUHT1414_MOSkriR/Clean/F*

zcat /home/ycoude01/F19FTSEUHT1414_MOSkriR/Clean/FG06_22_B2/V300042688_L3_AE47136387-610_1.fq.gz | head -80000 > ../ref/SAMPLE_V300042688_L3_AE47136387-610_1.fq.gz


Creation d'echantionn 
zcat file.gz | head -10000 > sample_file 

a recopier la doc https://nf-co.re/rnaseq/3.0/usage

Test avec les vrais echantillons :
```rmarin@cl6226comp1:~/Pipeline_variant_RDP$ ls ../ref/reads/
V300042688_L2_AE06084935-608_1.fq.gz  V300042688_L2_AE97758923-605_1.fq.gz  V300042688_L4_AE59776336-607_1.fq.gz  V300042688_L2_AE06084935-608_2.fq.gz  V300042688_L2_AE97758923-605_2.fq.gz  V300042688_L4_AE59776336-607_2.fq.gz
V300042688_L2_AE06354351-606_1.fq.gz  V300042688_L3_AE59776336-607_1.fq.gz
V300042688_L2_AE06354351-606_2.fq.gz  V300042688_L3_AE59776336-607_2.fq.gz
```
En sachant que l strating strain est SAMPLE_V300042688_L2_AE97758923-605_2 

## Determination du type recalibration pour GATK

Ce test est fait sur les fichiers de SNP globeau apres la separation en snp/indel avec les data ci dessus 

Lien de la methode avec BQSR https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/ 

| |Sans filtre | Normal | BQSR | VQSR |
| :---------------|:---------------: |:---------------:|:---------------:| -----:|
SNP | 913977 |812421 | 814712 | 
INDEL |59643| 59643 | 57523 |

VariantFiltration de GATK https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration
Certain sont "pass" et d'autre variant sont annotés avec le nom du filtre et n'ont pas été supprimés Pour cela ajour d'une etape en plus dans le process Extract_SNPIndel_Filtration
`gatk SelectVariants --exclude-filtered -V pre_filtered_indels.vcf -O indels.vcf`


## Determination du type de Bam a donner pour les programmes de VS

Certain bam sont filtés afin d'enlever tous les reads non mappé dans le process FilteringandIndexing_Bam grace a la commande samtools : 

`samtools view -hu -F4 ${bam_RD_MD} | samtools view -hu -F256 - | samtools view -hb -f3 - > ${pair_id}_ufilter.bam `
### Test du type de données d'entrées dans Breakdancer : 

| Sample | Add_ReadGroup_and_MarkDuplicates_bam | Filtering_and_indexing_bam |
| :--------------- |:---------------:| -----:|
 V300042688_L2_AE06084935-608 | 5606 | 1543 |
V300042688_L2_AE06354351-606 | 5415 | 1357
V300042688_L2_AE97758923-605 (SS) | 5250 | 1391
V300042688_L3_AE59776336-607 | 3011 | 601 | 
V300042688_L4_AE59776336-607 | 3087 | 637

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3661775/

Filtrer red non pairer pour BDancer
Pindel si pas sa paire pas grave parce que breakpoint donc on garde tout 


Fonctionnement de Breakdancer qui fait qu'il faut mieux les filtres ou pas ? 

### Test du type de données d'entrées dans Pindel : 

|Sample | Add_ReadGroup_and_MarkDuplicates_bam and insert lenght from breakdancer | Filtering_and_indexing_bam |
| :--------------- |:---------------:| -----:|
V300042688_L2_AE06084935-608 | |676709
V300042688_L2_AE06354351-606 | | 655720
V300042688_L2_AE97758923-605 (SS)  | | 63801
V300042688_L3_AE59776336-607 | 536783 | 470922
V300042688_L4_AE59776336-607 | 539882 | 473690

A prendre en compte les données obtenu avec le process Filtering_and_indexing_bam avait une valeur de taille d'insert par defaut de 500. Le second test a été fait en corrigant cela grace au fichier de config generé par breakdancer en recuperant la taille d'insert moyen pour plus de précision

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2781750/
Pindel utiliser les breakpoint comme point d'ancrage pour detecter les indel de plus ou moins grande taille donc a priori mieux vaut enever les unampped et les unpair



Appliquer VQSR https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-


https://www.nextflow.io/blog/2016/error-recovery-and-automatic-resources-management.html
Trop de ram utilisé -> error 140 
chnegemebt de queues dans la config