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

installation de nextflow 
https://www.nextflow.io/docs/latest/getstarted.html

Récuperation d'un pipeline existant avec fichier de config 

Modicication du pipeline afin de faire comme le script  GATK_fq-to-gVCF.v4_PSMN.sh

Possibilité de faire un index

ou d'utiliser un  index deja crée

`rmarin@cl6242comp1:~/Pipeline_variant_RDP$ ./nextflow run script/GATK_to_gVCF.nf -c script/GATK_to_gVCF.config --reads "/home/rmarin/V300042688_L2_AE06084935-608_{1,2}.fq.gz" --genomeindex "../ref/index/
Physcomitrella_patens_Phypa_V3_dna_rm_toplevel.fa.*" -profile psmn -resume`


ADD read group 
https://gatk.broadinstitute.org/hc/en-us/articles/360036713171-AddOrReplaceReadGroups-Picard-

markduplicate 
https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-