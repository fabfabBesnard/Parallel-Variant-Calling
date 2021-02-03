# Journal de Bord



## Objectif : 
Faire tourner le script avec un jeu de données echantionnée sur le PSMN

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