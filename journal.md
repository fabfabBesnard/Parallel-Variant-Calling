# Journal de Bord



## Objectif : 
Faire tourner le script avec un jeu de données echantionnée sur le PSMN

- Les données fastq paired end de test 
`/home/ycoude01/F19FTSEUHT1414_MOSkriR/Clean/FG06_13_C4/V300042688_L2_AE06084935-608_1.fq.gz /home/ycoude01/F19FTSEUHT1414_MOSkriR/Clean/FG06_13_C4/V300042688_L2_AE06084935-608_2.fq.gz`

- Le script fastq to gvcf  `/home/fbesnard/SCRIPTS_PSMN/GATK_fq-to-gVCF.v4_PSMN.sh`

Copie des fichiers en local 
```
romuald@RDP-ZBook15-YL:~$ scp -oProxyCommand="ssh rmarin@ssh.psmn.ens-lyon.fr netcat -w1 allo-psmn %p" rmarin@allo-psmn:~/V300042688_L2_AE06084935-608_1.fq.gz stage/.
rmarin@ssh.psmn.ens-lyon.fr's password: 
V300042688_L2_AE06084935-608_1.fq.gz          100% 3255MB   5.1MB/s   10:38 
```

Installation de htseq pour echantilloner les read en pair end avec un script trouvé (https://bioinfo-fr.net/rna-seq-plus-de-profondeur-ou-plus-dechantillons)

sudo apt install python3-pip
pip3 install HTSeq