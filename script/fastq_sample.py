import sys, random
import HTSeq
 
"""
Prélève aléatoirement une fraction de reads dans des fichiers FastQ
Exemple pour des RNA-seq paired-end reads
 
Usage:
python getRandomReads.py 0.5 input_1.fastq.gz input_2.fastq.gz output_1.fastq.gz output_2.fastq.gz
"""
 
 
# Pourcentage de read à prélever
fraction = float( sys.argv[1] )
 
# Lire fichier FastQ read-pair 1
in1 = iter( HTSeq.FastqReader( sys.argv[2] ) )
 
# Lire fichier FastQ read-pair 2
in2 = iter( HTSeq.FastqReader( sys.argv[3] ) )
 
# Fichier de sortie FastQ read-pair 1
out1 = open( sys.argv[4], "w" )
 
# Fichier de sortie FastQ read-pair 2
out2 = open( sys.argv[5], "w" )
 
# La magie se passe ici
while True:
   read1 = next( in1 )
   read2 = next( in2 )
   if random.random() < fraction:
      read1.write_to_fastq_file( out1 )
      read2.write_to_fastq_file( out2 )
 
# Fermer les fichier créés
out1.close()
out2.close()