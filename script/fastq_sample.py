import sys, random
import HTSeq
from itertools import takewhile
 
"""
Prélève aléatoirement une fraction de reads dans des fichiers FastQ
Exemple pour des RNA-seq paired-end reads
 
Usage:
python fast_sample.py 0.5 input_1.fastq.gz input_2.fastq.gz output_1.fastq output_2.fastq check

"""

# Compteur du nombre de read ajouté
nb_out = 0

# Initilisation du dict paired
# https://htseq.readthedocs.io/en/master/sequences.html#fastareader
print("Lecture du premier fichier")
sequences1 = dict((s.name, s) for s in HTSeq.FastqReader(sys.argv[2]))
print("lecture du second fichier")
sequences2 = dict((s.name, s) for s in HTSeq.FastqReader(sys.argv[3]))

# Verification de le taille des fichiers
try :
   if sys.argv[6] == "check" :
      if len(sequences1) != len(sequences2):
         print("Les fichiers fournis n'ont pas la même taille")
         sys.exit()
      else :
         print("Les fichiers ont la même taille")
except IndexError :
   pass

# Nombre de read à prélever
fraction = int(float( sys.argv[1] )*len(sequences1))

# Choix aléatoire de reads
# https://www.geeksforgeeks.org/python-random-sample-function/?ref=lbp
choice = random.sample(list(sequences1), fraction)

# Creation des fichiers de sortie
out1 = open( sys.argv[4], "w" )
out2 = open( sys.argv[5], "w" )

for i in choice:

   # On ajoute les reads choisis aux deux fichiers
   try :
      sequences2[i.split("/")[0]+"/2"].write_to_fastq_file( out2 )
      nb_out = nb_out + 1
      sequences1[i].write_to_fastq_file( out1 )

   # Si un read n'est pas présent dans sequences 2, on ne l'ajoute pas
   except KeyError :
      print("Le read complémentaire de %s n'a pas été trouvé dans le fichier 2"%(i.split("/")[0]))

print("Nombre de reads choisis aléatoirement : ", nb_out)

out1.close()
out2.close()