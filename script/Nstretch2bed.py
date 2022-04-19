#!/usr/bin/python3.6
import sys
import csv

#import the SeqIO module from Biopython
from Bio import SeqIO

# File author(s):
#           Malik Talbi <malik.talbi.2b@gmail.com>
#
#       File contributor(s):
#           Romuald Marin
#           Malik Talbi <malik.talbi.2b@gmail.com>
#           Fabrice Besnard <fabrice.besnard@ens-lyon.fr>
#
#       File maintainer(s) and contact :
#           Fabrice Besnard <fabrice.besnard@ens-lyon.fr>
#
#       RDP Lab, Signal Team, Lyon - INRAe
# ------------------------------------------------------------------------------

bed_out=open(sys.argv[2], mode="w")

dict_pos={}

with open(sys.argv[1], mode="r") as fasta_handle:
	for record in SeqIO.parse(fasta_handle, "fasta"):
		print("Sequence "+record.id+" is being read")
		dict_pos[record.id]=[]
		start_pos=0
		counter=0
		gap=False
		gap_length = 0
		for char in record.seq:
			if char == 'N':
				if gap_length == 0:
					start_pos=counter
					gap_length = 1
					gap = True
				else:
					gap_length += 1
			else:
				if gap:
					bed_out.write(record.id + "\t" + str(start_pos) + "\t" + str(start_pos + gap_length) + "\n")
					dict_pos[record.id].append([start_pos,start_pos + gap_length])
					gap_length = 0
					gap = False
			counter += 1

#tsv_handle=open(sys.argv[3], mode="r")
#reader = csv.DictReader(tsv_handle, delimiter = '\t')
bed_out.close()

'''for variant in reader:
   for inter in dict_pos[variant['CHROM']]:
	   in_repeted_part=False
	   if int(inter[0])<=int(variant['POS'])<=int(inter[1]):
		   in_repeted_part=True
		   print(in_repeted_part)'''