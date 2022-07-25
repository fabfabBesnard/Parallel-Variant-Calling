#!/usr/bin/python3
import sys

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

vcf_in=open(sys.argv[1], mode="r")
count=0

for line in vcf_in.readlines():

    if '#' not in line:

        line=line.split('\t') #Parse the line
        for othervar in line[7].split(';') :
            if "OTHERVAR" in othervar :
                if "./." in othervar or ".|." in othervar or "." in othervar:
                    print(str(othervar)+" "+str(line[9]))
            #othervar = str(line[7].split(';')[11].split("=")[1][1:-1])
            #print(othervar)
        allele=line[9].split(':')
        if "./." in allele[0] or ".|." in allele[0] or "." in allele[0]:
            count=count+1

print("./.="+str(count))
