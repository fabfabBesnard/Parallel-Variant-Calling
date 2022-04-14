#!/usr/bin/python3.6
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
vcf_out=open(sys.argv[2], mode="w")

for line in vcf_in.readlines():

    if '#' in line:
        vcf_out.write(line) #Add the '#' part of a VCF file

    else:

        line=line.split('\t') #Parse the line

        if 'ANN' in line[7] and line[10]=='1\n': #If Ann is in the line and the variant is in masked region with parse it
            
            info=line[7].split(';')#Parse the info part of vcf file
            keep_str=line #Variable useful to rebuild the correct line in the file
            
            for ann in info:

                if 'ANN' in ann: #We get the ANN part 
                    ann=ann.split('=') #Separate ANN= and the content of ANN
                    keep_str[7]=ann
                    ann=ann[1].split(',') #Parse the different annotation
                    keep_str[7][1]=[]

                    for variant_ann in ann:
                        variant_ann=variant_ann.split('|') #Parse the different ANN part
                        keep_str[7][1].append(variant_ann) #Fill in keep_str

                    for variant_ann_pos in range(0,len(keep_str[7][1])): #We add the WARNING_MASKED_REGION info

                        if keep_str[7][1][variant_ann_pos][15]:
                            keep_str[7][1][variant_ann_pos][15]=keep_str[7][1][variant_ann_pos][15]+'&WARNING_MASKED_REGION'

                        else:
                            keep_str[7][1][variant_ann_pos][15]='WARNING_MASKED_REGION'
                            
                        keep_str[7][1][variant_ann_pos]='|'.join(keep_str[7][1][variant_ann_pos])

            keep_str[7][1]=','.join(keep_str[7][1]) #We rebuild the line
            keep_str[7]='='.join(keep_str[7])
            del(keep_str[-1])
            keep_str='\t'.join(keep_str)
            vcf_out.write(keep_str+'\n')

        else:
            del(line[-1])
            vcf_out.write('\t'.join(line)+'\n')