#!/usr/bin/env python3
import glob
import sys

#en poucrentage
tolerance = 1
offset = sys.argv[1]
psize= sys.argv[2]

vcf = {}

listfichier = glob.glob('*.vcf')

for fichiervcf in listfichier:
    var = {}
    header = ''
    for ligne in open(fichiervcf, 'r'):
        if ligne.startswith('#'):
            header = header+ligne
        else :
            #1	471760	1	N	<DEL>	.	.	SVTYPE=DEL;STRANDS=+-:8;SVLEN=-6373;END=478133;CIPOS=-10,271;CIEND=-253,9;CIPOS95=0,82;CIEND95=-69,2;IMPRECISE;SU=8;PE=8;SR=0	GT:SU:PE:SR	./.:8:8:0
            li = ligne.split('\t')
            chrom = li[0]
            pos = int(li[1])
            for i in li[7].split(';'):
                if i.startswith('SVLEN='):
                    svlen = int(i[6:])
                    break
            if chrom not in var : 
                var[chrom] = [ [pos , svlen , ligne , True] ]
            else : 
                var[ chrom ].append( [pos , svlen , ligne , True])
        vcf[fichiervcf] = {'header' : header , 'dichrom' : var}


for fichierencours in vcf:
    listfichieracomparer = listfichier.copy()
    listfichieracomparer.remove(fichierencours)
    print( fichierencours)
    for chromencours in vcf[fichierencours]['dichrom']:
        for variant in vcf[fichierencours]['dichrom'][chromencours] :
            for fichiercomparer in listfichieracomparer:
                try :
                    for variant_bis in vcf[fichiercomparer]['dichrom'][chromencours] :
                        position1 = variant[0]-int(offset)
                        position2 = variant[0]+int(offset)
                        #si la position du variant est dans l'interval 
                        if (variant_bis[0] >= position1 ) & (variant_bis[0] <=position2 ):
                            taille1 = abs(variant[1])*(1-float(psize))
                            taille2 = abs(variant[1])*(1+float(psize))
                            if (abs(variant_bis[1]) >= taille1 ) & (abs(variant_bis[1]) <= taille2 ):
                                variant[3] = False
                                variant_bis[3] = False
                                break
                except KeyError:
                    pass
                            

for fichierencours in vcf:
    c = 0 
    n = 0
    filename = fichierencours[:-4]+'_filtered_SV.vcf'
    print(filename)
    f = open( filename, "w")
    f.write(vcf[fichierencours]['header'])
    for chromencours in vcf[fichierencours]['dichrom']:
        for variant in vcf[fichierencours]['dichrom'][chromencours] :
            if variant[3] == True:
                c+=1
                f.write(variant[2])
            else : 
                n +=1

    print( c , n)
    f.close()


