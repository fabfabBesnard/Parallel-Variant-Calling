#!/usr/bin/env python3


import sys

vcfname = sys.argv[1]


#list de variant [ chrom, pos, qual, istrue]


print( "CHROM\tPOS\tSVLEN")

li = []
# on parcours le contenu du fichier ligne par ligne
for ligneN in open(vcfname, 'r'):
    ligne = ligneN.strip('\n')
    if ligne.startswith('#'):
        pass
    else :
        #1	471760	1	N	<DEL>	.	.	SVTYPE=DEL;STRANDS=+-:8;SVLEN=-6373;END=478133;CIPOS=-10,271;CIEND=-253,9;CIPOS95=0,82;CIEND95=-69,2;IMPRECISE;SU=8;PE=8;SR=0	GT:SU:PE:SR	./.:8:8:0
        li = ligne.split('\t')
        chrom = li[0]
        pos = int(li[1])
        if li[-1][:3] == "0/0":
            continue
        for i in li[7].split(';'):
            if i.startswith('SVLEN='):
                svlen = int(i[6:])
                print('{0}\t{1}\t{2}'.format(chrom,pos , svlen) )
                break
        
    