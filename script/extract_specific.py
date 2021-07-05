#!/usr/bin/env python3

# Garder les lignes ## info 
# Les stokcker dans une variable pour editer les prochains VCF 
# faire de la comparaison ligne par ligne et si un sample ressort en etant 1/1 ou 2/2 ou 0/0 tout seul parmis les autres 
# Editer un fichier vcf par echantillon 

#definir le variant avec la list des 
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Mutant1	Mutant2	Mutant3	Mutant4	Mutant5	StartingStrain
#1	7217	.	T	C	2782.01	PASS	AC=12;AF=1.00;AN=12;DP=73;ExcessHet=3.0103;FS=0.000;MLEAC=12;MLEAF=1.00;MQ=60.00;QD=30.02;SOR=0.980	GT:AD:DP:GQ:PL	1/1:0,19:19:57:758,57,0	1/1:0,5:5:15:206,15,0	1/1:0,9:9:27:358,27,0	1/1:0,16:16:48:622,48,0	1/1:0,5:5:15:210,15,0	1/1:0,15:15:45:619,45,0
#1	346347	.	C	T,*	241.56	PASS	AC=2,8;AF=0.200,0.800;AN=10;DP=39;ExcessHet=3.0103;FS=0.000;MLEAC=2,9;MLEAF=0.200,0.900;MQ=60.00;QD=8.63;SOR=1.329	GT:AD:DP:GQ:PGT:PID:PL:PS	2|2:0,0,8:8:24:1|1:346312_A_ACC:357,357,357,24,24,0:346312	1|1:0,4,0:4:12:1|1:346312_A_ACC:180,12,0,180,12,180:346312	./.:4,0,0:4:.:.:.:0,0,0,0,0,0	2/2:0,0,5:7:35:.:.:482,344,312,36,35,0	2|2:0,0,5:5:15:1|1:346312_A_ACC:221,221,221,15,15,0:346312	2|2:0,0,6:6:18:1|1:346312_A_ACC:270,270,270,18,18,0:346312


import sys

def extractgoodvariant( ligne ):
    toadd = False
    samplerank = 0
    variantline = ''
    variant_in_liste = ligne.split('\t')
    #Keep sample information only 
    variant_in_liste = variant_in_liste[9:]
    
    #identifier variant stratingstrain !!
    # variantSS = variant_in_liste[startingstrainrank]
    # GTSS = variantSS.split(":")[0]
    # del variant_in_liste[startingstrainrank]

    dicoGT = {}
    for i in variant_in_liste:
        GT = i.split(":")[0]
        #Si la taille du variant est 1 -> polidy 1
        if len(GT) == 1 :
            if GT == '0':
                continue
            elif GT not in dicoGT.keys():
                #Ajout d'une clef GT avec comme valeur une liste des rang du variant
                dicoGT[ GT ] = [ samplerank ]
            else:
                dicoGT[ GT ].append(samplerank)
            samplerank+=1
        #Sinon diploide 
        else :
            if GT[0] == '0':
                continue
            elif GT not in dicoGT.keys():
                #Creation d'une clef GT avec comme valeur une liste des rang du variant
                dicoGT[ GT ]= [samplerank]
            else:
                dicoGT[ GT ].append(samplerank)
            samplerank+=1
    #Parcours le dico pour extraire un variant unique 
    for GT in dicoGT:
        #if GT != GTSS:
        if len( GT )!= 1:
            #if homozyogous 
            if (GT[0] == GT[2]) and ('.' not in GT) :
                #if present only one time 
                if len(dicoGT[GT]) == 1 :
                    toadd = True
                    samplerank = dicoGT[GT]
        else:
            if ('.' not in GT):
                if len(dicoGT[GT]) == 1 :
                    toadd = True
                    samplerank = dicoGT[GT]

    return toadd, samplerank

def create_new_variant_line( variant , rank):
    #Creer la ligne a ajouter en gardant le bon variant et en ajoutant les autre dans INFO
    newvarantlineinlist = variant.split('\t')[:9]
    #recupere les informations des autres variants
    othersampleinfo = variant.split('\t')[9:]
    del( othersampleinfo[rank] )
    newvarantlineinlist[7] = newvarantlineinlist[7]+";"+"OTHERVAR="+str(othersampleinfo)+";"
    scorevariant = variant.split('\t')[9+rank]
    newvarantlineinlist.append(scorevariant)
    newvarantline = '\t'.join( newvarantlineinlist)
    return newvarantline

vcfname = sys.argv[1]
qualseuil = sys.argv[2]
dp = sys.argv[3]

#vcfname = "/home/romuald/Documents/stage_rdp/filtered_snps.vcf"

# init vcf information
vcfheader = ''
dicovariant = {}


number_removed_qual = 0
number_removed_dp = 0 

# on parcours le contenu du fichier ligne par ligne
for ligneN in open(vcfname, 'r'):
    ligne = ligneN.strip('\n')
    if ligne.startswith('##'):
        vcfheader = vcfheader + ligne + '\n'
    elif ligne.startswith('#CHROM'):
        vcfheader = vcfheader + '##INFO=<ID=OTHERVAR,Number=.,Type=String,Description="Othersample information">\n##source=ExtractGoodvariantnextflowprocess'+ '\n'
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Mutant1	Mutant2	Mutant3	Mutant4	Mutant5	StartingStrain
        header=ligne
        info = ligne.split('\t')
        samplelist = info[9:]
        #check rank of stratin strain 
        #indexSS = samplelist.index(startingstrain)
         # Variant line 
        #CASE ploidy 2
        #1	7217	.	T	C	2782.01	PASS	AC=12;AF=1.00;AN=12;DP=73;ExcessHet=3.0103;FS=0.000;MLEAC=12;MLEAF=1.00;MQ=60.00;QD=30.02;SOR=0.980	GT:AD:DP:GQ:PL	1/1:0,19:19:57:758,57,0	1/1:0,5:5:15:206,15,0	1/1:0,9:9:27:358,27,0	1/1:0,16:16:48:622,48,0	1/1:0,5:5:15:210,15,0	1/1:0,15:15:45:619,45,0
        #CASE ploidy 1
        #1	25062	.	T	C	2805.73	.	AC=6;AF=1.00;AN=6;DP=75;FS=0.000;MLEAC=6;MLEAF=1.00;MQ=60.00;QD=32.49;SOR=0.804	GT:AD:DP:GQ:PL	1:0,22:22:99:857,0	1:0,9:9:99:309,0	1:0,12:12:99:475,0	1:0,16:16:99:652,0	1:0,7:7:99:239,0	1:0,8:8:99:287,0
    else :
        qualvariant = ligne.split('\t')[5]
        if float(qualvariant) < float(qualseuil):
            number_removed_qual += 1
            #go to next line
            continue
        else :
            add , samplerank = extractgoodvariant(ligne )
            if add :
                for i in samplerank:
                    newline = create_new_variant_line( ligne , i)
                    dpsample = newline.split('\t')[-1].split(':')[2]
                    if int(dpsample) < int(dp):
                        number_removed_dp += 1
                        continue
                    if samplelist[i] not in dicovariant.keys():
                        #Creation d'une clef avec comme valeur le rang du variant
                        #Creation de la string a ajouter au fichier vcf en ne selcetionnant que le variant d'interet
                        firstline =  "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO    "+samplelist[i]
                        dicovariant[ samplelist[i] ] = [ firstline , newline ]
                    else:
                        dicovariant[ samplelist[i] ].append( newline )

for filename in dicovariant:
    fichier = open(filename+"_"+vcfname, "w")
    fichier.write(vcfheader)
    for i in dicovariant[filename]:
        if i.endswith("\n"):
            fichier.write(i)
        else:
            fichier.write(i + '\n')
    fichier.close()

fstat = open( "nb_removed",'w')
fstat.write("number_removed_qual\tnumber_removed_dp\n")
fstat.write(str(number_removed_qual)+"\t"+str(number_removed_dp))