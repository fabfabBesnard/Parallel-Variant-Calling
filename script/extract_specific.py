#!/usr/bin/env python3

# Garder les lignes ## info 
# Les stocker dans une variable pour editer les prochains VCF 
# faire de la comparaison ligne par ligne et si un sample ressort en etant 1/1 ou 2/2 ou 0/0 tout seul parmis les autres 
# Editer un fichier vcf par echantillon 

#definir le variant avec la list des 
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Mutant1	Mutant2	Mutant3	Mutant4	Mutant5	StartingStrain
#1	7217	.	T	C	2782.01	PASS	AC=12;AF=1.00;AN=12;DP=73;ExcessHet=3.0103;FS=0.000;MLEAC=12;MLEAF=1.00;MQ=60.00;QD=30.02;SOR=0.980	GT:AD:DP:GQ:PL	1/1:0,19:19:57:758,57,0	1/1:0,5:5:15:206,15,0	1/1:0,9:9:27:358,27,0	1/1:0,16:16:48:622,48,0	1/1:0,5:5:15:210,15,0	1/1:0,15:15:45:619,45,0
#1	346347	.	C	T,*	241.56	PASS	AC=2,8;AF=0.200,0.800;AN=10;DP=39;ExcessHet=3.0103;FS=0.000;MLEAC=2,9;MLEAF=0.200,0.900;MQ=60.00;QD=8.63;SOR=1.329	GT:AD:DP:GQ:PGT:PID:PL:PS	2|2:0,0,8:8:24:1|1:346312_A_ACC:357,357,357,24,24,0:346312	1|1:0,4,0:4:12:1|1:346312_A_ACC:180,12,0,180,12,180:346312	./.:4,0,0:4:.:.:.:0,0,0,0,0,0	2/2:0,0,5:7:35:.:.:482,344,312,36,35,0	2|2:0,0,5:5:15:1|1:346312_A_ACC:221,221,221,15,15,0:346312	2|2:0,0,6:6:18:1|1:346312_A_ACC:270,270,270,18,18,0:346312


import sys
import numpy as np
import pandas as pd
import plotly.express as px

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
        #Si la taille du variant est 1 -> ploidy 1
        if len(GT) == 1 :
            if GT not in dicoGT.keys():
                #Ajout d'une clef GT avec comme valeur une liste des rang du variant
                dicoGT[ GT ] = [ samplerank ]
            else:
                dicoGT[ GT ].append(samplerank)
            samplerank+=1
        #Sinon diploide 
        else :
            if GT not in dicoGT.keys():
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
            if (GT[0] == GT[2]) and ('.' not in GT) and ('0' not in GT) :
                #if present only one time 
                if len(dicoGT[GT]) == 1 :
                    toadd = True
                    samplerank = dicoGT[GT]
        else:
            if ('.' not in GT) and ('0' not in GT):
                if len(dicoGT[GT]) == 1 :
                    toadd = True
                    samplerank = dicoGT[GT]
    if toadd == True:
        for GT in dicoGT:
            #if GT != GTSS:
            if '.' in GT :
                nbundefini = len(dicoGT[GT])
                if nbundefini in dico_undefined:
                    dico_undefined[ nbundefini ] = dico_undefined[ nbundefini ] + 1
                else :
                    dico_undefined[ nbundefini ] = 1
            if "." not in dicoGT:
                if 0 in dico_undefined:
                    dico_undefined[ 0 ] = dico_undefined[ 0 ] + 1
                else :
                    dico_undefined[ 0 ] = 1

    return toadd, samplerank

def good( ligne ):
    samplerank = 0
    toadd = False
    variant_in_liste = ligne.split('\t')
    #Keep sample information only 
    variant_in_liste = variant_in_liste[9:]

    dicoGT = {}
    for i in variant_in_liste:
        GT = i.split(":")[0]
        #Si la taille du variant est 1 -> ploidy 1
        if len(GT) == 1 :
            if GT not in dicoGT.keys():
                #Ajout d'une clef GT avec comme valeur une liste des rang du variant
                dicoGT[ GT ] = [ samplerank ]
            else:
                dicoGT[ GT ].append(samplerank)
            samplerank+=1
        #Sinon diploide 
        else :
            if GT not in dicoGT.keys():
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
            if (GT[0] == GT[2]) and ('.' not in GT) and ('0' not in GT) :
                #if present only one time 
                if len(dicoGT[GT]) == 1 :
                    toadd = True
                    samplerank = dicoGT[GT]
        else:
            if ('.' not in GT) and ('0' not in GT):
                if len(dicoGT[GT]) == 1 :
                    toadd = True
                    samplerank = dicoGT[GT]
    return toadd

def check_othervar_unknown( othersampleinfo ): #Verifie si les autres échantillons ne possède aucun read aligné et donc sont tous undiefined
    nb_read_tot=[]
    for sample in othersampleinfo:
        read=sample.split(":")[1].split(",")
        nb_read=0
        for i in read:
            nb_read=nb_read+int(i)
        nb_read_tot.append(nb_read)
    if sum(nb_read_tot)==0:
        return True
    else :
        return False

def check_othervar_false_positive ( othersampleinfo, sampleinfo ):
    for sample in othersampleinfo:
        allele=sample.split(":")[0]
        read=sample.split(":")[1].split(",")
        if allele == "." or allele =="./." or allele == ".|.":
            allele_count=read[int(sampleinfo[0])]
            total_read_count=0
            for i in read:
                total_read_count=total_read_count+int(i)
            if total_read_count != 0:
                prop_allele = int(allele_count)//int(total_read_count)
                #print(prop_allele)
                if prop_allele >= 0.25 :
                    return True
    return False

def create_new_variant_line( variant , rank, samplelist):
    #Creer la ligne a ajouter en gardant le bon variant et en ajoutant les autre dans INFO
    newvarantlineinlist = variant.split('\t')[:9]
    #recupere les informations des autres variants
    othersampleinfo = variant.split('\t')[9:]
    for i in range(0, len(othersampleinfo)) :
        othersampleinfo[i]=str(othersampleinfo[i])+":"+str(samplelist[i])
    sampleinfo=othersampleinfo[rank]
    del( othersampleinfo[rank] ) #on suprimme les infos concernant l'échantillon en cours de modification
    if check_othervar_unknown(othersampleinfo): #
        newvarantlineinlist[7] = newvarantlineinlist[7]+";"+"OTHERVAR="+str(othersampleinfo)+";WARNING_SPECIFIC_ALLELE=unread_position_in_other_samples"
    elif check_othervar_false_positive (othersampleinfo, sampleinfo):
        newvarantlineinlist[7] = newvarantlineinlist[7]+";"+"OTHERVAR="+str(othersampleinfo)+";WARNING_SPECIFIC_ALLELE=allele_detected_in_other_sample(s)"
    else :
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

dico_undefined = {}

number_removed = 0 

quallist = []
# On parcours le contenu du fichier ligne par ligne
for ligneN in open(vcfname, 'r'):
    ligne = ligneN.strip('\n')
    if ligne.startswith('##'):
        vcfheader = vcfheader + ligne + '\n'
    elif ligne.startswith('#CHROM'):
        vcfheader = vcfheader + '##INFO=<ID=OTHERVAR,Number=.,Type=String,Description="INFO field from other variants grouped in gVCF WRITTEN AS INFO;SAMPLENAME">\n##INFO=<ID=WARNING_SPECIFIC_ALLELE,Type=String,Description="indicates potential warning cases where the sample allele may not be unique to this sample. Check other samples of the cohort."\n##source=ExtractGoodvariantnextflowprocess'+ '\n'
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
        if good(ligne):
            quallist.append(float(qualvariant))
        if float(qualvariant) < float(qualseuil):
            number_removed += 1
            #go to next line
            continue
        else :
            #print("vcf accepte")
            add , samplerank = extractgoodvariant(ligne )
            if add :
                for i in samplerank:
                    newline = create_new_variant_line( ligne , i, samplelist)
                    dpsample = newline.split('\t')[-1].split(':')[2]
                    #print(dpsample)
                    if int(dpsample) < int(dp):
                        number_removed += 1
                        continue
                    if samplelist[i] not in dicovariant.keys():
                        #Creation d'une clef avec comme valeur le rang du variant
                        #Creation de la string a ajouter au fichier vcf en ne selcetionnant que le variant d'interet
                        firstline =  "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO    "+samplelist[i]
                        dicovariant[ samplelist[i] ] = [ firstline , newline ]
                    else:
                        dicovariant[ samplelist[i] ].append( newline )
if not dicovariant:
    fichier = open("empty_"+vcfname.split("/")[-1], "w")
    fichier.close()
else:
    for filename in dicovariant:
        fichier = open(filename+"_"+vcfname.split("/")[-1], "w")
        fichier.write(vcfheader)
        for i in dicovariant[filename]:
            if i.endswith("\n"):
                fichier.write(i)
            else:
                fichier.write(i + '\n')
        fichier.close()

fstat = open( "nb_removed",'w')
fstat.write("number_removed\n")
fstat.write( str(number_removed))

print(dico_undefined)

k = list(dico_undefined.keys()) 
v = list(dico_undefined.values())

data_undefined = {'Variant_undefined':k, 'Count': v}

df_undefined = pd.DataFrame.from_dict(data_undefined)

bar = px.bar(df_undefined, x='Variant_undefined', y='Count')
bar.write_html("Undefined_variants_number_for_"+vcfname.split('/')[-1]+"_mqc.html")

print(vcfname)

dfqual = pd.DataFrame(quallist, columns =['Quality'])
dfqual["filtered"] = np.where(dfqual["Quality"] < float(qualseuil), True, False)

print(dfqual.describe())
histqual = px.histogram(dfqual, x="Quality", color="filtered", title="Quality threshold = "+str(qualseuil), log_y = True )
histqual.write_html( "Quality_of_specific_variation_for_"+vcfname.split('/')[-1]+"_mqc.html")