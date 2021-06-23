# Variant calling and prediction of functional impact

## Intro

The aim of this pipeline is to search genetic variants within different mutated plants in order to highlight variants with a phenotypic impact and shed light on the function of certain genes. Thanks to NGS data and by comparing the variants present in the plants, we are able to identify and compare the variants that have a phenotypic impact.

For this project, I choose Nextflow, because this pipeline framework has a lot of advantage and it's very easy to install ( [How to install nextflow](https://www.nextflow.io/docs/latest/getstarted.html) ).

This pipeline has 3 steps:

- Mapping and processing reads : Mapping reads from differents mutated sample against reference genome with bwa mem. After that, sam file are filtering and annotate with samtools and picardtools.

- Variant Calling: A variant calling of short variant, snp and short indel, with [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller). And a structural variant calling, with pindel, Breakdancer, CNVnator, Lumpy group with [metaSV](https://github.com/bioinform/metasv).


- Effet prediction: A prediction of variant effect with [snpeff](http://pcingola.github.io/SnpEff/).


<img src="img/Tech-Flowchart.jpg" alt="Flowchart" width="600"/>

## requirements

## install

You need to have :
- Nextflow
- Plotly (python3 -m pip install plotly --user)


## General Usage

A workflow is run in the following way:

```
nextflow run -profile [psmn,singularity] VariantCaller.nf
-c file.config
--reads "reads/*_{1,2}.fq.gz"
--genomefasta genome.fa
--annotation genome_annotation.gff
--outdir My_analyse
```
The best way to do this is to run nextflow using a screen or tmux terminal.


Example with tmux :
```
# Open a new terminal
tmux
# run nextflow
nextflow run -c -profile
# "Detach" terminal
+ d
# "Attach" terminal
tmux attach
```

## Input files and parametres

- `-c` : Path to configuration file. This file contain all defaut values and configuttion for psmn. You can also change value of option directyle in configuration file of in command line with the following options.

- `--scriptdir` : Path to directory who contain script. For example /Pipeline_variant_RDP/script.nf'

- `--genomefasta` : Full path to file of reference genome (.fa or .fasta or .gz)

- `--reads` : Full path to directory and name of reads in fastq.gz. You can use a pattern to select several files. Example : Sequence*_{1,2}.fastq ( `{1,2}` for paired reads )

You can provide your own annotation file with `--annotationgff` or you can use snpeff database. Only for Physcomitrella patens or Arabidopsis thaliana ! In this case use `--annotationname`

- `--annotationgff` : Full path to file of annotation genome (.gff)

- `--annotationname` : Name of organism : 'Physcomitrella_patens' or 'Arabidopsis_thaliana'


### Optional:

- `--vqsr_file` : You can provide a reference variant file (.vcf) in order to apply a variant recalibartion score. For exemple [Arabidopsis_thaliana]( https://1001genomes.org/data/GMI-MPI/releases/v3.1/)

- `--sampletable` In order to simplify results, you can provide a table of sample name for your reads in csv format. Example:
```
V300042688_L2_AE06084935-608,Mutant1
V300042688_L4_AE49584879-612,Mutant2
V300042688_L4_AE47136387-610,Mutant3
```

- `--outdir` : Name of directory generate. Default: 'results'

- `--ploidy` : Number of chromosome copy. Default: 1 (1 or 2)

- `--minglobalqual` : (Only for short indel and snp) Threeshold of global quality per variant (for all sample). Default: 200

- `--mindepth` : (Only for short indel and snp) Threeshold of deepth (number of reads) for a variant per sample. Default: 4

- `-resume` : Add resume to command line and files from same analyse are retrieved from the cache.


## Example 

Command line : 

```
./nextflow run script/VariantCaller.nf -c script/VariantCaller.config 
--reads "/home/rmarin/Mydata/V30001743*{1,2}.fq.gz" 
--genomefasta ../Mydata/Arabidopsis_thaliana.TAIR10.31.dna.toplevel.fa 
--vqsrfile ../Mydata/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
-profile psmn 
--annotationname 'Arabidopsis_thaliana' 
--outdir My_analyse
```

```
.
├── Mydata
│   ├── 1001genomes_snp-short-indel_only_ACGTN.vcf
│   ├── Arabidopsis_thaliana.TAIR10.31.dna.toplevel.fa
│   ├── V300017433_L4_B5GARAwyvRAAAAABA-501_1.fq.gz
│   ├── V300017433_L4_B5GARAwyvRAAAAABA-501_2.fq.gz
│   ├── V300017433_L4_B5GARAwyvRAAABABA-502_1.fq.gz
│   └── V300017433_L4_B5GARAwyvRAAABABA-502_2.fq.gz
├── script
│   ├── breakdancer2vcf.py
│   ├── extract_specific.py
│   ├── extract_specific_SV.py
│   ├── extractSplitReads_BwaMem
│   ├── extractsvlen.py
│   ├── fastq_sample.py
│   ├── graph_qual_fromraw copy.py
│   ├── graph_qual_fromraw.py
│   ├── splitfa.py
│   ├── testmappingbowtie2.nf
│   ├── testmapping.nf
│   ├── VariantCaller.config
│   └── VariantCaller.nf
└── nextflow
```
