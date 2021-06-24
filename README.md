# Parallel genomic variant calling and prediction of functional impacts

Improve and speed-up your forward genetics by listing all mutations in several samples at a time !

### Content
  1. [Introduction]()
  2. [Pipeline overview]()
  3. [Requirements]()
  4. [General Usage]()
  5. [Input files and parameters]()
  6. [Example]

## Introduction

This pipeline is designed to automatically provide the most exhaustive and accurate possible list of genes affected by genomic variations (e.g. natural polymorphisms, mutations) in a sample's DNA, using Illumina paired-end sequencing data. 

The output is a user-friendly tsv table that can be parsed and filter with classical spreadsheet software (LibreOffice, Excel, ...). This table is sorted by gene and predict functional impact(s) of the identified genomic variation(s) to help biologists find the best candidates genes modified in the sample(s) provided.
http://www.ens-lyon.fr/PSMN/doku.php
The two main strengths of this pipeline are:
- **Automatic parallel analysis of a cohort of samples**: several samples' sequencing data can be provided at once and the pipeline automatically select genomic variations that are specific to each sample of the cohort. For example, in a mutagenesis experiment with a starting strain and several derived mutant strains, all inherited mutations from the starting strain will be discarded. This parallel analysis and multiple pairwise comparisons significantly **improve the specificity of the mutation search** by reducing false positive rate, while the automation of the workflow makes it easily scalable to large cohorts.
- **Exhaustive variant calling**: the pipeline automatically combines several variant callers to cover a **large spectrum of possible genomic variations**, from single nucleotide polymorphisms (**SNP**) up to **structural variations** (SV) of several kbps (deletions, insertions, inversions, translocations, copy number variations, etc...). This improves variant calling accuracy and resolution, especially for SV, while again pipeline automation ensure a simple workflow for biologist end users.

The pipeline is implemented in [Nextflow](https://www.nextflow.io/): it's very easy to install and allows to monitor the completion of all processes of the pipeline, can be deployed in clusters/clouds for parallel computing, it ensures reproducible analysis (simple configuration, supports Docker technology, keeps track of command lines and parameters), promotes efficient re-run and debugs, generates reports.

In theory, it can be used on every organisms for which a reference genome and annotation files are available (flexible input data provided by the user).
Organisms in which the pipeline has been tested: *Arabidopsis thaliana*, *Physcomitrium patens*.

## Pipeline overview (technical description)
This pipeline has 3 steps:

- Mapping and processing reads : Mapping reads from different mutated samples against reference genome with bwa mem. After that, sam file are filtering and annotate with samtools and picardtools.

- Variant Calling: A variant calling of short variants, snp and short indels, with [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller). And a structural variant calling, with pindel, Breakdancer, CNVnator, Lumpy combined with [metaSV](https://github.com/bioinform/metasv). If several samples are provided together, only variants unique to each strain will be selected, increasing specificity of the genomic variation profiling.


- Effet prediction: A prediction of variant effects with [snpeff](http://pcingola.github.io/SnpEff/).


<img src="img/Tech-Flowchart.jpg" alt="Flowchart" width="600"/>

## requirements

### Software
You will need to install first:
- Nextflow ([How to install nextflow](https://www.nextflow.io/docs/latest/getstarted.html) in three little steps ! ).

> note
> 
> As explained in the link, we recommand to place the nextflow binary executable file at a location accessible to your $PATH. Alternatively, you can permanently edit your $PATH. For example, in a bash shell, execute: `echo 'export PATH=$PATH:/path/to/nextflow'  >> ~/.bash_profile` (in a sh and ksh shell, `echo 'export PATH=$PATH:/path/to/nextflow' >> ~/.profile`)

- Plotly (python3 -m pip install plotly --user)

- just git clone or copy this repository (not install or compilation required !)

### input data
The sequencing data are typically **Illumina paired-end** sequencing fastq files generated from the genomic DNA of a unique strain.
Retained genomic variations are only **homozygous positions**. Both haploid and diploid organisms can be studied.

The sequence of the reference genome of the organism is required at the *fasta* format, annotation file is required at the *gff* format.

## General Usage

We have implemented a single workflow called 'VariantCaller.nf'. This pipeline can be called using the following command:

```
nextflow run -profile [psmn,singularity] VariantCaller.nf
-c file.config
--reads "reads/*_{1,2}.fq.gz"
--genomefasta genome.fa
--annotation genome_annotation.gff
--outdir My_analyse
```
We recommand to run nextflow using a *screen* or *tmux* terminal to avoid killing the workflow if the terminal closes by accident.

Example with tmux :
```
# Open a new terminal
tmux
# run nextflow
nextflow run -c -profile
# "Detach" terminal (the terminal vanishes from the screen but do not stop, so the running process is not killed)
ctrl+B (simultaneous keys) then type 'd'
# "Attach" terminal (re-open the remote terminal session which is still active)
tmux attach
```

## Input files and parameters

### configuration file
First define the config file of your analysis workflow.
Copy and paste the default config file provided in `Pipeline_variant_RDP/script/VariantCaller.config` to your analysis folder.

> note
> 
>  Only two profiles are defined in the second part of the config file. The pipeline has been optimized to run in our lab cluster. Running in an other environment may require to create a new profile. Expertise in nextflow is preferable to edit the profile configurations.

The configuration file must be given to the command line with the -c option:
- `-c` : Path to configuration file. This file contains all defaut values and configurations for run at our lab cluster ([PSMN](http://www.ens-lyon.fr/PSMN/doku.php)). 

Most of the main paramters can be either directly configured at the beginning of this config file. Alternatively, they can be specified in the command line as explained below.

### main parameters
- `-profile` : the profile adapted to your computig environment, deined in the config file (available: psmn or singularity)

- `VariantCaller.nf`: the typical pipeline to be executed (only one available)

- `--scriptdir` : Path to the directory that contains the pipeline scripts. For example /Pipeline_variant_RDP/script'

- `--genomefasta` : Full path to file of reference genome (.fa or .fasta or .gz)

- `--reads` : Full path to directory and name of reads in fastq.gz. You can use a pattern to select several files. 
Example : Sequence\* _{1,2}.fastq ( `{1,2}` for paired reads ). 

Symbolic links are accepted.

Provide an annotation file at gff format with `--annotationgff`. For *Physcomitrium patens* and *Arabidopsis thaliana* only, you do not need to provide an annotation file as the data are incorporated in our snpEff docker image. Use the option `--annotationname` instead.

- `--annotationgff` : Full path to file of annotation genome (.gff)

- `--annotationname` : Name of the organism  (either 'Physcomitrella_patens' or 'Arabidopsis_thaliana')


### optional parameters:

- `--vqsr_file` : You can provide a reference variant file (.vcf) in order to apply a variant recalibartion score. For exemple [Arabidopsis_thaliana]( https://1001genomes.org/data/GMI-MPI/releases/v3.1/)

- `--sampletable` In order to simplify result file naming, you can provide a table of correspondance between your regular sample names and the (often ) for your reads in csv format. Example:
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
