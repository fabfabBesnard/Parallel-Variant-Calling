# Variant calling and prediction of functional impact

## Intro

The aim of this pipeline is to search genetic variants within different mutated plants in order to highlight variants with a phenotypic impact and shed light on the function of certain genes. Thanks to NGS data and by comparing the variants present in the plants, we are able to identify and compare the variants that have a phenotypic impact. 

For this project, I choose Nextflow, because this pipeline framework have a lot of avantage and it's very easy to install (  [How to install nextflow](https://www.nextflow.io/docs/latest/getstarted.html) ).

This pipeline have 3 steps :

- Mapping and processing reads : Mapping reads from diffeentas mutated sample against reference genome with bwa mem. After that, sam file are filtering and annotate with samtools and picardtools.

- Variant Calling : ...

- Effet prediction : ...


<img src="img/Tech-Flowchart.jpg" alt="Flowchart" width="800"/>


