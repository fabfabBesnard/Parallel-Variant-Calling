# READ ME 
## Pipeline de recherche de variant 

The aim of this pipeline is to search genetic variants within different mutated plants in order to highlight variants with a phenotypic impact and shed light on the function of certain genes. Thanks to NGS data and by comparing the variants present in the plants, we are able to identify and compare the variants that have a phenotypic impact. 

For this project, I choose Nextflow, because this pipeline framework have a lot of avantage : 
 - Very easy to install and it only requires Bash 3.2 (or later) and Java 8 (or later). It does not require any special installation procedure, follow this 2 steps tutorial : [How to install nextflow](https://www.nextflow.io/docs/latest/getstarted.html)


 - Allows you to work locally on your personal computer or on a cluster (PSMN) or a cloud ( Amazon cloud, google ...). 

 - Nextflow use a system of processes and channels. A process is a simple task, with a input, an output, a condition (optionnal) and a script. A channel is a link between processes : a colection of files, a path ... All the processes and chanel are stocked in a cache so you dont have to re-execute it from scratch. 

 - One of the bigget advantage is the containers. You can Processes works with a specific docker 
 https://www.nextflow.io/docs/latest/docker.html
 et aaussi un groupe d'image docker https://www.nextflow.io/docs/latest/docker.html#multiple-containers
 blabla docker hub 
 


