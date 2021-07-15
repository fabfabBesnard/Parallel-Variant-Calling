# How to update snpeff docker
Snpeff provides a large number of annotations directly included in the database. To save memory space, it would not be wise to include to much of them in the docker, so we arbitrarily chose a few genomes of model organisms. In case you do not have your favorite organism annotation in the genomes already included in our snpeff docker, you can try to create your own docker containing a snpeff buil.
>Note: 
>
>the pipeline offers anoter alternative: you can try to build the annotation prediction from a gff with snpeff using the main option `--annotationgff`.


The first step is to find if your ogranism are included in the snpeff database. For this download [snpeff](https://pcingola.github.io/SnpEff/download/) and install it.

One you have installed snpeff, do this command line to find out if snpeff provides an annotation build of your organism : 

`java -jar snpEff.jar databases | grep "Myspecies" `

For example : you want to run the pipeline with the genome of *Nicotiana attenuata*.

```sh
Linux:~/snpEff$ java -jar snpEff.jar databases | grep "Nicotina"
Nicotiana_attenuata                                         	Nicotiana_attenuata                                         	          	                              	https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_Nicotiana_attenuata.zip

```

Once you have found an annotation build for your organism, carefully note the name of this organism in the snpeff database (warning: name could be different, you should always take the first name).

```sh
Database_name                                         	Species_name                                         	          	                              	https://snpeff.blob.core.windows.net/databases/v5_0/.....zip

```

Once you have the database name of your organism, use it to edit the dockerfile by changing this part:
```docker
RUN cd snpEff && \
	java -jar snpEff.jar download Arabidopsis_thaliana && \
	java -jar snpEff.jar download Physcomitrella_patens && \
	java -jar snpEff.jar download Caenorhabditis_elegans && \
	java -jar snpEff.jar download Caenorhabditis_briggsae &&\
	java -jar snpEff.jar download Populus_trichocarpa && \
	java -jar snpEff.jar download Saccharomyces_cerevisiae && \
	java -jar snpEff.jar download Zea_mays && \
	java -jar snpEff.jar download Drosophila_melanogaster &&\
	java -jar snpEff.jar download Schizosaccharomyces_pombe
```
to this ( add or replace depending on the final size of the snpeff docker you would like to have): 

```docker
RUN cd snpEff && \
	java -jar snpEff.jar download Nicotiana_attenuata
```

Save and create the docker with this command line (you should be in the docker directory) `docker build. -T snpeff: latest `


Finally, you can create your dockerhub account and push your container to dockerhub. (Use Docker hub [documentation] (https://docs.docker.com/docker-hub/))
Modify configuration file of your pipeline:
Change this `container = "romudock/snpeff: latest"` to this `container = "YOURACCOUNTNAME/snpeff: latest"` (use find ctrl+F and replace)
