# How to update snpeff docker
Snpeff provides a large number of annotations directly included in the database. Unfortunately, it is not possible to include them all in the docker. In case of you do not have your model organism annotation in the genomes already included in the docker and your annotation file is badly formatted and causes an error, you can try to create your own docker.
The first step is to find if your organism is included in the
The first step is to find if your ogranism are included in the snpeff database. For this download [snpeff](https://pcingola.github.io/SnpEff/download/) 

One you have installed snpeff, do this command line to find out if snpeff contain your model organism : 

`java -jar snpEff.jar databases | grep "Myspecies" `

For exemple : 
I want to run the pipeline with Nicotiana attenuata as model organism

```sh
Linux:~/snpEff$ java -jar snpEff.jar databases | grep "Nicotina"
Nicotiana_attenuata                                         	Nicotiana_attenuata                                         	          	                              	https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_Nicotiana_attenuata.zip

```

Once you have found your model organism, take the name of your organism in the database (warning name could be different, you should always take the first name).

```sh
Database_name                                         	Species_name                                         	          	                              	https://snpeff.blob.core.windows.net/databases/v5_0/.....zip

```

When you have the databse name of your organism : you can write dockerfile. Go to the dockerfile and change this part : 
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
to this 

```docker
RUN cd snpEff && \
	java -jar snpEff.jar download Nicotiana_attenuata
```

Save and create the docker with this command line (you should be in the docker directory) `docker build. -T snpeff: latest `


And then you can create your dockerhub account and push your container to dockerhub. (Use Docker hub [documentation] (https://docs.docker.com/docker-hub/))
Modify configuration file of your pipeline:
Change this `container = "romudock/snpeff: latest"` to this `container = "YOURACCOUNTNAME/snpeff: latest"` (use find ctrl+F and replace)
