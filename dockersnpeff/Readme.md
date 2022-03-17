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

Once you have found an annotation build for your organism, carefully note the name of this organism in the snpEff database (warning: name could be different, you should always take the first name).

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

Provided that Docker (tested for >v20.10.3) is installed on your system, build and save the docker image with this command line:
```
docker build -t snpeff:latest . # should be run in this 'dockersnpeff' directory, so that by default the Dockerfile is sourced here
```
	! *note 1*: do not worry if the terminal prints "ERROR while connecting to ..." the species database: the database is correctly downloaded however.
	! *note 2*:* verify that the image has the right content by mounting the container and exploring its content interactively:
		`docker image list` #your new image shoulde appear in the list
		`docker run -it snpeff:latest bash` # if needed, replace snpeff:latest by correct image_name:tag. Your prompt should have changed with something like "root@9fa82673eb90:/myapp#", indciating that you are now in the image container
		`java -jar snpEff/snpEff.jar -version` # verify that snpEff command is working inside the container
		`ls snpEff/data` # With the database downloads listed above, each species should have a dedicated folder, so the result should be: Arabidopsis_thaliana     Caenorhabditis_elegans   Physcomitrella_patens  Saccharomyces_cerevisiae   Zea_mays
Caenorhabditis_briggsae  Drosophila_melanogaster  Populus_trichocarpa    Schizosaccharomyces_pombe. Each folder contains a unique file called "snpEffectPredictor.bin"

To share and distribute you docker image, push your new docker image to dockerhub. You need a dockerhub account (Use Docker hub [documentation] (https://docs.docker.com/docker-hub/)). Then simply push your image on the distant dockerhub server, e.g.:
	`docker push yourdockerhub_userid/snpeff:yourtag`
	! *note 1*: the names of your local and your distant repositories should match. You can change the name you give to the local repository at the build step with `docker tag local-image:tagname yourdockerhub_userid/snpeff:tagname`
	! *note 2*: to be able to push, ensure you are connected to dockerhub: `docker login -u user_id -p <password>` (Use --password-stdin instead of -p to be more secure)


Finally, edit the relevant profiles in nextflow configuration file of your 'Parallel-Variant-Calling' pipeline (in script/nextflow.config), to specify when you would like to use this new docker for snpEff:
Change this `container = "dockerhub_previous_userid/snpeff:latest"` to this `container = "YOURACCOUNTNAME/snpeff:latest"` (use find ctrl+F and replace)
