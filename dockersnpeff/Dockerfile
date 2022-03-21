FROM ubuntu:20.04
LABEL maintener="Fabrice Besnard <fabrice.besnard@ens-lyon.fr>"
LABEL repository="https://github.com/fabfabBesnard/Parallel-Variant-Calling/"
LABEL description="Latest snpEff program coming with a few pre-built databases for popular model organisms in biology"

# Build a context directory (e.g.: https://newbedev.com/what-is-app-working-directory-for-a-dockerfile)
WORKDIR /myapp
# Change Shell to 'bash', default is 'sh'
SHELL [ "/bin/bash", "-c" ]

RUN apt-get -y update && \
	apt-get -y install wget unzip openjdk-11-jre-headless && \ 
    apt-get -y clean

RUN	wget "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip" && \
    unzip snpEff_latest_core.zip

# Install pre-built databases    
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

# Add snpeff (or snpsift and all scripts) executable commands to the $PATH when the docker container is running
RUN cp -r /myapp/snpEff/scripts/* /myapp/snpEff/exec/. 
ENV PATH=$PATH:/myapp/snpEff/exec/