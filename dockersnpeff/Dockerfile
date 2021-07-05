FROM ubuntu:20.04
LABEL maintener="romudock"

WORKDIR .

SHELL [ "/bin/bash", "-c" ]

RUN apt-get -y update && \
	apt-get -y install wget unzip openjdk-11-jre-headless && \ 
    apt-get -y clean

RUN	wget "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip" && \
    unzip snpEff_latest_core.zip
    
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

RUN cp -r snpEff/scripts/* snpEff/exec/.
 
ENV PATH=$PATH:/snpEff/exec/