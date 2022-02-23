#!/bin/sh

# File author(s):
#           Malik Talbi <malik.talbi.2b@gmail.com>
#
#       File contributor(s) :
#           Malik Talbi <malik.talbi.2b@gmail.com>
#           Fabrice Besnard <fabrice.besnard@ens-lyon.fr>
#
#       File maintainer(s) and contact :
#           Fabrice Besnard <fabrice.besnard@ens-lyon.fr>
#
#       RDP Lab, Signal Team, Lyon - INRAe
# ------------------------------------------------------------------------------

# Anaconda is required to create the environnement

#Execute the pipeline using bash -i 

#Environement creation with Python3.8

conda create -n var_call python=3.8

#Activation of var_call

conda activate var_call

#Packages Installation

pip install htseq #HTSeq
conda install -c bioconda picard #picard