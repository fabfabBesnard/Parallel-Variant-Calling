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

# Anaconda is required to install the environnement

# Creation of the environnement using the requirements.yml file

conda env create -f requirements.yml

# Activation of the environnement

conda init bash 
conda activate var_call

# Check if HTSeq and numpy versions are correctly installed

python3 -m pip install htseq
python3 -m pip install --upgrade numpy