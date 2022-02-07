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

# On charge python/3.8.3

module load Python/3.8.3

# Creation de l'environnement virtuel var_call

python3.8 -m venv ~/var_call

# Activation of the environnement

source ~/var_call/bin/activate

# Installation de htseq

python3.8 -m pip install -U HTSeq
