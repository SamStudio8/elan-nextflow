#!/usr/bin/bash
source ~/.bootstrap.sh
source "$EAGLEOWL_CONF/ocarina.env"
source "$EAGLEOWL_CONF/envs.env"

eval "$(conda shell.bash hook)"
conda activate $CONDA_OCARINA

ocarina --quiet --env oauth refresh
