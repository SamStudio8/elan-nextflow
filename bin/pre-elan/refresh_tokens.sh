#!/usr/bin/bash
source ~/.bootstrap.sh
source "$EAGLEOWL_CONF/ocarina.env"

eval "$(conda shell.bash hook)"
conda activate samstudio8

ocarina --quiet --env oauth refresh
