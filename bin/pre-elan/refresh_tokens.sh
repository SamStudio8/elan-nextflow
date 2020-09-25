#!/usr/bin/bash
eval "$(conda shell.bash hook)"

source ~/.ocarina
conda activate samstudio8

ocarina --env oauth refresh
