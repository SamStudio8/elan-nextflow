#!/usr/bin/bash
source ~/.path
source ~/.ocarina

eval "$(conda shell.bash hook)"

conda activate samstudio8

ocarina --env oauth refresh
