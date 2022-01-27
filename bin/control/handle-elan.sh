#!/usr/bin/bash

DATESTAMP=$1
while read var; do
      [ -z "${!var}" ] && { echo 'Global Elan variable '$var' is empty or not set. Environment has not been initialised from Eagle Owl, or you are trying to call handle-elan manually which is a bad idea. Aborting.'; exit 64; }
done << EOF
DATESTAMP
ELAN_STEP1_STDOUTERR
EOF

RAISE=0
if grep --quiet 'Batch job submission failed' $ELAN_STEP1_STDOUTERR; then
    MSG='{"text":"*COG-UK inbound pipeline* Automatically re-raising Elan after a BLURM error (Batch job submission failed)"}'
    RAISE=1
elif grep --quiet 'terminated for an unknown reason' $ELAN_STEP1_STDOUTERR; then
    MSG='{"text":"*COG-UK inbound pipeline* Automatically re-raising Elan after a BLURM error (Terminated for an unknown reason)"}'
    RAISE=1
else
    MSG='{"text":"<!channel> *COG-UK inbound pipeline failed* and could not be automatically re-raised."}'
fi

sleep 30 # breathe

# kill any remaining jobs that might be stranded
#squeue | awk '$3 ~ /^elan-/ { print $1 }' | xargs scancel

curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK

if [ $RAISE -ne 0 ]; then
    sleep 300 # give java some more time to calm down
    $ELAN_SOFTWARE_DIR/bin/control/go-full-elan.sh $DATESTAMP
fi
