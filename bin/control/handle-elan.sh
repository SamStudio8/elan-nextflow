#!/usr/bin/bash
#source ~/.bootstrap.sh
#source "$EAGLEOWL_CONF/paths.env"
#source "$EAGLEOWL_CONF/slack.env"

DATESTAMP=$1

RAISE=0
if grep --quiet 'Batch job submission failed' nf.elan.$DATESTAMP.log; then
    MSG='{"text":"*COG-UK inbound pipeline* Automatically re-raising Elan after a BLURM error (Batch job submission failed)"}'
    RAISE=1
elif grep --quiet 'terminated for an unknown reason' nf.elan.$DATESTAMP.log; then
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
