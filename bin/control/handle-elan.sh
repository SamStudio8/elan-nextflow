#!/usr/bin/bash
source ~/.path
source ~/.ocarina
cd $ELAN_SOFTWARE_DIR

DATESTAMP=$1

RAISE=0
if grep --quiet 'Batch job submission failed' nf.elan.$DATESTAMP.log; then
    MSG='{"text":"*COG-UK inbound pipeline* Automatically re-raising Elan after a BLURM error (Batch job submission failed)"}'
    RAISE=1
elif grep --quiet 'terminated for an unknown reason' nf.elan.$DATESTAMP.log; then
    MSG='{"text":"*COG-UK inbound pipeline* Automatically re-raising Elan after a BLURM error (Terminated for an unknown reason)"}'
    RAISE=1
else
    MSG='{"text":"<!channel> *COG-UK inbound pipeline failed* and could not be automatically re-raised"}'
fi

curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK
if [ $RAISE -ne 0 ]; then
    sleep 300 # give java some time to calm down
    $ELAN_SOFTWARE_DIR/bin/control/go-full-elan.sh $DATESTAMP
fi
