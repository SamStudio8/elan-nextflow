#!/usr/bin/bash
source ~/.path
source ~/.ocarina
cd $ELAN_SOFTWARE_DIR

DATESTAMP=$1

if grep --quiet 'Batch job submission failed' nf.elan.$DATESTAMP.log; then
    MSG='{"text":"*COG-UK inbound pipeline* Automatically re-raising Elan after a BLURM error (Batch job submission failed)"}'
    curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK
    $ELAN_SOFTWARE_DIR/bin/control/go-full-elan.sh
elif grep --quiet 'terminated for an unknown reason' nf.elan.$DATESTAMP.log; then
    MSG='{"text":"*COG-UK inbound pipeline* Automatically re-raising Elan after a BLURM error (Terminated for an unknown reason)"}'
    curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK
    $ELAN_SOFTWARE_DIR/bin/control/go-full-elan.sh
fi
