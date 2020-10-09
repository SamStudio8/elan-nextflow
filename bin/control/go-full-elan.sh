#!/usr/bin/bash
source ~/.path
source ~/.ocarina
cd $ELAN_SOFTWARE_DIR

MSG='{"text":"*COG-UK inbound pipeline begins...*
_HERE WE GO!_"}'
curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK

DATESTAMP=`date '+%Y%m%d'`

# OCARINA_FILE only written if elan processed at least one sample
OCARINA_FILE="$ELAN_DIR/staging/summary/$DATESTAMP/ocarina.files.ls"

if [ ! -f "$OCARINA_FILE" ]; then
    $NEXTFLOW_BIN run elan.nf -c elan.config --dump $PRE_ELAN_DIR/latest.tsv --publish $ELAN_DIR --schemegit /cephfs/covid/software/sam/artic-ncov2019 --datestamp $DATESTAMP > nf.elan.$DATESTAMP.log 2>&1;
    ret=$?
    mv .nextflow.log elan.nextflow.log

    if [ $ret -ne 0 ]; then
        lines=`tail -n 25 nf.elan.$DATESTAMP.log`
    else
        lines=`awk -vRS= 'END{print}' nf.elan.$DATESTAMP.log`
    fi

    MSG='{"text":"*COG-UK inbound pipeline finished...*
...with exit status '"$ret"'
'"\`\`\`${lines}\`\`\`"'
_Have a nice day!_"}'
    curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK

    if [ $ret -ne 0 ]; then
        MSG='{"text":"<!channel> *COG-UK inbound pipeline failed...*"}'
        curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK
        exit $ret
    fi
else
    MSG='{"text":"*COG-UK inbound pipeline* Cowardly skipping Elan as the ocarina file already exists for today"}'
    curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK
fi

# If the OCARINA_FILE has still not been written at this point, it means Elan ran but today's pipeline is empty - abort early but successfully, dont send a tael
if [ ! -f "$OCARINA_FILE" ]; then
    curl -X POST -H 'Content-type: application/json' --data '{"text":"\n*COG-UK inbound pipeline empty*\nNo new valid files today, try again tomorrow."}' $SLACK_REAL_HOOK
    exit 0
fi

$NEXTFLOW_BIN run ocarina.nf -c elan.config --manifest $OCARINA_FILE > nf.ocarina.$DATESTAMP.log 2>&1;
ret=$?
lines=`awk -vRS= 'END{print}' nf.ocarina.$DATESTAMP.log`
MSG='{"text":"*COG-UK QC pipeline finished...*
...with exit status '"$ret"'
'"\`\`\`${lines}\`\`\`"'"
}'
curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK
cat elan.nextflow.log .nextflow.log > inbound.nextflow.log

if [ $ret -ne 0 ]; then
    MSG='{"text":"<!channel> *COG-UK inbound pipeline failed...*"}'
    curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK
    exit $ret
fi

bash $ELAN_SOFTWARE_DIR/bin/control/cog-publish.sh $DATESTAMP
ret=$?
MSG='{"text":"*COG-UK publishing pipeline finished...*
...with exit status '"$ret"'"}'
curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK
if [ $ret -ne 0 ]; then
    MSG='{"text":"<!channel> *COG-UK inbound pipeline failed...*"}'
    curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK
    exit $ret
fi

mv inbound.nextflow.log $ELAN_DIR/staging/summary/$DATESTAMP/nf.elan.$DATESTAMP.log

# Scream into the COGUK/ether
eval "$(conda shell.bash hook)"
conda activate sam-ipc
python $ELAN_SOFTWARE_DIR/bin/ipc/mqtt-message.py -t 'COGUK/infrastructure/pipelines/elan/status' --attr status finished --attr date $DATESTAMP
