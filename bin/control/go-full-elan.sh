#!/usr/bin/bash
source ~/.path
source ~/.ocarina
cd /cephfs/covid/software/sam/elan

MSG='{"text":"<!channel> *COG-UK inbound pipeline begins...*
_HERE WE GO!_"}'
curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK

DATESTAMP=`date '+%Y%m%d'`
$NEXTFLOW_BIN run elan.nf -c elan.config --dump $PRE_ELAN_DIR/latest.tsv --publish $ELAN_DIR --schemegit /cephfs/covid/software/sam/artic-ncov2019 --datestamp $DATESTAMP > nf.elan.$DATESTAMP.log 2>&1;
ret=$?
lines=`awk -vRS= 'END{print}' nf.elan.$DATESTAMP.log`

MSG='{"text":"<!channel> *COG-UK inbound pipeline finished...*
...with exit status '"$ret"'
'"\`\`\`${lines}\`\`\`"'
 _Happy Friday!_"}'
curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK

if [ $ret -ne 0 ]; then
    exit $ret
fi

$NEXTFLOW_BIN run ocarina.nf -c elan.config --manifest $ELAN_DIR/staging/summary/$DATESTAMP/ocarina.files.ls > nf.ocarina.$DATESTAMP.log 2>&1;
ret=$?
lines=`awk -vRS= 'END{print}' nf.ocarina.$DATESTAMP.log`
MSG='{"text":"<!channel> *COG-UK QC pipeline finished...*
...with exit status '"$ret"'
'"\`\`\`${lines}\`\`\`"'"
}'
curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK

if [ $ret -ne 0 ]; then
    exit $ret
fi

bash $PRE_ELAN_DIR/new-publish.sh $DATESTAMP
ret=$?
MSG='{"text":"<!channel> *COG-UK publishing pipeline finished...*
...with exit status '"$ret"'"}'
curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK
if [ $ret -ne 0 ]; then
    exit $ret
fi

mv .nextflow.log $ELAN_DIR/staging/summary/$DATESTAMP/nf.elan.$DATESTAMP.log
