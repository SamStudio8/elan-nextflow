#!/usr/bin/bash
# go-full-elan is composed of three steps
# 	* elan-nextflow (nextflow)  processes daily uploads for QC and kicks out anything very terrible
#	* elan-ocarina  (nextflow)  updates Majora en masse for all files from step 1
#	* cog-publish   (bash)      publishes QC pass files to users
#
# go-full-elan takes the following inputs
# 	$1    DATESTAMP    elan datestamp, usually today (format YYYYMMDD)

while read var; do
      [ -z "${!var}" ] && { echo 'Global Eagle Owl variable '$var' is empty or not set. Environment likely uninitialised. Aborting.'; exit 64; }
done << EOF
EAGLEOWL_CONF
EOF

# Init the Elan environment from the Eagle Owl config dir
#source "$EAGLEOWL_CONF/elan/dev.env"
source "$EAGLEOWL_CONF/elan/prod.env"

DATESTAMP=$1
while read var; do
      [ -z "${!var}" ] && { echo 'Global Elan variable '$var' is empty or not set. Environment has not been initialised from Eagle Owl. Aborting.'; exit 64; }
done << EOF
DATESTAMP
UPLOADS_DIR_GLOB
ELAN_CONFIG
ELAN_DIR
ELAN_SOFTWARE_DIR
ELAN_RUN_DIR
NEXTFLOW_BIN
ELAN_SLACK_MGMT_HOOK
ELAN_SLACK_INBOUND_HOOK
MQTT_HOST
COG_PUBLISH_MODE
ELAN_LOG_DIR
OCARINA_CONF_FILE
OCARINA_PROFILE
ARTIFACTS_ROOT
CONDA_OCARINA
CONDA_IPC
NXF_WORK
NXF_CONDA_CACHEDIR
MQTT_ENV
EOF

INBOUND_MANIFEST=${INBOUND_MANIFEST:-}

if ! command -v conda &> /dev/null
then
    echo "conda not in PATH"
    exit 78 # EX_CONFIG
fi

if ! command -v 'mqtt-message.py' &> /dev/null
then
    echo "mqtt-message.py (tael) not in PATH"
    exit 78 # EX_CONFIG
fi

if ! command -v 'ocarina_resolve.py' &> /dev/null
then
    echo "ocarina_resolve.py (elan bin dir) not in PATH"
    exit 78 # EX_CONFIG
fi

cd $ELAN_RUN_DIR
echo $DATESTAMP

# Centralise .nextflow.log location
export ELAN_DAY_LOG_DIR="$ELAN_LOG_DIR/$DATESTAMP" # export for cog-publish
mkdir -p $ELAN_DAY_LOG_DIR

export ELAN_STEP1_NFLOG="$ELAN_DAY_LOG_DIR/nf.elan.log" # export for handle-elan
ELAN_STEP2_NFLOG="$ELAN_DAY_LOG_DIR/nf.ocarina.log"
ELAN_STEP3_LOG="$ELAN_DAY_LOG_DIR/publish.log"

MSG='{"text":"*COG-UK inbound pipeline begins...*
*Majora user* `'$OCARINA_PROFILE'`
*Uploads dir* `'$UPLOADS_DIR_GLOB'`
*Artifacts dir* `'$ARTIFACTS_ROOT'`
*Staging dir* `'$ELAN_DIR'`
*Log dir* `'$ELAN_LOG_DIR'/'$DATESTAMP'`"}'
curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK

export ELAN_STEP1_STDOUTERR="$ELAN_LOG_DIR/$DATESTAMP/nf.elan.$DATESTAMP.log" # export for handle-elan
ELAN_STEP2_STDOUTERR="$ELAN_LOG_DIR/$DATESTAMP/nf.ocarina.$DATESTAMP.log"

# OCARINA_FILE only written if elan processed at least one sample
OCARINA_FILE="$ARTIFACTS_ROOT/elan/$DATESTAMP/ocarina.files.ls"
ELAN_OK_FLAG="$ELAN_DAY_LOG_DIR/elan.ok.flag"
OCARINA_OK_FLAG="$ELAN_DAY_LOG_DIR/ocarina.ok.flag"
PUBLISH_OK_FLAG="$ELAN_DAY_LOG_DIR/publish.ok.flag"

if [ ! -f "$ELAN_OK_FLAG" ]; then
    # If a log already exists, then the pipeline needs to be resumed
    # NOTE --uploads MUST be quoted to prevent a premature shell expansion of the uploads glob
    RESUME_FLAG=""
    if [ -f "$ELAN_STEP1_STDOUTERR" ]; then
        RESUME_FLAG="-resume"
	MSG='{"text":"*COG-UK inbound pipeline* Using -resume to re-raise Elan without trashing everything. Delete today'\''s log (`rm '$ELAN_STEP1_STDOUTERR'`) to force a full restart."}'
        curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK
    fi

    INBOUND_MANIFEST_FLAG=""
    if  [ ! -z "$INBOUND_MANIFEST" ]; then
        INBOUND_MANIFEST_FLAG="--inbound_manifest $INBOUND_MANIFEST"
        MSG='{"text":"*COG-UK inbound pipeline* Using `'$INBOUND_MANIFEST_FLAG'`. Unset `INBOUND_MANIFEST` to prevent this from happening and force a full restart."}'
        curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK
    fi
    /usr/bin/flock -w 1 /dev/shm/.sam_elan -c "$NEXTFLOW_BIN -log $ELAN_STEP1_NFLOG run $ELAN_SOFTWARE_DIR -c $ELAN_CONFIG --mode inbound --ocarina_profile $OCARINA_PROFILE --artifacts_root $ARTIFACTS_ROOT --publish $ELAN_DIR --uploads \"$UPLOADS_DIR_GLOB\" --datestamp $DATESTAMP $RESUME_FLAG $INBOUND_MANIFEST_FLAG > $ELAN_STEP1_STDOUTERR 2>&1;"
    ret=$?

    if [ $ret -ne 0 ]; then
        lines=`tail -n 25 $ELAN_STEP1_STDOUTERR`
    else
        lines=`awk -vRS= 'END{print}' $ELAN_STEP1_STDOUTERR`
    fi

    MSG='{"text":"*COG-UK Elan-Inbound pipeline finished...*
...with exit status '"$ret"'
'"\`\`\`${lines}\`\`\`"'
_Have a nice day!_"}'
    curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK

    if [ $ret -ne 0 ]; then
        $ELAN_SOFTWARE_DIR/bin/control/handle-elan.sh $DATESTAMP
        exit $ret # get out of here before we loop ourselves into infinity
    else
        touch $ELAN_OK_FLAG
    fi
else
	MSG='{"text":"*COG-UK inbound pipeline* Cowardly skipping Elan as the OK flag (`'$ELAN_OK_FLAG'`) already exists for today. If elan-nextflow ran correctly (which it looks like it did), do not delete the OK flag!"}'
    curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK
fi

# If the OCARINA_FILE has still not been written at this point, it means Elan ran but today's pipeline is empty - abort early but successfully, dont send a tael
if [ ! -f "$OCARINA_FILE" ]; then
    curl -X POST -H 'Content-type: application/json' --data '{"text":"\n*COG-UK inbound pipeline empty*\nNo new valid files today, try again tomorrow."}' $ELAN_SLACK_INBOUND_HOOK
    exit 0
fi


# NOTE samstudio8 2022-01-28
# Added crappy shim to stop Ocarina blowing up until we can patch Majora properly
# See https://github.com/COG-UK/dipi-group/issues/183
$ELAN_SOFTWARE_DIR/bin/control/d183_fix_mag.sh $DATESTAMP
if [ ! -f "$ELAN_DAY_LOG_DIR/publish.d183mag.ok" ]; then
    MSG='{"text":"<!channel> *COG-UK inbound pipeline failed...* Could not create daily BAM and FASTA directory MAG in Majora. This is most likely due to a unfortunately timed OAuth rotation, try running Elan again."}'
    curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK
    exit 70
fi


if [ ! -f "$OCARINA_OK_FLAG" ]; then
    $NEXTFLOW_BIN -log $ELAN_STEP2_NFLOG run $ELAN_SOFTWARE_DIR -c $ELAN_CONFIG --mode ocarina --ocarina_profile $OCARINA_PROFILE --manifest $OCARINA_FILE > $ELAN_STEP2_STDOUTERR 2>&1;
    ret=$?
    lines=`awk -vRS= 'END{print}' $ELAN_STEP2_STDOUTERR`
    MSG='{"text":"*COG-UK Elan-Ocarina pipeline finished...*
...with exit status '"$ret"'
'"\`\`\`${lines}\`\`\`"'"
}'
    curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK

    if [ $ret -ne 0 ]; then
        MSG='{"text":"<!channel> *COG-UK inbound pipeline failed (Ocarina)*"}'
        curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK
        exit $ret
    else
        touch $OCARINA_OK_FLAG
    fi
else
    MSG='{"text":"*COG-UK inbound pipeline* Cowardly skipping Ocarina as the OK flag (`'$OCARINA_OK_FLAG'`) already exists for today. You can safely re-run elan-ocarina by deleting the OK flag, but it will re-emit all requests to Majora!"}'
    curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK
fi


if [ ! -f "$PUBLISH_OK_FLAG" ]; then
    SECONDS=0
    bash $ELAN_SOFTWARE_DIR/bin/control/cog-publish.sh $DATESTAMP > $ELAN_STEP3_LOG
    ret=$?
    TIMER=$(python -c "import datetime; print('(publish)', str(datetime.timedelta(seconds=$SECONDS)))")
    MSG='{"text":"*COG-UK Elan-Publish pipeline finished...*
...with exit status '"$ret"' in '"$TIMER"'"}'
    curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK
    if [ $ret -ne 0 ]; then
        MSG='{"text":"<!channel> *COG-UK inbound pipeline failed...*"}'
        curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK
        exit $ret
    else
        touch $PUBLISH_OK_FLAG
    fi
else
    MSG='{"text":"*COG-UK inbound pipeline* Cowardly skipping Publish as the OK flag (`'$PUBLISH_OK_FLAG'`) already exists for today. You cannot safely re-run cog-publish by deleting the OK flag, you must also use cog-publish-link to repoint the artifacts head symlink."}'
    curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK
fi


# Scream into the COGUK/ether
eval "$(conda shell.bash hook)"
conda activate $CONDA_IPC
mqtt-message.py -t "$MQTT_ENV/infrastructure/pipelines/elan/status" --host $MQTT_HOST --attr status finished --attr date $DATESTAMP
