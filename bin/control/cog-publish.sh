#!/usr/bin/bash
while read var; do
      [ -z "${!var}" ] && { echo 'Global Eagle Owl variable '$var' is empty or not set. Environment likely uninitialised. Aborting.'; exit 64; }
done << EOF
EAGLEOWL_CONF
ELAN_DAY_LOG_DIR
EOF

# Init the Elan environment from the Eagle Owl config dir
source "$EAGLEOWL_CONF/elan/dev.env"
#source "$EAGLEOWL_CONF/elan/prod.env"

# Activate env
eval "$(conda shell.bash hook)"
conda activate $CONDA_OCARINA

set -euo pipefail

echo $1

# Get last successful pipe date based on latest symlink
LAST_DIR_NAME=`readlink $ARTIFACTS_ROOT/elan/head`
LAST_DIR_DATE=`basename $LAST_DIR_NAME`

if [ "$LAST_DIR_DATE" = "$1" ]; then
    MSG='{"text":"<!channel>*COG-UK inbound pipeline* Cowardly refusing to publish, as the `LAST_DIR_DATE` is set to the same date you are trying to publish (`'$LAST_DIR_DATE'`). My guess is you are trying to re-publish the data set for today after it has already been released. I cannot let you do this as bad things will happen (https://github.com/COG-UK/dipi-group/issues/179). If you really need to re-publish the data set for today, use the cog-publish-link script to push the artifacts head back to yesterday, remove all `publish.*.ok` files and run this command again."}'
curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK
    exit 78;
fi


LAST_DATE=`date -d $LAST_DIR_DATE '+%Y-%m-%d'`
echo "[CPUB] LAST_DATE=$LAST_DATE"

PAGS_OK_FLAG="$ELAN_DAY_LOG_DIR/publish.pags.ok"

if [ ! -f "$PAGS_OK_FLAG" ]; then

    # Get files that pass QC since last pipe for reconcile
    ocarina --oauth --quiet --env get pag --mode pagfiles --test-name 'cog-uk-elan-minimal-qc' --pass --published-after $LAST_DATE --task-wait --task-wait-attempts 15 --task-wait-minutes 1 > latest_elan.pass_pag_lookup.tsv
    cp latest_elan.pass_pag_lookup.tsv $ARTIFACTS_ROOT/elan/$1/

    # Get files that were suppressed and withdrawn since last pipe for reconcile
    ocarina --oauth --quiet --env get pag --mode pagfiles --test-name 'cog-uk-elan-minimal-qc' --pass --suppressed-after $LAST_DATE --task-wait --task-wait-attempts 15 --task-wait-minutes 1 > latest_elan.kill_pag_lookup.tsv
    cp latest_elan.kill_pag_lookup.tsv $ARTIFACTS_ROOT/elan/$1/

    # NOTE 2022-01-25
    #   This is the intermediate future of presenting artifacts to users,
    #   we are taking away the assumption they will be stored in a particular
    #   place and providing a lookup for users to be able to find them. Long term
    #   I would like to see all these artifacts stored in our S3 gateway instead.
    # Make big super lookup table
    ocarina --oauth --quiet --env get pag --mode pagfiles --test-name 'cog-uk-elan-minimal-qc' --task-wait --task-wait-attempts 15 --output-header > majora.pag_lookup.tsv
    # append all suppressed files to lookup table
    ocarina --oauth --quiet --env get pag --mode pagfiles --test-name 'cog-uk-elan-minimal-qc' --task-wait --task-wait-attempts 15 --suppressed-after 1970-01-01 >> majora.pag_lookup.tsv
    # and unpack the pag for a slightly more end-user friendly lookup table
    python $ELAN_SOFTWARE_DIR/bin/control/cog-publish-unpag.py --tsv majora.pag_lookup.tsv > $ARTIFACTS_ROOT/elan/$1/majora.pag_lookup.tsv

    touch $PAGS_OK_FLAG
else
    echo "[CPUB] Skipping PAG table, delete $PAGS_OK_FLAG to repeat"
fi

echo "[CPUB]" `date` " - Parse quickcheck file"
set +e
# Guard against grep returning 1 in the case of no failed files
grep -v '^0' $ARTIFACTS_ROOT/elan/$1/elan.quickcheck.ls > $ARTIFACTS_ROOT/elan/$1/elan.quickcheck.bad.ls
ret=$?
set -e
BAD_EGGS=''
if [ $ret -eq 0 ]; then
    BAD_EGGS=`grep -v '^0' $ARTIFACTS_ROOT/elan/$1/elan.quickcheck.ls | cut -f2,3 -d' ' | sort | uniq -c | column -t -o$'\t' | sed 's,bam,bam failed samtools quickcheck,' | sed 's,fasta,fasta was truncated or contained non IUPAC characters,' | sed 's,swell,bam was aligned to wrong reference or had no alignments,' | column -t -s$'\t' | sort -nr`
elif [ $ret -eq 1 ]; then
    # No bad eggs
    BAD_EGGS=''
elif [ $ret -gt 1 ]; then
    # POSIX says grep will exit >1 if there is a true error
    exit $ret
fi

# An easier to use consensus and metadata table (samstudio8/majora/27)
# NOTE samstudio8/2021-01-28
#      Send the reconcile job to SLURM where the I/O is faster to save some time
#      and improve consistency. --wait will block until complete and the script is
#      set to fail on failure.
# NOTE samstudio8/2021-01-30
#      `until` will resubmit the reconcile job until it exits 0
#      Hopefully pizza night will not be ruined by NODE_FAIL bullshit again
RECONCILE_OK_FLAG="$ELAN_DAY_LOG_DIR/publish.reconcile.ok"
if [ ! -f "$RECONCILE_OK_FLAG" ]; then

    echo "[CPUB]" `date` " - Reconciling consensus"

    if [ "$COG_PUBLISH_MODE" == "slurm" ]; then
        until sbatch --export=ELAN_SOFTWARE_DIR=$ELAN_SOFTWARE_DIR,DATESTAMP=$1,ELAN_DIR=$ELAN_DIR,ARTIFACTS_ROOT=$ARTIFACTS_ROOT -o $ARTIFACTS_ROOT/elan/$1/epubrcn-slurm-%j.out --wait $ELAN_SOFTWARE_DIR/bin/control/reconcile_downstream.sjob
        do
            ret=$?
            echo "[CPUB]" `date` " - Reconciling consensus (SLURM) - Last exit $ret"
            MSG='{"text":"*COG-UK inbound pipeline* Restarting publish reconcile"}'
            curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK
            sleep 60
        done
    else
        export DATESTAMP=$1
        until bash $ELAN_SOFTWARE_DIR/bin/control/reconcile_downstream.sjob 2> $ELAN_DAY_LOG_DIR/reconcile.log
        do
            ret=$?
            echo "[CPUB]" `date` " - Reconciling consensus (LOCAL) - Last exit $ret"
            MSG='{"text":"*COG-UK inbound pipeline* Restarting publish reconcile"}'
            curl -X POST -H 'Content-type: application/json' --data "$MSG" $ELAN_SLACK_MGMT_HOOK
            sleep 60
        done
    fi
    touch $RECONCILE_OK_FLAG
else
    echo "[CPUB] Skipping reconcile, delete $RECONCILE_OK_FLAG to repeat"
fi

# NOTE samstudio8/2021-11-18
# Start indexing the daily consensus so we can use it with the caffeine cat reconciler
echo "[CPUB]" `date` " - Indexing latest"
samtools faidx $ARTIFACTS_ROOT/elan/$1/elan.consensus.fasta

# Repoint latest
echo "[CPUB]" `date` " - Linking latest"
bash $ELAN_SOFTWARE_DIR/bin/control/cog-publish-link.sh $1

# Announce summary table and failures
DASH_DATE=`date -d $1 +%Y-%m-%d`
TABLE=`ocarina --oauth -q --env get summary --md --gte-date $DASH_DATE | column -t -s'|'`


POST='{"text":"
*COG-UK inbound pipeline QC summary* '"\`\`\`${TABLE}\`\`\`"'

_Note that these figures are only intended as feedback to individual sites - quality control is not a contest!_
_Due to the many different extraction and sequencing techniques, as well as local policies with regard to uploading sequences that failed local QC, these QC metrics should not be used to compare performance between sequencing sites._

**
*COG-UK inbound pipeline failure summary* '"\`\`\`${BAD_EGGS}\`\`\`"'

_These sequences have failed fatally and cannot be processed by Elan._
_Please refer to '"\`$ARTIFACTS_ROOT/elan/$1/elan.quickcheck.bad.ls\`"' to identify the specific sequences for your organisation._
_These errors will appear every day, forever, until the data in question has been corrected or removed._
"}'
curl -X POST -H 'Content-type: application/json' --data "$POST" $SLACK_REAL_HOOK

# Final summary
COUNT_PASS=`ocarina --oauth --env get summary --md | awk '{sum+=$8} END {print sum}'`
COUNT_NEW=`wc -l $ARTIFACTS_ROOT/elan/$1/swell.qc.tsv | cut -f1 -d' '`
POST='{
    "attachments": [
        {
            "blocks": [
                {
                    "type": "header",
                    "text": {
                        "type": "plain_text",
                        "text": ":tada: COG-UK inbound-distribution pipeline finished",
                        "emoji": true
                    }
                },
                {
                    "type": "section",
                    "text": {
                        "type": "mrkdwn",
                        "text": "
*'$COUNT_NEW'* new sequences matched to Majora metadata today
*'$COUNT_PASS'* sequences passed basic quality control to date
Artifacts successfully published by elan-nextflow to `'$ARTIFACTS_ROOT'/elan/latest/`"
                    },
                    "accessory": {
                        "type": "image",
                        "image_url": "https://avatars.slack-edge.com/2019-05-03/627972616934_a621b7d3a28c2b6a7bd1_512.jpg",
                        "alt_text": "Majora is watching."
                    }
                },
                {
                    "type": "context",
                    "elements": [
                        {
                            "type": "mrkdwn",
                            "text":"
_QC Reports have been calculated and can be reached from your Majora profile._
_The outbound distribution pipeline will run next Monday._
_Thanks for your patience, have a nice day!_"
                        }
                    ]
                }
            ],
            "color": "#36C5F0",
        }
    ]
}'
curl -X POST -H 'Content-type: application/json' --data "$POST" $SLACK_REAL_HOOK
