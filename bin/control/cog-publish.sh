#!/usr/bin/bash

source ~/.bootstrap.sh
source "$EAGLEOWL_CONF/envs.env"
source "$EAGLEOWL_CONF/paths.env"
source "$EAGLEOWL_CONF/slack.env"
source "$EAGLEOWL_CONF/service_elan.env"

# Activate env
eval "$(conda shell.bash hook)"
conda activate $CONDA_OCARINA

set -euo pipefail

echo $1

# Get last successful pipe date based on latest symlink
LAST_DIR_NAME=`readlink $COG_PUBLISHED_DIR/head`
LAST_DIR_DATE=`basename $LAST_DIR_NAME`
LAST_DATE=`date -d $LAST_DIR_DATE '+%Y-%m-%d'`
echo "[CPUB] LAST_DATE=$LAST_DATE"

LINKS_OK_FLAG="$ELAN_DIR/staging/summary/$1/publish.links.ok"

if [ ! -f "$LINKS_OK_FLAG" ]; then

    # Get files that pass QC since last pipe
    ocarina --oauth --quiet --env get pag --mode pagfiles --test-name 'cog-uk-elan-minimal-qc' --pass --published-after $LAST_DATE --task-wait --task-wait-attempts 15 --task-wait-minutes 1 > elan.pass.latest
    cp elan.pass.latest $ELAN_DIR/staging/summary/$1/

    # Get files that were suppressed and withdrawn since last pipe
    ocarina --oauth --quiet --env get pag --mode pagfiles --test-name 'cog-uk-elan-minimal-qc' --pass --suppressed-after $LAST_DATE --task-wait --task-wait-attempts 15 --task-wait-minutes 1 > elan.kill.latest
    cp elan.kill.latest $ELAN_DIR/staging/summary/$1/


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


    # Files to add (there may be none in the strange case where we are rerunning after a successful day - so protect the grep from failure)
    set +e
    grep 'consensus' elan.pass.latest > elan.pass.latest.consensus.ls
    set -e
    cut -f3 elan.pass.latest.consensus.ls > pass.fasta.ls
    wc -l pass.fasta.ls

    set +e
    grep 'alignment' elan.pass.latest > elan.pass.latest.alignment.ls
    set -e
    cut -f3 elan.pass.latest.alignment.ls > pass.bam.ls
    wc -l pass.bam.ls


    # Files to kill (there may be none - so protect the grep from failure)
    set +e
    grep 'consensus' elan.kill.latest > elan.kill.latest.consensus.ls
    set -e
    cut -f3 elan.kill.latest.consensus.ls > kill.fasta.ls
    wc -l kill.fasta.ls

    set +e
    grep 'alignment' elan.kill.latest > elan.kill.latest.alignment.ls
    set -e
    cut -f3 elan.kill.latest.alignment.ls > kill.bam.ls
    wc -l kill.bam.ls

    # mkdirs
    mkdir -p $COG_PUBLISHED_DIR/$1/
    mkdir -p $COG_PUBLISHED_DIR/$1/summary
    chmod 700 $COG_PUBLISHED_DIR/$1/ # Initially prevent users accessing this directory until complete
    chmod 755 $COG_PUBLISHED_DIR/$1/summary

    # Linky
    # Use -f force in case a late publishing pipeline from the previous day leaves
    # some PAGs published today (after midnight)
    echo "[CPUB]" `date` " - Linking new FASTA"
    for fas in `cat pass.fasta.ls`;
    do
        ln -sf $fas $COG_PUBLISHED_DIR/latest/fasta/
    done

    echo "[CPUB]" `date` " - Linking new BAM"
    for bam in `cat pass.bam.ls`;
    do
        ln -sf $bam $COG_PUBLISHED_DIR/latest/alignment/
        ln -sf $bam.bai $COG_PUBLISHED_DIR/latest/alignment/
    done


    # Unlinky
    # Ignore unlinking errors as yesterday's unlinks may still be included
    echo "[CPUB]" `date` " - Unlinking suppressed FASTA"
    for fas in `cat kill.fasta.ls`;
    do
        base=`basename $fas`
        set +e
        unlink $COG_PUBLISHED_DIR/latest/fasta/$base
        ret=$?
        set -e

        if [ $ret -eq 0 ]; then
            echo "[KILL][FAS] $base"
        fi
    done

    echo "[CPUB]" `date` " - Unlinking suppressed BAM"
    for bam in `cat kill.bam.ls`;
    do
        base=`basename $bam`
        set +e
        unlink $COG_PUBLISHED_DIR/latest/alignment/$base
        ret=$?
        unlink $COG_PUBLISHED_DIR/latest/alignment/$base.bai
        set -e

        if [ $ret -eq 0 ]; then
            echo "[KILL][BAM] $base"
        fi
    done

    touch $LINKS_OK_FLAG
else
    echo "[CPUB] Skipping linking, delete $LINKS_OK_FLAG to repeat"
fi

# Make metadata available
echo "[CPUB]" `date` " - Make metadata available"
cp $ELAN_DIR/staging/summary/$1/majora.metadata.tsv $COG_PUBLISHED_DIR/$1/majora.$1.metadata.tsv
chmod 644 $COG_PUBLISHED_DIR/$1/majora.$1.metadata.tsv
ls -lah $COG_PUBLISHED_DIR/$1/majora.$1.metadata.tsv

# Make other summary information available
echo "[CPUB]" `date` " - Make summary available"
cp $ELAN_DIR/staging/summary/$1/elan.quickcheck.ls $COG_PUBLISHED_DIR/$1/summary/elan.quickcheck.all.ls

echo "[CPUB]" `date` " - Parse quickcheck file"
set +e
# Guard against grep returning 1 in the case of no failed files
grep -v '^0' $ELAN_DIR/staging/summary/$1/elan.quickcheck.ls > $COG_PUBLISHED_DIR/$1/summary/elan.quickcheck.bad.ls
ret=$?
set -e
BAD_EGGS=''
if [ $ret -eq 0 ]; then
    BAD_EGGS=`grep -v '^0' $ELAN_DIR/staging/summary/$1/elan.quickcheck.ls | cut -f2,3 -d' ' | sort | uniq -c | column -t -o$'\t' | sed 's,bam,bam failed samtools quickcheck,' | sed 's,fasta,fasta was truncated or contained non IUPAC characters,' | sed 's,swell,bam was aligned to wrong reference or had no alignments,' | column -t -s$'\t' | sort -nr`
elif [ $ret -eq 1 ]; then
    # No bad eggs
    BAD_EGGS=''
elif [ $ret -gt 1 ]; then
    # POSIX says grep will exit >1 if there is a true error
    exit $ret
fi

chmod 644 $COG_PUBLISHED_DIR/$1/summary/*

# An easier to use consensus and metadata table (samstudio8/majora/27)
# NOTE samstudio8/2021-01-28
#      Send the reconcile job to SLURM where the I/O is faster to save some time
#      and improve consistency. --wait will block until complete and the script is
#      set to fail on failure.
# NOTE samstudio8/2021-01-30
#      `until` will resubmit the reconcile job until it exits 0
#      Hopefully pizza night will not be ruined by NODE_FAIL bullshit again
RECONCILE_OK_FLAG="$ELAN_DIR/staging/summary/$1/publish.reconcile.ok"
if [ ! -f "$RECONCILE_OK_FLAG" ]; then

    echo "[CPUB]" `date` " - Reconciling consensus"

    if [ "$COG_PUBLISH_MODE" == "slurm" ]; then
        until sbatch --export=ELAN_SOFTWARE_DIR=$ELAN_SOFTWARE_DIR,COG_PUBLISHED_DIR=$COG_PUBLISHED_DIR,DATESTAMP=$1,ELAN_DIR=$ELAN_DIR -o $COG_PUBLISHED_DIR/$1/summary/epubrcn-slurm-%j.out --wait $ELAN_SOFTWARE_DIR/bin/control/reconcile_downstream.sjob
        do
            ret=$?
            echo "[CPUB]" `date` " - Reconciling consensus (SLURM) - Last exit $ret"
            MSG='{"text":"*COG-UK inbound pipeline* Restarting publish reconcile"}'
            curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK
            sleep 60
        done
    else
        export DATESTAMP=$1
        until bash $ELAN_SOFTWARE_DIR/bin/control/reconcile_downstream.sjob 2> $ELAN_DIR/staging/summary/$DATESTAMP/reconcile.log
        do
            ret=$?
            echo "[CPUB]" `date` " - Reconciling consensus (LOCAL) - Last exit $ret"
            MSG='{"text":"*COG-UK inbound pipeline* Restarting publish reconcile"}'
            curl -X POST -H 'Content-type: application/json' --data "$MSG" $SLACK_MGMT_HOOK
            sleep 60
        done
    fi
    touch $RECONCILE_OK_FLAG
else
    echo "[CPUB] Skipping reconcile, delete $RECONCILE_OK_FLAG to repeat"
fi

chmod 644 $COG_PUBLISHED_DIR/$1/majora.$1.metadata.matched.tsv
chmod 644 $COG_PUBLISHED_DIR/$1/elan.$1.consensus.matched.fasta

# NOTE samstudio8/2021-01-28
#      I've scrapped the consensus merging (cat) step as we're basically double
#      handling data to generate the FASTA. This symlinks the matched FASTA to
#      replace the "unmatched" FASTA. The files will be the same going forward.
ln -fn -s $COG_PUBLISHED_DIR/$1/elan.$1.consensus.matched.fasta $COG_PUBLISHED_DIR/$1/elan.$1.consensus.fasta

# Make QC tables available
# NOTE(samstudio8 20210107) - As all QC data is shared, just softlink the qc dir root instead of wasting time linking each qc table
# NOTE(samstudio8 20210128) - No need to do this anymore as we have canonical latest/ dir

# NOTE samstudio8/2021-11-18
# Start indexing the daily consensus so we can use it with the caffeine cat reconciler
echo "[CPUB]" `date` " - Indexing latest"
samtools faidx $COG_PUBLISHED_DIR/$1/elan.$1.consensus.fasta

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
_Please refer to '"\`$COG_PUBLISHED_DIR/$1/summary/elan.quickcheck.bad.ls\`"' to identify the specific sequences for your organisation._
_These errors will appear every day, forever, until the data in question has been corrected or removed._
"}'
curl -X POST -H 'Content-type: application/json' --data "$POST" $SLACK_REAL_HOOK

# Final summary
COUNT_PASS=`ocarina --oauth --env get summary --md | awk '{sum+=$8} END {print sum}'`
COUNT_NEW=`wc -l $ELAN_DIR/staging/summary/$1/swell.qc.tsv | cut -f1 -d' '`
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
Artifacts successfully published by elan-nextflow to `'$COG_PUBLISHED_DIR'/latest`"
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
