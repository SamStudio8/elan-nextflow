#!/usr/bin/bash

eval "$(conda shell.bash hook)"
conda activate samstudio8

set -euo pipefail

echo $1

# Get last successful pipe date based on latest symlink
LAST_DIR_NAME=`readlink $COG_PUBLISHED_DIR/latest`
LAST_DIR_DATE=`basename $LAST_DIR_NAME`
LAST_DATE=`date -d $LAST_DIR_DATE '+%Y-%m-%d'`
echo "[CPUB] LAST_DATE=$LAST_DATE"

source ~/.ocarina

# Get files that pass QC since last pipe
ocarina --quiet --env get pagfiles --test-name 'cog-uk-elan-minimal-qc' --pass --published-after $LAST_DATE --task-wait --task-wait-attempts 15 --task-wait-minutes 1 > elan.pass.latest
cp elan.pass.latest $ELAN_DIR/staging/summary/$1/

# Get files that were suppressed and withdrawn since last pipe
ocarina --quiet --env get pagfiles --test-name 'cog-uk-elan-minimal-qc' --pass --suppressed-after $LAST_DATE --task-wait --task-wait-attempts 15 --task-wait-minutes 1 > elan.kill.latest
cp elan.kill.latest $ELAN_DIR/staging/summary/$1/


# Files to add
grep 'consensus' elan.pass.latest > elan.pass.latest.consensus.ls
cut -f3 elan.pass.latest.consensus.ls > pass.fasta.ls
wc -l pass.fasta.ls

grep 'alignment' elan.pass.latest > elan.pass.latest.alignment.ls
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
mkdir $COG_PUBLISHED_DIR/$1/
mkdir $COG_PUBLISHED_DIR/$1/fasta
mkdir $COG_PUBLISHED_DIR/$1/alignment
mkdir $COG_PUBLISHED_DIR/$1/summary
chmod 700 $COG_PUBLISHED_DIR/$1/ # Initially prevent users accessing this directory until complete
chmod 755 $COG_PUBLISHED_DIR/$1/fasta
chmod 755 $COG_PUBLISHED_DIR/$1/alignment
chmod 755 $COG_PUBLISHED_DIR/$1/summary

# Copy the links of the last alignment and fasta to today (faster than the previous individual linking)
cp -r $COG_PUBLISHED_DIR/latest/alignment/ $COG_PUBLISHED_DIR/$1/
cp -r $COG_PUBLISHED_DIR/latest/fasta/ $COG_PUBLISHED_DIR/$1/


# Linky
# Use -f force in case a late publishing pipeline from the previous day leaves
# some PAGs published today (after midnight)
echo "[CPUB]" `date` " - Linking new FASTA"
for fas in `cat pass.fasta.ls`;
do
    ln -sf $fas $COG_PUBLISHED_DIR/$1/fasta/
done

echo "[CPUB]" `date` " - Linking new BAM"
for bam in `cat pass.bam.ls`;
do
    ln -sf $bam $COG_PUBLISHED_DIR/$1/alignment/
    ln -sf $bam.bai $COG_PUBLISHED_DIR/$1/alignment/
done


# Unlinky
# Ignore unlinking errors as yesterday's unlinks may still be included
echo "[CPUB]" `date` " - Unlinking suppressed FASTA"
for fas in `cat kill.fasta.ls`;
do
    base=`basename $fas`
    set +e
    unlink $COG_PUBLISHED_DIR/$1/fasta/$base
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
    unlink $COG_PUBLISHED_DIR/$1/alignment/$bam
    ret=$?
    unlink $COG_PUBLISHED_DIR/$1/alignment/$bam.bai
    set -e

    if [ $ret -eq 0 ]; then
        echo "[KILL][BAM] $base"
    fi
done


# Merge FASTA into one location
echo "[CPUB]" `date` " - Merging consensus"
find $COG_PUBLISHED_DIR/$1/fasta/ -name '*.fasta'| xargs cat > $COG_PUBLISHED_DIR/$1/elan.$1.consensus.fasta
chmod 644 $COG_PUBLISHED_DIR/$1/elan.$1.consensus.fasta
ls -lah $COG_PUBLISHED_DIR/$1/elan.$1.consensus.fasta

# Make metadata available
cp $ELAN_DIR/staging/summary/$1/majora.metadata.tsv $COG_PUBLISHED_DIR/$1/majora.$1.metadata.tsv
chmod 644 $COG_PUBLISHED_DIR/$1/majora.$1.metadata.tsv
ls -lah $COG_PUBLISHED_DIR/$1/majora.$1.metadata.tsv

# Make other summary information available
cp $ELAN_DIR/staging/summary/$1/elan.quickcheck.ls $COG_PUBLISHED_DIR/$1/summary/elan.quickcheck.all.ls
grep -v '^0' $ELAN_DIR/staging/summary/$1/elan.quickcheck.ls > $COG_PUBLISHED_DIR/$1/summary/elan.quickcheck.bad.ls
chmod 644 $COG_PUBLISHED_DIR/$1/summary/*

# An easier to use consensus and metadata table (samstudio8/majora/27)
# NOTE samstudio8/2021-01-28
#      Send the reconcile job to SLURM where the I/O is faster to save some time
#      and improve consistency. --wait will block until complete and the script is
#      set to fail on failure. We'll need to wrap this up to ensure it runs in future.
echo "[CPUB]" `date` " - Reconciling consensus (SLURM)"
DATESTAMP=$1
sbatch --export=ALL --wait test_wait.sjob $ELAN_SOFTWARE_DIR/bin/control/reconcile_downstream.sjob
chmod 644 $COG_PUBLISHED_DIR/$1/majora.$1.metadata.matched.tsv
chmod 644 $COG_PUBLISHED_DIR/$1/elan.$1.consensus.matched.fasta

# Make QC tables available
# NOTE(samstudio8 20210107) - As all QC data is shared, just softlink the qc dir root instead of wasting time linking each qc table
echo "[CPUB]" `date` " - Linking QC"
ln -fn -s $ELAN_DIR/staging/qc $COG_PUBLISHED_DIR/$1/qc

# Repoint latest
echo "[CPUB]" `date` " - Linking latest"
chmod 755 $COG_PUBLISHED_DIR/$1/
ln -fn -s $COG_PUBLISHED_DIR/$1/ $COG_PUBLISHED_DIR/latest
ln -fn -s $COG_PUBLISHED_DIR/$1/elan.$1.consensus.fasta $COG_PUBLISHED_DIR/elan.latest.consensus.fasta
ln -fn -s $COG_PUBLISHED_DIR/$1/majora.$1.metadata.tsv $COG_PUBLISHED_DIR/majora.latest.metadata.tsv
ln -fn -s $COG_PUBLISHED_DIR/$1/elan.$1.consensus.matched.fasta $COG_PUBLISHED_DIR/elan.latest.consensus.matched.fasta
ln -fn -s $COG_PUBLISHED_DIR/$1/majora.$1.metadata.matched.tsv $COG_PUBLISHED_DIR/majora.latest.metadata.matched.tsv

# Init phylo directory
mkdir $COG_RESULTS_DIR/phylogenetics/$1
#mkdir $COG_PUBLISHED_DIR/$1/phylogenetics
#sudo setfacl -Rm g::rx $COG_PUBLISHED_DIR/$1/phylogenetics
#sudo setfacl -d -Rm g::rx $COG_PUBLISHED_DIR/$1/phylogenetics
#sudo setfacl -Rm g:phylogenetics:rwx $COG_PUBLISHED_DIR/$1/phylogenetics
#sudo setfacl -d -Rm g:phylogenetics:rwx $COG_PUBLISHED_DIR/$1/phylogenetics

# Announce summary table and failures
DASH_DATE=`date -d $1 +%Y-%m-%d`
TABLE=`ocarina -q --env get summary --md --gte-date $DASH_DATE | column -t -s'|'`

BAD_EGGS=`grep -v '^0' $ELAN_DIR/staging/summary/$1/elan.quickcheck.ls | cut -f2,3 -d' ' | sort | uniq -c | column -t -o$'\t' | sed 's,bam,bam failed samtools quickcheck,' | sed 's,fasta,fasta had short or no sequence,' | sed 's,swell,bam was aligned to wrong reference or had no alignments,' | column -t -s$'\t' | sort -nr`

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
COUNT_PASS=`ocarina --env get summary --md | awk '{sum+=$8} END {print sum}'`
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
