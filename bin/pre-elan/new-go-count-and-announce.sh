#!/usr/bin/bash

echo "[ELAN]" `date` " - LETS ROLL"

set -euo pipefail

# Credentials
source ~/.ocarina

# Pull down the entire sequencing run manifest
/rds/homes/n/nicholsz/.conda/envs/samstudio8/bin/ocarina --quiet --env get sequencing --run-name '*' --tsv --task-wait-minutes 2 --task-wait > latest.tsv
# and link it to the file system

set +o pipefail
find /cephfs/covid/bham/*/upload -type f -name "*fa*" | grep -v '\.fai$' | python /cephfs/covid/software/sam/elan/bin/ocarina_resolve.py latest.tsv > q 2> t
set -o pipefail

COUNT_MAJORA=`wc -l latest.tsv | cut -f1 -d' '`
COUNT_ELAN=`wc -l q | cut -f1 -d' '`
COUNT_ELAN_NEW=`grep -c '^1' q`
SITE_COUNTS=`awk '$14=="SANG" {print $14 " ("$13")"; next}; {print $14}' q | sort | uniq -c | sort -nr`
SITE_COUNTS_NEW=`grep '^1' q | awk '$14=="SANG" {print $14 " ("$13")"; next}; {print $14}' | sort | uniq -c | sort -nr`

SITE_MISSING_FILE=`grep 'ORPHAN-SITE' t | awk '{print $6 " " $2}' | sort -k2nr`
FILE_MISSING_META=`grep 'ORPHAN-DIRX' t | awk '$2 > 2 {print $2,$8}' | sort -nr | column -t`

###############################################################################
PRE='{"text":"<!channel>

*COG-UK inbound-distribution pre-pipeline report*
'$COUNT_ELAN_NEW' new sequences this week

***
*Samples with metadata but missing uploaded sequences on CLIMB, by sequencing centre*
Please check your upload directories...'"\`\`\`${SITE_MISSING_FILE}\`\`\`"'

*Uploaded sequences missing metadata by secondary directory*
These directories contain one or more directories with samples that do not have metadata.
Please check you have uploaded all your metadata this week...'"\`\`\`${FILE_MISSING_META}\`\`\`"'
***

*New sequences by centre*'"\`\`\`${SITE_COUNTS_NEW}\`\`\`"'

***
_Files and metadata uploaded in the past 10-15 minutes may not have been seen in time to be counted here._
_The inbound pipeline will be run autonomously at ten minutes past the next hour. Files and metadata must be uploaded before one minute past the next hour._
_Not even Sam can stop the pipeline now..._"}'
###############################################################################
POST='{"text":"<!channel>

*COG-UK inbound pipeline ready*
'$COUNT_MAJORA' sample sequencing experiments in Majora
'$COUNT_ELAN_NEW' new sequences this week
'$COUNT_ELAN' sequences matched to Majora metadata

***

*New sequences by centre*'"\`\`\`${SITE_COUNTS_NEW}\`\`\`"'

*Cumulative uploaded sequences by centre*'"\`\`\`${SITE_COUNTS}\`\`\`"'

_The pipeline will start in a few minutes. Happy Friday!_"}'
###############################################################################

# Announce
curl -X POST -H 'Content-type: application/json' --data "${!1}" "${!2}"
