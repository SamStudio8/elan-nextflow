#!/usr/bin/bash

echo "[ELAN]" `date` " - LETS ROLL"

set -euo pipefail
DATESTAMP=`date '+%Y%m%d'`

# Credentials
source ~/.ocarina

# Pull down the entire sequencing run manifest
/rds/homes/n/nicholsz/.conda/envs/samstudio8/bin/ocarina --quiet --env get sequencing --run-name '*' --tsv --task-wait-minutes 3 --task-wait > latest.tsv
# and link it to the file system

set +o pipefail
find /cephfs/covid/bham/*/upload -type f -name "*fa*" | grep -v '\.fai$' | python /cephfs/covid/software/sam/elan/bin/ocarina_resolve.py latest.tsv > q 2> t
set -o pipefail
cp t $COG_PUBLISHED_DIR/elan/$DATESTAMP.missing.ls
chmod 644 $COG_PUBLISHED_DIR/elan/$DATESTAMP.missing.ls

COUNT_MAJORA=`wc -l latest.tsv | cut -f1 -d' '`
COUNT_ELAN_NEW=`grep -c '^1' q`
COUNT_ELAN_OLD=`/rds/homes/n/nicholsz/.conda/envs/samstudio8/bin/ocarina --env get summary --md | awk '{sum+=$6} END {print sum}'`
COUNT_ELAN_OLDANDNEW=`expr $COUNT_ELAN_NEW + $COUNT_ELAN_OLD`
SITE_COUNTS=`awk '$14=="SANG" {print $14 " ("$13")"; next}; {print $14}' q | sort | uniq -c | sort -nr`
SITE_COUNTS_NEW=`grep '^1' q | awk '$14=="SANG" {print $14 " ("$13")"; next}; {print $14}' | sort | uniq -c | sort -nr`

SITE_MISSING_FILE=`grep 'ORPHAN-SITE' t | awk '{print $6 " " $2}' | sort -k2nr`
FILE_MISSING_META=`grep 'ORPHAN-DIRX' t | grep -v jacksond | awk '$2 > 2 {print $2,$8}' | sort -nr | column -t`
SANG_MISSING_META=`grep 'ORPHAN-FILE' t | grep jacksond | cut -f2 -d' ' | awk -F'/' '{print $(NF-1)}' | cut -c1-4 | tr -d '[0-9]_' | sort | uniq -c | sort -nr`

###############################################################################
PRE='{"text":"<!channel>

*COG-UK inbound-distribution pre-pipeline report*
'$COUNT_ELAN_NEW' new sequences today

***
*Samples with metadata but missing uploaded sequences on CLIMB, by sequencing centre*
Please check your upload directories...'"\`\`\`${SITE_MISSING_FILE}\`\`\`"'
To inspect the barcodes with metadata but missing an uploaded sequence from your site, execute:
```grep ORPHAN-COGX '$COG_PUBLISHED_DIR'/elan/'$DATESTAMP'.missing.ls | grep '"'"'\\[BIRM\\]'"'"'```
_Replace BIRM with your site code from the table above. Ensure to maintain the brackets and quotes._

***
*Uploaded sequences missing metadata by secondary directory*
These directories contain one or more directories with samples that do not have metadata.
Please check you have uploaded all your metadata...'"\`\`\`${FILE_MISSING_META}\`\`\`"'
To inspect the uploaded sequences that are missing metadata from your site, execute:
```grep ORPHAN-FILE '$COG_PUBLISHED_DIR'/elan/'$DATESTAMP'.missing.ls | grep climb-covid19-nichollss```
_Replace climb-covid19-nichollss with the username of the uploader from the table above_

***
*Samples sequenced by Sanger missing biosample metadata*
Metadata is missing for samples submitted from these sites, for sequencing at Sanger.
Local sites should ensure they have uploaded the biosample-only metadata for samples not sequenced locally.'"\`\`\`${SANG_MISSING_META}\`\`\`"'
To inspect the sequences uploaded from Sanger that are missing metadata from your site, execute:
```grep ORPHAN-FILE '$COG_PUBLISHED_DIR'/elan/'$DATESTAMP'.missing.ls | grep climb-covid19-jacksond | grep BIRM```
_Replace BIRM with your barcode prefix from the table above, without changing the username._
***

*New sequences by centre*'"\`\`\`${SITE_COUNTS_NEW}\`\`\`"'

***
_Files and metadata uploaded in the past 10-15 minutes may not have been seen in time to be counted here._
_Please try and help yourself with the commands above to self-diagnose basic missing data issues before contacting Sam._
_Please note that Sam will not be able to assist in queries to disambiguate barcode prefixes that are shared across sites (e.g. LOND)._
_The inbound pipeline will be run autonomously at ten minutes past the next hour._
_Files and metadata must be uploaded before one minute past the next hour. Not even Sam can stop the pipeline now..._"}'
###############################################################################
POST='{"text":"<!channel>

*COG-UK inbound pipeline ready*
'$COUNT_MAJORA' sample sequencing experiments in Majora
'$COUNT_ELAN_NEW' new sequences today
'$COUNT_ELAN_OLDANDNEW' sequences matched to Majora metadata

***

*New sequences by centre*'"\`\`\`${SITE_COUNTS_NEW}\`\`\`"'

_The pipeline will start in a few minutes. Happy Friday!_"}'
###############################################################################

# Announce
curl -X POST -H 'Content-type: application/json' --data "${!1}" "${!2}"
