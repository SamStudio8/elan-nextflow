#!/usr/bin/bash
source ~/.path
source ~/.ocarina

# conda env and credentials
eval "$(conda shell.bash hook)"
conda activate samstudio8

set -euo pipefail
DATESTAMP=`date '+%Y%m%d'`
echo "[ELAN]" `date` " - LETS ROLL"

# Pull down the entire sequencing run manifest as JSON over the new v3 API
#ocarina --quiet --env --oauth get dataview --mdv ELAN1 -o elan1.mdv.json --task-wait --task-wait-attempts 60
#set +o pipefail
#find /cephfs/covid/bham/*/upload -type f -name "*fa*" | grep -v '\.fai$' | python /cephfs/covid/software/sam/elan/bin/ocarina_resolve_mdv.py elan1.mdv.json latest.tsv > q 2> t
#set -o pipefail

# Keep the v2 requests for now
ocarina --oauth --quiet --env get sequencing --run-name '*' --faster --tsv --task-wait-attempts 75 --task-wait > latest.$1.tsv

# Map filesystem to metadata to find work to do
set +o pipefail
echo "[ELAN]" `date` " - Find files"
find /cephfs/covid/bham/*/upload -type f -name "*fa*" | grep -v '\.fai$' | python /cephfs/covid/software/sam/elan/bin/ocarina_resolve.py latest.$1.tsv > q.$1.tsv 2> t.$1.txt
set -o pipefail

# Save manifest at the end of the day (SHORT message fires last)
if [ "$1" = "SHORT" ]; then
    set +o pipefail
    grep '^1' q.$1.tsv | awk '{print $8"/"$10"\n"$8"/"$12}' > synclist/sync.$DATESTAMP.ls
    set -o pipefail
    SYNCPATH=$(readlink -f synclist/sync.$DATESTAMP.ls)
    python $ELAN_SOFTWARE_DIR/bin/ipc/mqtt-message.py -t 'COGUK/infrastructure/housekeeping/synclisting/status' --host $MQTT_HOST --attr status finished --attr date $DATESTAMP --attr synclist $SYNCPATH
fi

cp t $COG_PUBLISHED_DIR/elan/$DATESTAMP.missing.ls
chmod 644 $COG_PUBLISHED_DIR/elan/$DATESTAMP.missing.ls

ocarina --quiet --oauth --env get summary --md > summary.$1.md
$ELAN_SOFTWARE_DIR/bin/message_uploads.sh latest.$1.tsv q.$1.tsv t.$1.txt $1 "${!2}" summary.$1.md
