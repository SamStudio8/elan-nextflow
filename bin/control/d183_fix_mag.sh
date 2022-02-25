#!/usr/bin/env bash
# Garbage script for a garbage plan to touch the new-style dated ARTIFACTS_ROOT
# directories to avoid a potential race condition for the creation of [D|M]AG
# (DigitalArtifactGroup or MajoraArtifactGroups) in Majora.
# See https://github.com/COG-UK/dipi-group/issues/183

# $1 DATESTAMP
DATESTAMP=$1

while read var; do
      [ -z "${!var}" ] && { echo 'Global Elan variable '$var' is empty or not set. Environment likely uninitialised. Aborting.'; exit 64; }
done << EOF
DATESTAMP
ARTIFACTS_ROOT
CONDA_OCARINA
ELAN_DAY_LOG_DIR
OCARINA_CONF_FILE
EOF

eval "$(conda shell.bash hook)"
conda activate $CONDA_OCARINA

set -euo pipefail

# Touch some pointless files
# NOTE relies on Elan creating these dirs first
touch $ARTIFACTS_ROOT/fasta/$DATESTAMP/.hoot
touch $ARTIFACTS_ROOT/bam/$DATESTAMP/.hoot

# Use Ocarina to tell Majora about them
# This will ensure a single call is emitted to add $ARTIFACTS_ROOT/{fasta,bam}/$DATESTAMP to Majora
# and avoid ocarina.nf falling over itself to add as many versions of the MAG as many times as it can
# in the opening seconds of Saturday
ocarina --oauth --angry --profile $OCARINA_PROFILE put file --path "$ARTIFACTS_ROOT/fasta/$DATESTAMP/.hoot" --type file --i-have-bad-files --full-path --node climb --pipeline-hook "cog-publish-$DATESTAMP";
ocarina --oauth --angry --profile $OCARINA_PROFILE put file --path "$ARTIFACTS_ROOT/bam/$DATESTAMP/.hoot" --type file --i-have-bad-files --full-path --node climb --pipeline-hook "cog-publish-$DATESTAMP";
touch "$ELAN_DAY_LOG_DIR/publish.d183mag.ok"
