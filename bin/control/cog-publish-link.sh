set -euo pipefail

# Repoint latest
# This script is usually run automatically by cog-publish
# Under exceptional circumstances one may need to repoint the latest data set manually
#
# Usage:
#   cog-publish-link.sh <DATESTAMP>
#
# $1 DATESTAMP that latest should link to

chmod 755 $COG_PUBLISHED_DIR/$1/
ln -fn -s $COG_PUBLISHED_DIR/$1/elan.$1.consensus.fasta $COG_PUBLISHED_DIR/elan.latest.consensus.fasta
ln -fn -s $COG_PUBLISHED_DIR/$1/elan.$1.consensus.fasta.fai $COG_PUBLISHED_DIR/elan.latest.consensus.fasta.fai
ln -fn -s $COG_PUBLISHED_DIR/$1/majora.$1.metadata.tsv $COG_PUBLISHED_DIR/majora.latest.metadata.tsv
ln -fn -s $COG_PUBLISHED_DIR/$1/elan.$1.consensus.matched.fasta $COG_PUBLISHED_DIR/elan.latest.consensus.matched.fasta
ln -fn -s $COG_PUBLISHED_DIR/$1/majora.$1.metadata.matched.tsv $COG_PUBLISHED_DIR/majora.latest.metadata.matched.tsv

# NOTE samstudio8 2021-01-28
#      https://github.com/SamStudio8/elan-nextflow/issues/12
#      We've now moved away from providing symlinks at a daily resolution because
#      its hilariously expensive to maintain them and nobody actually uses them.
#      latest/ is now a real dir rather than a symlink. To match the expected
#      behaviour of the latest files being accessible via latest/ we'll additionally
#      now symlink the FASTA and TSV. Eventually we can deprecate providing them
#      above this dir.
#      I now use a `head` symlink to keep track of the last successful run.
ln -fn -s $COG_PUBLISHED_DIR/$1/elan.$1.consensus.fasta $COG_PUBLISHED_DIR/latest/elan.consensus.fasta
ln -fn -s $COG_PUBLISHED_DIR/$1/elan.$1.consensus.fasta.fai $COG_PUBLISHED_DIR/latest/elan.consensus.fasta.fai
ln -fn -s $COG_PUBLISHED_DIR/$1/majora.$1.metadata.tsv $COG_PUBLISHED_DIR/latest/majora.metadata.tsv
ln -fn -s $COG_PUBLISHED_DIR/$1/elan.$1.consensus.matched.fasta $COG_PUBLISHED_DIR/latest/elan.consensus.matched.fasta
ln -fn -s $COG_PUBLISHED_DIR/$1/majora.$1.metadata.matched.tsv $COG_PUBLISHED_DIR/latest/majora.metadata.matched.tsv
ln -fn -s $COG_PUBLISHED_DIR/$1/summary $COG_PUBLISHED_DIR/latest/summary
ln -fn -s $COG_PUBLISHED_DIR/$1/ $COG_PUBLISHED_DIR/head

# NOTE samstudio8 2022-01-25
#   ARTIFACTS_ROOT is the future, first map the elan status information to latest
#   Eventually this will be the home of the daily consensus data too
ln -fn -s $ARTIFACTS_ROOT/elan/$1 $ARTIFACTS_ROOT/elan/latest
