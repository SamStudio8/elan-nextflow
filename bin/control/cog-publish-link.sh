set -euo pipefail

# Repoint latest
# This script is usually run automatically by cog-publish
# Under exceptional circumstances one may need to repoint the latest data set manually
#
# Usage:
#   cog-publish-link.sh <DATESTAMP>
#
# $1 DATESTAMP that latest should link to

# NOTE samstudio8 2022-01-25
#   ARTIFACTS_ROOT is the future, first map the elan status information to latest
#   Eventually this will be the home of the daily consensus data too
#   We used to protect access to the datestamp directory but users should just use latest
#   and expect unexpected behaviour if they don't
ln -fn -s $ARTIFACTS_ROOT/elan/$1 $ARTIFACTS_ROOT/elan/latest
ln -fn -s $ARTIFACTS_ROOT/elan/$1 $ARTIFACTS_ROOT/elan/head

