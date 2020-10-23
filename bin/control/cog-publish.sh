#!/usr/bin/bash

eval "$(conda shell.bash hook)"
conda activate samstudio8

set -euo pipefail

echo $1

source ~/.ocarina

# Get files that pass QC
ocarina --quiet --env get pag --test-name 'cog-uk-elan-minimal-qc' --ls-files --pass --task-wait --task-wait-attempts 75 --task-wait-minutes 1 > elan.pass.latest

grep 'consensus' elan.pass.latest > elan.pass.latest.consensus.ls
cut -f3 elan.pass.latest.consensus.ls > pass.fasta.ls
wc -l pass.fasta.ls

grep 'alignment' elan.pass.latest > elan.pass.latest.alignment.ls
cut -f3 elan.pass.latest.alignment.ls > pass.bam.ls
wc -l pass.bam.ls

# mkdirs
mkdir $COG_PUBLISHED_DIR/$1/
mkdir $COG_PUBLISHED_DIR/$1/fasta
mkdir $COG_PUBLISHED_DIR/$1/alignment
mkdir $COG_PUBLISHED_DIR/$1/qc
mkdir $COG_PUBLISHED_DIR/$1/summary
chmod 755 $COG_PUBLISHED_DIR/$1/
chmod 755 $COG_PUBLISHED_DIR/$1/fasta
chmod 755 $COG_PUBLISHED_DIR/$1/alignment
chmod 755 $COG_PUBLISHED_DIR/$1/qc
chmod 755 $COG_PUBLISHED_DIR/$1/summary

# Linky
echo "[CPUB]" `date` " - Linking FASTA"
for fas in `cat pass.fasta.ls`;
do
    #echo "Linking $fas"
    ln -s $fas $COG_PUBLISHED_DIR/$1/fasta/
done

echo "[CPUB]" `date` " - Linking BAM"
echo 'Linking BAMs'
for bam in `cat pass.bam.ls`;
do
    #echo "Linking $bam"
    ln -s $bam $COG_PUBLISHED_DIR/$1/alignment/
    ln -s $bam.bai $COG_PUBLISHED_DIR/$1/alignment/
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
python $ELAN_SOFTWARE_DIR/bin/control/reconcile_downstream.py $COG_PUBLISHED_DIR/$1/majora.$1.metadata.tsv $COG_PUBLISHED_DIR/$1/fasta/ > $COG_PUBLISHED_DIR/$1/elan.$1.consensus.matched.fasta
chmod 644 $COG_PUBLISHED_DIR/$1/majora.$1.metadata.matched.tsv
chmod 644 $COG_PUBLISHED_DIR/$1/elan.$1.consensus.matched.fasta

# Make QC tables available
echo "[CPUB]" `date` " - Linking QC"
for qcc in `find $ELAN_DIR/staging/qc/ -name '*.qc'`;
do
    #echo "Linking $qcc"
    ln -s $qcc $COG_PUBLISHED_DIR/$1/qc/
done

# Fix perms
#echo "Fixing perms..."
#find $ELAN_DIR -type d -exec chmod 755 {} \;
#find $ELAN_DIR -type f -exec chmod 644 {} \;
#find $COG_PUBLISHED_DIR/$1/ -type l -exec chmod 644 {} \;

# Repoint latest
echo "[CPUB]" `date` " - Linking latest"
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
COUNT_PASS=`wc -l pass.fasta.ls | cut -f1 -d' '`
COUNT_NEW=`wc -l $ELAN_DIR/staging/summary/$1/swell.qc.tsv | cut -f1 -d' '`
curl -X POST -H 'Content-type: application/json' --data '{"text":"\n*COG-UK inbound pipeline complete*\n'$COUNT_NEW' new sequences matched to Majora metadata\n'$COUNT_PASS' sequences passed basic quality control to date!\nArtifacts successfully published by elan-nextflow to `'$COG_PUBLISHED_DIR'/latest`\n***\n_QC Reports have been calculated and can be reached from your Majora profile._\n_The outbound distribution pipeline will run next Monday._\n_Thanks for your patience, have a nice day!_"}' $SLACK_REAL_HOOK
