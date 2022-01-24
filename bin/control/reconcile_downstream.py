# Thanks Heng
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        #name, seqs, last = last[1:].partition(" ")[0], [], None
        name, seqs, last = last[1:].replace(" ", "_"), [], None #NOTE/samstudio8, we change _ to spaces here
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

import csv
import os
import sys
import datetime
import time
import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--imeta", required=True)
parser.add_argument("--ometa", required=False) # reconcile_downstream does not pass ometa by default
parser.add_argument("--ifasta", required=True)
parser.add_argument("--ofasta", required=False) # reconcile_downstream uses stdout for FASTA by default
parser.add_argument("--kill-pagfiles", required=True) # required as you'll almost certainly mess it up without it, COG-UK PAG converted internally to COGUK PAG
parser.add_argument("--pass-pagfiles", required=True) # required as you won't get any new sequences without it, COG-UK PAG converted internally to COGUK PAG
args = parser.parse_args()

METADATA_FP = args.imeta
if args.ometa:
    NEW_METADATA_FP = args.ometa
else:
    NEW_METADATA_FP = METADATA_FP.replace(".tsv", ".matched.tsv")

if args.ofasta:
    NEW_FASTA_FP = open(args.ofasta, 'w')
else:
    NEW_FASTA_FP = sys.stdout

# 20211118/samstudio8 The caffeine-cat update
# This step takes nearly six hours on a bad day (and an hour on a good day),
# so let's kick it up a notch with a "cache". We'll grab the most recent published
# reconciled consensus set and use it as a cache to reference existing sequences
# rather than hitting the individual FASTA objects a million times
PUBLISHED_FASTA = pysam.FastaFile(args.ifasta)

# The published FASTA header contains the PAG (COGUK/SAMPLE/SITE:RUN) and other pipe delimited information
# We'll need to keep a PAG to FASTA header map, so we can translate the metadata PAG to a header for fetching from the FASTA
# We can assume the PAG is always the first key in the header because we control how Elan produces it
# NOTE You must be aware the PAG name contains COG-UK but the FASTA headers use COGUK as standard, and require translation for kill.ls and pass.ls
PUBLISHED_FASTA_PAG_TO_HEADER = {ref.split('|')[0]:ref for ref in PUBLISHED_FASTA.references}

if len(PUBLISHED_FASTA_PAG_TO_HEADER) != len(PUBLISHED_FASTA.references):
    # Should not happen but we'll guard the case anyway
    sys.stderr.write("[FAIL] Published FASTA appears to have a PAG map to more than one sequence, which should be impossible...\n")
    sys.exit(65) # EX_DATAERR

# The kill ls file will list all PAGs suppressed since Elan last ran successfully
# We'll use this to remove any existing sequences that should no longer appear in the consensus
# Previously we relied on unlinking and hitting a file that no longer exists, this is faster and more explicit
# This is required because the metadata file has no idea if a sequence has been suppressed and will blindly copy it from the latest FASTA if we don't tell it otherwise
suppressed_pags = set([])
if args.kill_pagfiles:
    with open(args.kill_pagfiles) as suppressed_fh:
        for line in suppressed_fh:
            fields = line.replace("COG-UK", "COGUK").strip().split()
            suppressed_pags.add(fields[0])

# The pass ls file will list all PAGs that have passed QC since Elan last ran successfully
# We can use this to determine whether a latest.fasta cache miss is really a new valid sequence
# or whether it was previously suppressed, failed QC or is just outright invalid
# This will save us hitting the file system for files that don't exist at all (pretty bad case for a full ceph MDS)
passed_pags = {}
if args.pass_pagfiles:
    with open(args.pass_pagfiles) as passed_fh:
        for line in passed_fh:
            # NOTE If the order of ocarina pagfiles changes this will need updating...
            pag, file_type, file_path, *_ = line.replace("COG-UK", "COGUK").strip().split('\t')
            if pag not in passed_pags:
                passed_pags[pag] = {"consensus": None, "alignment": None}
            passed_pags[pag][file_type] = file_path

# Open the metadata file
metadata_fh = csv.DictReader(open(METADATA_FP), delimiter="\t")
new_metadata_fh = csv.DictWriter(open(NEW_METADATA_FP, 'w'), metadata_fh.fieldnames + ["fasta_header"], delimiter="\t")
new_metadata_fh.writeheader()

seen_pags = set([])

start = time.time()
end = None

valid_rn = 1
n_hit = n_miss = n_skip = n_kill = 0
for row in metadata_fh:

    if valid_rn % 50_000 == 0:
        end = time.time()
        delta = datetime.timedelta(seconds= (end - start) )
        sys.stderr.write("[TIME] Processed block (%d) in %s\n" % (valid_rn, str(delta)))
        start = time.time()

    # Construct the PAG name (should probably be passed in from Majora somehow rather than manually constructed)
    # NOTE Reminder: for accidental historical reasons the PAG name contains COG-UK but the FASTA headers use COGUK as standard
    c_pag = "COGUK/%s/%s:%s" % (row["central_sample_id"], row["sequencing_org_code"], row["run_name"])

    # Attempt to fetch the key from the published FASTA
    # It won't be here if the FASTA is new so we'll still need to hit the file system in those cases
    c_header = PUBLISHED_FASTA_PAG_TO_HEADER.get(c_pag)

    if c_pag in seen_pags:
        # Cover the special case where users have managed to upload a sample on multiple libraries for the same run...
        # In these cases we can just write out the first line of the metadata file with no loss of information on the biosample and run
        # We risk losing library-specific properties (eg. primer set etc), fixing that is left as an exercise to an avid downstream analyst
        #sys.stderr.write("[SEEN] %s\n" % c_pag)
        continue
    elif c_pag in suppressed_pags:
        n_kill += 1
        sys.stderr.write("[KILL] %s\n" % c_pag)
        continue
    else:
        seen_pags.add(c_pag)


    if c_header:
        # The current PAG matches to a header in the latest published FASTA which
        # means we can short circuit fetching the individual FASTA from ceph and
        # grab it from the open file handle
        #sys.stderr.write("[_HIT] %s\n" % c_pag)
        c_seq = PUBLISHED_FASTA.fetch(c_header)
        n_hit += 1
    else:
        # If we're hitting the file system the file must be
        # * new today
        # * previously suppressed, failed QC, or doesnt even exist in the first place
        # So we can determine the former from the latter by checking if the PAG was
        # expected to be added today and short circuit checking the FS entirely
        # ceph will be very grateful
        if c_pag not in passed_pags:
            #sys.stderr.write("[SKIP] %s\n" % c_pag) # uncomment for big pointless log
            n_skip += 1
            continue

        # Elsewise this file really is new and we'll go and open it
        c_fn = passed_pags[c_pag].get("consensus")
        if not c_fn:
            sys.stderr.write("[WARN] %s does not have consensus entry in --pass-pagfiles which should be impossible\n" % c_pag)
            continue

        try:
            # NOTE(samstudio8, 20210106) - Use a try here to open the file straight away,
            # this allows us to skip a `stat` on the file system from os.path.exists
            sys.stderr.write("[MISS] %s\n" % c_pag)
            n_miss += 1
            c_fh = open(c_fn)
        except FileNotFoundError:
            sys.stderr.write("[WARN] %s had consensus entry at %s that could not be opened which should be impossible\n" % (c_pag, c_fn))
            continue

        for name, seq, qual in readfq(c_fh):
            c_seq = seq
            break # Take the first sequence (there should only ever be one anyway)
        c_fh.close()

    if not c_seq:
        continue

    # Write the FASTA to stdout, update the valid_rn
    # For legacy reasons we modify the header slightly to include the current
    # valid row number for downstream tools to map the FASTA to a metadata row
    # We can actually do whatever we like to the header here, as all of the
    # metadata is in scope, so knock yourself out if you want to pipe delimit all
    # sorts of fun stuff
    NEW_FASTA_FP.write(">%s|%d\n" % (c_pag, valid_rn))
    NEW_FASTA_FP.write(c_seq)
    NEW_FASTA_FP.write('\n')

    # Write the metadata out for perfect pairing
    row["fasta_header"] = c_pag
    new_metadata_fh.writerow(row)
    valid_rn += 1

if args.ofasta:
    NEW_FASTA_FP.close()

# TODO Check refd size == input
sys.stderr.write("[DONE] %d hits, %d misses, %d skips, %d kills\n" % (n_hit, n_miss, n_skip, n_kill))

