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
        name, seqs, last = last[1:].replace(" ", "_"), [], None #NOTE/samstudio8, add _ to spaces here
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

METADATA_FP = sys.argv[1]
PUBLISHED_DIR = sys.argv[2]
if len(sys.argv) == 4:
    NEW_METADATA_FP = sys.argv[3]
else:
    NEW_METADATA_FP = METADATA_FP.replace(".tsv", ".matched.tsv")

# Open the metadata file
metadata_fh = csv.DictReader(open(METADATA_FP), delimiter="\t")
new_metadata_fh = csv.DictWriter(open(NEW_METADATA_FP, 'w'), metadata_fh.fieldnames + ["fasta_header"], delimiter="\t")
new_metadata_fh.writeheader()

seen_files = set([])

valid_rn = 1
for row in metadata_fh:

    # Build the file name to the artifact from the central_sample_id and run_name
    c_fn = "%s/%s.%s.climb.fasta" % (PUBLISHED_DIR, row["central_sample_id"], row["run_name"])

    if c_fn in seen_files:
        # Cover the special case where users have managed to upload a sample on multiple libraries for the same run...
        # In these cases we can just write out the first line of the metadata file with no loss of information on the biosample and run
        # We risk losing library-specific properties (eg. primer set etc), fixing that is left as an exercise to an avid downstream analyst
        sys.stderr.write("[SKIP] %s\n" % c_fn)
        continue
    else:
        seen_files.add(c_fn)

    # If the path exists, open it and dump out the FASTA.
    try:
        # NOTE(samstudio8, 20210106) - Use a try here to open the file straight away,
        # this allows us to skip a `stat` on the file system from os.path.exists
        c_fh = open(c_fn)
    except FileNotFoundError:
        # Skip to the next record in the hope it actually exists, this is a dud
        continue

    for name, seq, qual in readfq(c_fh):
        # Modify the header slightly to include the current valid row number for mapping later
        # (although this is unnessary now the table will be perfectly paired)
        # We can actually do whatever we like to the header here, as all of the metadata is in scope, so knock yourself out.
        name = "COGUK/%s/%s:%s" % (row["central_sample_id"], row["sequencing_org_code"], row["run_name"])
        row["fasta_header"] = name
        print(">%s|%d" % (name, valid_rn))
        print(seq)
    c_fh.close()
    valid_rn += 1

    # Write the metadata out for perfect pairing
    new_metadata_fh.writerow(row)
