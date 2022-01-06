#!/usr/bin/env python
import sys

# thanks heng
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
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

iupac_chars = set("ACGTRYKMSWBDHVN")

if len(sys.argv) != 3:
    sys.exit(1)
min_len = int(sys.argv[2])

# Handle rename of FASTA sequence
try:
    fasta_fh = open(sys.argv[1])
except IOError:
    sys.exit(1)

heng_iter = readfq(fasta_fh)
for name, seq, qual in heng_iter:
    break # first record only?
if not name:
    sys.exit(1)
if not seq:
    sys.exit(2)
if any(base.upper() not in iupac_chars for base in seq):
    sys.exit(65)
if len(seq) < min_len:
    sys.exit(3)
fasta_fh.close()
