import sys

def get_by_contains_ext(s, exts):
    fns = [(i, x) for i, x in enumerate(s.split(' ')) for ext in exts if x.lower().endswith(ext)]
    return len(fns), fns

FH = sys.argv[1]
COG = sys.argv[2]

for line in open(FH):
    if line.startswith("@PG"):
        fields = line.strip().split('\t')
        for i, field in enumerate(fields):
            if field.startswith("CL"):
                cl_fields = field.split(' ')

                nfq, fqs = get_by_contains_ext(field, ["fastq", "fq", "fq.gz", "fastq.gz"])
                for fqi, curr_fq in enumerate(fqs):
                    field_j, fn = curr_fq
                    if nfq == 1:
                        # SINGLE
                        new_fn = "%s.fq.gz" % COG
                    elif nfq == 2:
                        # PAIRED
                        new_fn = "%s_%d.fq.gz" % (COG, fqi+1)
                    else:
                        sys.stderr.write("Cowardly refusing to update CL line with %d apparent FQ\n" % nfq)
                        break

                    sys.stderr.write("Replaced CL:%d <%s> -> <%s>\n" % (field_j, fn, new_fn))
                    cl_fields[field_j] = new_fn

                nsam, sams = get_by_contains_ext(field, ["sam"])
                if nsam == 1:
                    field_j, fn = sams[0]
                    new_fn = "%s.sam" % COG
                    sys.stderr.write("Replaced CL:%d <%s> -> <%s>\n" % (field_j, fn, new_fn))
                    cl_fields[field_j] = new_fn
                elif nsam > 1:
                    sys.stderr.write("Cowardly refusing to update CL line with more than one apparent SAM\n")


                nbam, bams = get_by_contains_ext(field, ["bam"])
                if nbam == 1:
                    field_j, fn = bams[0]
                    new_fn = "%s.bam" % COG
                    sys.stderr.write("Replaced CL:%d <%s> -> <%s>\n" % (field_j, fn, new_fn))
                    cl_fields[field_j] = new_fn
                elif nbam > 1:
                    sys.stderr.write("Cowardly refusing to update CL line with more than one apparent BAM\n")


                field = ' '.join(cl_fields)
                fields[i] = field
        print('\t'.join(fields))
    else:
        # This does not concern us
        print(line.strip())
