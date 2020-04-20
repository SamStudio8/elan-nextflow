#!/usr/bin/python
import sys
import csv
import os
import uuid
import re

def best_effort_path(exts, dir_path, coguk_id, use_id=False):
    candidates = []
    if not os.path.exists(dir_path):
        sys.stderr.write("[FAIL] Dir %s for %s no longer exists\n" % (dir_path, coguk_id))
        return None

    for f in os.listdir(dir_path):
        if f.startswith('.'):
            continue
        for ext in exts:
            #if ext == "bam" and "primertrimmed" in f:
            #    continue
            if f.endswith(ext):
                if use_id:
                    if coguk_id in f:
                        candidates.append(f)          
                        break
                else:
                    candidates.append(f)
                    break

    if len(candidates) == 0:
        sys.stderr.write("[WARN] No %s found for %s in %s\n" % (str(exts), coguk_id, dir_path))
        return None
    elif len(candidates) > 1:        
        if use_id:
            sys.stderr.write("[WARN] Multiple %s found for %s in %s\n" % (str(exts), str(candidates), dir_path))
            return None
        else:
            return best_effort_path(exts, dir_path, coguk_id, use_id=True)
    else:
        return candidates[0]

def lookup_instrument(make):
    makes = {
        "ILLUMINA": "ILL",
        "OXFORD_NANOPORE": "ONT",
    }
    if make.upper() in makes:
        return makes[make.upper()]
    else:
        return "RUN"

manifest = csv.DictReader(open(sys.argv[1]), delimiter='\t')
runs_by_sample = {}
for row in manifest:
    run_name = row.get("run_name")
    central_sample_id = row.get("central_sample_id")

    if run_name and central_sample_id:
        if central_sample_id not in runs_by_sample:
            runs_by_sample[central_sample_id] = {}
        runs_by_sample[central_sample_id][run_name] = {
            "path": None,
            "user": None,
            "site": row.get("sequencing_site_code"),
            "date": row.get("sequencing_submission_date"),
            "pipe": row.get("sequencing_uuid"),
            "tech": lookup_instrument(row.get("instrument_make")),
        }
sys.stderr.write("[NOTE] Detected %d samples\n" % len(runs_by_sample))

for line in sys.stdin.readlines():
    fields = line.strip().split(os.path.sep)

    # Assume base layout for now
    username = fields[4]

    current_sample = None
    for sample_name in runs_by_sample:
        if sample_name.upper() == fields[-1].upper() or sample_name.upper() == fields[-2].upper():
            current_sample = sample_name
            break

        for field in re.split('[^A-Za-z0-9-]', fields[-1]):
            if field.upper() == sample_name.upper():
                current_sample = sample_name
                break

        #if sample_name in fields[-1] or sample_name in fields[-2]:
        #    current_sample = sample_name
        #    break

    if not current_sample:
        sys.stderr.write("[WARN] %s not in Majora\n" % str(fields))
        continue
    elif len(runs_by_sample[current_sample]) == 1:
        # Assume if there is only one run for this sample, this is the right one...?
        target_run_name = list(runs_by_sample[current_sample].keys())[0]
        runs_by_sample[current_sample][target_run_name] = {"user": username, "path": os.path.sep.join(fields[:-1])}
    else:
        # Attempt to disambiguate by looking for a directory named after the run
        for target_run_name in runs_by_sample[current_sample]:
            for f in fields:
                if target_run_name in f:
                    runs_by_sample[current_sample][target_run_name] = {"user": username, "path": os.path.sep.join(fields[:-1])}

matched_samples = {}
for sample_name in runs_by_sample:
    if sample_name not in matched_samples:
        matched_samples[sample_name] = {}
    for run_name in runs_by_sample[sample_name]:
        if run_name not in matched_samples[sample_name]:
            matched_samples[sample_name][run_name] = 0

        d = runs_by_sample[sample_name][run_name]["path"]
        best_fa = best_bam = ""

        if d:
            best_fa = best_effort_path(["fa", "fas", "fasta", "fna"], d, sample_name) or ""
            best_bam = best_effort_path(["bam"], d, sample_name) or ""

        if not d or (not best_fa or not best_bam):
            sys.stderr.write("[FAIL] Directory not resolved for %s %s\n" % (sample_name, run_name))
            continue
        else:
            matched_samples[sample_name][run_name] += 1

        username = runs_by_sample[sample_name][run_name]["user"]
        pipeuuid = runs_by_sample[sample_name][run_name]["pipe"]

        print("\t".join([
            sample_name,
            run_name,
            username,
            pipeuuid,
            "%s-%s-%s-%s" % (runs_by_sample[sample_name][run_name]["site"], runs_by_sample[sample_name][run_name]["tech"], runs_by_sample[sample_name][run_name]["date"].replace('-', ''), pipeuuid[:5]),
            d or '-',
            'consensus',
            best_fa,
            'alignment',
            best_bam,
        ]))

for sample_name in matched_samples:
    if sum(matched_samples[sample_name][run_name] for run_name in matched_samples[sample_name]) == 0:
        sys.stderr.write("[FAIL] No matches found for %s: [%s]\n" % (sample_name, str(list(matched_samples[sample_name].keys()))))
