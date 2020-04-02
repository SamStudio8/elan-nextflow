#!/usr/bin/python
import sys
import csv
import os
import uuid

def best_effort_path(exts, dir_path, coguk_id, use_id=False):
    candidates = []
    if not os.path.exists(dir_path):
        sys.stderr.write("[FAIL] Dir %s for %s no longer exists\n" % (dir_path, coguk_id))
        return None

    for f in os.listdir(dir_path):
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
            sys.stderr.write("[WARN] Multiple %s found for %s in %s\n" % (str(exts), coguk_id, dir_path))
            return None
        else:
            return best_effort_path(exts, dir_path, coguk_id, use_id=True)
    else:
        return candidates[0]

manifest = csv.DictReader(open(sys.argv[1]), delimiter='\t')
runs_by_sample = {}
for row in manifest:
    run_name = row.get("run_name")
    central_sample_id = row.get("central_sample_id")

    if run_name and central_sample_id:
        if central_sample_id not in runs_by_sample:
            runs_by_sample[central_sample_id] = {}
        runs_by_sample[central_sample_id][run_name] = {"path": None, "user": None}
sys.stderr.write("[NOTE] Detected %d samples\n" % len(runs_by_sample))

for line in sys.stdin.readlines():
    fields = line.strip().split(os.path.sep)

    # Assume base layout for now
    username = fields[4]

    current_sample = None
    for sample_name in runs_by_sample:
        if sample_name in fields[-1] or sample_name in fields[-2]:
            current_sample = sample_name
            break
    if not current_sample:
        sys.stderr.write("[WARN] %s not in Majora\n" % str(fields))

    elif len(runs_by_sample[sample_name]) == 1:
        # Assume if there is only one run for this sample, this is the right one...?
        target_run_name = list(runs_by_sample[sample_name].keys())[0]
        runs_by_sample[sample_name][target_run_name] = {"user": username, "path": os.path.sep.join(fields[:-1])}
    else:
        # Attempt to disambiguate by looking for a directory named after the run
        for target_run_name in runs_by_sample[sample_name]:
            for f in fields:
                if target_run_name in f:
                    runs_by_sample[sample_name][target_run_name] = {"user": username, "path": os.path.sep.join(fields[:-1])}

uuids_for_runs = {}
for sample_name in runs_by_sample:
    for run_name in runs_by_sample[sample_name]:
        if run_name not in uuids_for_runs:
            uuids_for_runs[run_name] = uuid.uuid4()

        d = runs_by_sample[sample_name][run_name]["path"]
        best_fa = best_bam = ""

        if d:
            best_fa = best_effort_path(["fa", "fas", "fasta", "fna"], d, sample_name) or ""
            best_bam = best_effort_path(["bam"], d, sample_name) or ""
        if not best_fa or not best_bam:
            continue

        username = runs_by_sample[sample_name][run_name]["user"]
        if username == "nicholsz":
            username = "test-s.nicholls"
        else:
            username = "test-climb-covid19-nichollsn"

        print("\t".join([
            sample_name,
            run_name,
            username,
            str(uuids_for_runs[run_name]),
            d or '-',
            'consensus',
            str(uuid.uuid4()),
            best_fa,
            'alignment',
            str(uuid.uuid4()),
            best_bam,
        ]))

