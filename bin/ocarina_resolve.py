#!/usr/bin/env python3
import sys
import csv
import os
import uuid
import re
import datetime

RECENT_DAYS_NO = 14
RECENT_DAYS_DESC = "two weeks"

def map_foel_user(username):
    # Map a user directory (eg. private producer name) to a Majora user
    # We only need the external-dhsc user for now, but if we onboard other
    # types of organisation we will need to load a mapping of ingest directories
    # (or regexes) to majora usernames
    # TODO Will need to handle the case where a user is unmapped
    return "external-dhsc"

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

# Read all the samples into struct
manifest = csv.DictReader(open(sys.argv[1]), delimiter='\t')
runs_by_sample = {}
new_count = 0
tot_count = 0
invalid_count = 0


countries = {"UK-ENG": "England", "UK-WLS": "Wales", "UK-NIR": "Northern_Ireland", "UK-SCT": "Scotland", "": "UNKNOWN"}
for row in manifest:
    if None in row.keys():
        sys.stderr.write("[SKIP] Ignoring %s:%s as it has too many columns (there is probably a tab character in the source data...)\n" % (row.get("central_sample_id"), row.get("run_name")))
        invalid_count += 1
        continue
    run_name = row.get("run_name")
    central_sample_id = row.get("central_sample_id")

    if run_name and central_sample_id:
        if central_sample_id not in runs_by_sample:
            runs_by_sample[central_sample_id] = {}

        pags = [x for x in row.get("published_as", "").split(",") if len(x) > 0] # shared across runs so just do it here at a minor penalty
        n_pags = len(pags)
        is_new = True
        has_pag = ''
        for pag in pags:
            if run_name in pag:
                has_pag = pag
                is_new = False

        source_site = row.get("submission_org_code")
        seq_site = row.get("sequencing_org_code")
        if seq_site == "SANG":
            if source_site == "SANG" or source_site == "PHEC":
                source_site = central_sample_id[:4] # SANG submits on behalf of some sites
            else:
                pass # skip and pass through the submission_org
        elif seq_site == "PHEC":
            if source_site == "PHEC":
                source_site = central_sample_id[:4] # PHEC is using prefix hackery like SANG
            else:
                pass # skip and pass through the submission_org

        lp = row.get("library_primers")
        if lp and len(lp) > 0 and lp != "None":
            tile_str = lp
            tiles_t = re.findall(r'\d+', tile_str)
        else:
            tile_str = row.get("meta.artic.primers", "0")
            tiles_t = re.findall(r'\d+', tile_str)

        tiles = 0
        if len(tiles_t) > 0:
            tiles = tiles_t[0]

        tot_count += 1
        if is_new:
            new_count += 1 # run must be new
        
        cor_date = row.get('collection_date')
        if not cor_date or cor_date == "None":
            cor_date = "R%s" % row.get('received_date')

        runs_by_sample[central_sample_id][run_name] = {
            "path": None,
            "user": None,
            "source_site": source_site,
            "site": seq_site,
            "date": row.get("sequencing_submission_date"),
            "pipe": row.get("sequencing_uuid", ""),
            "tech": row.get("instrument_make", ""),
            "n_pags": n_pags,
            "pag": has_pag,
            "is_new": is_new,
            "run_name": run_name,
            "tiles": tiles,
            "adm0": row.get("adm0", "UNKNOWN"),
            "adm1_mapped": countries[row.get("adm1", "")],
            "cor_date": cor_date,
            "seq_date": row.get("sequencing_submission_date", ""),
        }
        #print(runs_by_sample[central_sample_id][run_name])

new_samples = 0
for sample in runs_by_sample:
    for run in runs_by_sample[sample]:
        if runs_by_sample[sample][run]["n_pags"] == 0:
            new_samples += 1
#sys.stderr.write("[NOTE] Detected %d total samples, %d currently unpublished samples\n" % (len(runs_by_sample), new_samples))
sys.stderr.write("[NOTE] Detected %d total consensus sequences, %d currently unpublished\n" % (tot_count, new_count))
sys.stderr.write("[NOTE][META][BAD-LINES]|%d\n" % invalid_count)
sys.stderr.write("[NOTE][META][RECENT-DAYS]|%d\n" % RECENT_DAYS_NO)

orphaned_dirs = {}
for line in sys.stdin.readlines():
    fields = line.strip().split(os.path.sep)

    # Assume base layout for now
    username = fields[4]

    #NOTE samstudio8/2021-10-13
    # Slowly moving away from our dependency on the file system, if the "username"
    # from the file path is from the new ingest solution directory, we must explicitly
    # map the producer name from the subdirectory after the upload/ subdir
    if username == "foel-ingest":
        username = map_foel_user(fields[6])

    current_sample = None
    if fields[-1].upper() in runs_by_sample:
        current_sample = fields[-1].upper()
    elif fields[-1] in runs_by_sample:
        current_sample = fields[-1]

    elif fields[-2].upper() in runs_by_sample:
        current_sample = fields[-2].upper()
    elif fields[-2] in runs_by_sample:
        current_sample = fields[-2]

    else:
        for field in re.split('[^A-Za-z0-9-]', fields[-1]):
            if field.upper() in runs_by_sample:
                current_sample = field.upper()
            elif field in runs_by_sample:
                current_sample = field

    if not current_sample:
        try:
            mod_time = datetime.datetime.fromtimestamp( os.path.getmtime(line.strip()) )
        except OSError:
            # Assume if the file can't be reached it is new for the sake of argument
            mod_time = datetime.datetime.now()

        if mod_time >= (datetime.datetime.now() - datetime.timedelta(days=RECENT_DAYS_NO)):
            recent = True
        else:
            recent = False

        d = os.path.sep.join(fields[:-2])
        if d not in orphaned_dirs:
            orphaned_dirs[d] = {"username": username, "count": 0, "files": [], "recent-count":0, "unrecent-count": 0}

        if recent:
            orphaned_dirs[d]["recent-count"] += 1
        else:
            orphaned_dirs[d]["unrecent-count"] += 1

        orphaned_dirs[d]["count"] += 1
        orphaned_dirs[d]["files"].append(fields[-2:])

        #sys.stderr.write("[WARN] %s not in Majora\n" % str(fields))
        continue
    
    #NOTE samstudio8/2020-07-02
    # I should have done this sooner, given the mistakes users are making when submitting samples
    # to SANG, there is a potential for a race condition where one site can accidentally pick up
    # the "wrong" run metadata. What with sites now going back to upload old failures we can
    # no longer give people the benefit of the doubt when one run is being uploaded,
    # Going forward the path must contain the run_name exactly to be correctly disambiguated.
    #
    #elif len(runs_by_sample[current_sample]) == 1:
    #    # Assume if there is only one run for this sample, this is the right one...?
    #    target_run_name = list(runs_by_sample[current_sample].keys())[0]
    #    runs_by_sample[current_sample][target_run_name].update({"user": username, "path": os.path.sep.join(fields[:-1])})
    #
    else:
        # Attempt to disambiguate by looking for a directory named after the run
        for target_run_name in runs_by_sample[current_sample]:
            for f in fields:
                if target_run_name in f:
                    runs_by_sample[current_sample][target_run_name].update({"user": username, "path": os.path.sep.join(fields[:-1])})

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
            #sys.stderr.write("[FAIL] Directory not resolved for %s %s\n" % (sample_name, run_name))
            continue
        else:
            matched_samples[sample_name][run_name] += 1

        username = runs_by_sample[sample_name][run_name]["user"]
        pipeuuid = runs_by_sample[sample_name][run_name]["pipe"]

        print("\t".join([
            '1' if runs_by_sample[sample_name][run_name]["is_new"] else '0',
            sample_name,
            run_name,
            username,
            pipeuuid,
            "%s-%s-%s-%s" % (runs_by_sample[sample_name][run_name]["site"], runs_by_sample[sample_name][run_name]["tech"], runs_by_sample[sample_name][run_name]["date"].replace('-', ''), pipeuuid[:5]),
            runs_by_sample[sample_name][run_name]["tech"],
            d or '-',
            'consensus',
            best_fa,
            'alignment',
            best_bam,
            runs_by_sample[sample_name][run_name]["source_site"],
            runs_by_sample[sample_name][run_name]["site"],
            runs_by_sample[sample_name][run_name]["pag"],
            str(runs_by_sample[sample_name][run_name]["tiles"]),
            runs_by_sample[sample_name][run_name]["adm0"],
            runs_by_sample[sample_name][run_name]["adm1_mapped"],
            runs_by_sample[sample_name][run_name]["cor_date"],
            runs_by_sample[sample_name][run_name]["seq_date"],
        ]))

missing_samples_by_site = {}
for sample_name in matched_samples:
    if sum(matched_samples[sample_name][run_name] for run_name in matched_samples[sample_name]) == 0:
        # If not seen at all
        sites_seen = set([])
        most_recent_run_date = {}

        for run_name in matched_samples[sample_name]:
            is_missing = False
            site = runs_by_sample[sample_name][run_name]["site"]
            if site not in missing_samples_by_site and site not in sites_seen:
                missing_samples_by_site[site] = {
                    "count": 0,
                    "sample_count": 0,

                    "recent_sample_run_count": 0,
                    "unrecent_sample_run_count": 0,

                    "missing": {},
                    "deleted": {},
                    "deleted_count": 0,
                }
            missing_samples_by_site[site]["count"] += 1
            if runs_by_sample[sample_name][run_name]["pag"]:
                missing_samples_by_site[site]["deleted_count"] += 1
                missing_samples_by_site[site]["deleted"][sample_name] = matched_samples[sample_name].keys()
            else:
                is_missing = True
                missing_samples_by_site[site]["missing"][sample_name] = matched_samples[sample_name].keys()
                if site not in sites_seen:
                    missing_samples_by_site[site]["sample_count"] += 1
            sites_seen.add(site)

            if is_missing:
                # Get the most recent run date
                run_dt = datetime.datetime.strptime(runs_by_sample[sample_name][run_name]["seq_date"], "%Y-%m-%d")
                if site not in most_recent_run_date:
                    most_recent_run_date[site] = run_dt
                else:
                    if most_recent_run_date[site] < run_dt:
                        most_recent_run_date[site] = run_dt

        for site in most_recent_run_date:
            if most_recent_run_date[site] >= (datetime.datetime.now() - datetime.timedelta(days=RECENT_DAYS_NO)):
                missing_samples_by_site[site]["recent_sample_run_count"] += 1
            else:
                missing_samples_by_site[site]["unrecent_sample_run_count"] += 1

for site in missing_samples_by_site:
    if missing_samples_by_site[site]["sample_count"] > 0:
        sys.stderr.write("[NOFILE][ORPHAN-SITE][%s] %d biosamples sequenced at %s missing a matched file on CLIMB\n" % (site, missing_samples_by_site[site]["sample_count"], site))
        sys.stderr.write("[NOFILE][ORPHAN-NEW-SITE][%s] %d biosamples sequenced at %s (%s ago or newer) missing a matched file on CLIMB\n" % (site, missing_samples_by_site[site]["recent_sample_run_count"], site, RECENT_DAYS_DESC))
        sys.stderr.write("[NOFILE][ORPHAN-OLD-SITE][%s] %d biosamples sequenced at %s (%s ago or older) missing a matched file on CLIMB\n" % (site, missing_samples_by_site[site]["unrecent_sample_run_count"], site, RECENT_DAYS_DESC))
        for cogid in sorted(missing_samples_by_site[site]["missing"]):
            locations = missing_samples_by_site[site]["missing"][cogid]
            sys.stderr.write("[MAJORA][ORPHAN-COGX][%s] %s in Majora but not CLIMB. Possible run names: %s\n" % (site, cogid, str(list(locations))))
        for cogid in sorted(missing_samples_by_site[site]["deleted"]):
            locations = missing_samples_by_site[site]["deleted"][cogid]
            #sys.stderr.write("[MAJORA][DELETE-COGX][%s] %s in Majora, and released. but no longer on CLIMB. Possible run names: %s\n" % (site, cogid, str(list(locations))))

old_userfiles = {}
d_i = 0
for d_key, d in orphaned_dirs.items():
    d_i += 1
    sys.stderr.write("[NOMETA][ORPHAN-DIRX][OD%s] %d FASTA orphaned without metadata inside %s\n" % (d_i, d["count"], d_key))
    sys.stderr.write("[NOMETA][ORPHAN-NEW-DIRX][OD%s] %d FASTA orphaned without metadata inside %s (%s ago or newer)\n" % (d_i, d["recent-count"], d_key, RECENT_DAYS_DESC))
    sys.stderr.write("[NOMETA][ORPHAN-OLD-DIRX][OD%s] %d FASTA orphaned without metadata inside %s (%s ago or older))\n" % (d_i, d["unrecent-count"], d_key, RECENT_DAYS_DESC))
    for f in d["files"]:
        sys.stderr.write("[NOMETA][ORPHAN-FILE][OD%s] %s/%s #Delete this file to suppress this error.\n" % (d_i, d_key, os.path.sep.join(f)))

    if d["username"] not in old_userfiles:
        old_userfiles[d["username"]] = 0
    old_userfiles[d["username"]] += d["unrecent-count"]

for user, count in old_userfiles.items():
    if count > 0:
        sys.stderr.write("[NOMETA][ORPHAN-USER-DIRX]|%s|%d|FASTA orphaned without metadata by user\n" % (user, count))
