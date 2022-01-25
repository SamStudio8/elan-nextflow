#!/usr/bin/env python3

import argparse
import csv
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--tsv", required=True)
parser.add_argument("--col", default="pag_name")
parser.add_argument("--sep-in", default='\t')
parser.add_argument("--sep-out", default='\t')
args = parser.parse_args()

with open(args.tsv) as fh:
    csv_in = csv.DictReader(fh, delimiter=args.sep_in)

    new_fieldnames = ["central_sample_id", "run_name", "sequencing_org_code"] + csv_in.fieldnames
    csv_out = csv.DictWriter(sys.stdout, delimiter=args.sep_out, fieldnames=new_fieldnames)
    csv_out.writeheader()

    for row in csv_in:
        pag = row[args.col]
        cog_uk, central_sample_id, org_run = pag.split('/')
        org_code, run_name = org_run.split(':')

        row["central_sample_id"] = central_sample_id
        row["run_name"] = run_name
        row["sequencing_org_code"] = org_code
        csv_out.writerow(row)

