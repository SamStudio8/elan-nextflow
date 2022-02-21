import hashlib
import argparse
import pathlib
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--reference_list", type=argparse.FileType("r"))
parser.add_argument("--test_list", type=argparse.FileType("r"))
parser.add_argument("--print_missing_ids", action="store_true")
args = parser.parse_args()


def hash_file(filepath, blocksize=2 ** 20):
    m = hashlib.md5()
    with open(filepath, "rb") as f:
        while True:
            buf = f.read(blocksize)
            if not buf:
                break
            m.update(buf)
    return m.hexdigest()


def compare_files(ref_file, test_file):
    ref_hash = hash_file(ref_file)
    test_hash = hash_file(test_file)
    return (ref_hash, test_hash)


def parse_file(fh):
    splits = (line.split("\t") for line in fh)
    out_dict = {}
    for split in splits:
        cog_id = split[0]
        fasta_path = pathlib.Path(split[5])
        bam_path = pathlib.Path(split[6])
        out_dict[cog_id] = {"fasta": fasta_path, "bam": bam_path}
    return out_dict


reference_dict = parse_file(args.reference_list)
test_dict = parse_file(args.test_list)


def test_id_list(
    id_list, reference_dictionary=reference_dict, test_dictionary=test_dict
):
    header_printed = False
    true_count = 0
    false_count = 0
    for ref_id in id_list:
        ref_files = reference_dictionary[ref_id]
        test_files = test_dictionary[ref_id]
        hashes = {}
        for filetype in ("fasta", "bam"):
            hashes[filetype] = compare_files(ref_files[filetype], test_files[filetype])
        fasta_match = hashes["fasta"][0] == hashes["fasta"][1]
        bam_match = hashes["bam"][0] == hashes["bam"][1]
        if fasta_match and bam_match:
            true_count += 1
        else:
            if not header_printed:
                sys.stdout.write(
                    f"cog_id\tfasta_match\tref_fastapath\ttest_fastapath\tbam_match\tref_bampath\ttest_bampath\n"
                )
                header_printed = True
            false_count += 1
            sys.stdout.write(
                "\t".join(
                    [
                        str(x)
                        for x in [
                            ref_id,
                            fasta_match,
                            ref_files["fasta"],
                            test_files["fasta"],
                            bam_match,
                            ref_files["bam"],
                            test_files["bam"],
                        ]
                    ]
                )
                + "\n"
            )
    sys.stdout.write(
        f"TEST EXECUTION FINISHED:\nCOG-IDs which passed inspection - {true_count}\nCOG-IDs which failed inspection - {false_count}\n"
    )


ref_ids = set(ref_id for ref_id in reference_dict.keys())
test_ids = set(test_id for test_id in test_dict.keys())

if ref_ids == test_ids:
    test_id_list(ref_ids)
elif test_ids.issubset(ref_ids):
    sys.stdout.write(
        "IDs do not Intersect, running on intersection\n-----------------------------------------------\n"
    )
    if args.print_missing_ids:
        sys.stdout.write("List of test_ids not present within reference list:\n")
        for cog_id in ref_ids.difference(test_ids):
            sys.stdout.write(str(cog_id) + "\n")
    id_intersection = test_ids.intersection(ref_ids)
    test_id_list(id_intersection)
else:
    sys.stdout.write("No intersection between test lists; exiting")
