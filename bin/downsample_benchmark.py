#!/usr/bin/env python
# AUTHOR: William Stafford Noble
# CREATE DATE: 10 May 2024, Branford, CT
import sys
import argparse
import random
import pyteomics.mgf
import os
import glob
import shutil

DESCRIPTION = """Select from a collection of sets of MGF files so that
each set has approximately the same number of spectra. This is done by
randomly permuting the order of MGFs within each set and then
accepting files in order until at least a specified number of spectra
is obtained."""

def count_spectra(mgf_filename):
    """Count the number of spectra in an MGF file."""
    return(len(pyteomics.mgf.read(mgf_filename)))

###########################################################################
# MAIN
###########################################################################
def main():
    global DESCRIPTION

    # Parse the command line.
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('mgf_dirs', nargs='+', help="MGF directories")
    parser.add_argument('--num_spectra', type=int, required=False,
                        default=100000, help="Target number of spectra per set")
    parser.add_argument('--seed', type=int, required=False, default=7718,
                        help="Random number seed")
    parser.add_argument('--root', type=str, required=True, 
                        help="Output directory name")
    args = parser.parse_args()

    random.seed(args.seed)
    os.makedirs(args.root, exist_ok=True)

    for mgf_dir in args.mgf_dirs:

        out_dir = os.path.join(args.root, os.path.basename(mgf_dir))
        # N.B. Hack to skip existing files.
        if os.path.isfile(out_dir):
            print(f"Skipping {mgf_dir}.", file=sys.stderr)
            continue
        os.makedirs(out_dir, exist_ok=True)
        mgf_filenames = glob.glob(os.path.join(mgf_dir, "*.mgf"))
        random.shuffle(mgf_filenames)
        total_spectra = 0
        for mgf_filename in mgf_filenames:
            new_filename = os.path.join(args.root,
                                        os.path.basename(mgf_dir),
                                        os.path.basename(mgf_filename))
            shutil.copy(mgf_filename, new_filename)
            num_spectra = count_spectra(mgf_filename)
            print(f"Read {num_spectra} from {mgf_filename}.")
            total_spectra += num_spectra
            if (total_spectra > args.num_spectra):
                print(f"Limit reached ({total_spectra} > {args.num_spectra})"
                      + f" for {mgf_dir}.",
                      file=sys.stderr)
                break

if __name__ == "__main__":
    main()
    sys.exit(0)
