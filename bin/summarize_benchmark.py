#!/usr/bin/env python
# AUTHOR: William Stafford Noble
# CREATE DATE: 10 May 2024, Branford, CT
import sys
import argparse
import pandas
import os
import glob

DESCRIPTION = """Write text and HTML tables describing a given
benchmark dataset."""

def count_spectra(mgf_filenames):
    "Count number of spectra in a set of MGF files."
    return_value = 0

    for mgf_filename in mgf_filenames:
        with open(mgf_filename, "r") as mgf_file:
            for line in mgf_file:
                if line.rstrip() == "BEGIN IONS":
                    return_value += 1

    return return_value

def count_lines(filename):
    with open(filename, "r") as my_file:
        return(len(my_file.readlines()))

###########################################################################
# MAIN
###########################################################################
def main():
    global DESCRIPTION

    # Parse the command line.
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--driver_filename', type=str, required=True,
                        help="Name of driver file.  Columns are PXD ID, "
                        + "species name, precursor tolerance (ppm), fragment "
                        + "bin size (m/z)")
    parser.add_argument('--data_dir', type=str, required=True,
                        help="Directory containing original data")
    parser.add_argument('--benchmark_dir', type=str, required=True,
                        help="Directory containing benchmark")
    parser.add_argument('--root', type=str, required=True, 
                        help="Output directory and filename root")
    args = parser.parse_args()

    driver = pandas.read_csv(args.driver_filename, sep="\t", header=None)
    all_species = list(driver[1])

    column_names = ["#raw", "#mgf", "#spectra", "#PSMs", "#peptides"]
    output = {} # Key = column name, value = list of entries.
    output["species"] = all_species
    for stat in column_names:
        output[stat] = []

    # Gather up all the statistics
    for species in all_species:
        for stat in column_names:
            if stat == "#raw":
                output[stat].append(
                    len(glob.glob(os.path.join(args.data_dir,
                                               species, "*.mgf")))
                )
            elif stat == "#mgf":
                output[stat].append(
                    len(glob.glob(os.path.join(args.benchmark_dir,
                                               species, "*.mgf")))
                )
            elif stat == "#spectra":
                output[stat].append(
                    count_spectra(glob.glob(os.path.join(args.data_dir,
                                                         species, "*.mgf")))
                )
            elif stat == "#PSMs":
                output[stat].append(
                    count_spectra(glob.glob(os.path.join(args.benchmark_dir,
                                                         species, "*.mgf")))
                )
            elif stat == "#peptides":
                peptide_filename = os.path.join(args.benchmark_dir,
                                                species, "peptides.txt")
                try:
                    num_peptides = count_lines(peptide_filename)
                except:
                    num_peptides = 0
                output[stat].append(num_peptides)

    output["precursor"] = list(driver[2])
    output["fragment"] = list(driver[3])

    # Convert output to a dataframe.
    output_df = pandas.DataFrame(output)
    output_df.to_csv(f"{args.root}.txt", index=False, sep="\t")
    output_df.to_html(f"{args.root}.html", index=False)
                

if __name__ == "__main__":
    main()
    sys.exit(0)
