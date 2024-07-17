#!/usr/bin/env python
# AUTHOR: William Stafford Noble
# CREATE DATE: 10 May 2024, Branford, CT
import sys
import argparse
import pandas
import os
import glob
import re

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
    if os.path.exists(filename):
        with open(filename, "r") as my_file:
            return(len(my_file.readlines()))
    return 0

def i2l(peptide):
    "Leucines to isoleucines."
    return re.sub("I", "L", peptide)

def clean_peptide (peptide):
    "Remove PTMs."

    return i2l(re.sub(r"[0-9\.\[\]]+", "", peptide))
    

def count_peptides(mgf_filenames):
    """Count number of distinct peptides in a set of MGF files.  Eliminates
    modifications and converts isoleucines first."""

    peptides = {} # Key = peptide, value = True
    for mgf_filename in mgf_filenames:
        with open(mgf_filename, "r") as mgf_file:
            for line in mgf_file:
                if line.startswith("SEQ="):
                    peptides[clean_peptide(line.rstrip().split("=")[1])] = True

    return len(peptides)
                            
    
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

    output = {} # Key = column name, value = list of entries.
    output["species"] = all_species
    column_names = ["#raw", "#mgf", "#spectra", "#PSMs", "#peptides"]
    for stat in column_names:
        output[stat] = []

    # Gather up all the statistics
    for species in all_species:
        print(species, file=sys.stderr)

        data_mgfs = glob.glob(os.path.join(args.data_dir, species, "*.mgf"))
        benchmark_mgfs = glob.glob(
            os.path.join(args.benchmark_dir, species, "*.mgf")
        )
        
        output["#raw"].append(len(data_mgfs))
        output["#mgf"].append(len(benchmark_mgfs))
        output["#spectra"].append(count_spectra(data_mgfs))
        output["#PSMs"].append(count_spectra(benchmark_mgfs))
        output["#peptides"].append(count_peptides(benchmark_mgfs))

    output["precursor"] = list(driver[2])
    output["fragment"] = list(driver[3])

    # Add totals.
    output["species"].append("total")
    for stat in column_names:
        output[stat].append(sum(output[stat]))
    output["precursor"].append("")
    output["fragment"].append("")

    # Convert output to a dataframe.
    output_df = pandas.DataFrame(output)
    output_df.to_csv(f"{args.root}.txt", index=False, sep="\t")
    output_df.to_html(f"{args.root}.html", index=False)

if __name__ == "__main__":
    main()
    sys.exit(0)
