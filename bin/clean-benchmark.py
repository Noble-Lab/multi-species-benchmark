#!/usr/bin/env python
# AUTHOR: William Stafford Noble
# CREATE DATE: 23 Sep 2022
import sys
import argparse
import os
import random
import re
from pyteomics import mgf

DESCRIPTION = """
Clean up the 9-species benchmark to eliminate peptides that are shared
between species.  The old root directory is assumed to contain one
subdirectory per species, with annotated MGFs in each subdirectory.
The script creates a copy of that directory structure in the specified
new location.  For any peptide that appears in more than one species,
one species is selected at random, and all PSMs associated with that
peptide in other species are removed.

At the same time, convert PTMs from Tide format to Casanovo format.
"""

def convert_ptms (peptide):
    "Convert PTMs in a peptide from Tide to Casanovo format."

    # Handle internal modifications.
    mod_dict = {"C":"C+57.021",
                "M[15.9949]":"M+15.995",
                "[0.9840]":"+0.984"}
    for dn_mod in mod_dict.keys():
        peptide = peptide.replace(dn_mod, mod_dict[dn_mod])

    # Handle terminal modifications.
    n_term_mods = { "[-17.0265]":"-17.027",
                    "[42.0106]":"+42.011",
                    "[43.0058]":"+43.006",
                    "[25.9803]":"+43.006-17.027"}
    for n_mod in n_term_mods.keys():
        if n_mod in peptide:
            peptide = n_term_mods[n_mod]+peptide.replace(n_mod, "")

    return peptide

def i2l(peptide):
    "Leucines to isoleucines."
    return re.sub("I", "L", peptide)

def clean_peptide (peptide, do_i2l):
    "Remove flanking amino acids and PTMs."

    no_mods = re.sub(r"[0-9\.\[\]]+", "", peptide)

    if do_i2l:
        return i2l(no_mods)
    return no_mods

###########################################################################
# MAIN
###########################################################################
def main():
    global DESCRIPTION

    random.seed(7718)

    # Parse the command line.
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--old_root', type=str, required=True,
                        help="Directory of old benchmark")
    parser.add_argument('--new_root', type=str, required=True,
                        help="Directory of new benchmark")
    parser.add_argument('--i2l', action=argparse.BooleanOptionalAction,
                        help="Isoleucine to leucine")
    args = parser.parse_args()
    
    # Get the list of species.
    species_list = [ f.name for f in os.scandir(args.old_root) if f.is_dir() ]

    # Construct the peptide-to-species mapping.
    species_mapping = {} # key = peptide, value = list of species
    for species in species_list:
        print(f"Extracting peptides from {species}.", file=sys.stderr)
        mgf_list = [ f.name for f in os.scandir(os.path.join(args.old_root, 
                                                             species))
                     if f.is_file() ]
        for mgf_filename in mgf_list:
            with open(os.path.join(args.old_root, species, 
                                   mgf_filename), "r") as mgf_file:
                for line in mgf_file:
                    if (line[:4] == "SEQ="):
                        peptide = clean_peptide(line[4:-1], args.i2l)
                        if (peptide not in species_mapping):
                            species_mapping[peptide] = [species]
                        elif (species not in species_mapping[peptide]):
                            species_mapping[peptide].append(species)
        print(f"Found {len(species_mapping)} peptides.", file=sys.stderr)

    # If a peptide appears in more than one species, select one
    # species randomly.
    num_duplicates = 0
    eliminated_psms = 0
    for peptide in species_mapping.keys():
        if (len(species_mapping[peptide]) > 1):
            selected_species = random.choice(species_mapping[peptide])
            species_mapping[peptide] = selected_species
            num_duplicates += 1
        else:
            species_mapping[peptide] = species_mapping[peptide][0]
    print(f"Found {num_duplicates} duplicated peptides.", file=sys.stderr)

    # Create the new benchmark.
    for species in species_list:
        os.makedirs(os.path.join(args.new_root, species), exist_ok=True)

        # Print the peptides for that species.
        with open(os.path.join(args.new_root, species, "peptides.txt"), "w") \
             as peptide_filename:
            num_peptides = 0
            for peptide in species_mapping.keys():
                if species_mapping[peptide] == species:
                    peptide_filename.write(f"{peptide}\n")
                    num_peptides += 1
        print(f"{num_peptides} distinct peptides in {species}.", 
              file=sys.stderr)

        print(f"Creating cleaned MGFs for {species}.", file=sys.stderr)
        mgf_list = [ f.name for f in os.scandir(
            os.path.join(args.old_root, species)
        )
                     if f.is_file() ]
        for mgf_file in mgf_list:
            old_mgf_filename = os.path.join(args.old_root, species, mgf_file)
            new_mgf_filename = os.path.join(args.new_root, species, mgf_file)
            print(f"Creating {mgf_file}.", file=sys.stderr)
            num_printed = 0
            num_skipped = 0
            with mgf.read(source=old_mgf_filename, use_header=False,
                          read_charges=False, use_index=False) as reader, \
                 open(new_mgf_filename, "w") as new_mgf:
                for spectrum in reader:
                    peptide = spectrum['params']['seq']
                    unmod_peptide = clean_peptide(peptide, args.i2l)

                    try:
                        if (species_mapping[unmod_peptide] == species):
                            # Convert PTM format.
                            spectrum['params']['seq'] = convert_ptms(peptide)
                            mgf.write([spectrum], output=new_mgf)
                            num_printed += 1
                        else:
                            num_skipped += 1
                    except:
                        print(f"Cannot find {peptide} from {mgf_file}.",
                              file=sys.stderr)
                        sys.exit(1)
                        
            print(f"Printed {num_printed} spectra and skipped {num_skipped}.",
                  file=sys.stderr)

if __name__ == '__main__':
    main()

#############################################################################
# TESTING
#############################################################################
import pytest
import tempfile

def test_convert_ptms():

    assert("+43.006HHVLHHQTVDK" == convert_ptms("H[43.0058]HVLHHQTVDK"))
    assert("+43.006IIQ+0.984N+0.984AYK"
           == convert_ptms("I[43.0058]IQ[0.9840]N[0.9840]AYK"))
