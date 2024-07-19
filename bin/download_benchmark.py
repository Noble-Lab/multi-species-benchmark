#!/usr/bin/env python
# AUTHOR: William Stafford Noble
# CREATE DATE: 30 Aug 2022
import sys
import datetime
import os
import ppx
import subprocess

USAGE = """USAGE: construct-benchmark.py <listing> <output dir>

Download the raw files associated with the DeepNovo benchmark dataset
from Tran et al., PNAS, 2017.  The input listing is a two-column text
file containing PRIDE IDs and species names.  The raw files are
downloaded into the specified output directory in subdirectories named
by species.
"""

# The location of the Thermo RAW file parser.
THERMORFP = "/Users/wnoble/opt/miniconda3/bin/thermorawfileparser"
if not os.path.exists(THERMORFP):
    THERMORFP = "~/proj/thermorawfileparser/ThermoRawFileParser.exe"
if not os.path.exists(THERMORFP):
    THERMORFP = "~/software/ThermoRawFileParser.v1.3.4/ThermoRawFileParser.exe"

#############################################################################
def run_command(command, outputFileName):
    """Run a command with error checking."""

    # Skip the command if the output file already exists.
    if (outputFileName != "") and (os.path.exists(outputFileName)):
        sys.stderr.write("%s exists.\n" % outputFileName)
        return
  
    sys.stderr.write("RUN: %s\n" % command)
    try:
        returnCode = subprocess.call(command, shell=True)
        if (returnCode != 0):
            sys.stderr.write("Child was terminated by signal %d\n" 
                             % -returnCode)
#            sys.exit(1)
    except OSError as e:
        sys.stderr.write("Execution failed: %s\n" % e)
        sys.exit(1)

###########################################################################
# MAIN
###########################################################################
def main():
    global USAGE

    # Parse the command line.
    if (len(sys.argv) != 3):
        sys.stderr.write(USAGE)
        sys.exit(1)
    listing_filename = sys.argv[1]
    output_dir = sys.argv[2]

    # Read the PRIDE IDs.
    prides = {} # Key = species, value = PRIDE identifier
    with open(listing_filename, "r") as listing_file:
        for line in listing_file:
            (pride, species) = line.rstrip().split()[:2]
            prides[species] = pride
    print(f"Read {len(prides)} entries from {listing_filename}.", 
          file=sys.stderr)

    # Make a small README in the output directory.
    with open(os.path.join(output_dir, "README"), "w") as readme:
        readme.write("Created by construct-benchmark.py on " +
                     f"{datetime.datetime.now()}.\n")
        readme.write("https://github.com/Noble-Lab/2021_melih_ms-chimera/tree/main/bin\n")
              

    # Download all raw files from each project.
    for species in prides.keys():
        print(f"Downloading {species} from {prides[species]}.",
              file=sys.stderr)
        local_dir = os.path.join(output_dir, species)
        os.makedirs(local_dir, exist_ok=True)

        # Locate this project on PRIDE.
        proj = ppx.find_project(prides[species], local=local_dir)
        try:
            raw_files = proj.remote_files("*.raw")
        except:
            sys.stderr.write(f"Error accessing {prides[species]}.\n")
            continue
        print(f"Downloading {len(raw_files)} files from {prides[species]}.",
              file=sys.stderr)

        # Traverse all the raw files in this project.
        for raw_file in raw_files:
            local_raw_file = f"{local_dir}/{raw_file}"

            # Skip download if file already exists.
            local_raw_file = os.path.join(local_dir, raw_file)
            mgf_file = os.path.join(local_dir, raw_file[:-3] + "mgf")
            if ( os.path.exists(mgf_file) or os.path.exists(local_raw_file) ):
                sys.stderr.write(f"Skipping {raw_file}.\n")
            else:
                sys.stderr.write(f"Downloading {raw_file}.")
                try:
                    proj.download(raw_file)
                except:
                    sys.stderr.write(f"Can't download {raw_file}.\n")
                    continue

            # Convert raw file to mgf
            if os.path.exists(mgf_file):
                print(f"Skipping conversion of {raw_file}.", file=sys.stderr)
            else:
                command = f"mono {THERMORFP} --output_file={mgf_file} "
                command += f"--logging=0 --format=0 --input={local_raw_file}"
                run_command(command, mgf_file)
                if os.path.exists(mgf_file):
                    os.remove(local_raw_file)
                else:
                    print(f"Failed to convert {local_raw_file}.",
                          file=sys.stderr)


if __name__ == '__main__':
    main()
