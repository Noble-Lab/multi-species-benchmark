#!/usr/bin/env python
# AUTHOR: William Stafford Noble
# CREATE DATE: 30 Aug 2022
import sys

USAGE = """USAGE: annotate_mgf.py <threshold> <log> <percolator> <mgf>

Use Percolator PSM-level output to annotate a given MGF file,
including only PSMs with q-value better than the specified threshold.

Inputs are
 o the desired FDR threshold,
 o the log file produced by Crux percolator,
 o the PSM-level Percolator output file, and
 o the MGF file to be annotated.

The resulting annotated MGF is printed to stdout.

Note that spectra with no assigned charge state are skipped.

Searches the log file for the assigned file index in a line like this:

INFO: Assigning index <integer> to <file>.

Assumes that the Percolator file has PSM IDs of the form
<foo>_<file-index>_<scan-number>

"""

###########################################################################
# MAIN
###########################################################################
def main():
    global USAGE

    # Parse the command line.
    if (len(sys.argv) != 5):
        sys.stderr.write(USAGE)
        sys.exit(1)
    fdr_threshold = float(sys.argv[1])
    log_filename = sys.argv[2]
    percolator_filename = sys.argv[3]
    mgf_filename = sys.argv[4]

    # Get the file index from the log file.
    target_file_index = None
    with open(log_filename, "r") as log_file:
        for line in log_file:
            words = line.rstrip().split()
            if ( (len(words) == 6) and
                 (words[0] == "INFO:") and
                 (words[1] == "Assigning") and
                 (words[2] == "index") and
                 (words[5][:-1] == mgf_filename) ): # Remove trailing period.
                target_file_index = int(words[3])
                break
    if (target_file_index == None):
        print(f"Cannot find {mgf_filename} in {log_filename}.", file=sys.stderr)
        sys.exit(1)
    print(f"File index: {target_file_index}.", file=sys.stderr)

    # Read the Percolator PSMs into a dictionary.
    psms = {} # Key = scan number, value = peptide
    with open(percolator_filename, "r") as percolator_file:

        header = percolator_file.readline().rstrip().split("\t")
        psm_id_col = header.index("PSMId")
        qvalue_col = header.index("q-value")
        peptide_col = header.index("peptide")

        for line in percolator_file:
            words = line.rstrip().split("\t")
            psm_id = words[psm_id_col].split("_")
            file_index = int(psm_id[1])
            scan_number = int(psm_id[2])
            qvalue = float(words[qvalue_col])
            peptide = words[peptide_col]

            # Skip peptides containing pyrrolysine or selenocysteine.
            if "O" in peptide or "U" in peptide:
                print(f"Skipping {peptide} with non-canonical amino acid.",
                      file=sys.stderr)
                continue

            if ( (file_index == target_file_index) and 
                 (qvalue <= fdr_threshold) ):
                psms[scan_number] = peptide
    print(f"Read {len(psms)} PSMs from {percolator_filename} " +
          f"at q<{fdr_threshold}.", file=sys.stderr)

    # Write the annotated MGF.
    num_printed = 0
    with open(mgf_filename, "r") as mgf_file:
        this_spectrum = ""
        print_me = False
        has_charge = False
        has_mass = False
        for line in mgf_file:

            # Starting a new spectrum? Then print and reset.
            if (line == "BEGIN IONS\n"):
                if print_me and has_charge and has_mass:
                    print(this_spectrum)
                    num_printed += 1
                print_me = False
                has_charge = False
                has_mass = False
                this_spectrum = line
            else:

                # Check for charge state.
                if (line[:6] == "CHARGE"):
                    has_charge = True
                if (line[:7] == "PEPMASS"):
                    has_mass = True

                # Add to this spectrum, possibly including annotation.
                this_spectrum += line
                words = line.split("=")
                if (words[0] == "SCANS"):
                    scan_number = int(words[1])
                    if (scan_number in psms):
                        # Strip flanking amino acids.
                        this_spectrum += f"SEQ={psms[scan_number][2:-2]}\n"
                        print_me = True

        # Print the final spectrum.
        if print_me and has_charge and has_mass:
            print(this_spectrum, end='')
            num_printed += 1
    print(f"Printed {num_printed} PSMs.", file=sys.stderr)

if __name__ == '__main__':
    main()
