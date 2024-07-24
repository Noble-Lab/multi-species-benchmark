#!/usr/bin/env python
# AUTHOR: WSN
# CREATE DATE: 17 July 2024
import sys
import pandas
import re
import os.path
import matplotlib.pyplot as plt

###########################################################################
# MAIN
###########################################################################
def main():

    # Check the command line.
    if len(sys.argv) < 3:
        print("USAGE: make_fdr_curves.py <plot> <percolator>+",
              file=sys.stderr)
        sys.exit(1)

    fig = plt.figure(figsize=(4,8))
    ax = plt.gca()
    ax.set_xlim(0, 0.1)
    ax.set_ylim(0, 2e6)
    plt.xlabel("FDR threshold")
    plt.ylabel("Accepted PSMs")

    i = 0
    for percolator_filename in sys.argv[2:]:
        species = re.sub(
            "-", " ", os.path.split(os.path.dirname(percolator_filename))[-1]
        )
        print(f"Reading {species}.", file=sys.stderr)

        percolator_df = pandas.read_csv(percolator_filename, sep="\t")
        qvalues = percolator_df["q-value"].values.tolist()

        plt.plot(qvalues, range(len(qvalues)), label=species)
        i += 1

    plt.legend(loc='center right')
    plt.savefig(sys.argv[1])

if __name__ == "__main__":
    main()
