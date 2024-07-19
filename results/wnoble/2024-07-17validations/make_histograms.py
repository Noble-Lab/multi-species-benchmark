#!/usr/bin/env python
# AUTHOR: WSN
# CREATE DATE: 17 July 2024
import sys
import matplotlib.pyplot as plt

###########################################################################
# MAIN
###########################################################################
def main():

    # Check the command line.
    if len(sys.argv) < 3:
        print("USAGE: make_histograms <output plot> <mgf>+", file=sys.stderr)
        sys.exit(1)
    
    # Read all the match percentages from species-specific files.
    match_df = {} # Key = species, value = list of match percentages.
    for match_filename in sys.argv[2:]:
        species = match_filename.split(".")[1] # match_by.<species>.txt
        with open(match_filename, "r") as match_file:
            match_file.readline()
            match_df[species] = [float(line.rstrip().split("\t")[2])
                                 for line in match_file]
            print(f"Read {len(match_df[species])} match percentages from " +
                  f"{species}.", file=sys.stderr)

    # Make the histogram.
    fig, axs = plt.subplots(len(match_df), 1, sharex=True, tight_layout=True)
    i = 0
    for species in match_df.keys():
        axs[i].hist(match_df[species], density=True, bins=100)
        axs[i].set_title(species, loc='left')
        axs[i].spines['top'].set_visible(False)
        axs[i].spines['right'].set_visible(False)
        i += 1
    plt.savefig(sys.argv[1])

if __name__ == "__main__":
    main()
