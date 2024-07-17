#!/usr/bin/env python
# AUTHOR: WSN
# CREATE DATE: 17 July 2024
import sys

###########################################################################
# MAIN
###########################################################################
def main():

    # Read all the match percentages from species-specific files.
    match_df = {} # Key = species, value = list of match percentages.
    for match_filename in sys.argv[1:]:
        species = match_filename.split(".")[1] # match_by.<species>.txt
        with open(match_filename, "r") as match_file:
            match_file.readline()
            match_df[species] = [line.rstrip().split("\t")[2]
                                 for line in match_file]
            print(f"Read {len(match_df[species])} match percentages from " +
                  f"{species}.", file=sys.stderr)
            

    # Make the histogram.
    fig, axs = plt.subplots(len(match_df), 1, sharex=True, tight_layout=True)
    i = 0
    for species in match_df.keys():
        axs[i].hist(match_df[species], bins=100)
        i += 1
    plt.show()

if __name__ == "__main__":
    main()
