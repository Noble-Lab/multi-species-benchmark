#!/usr/bin/env python
# AUTHOR: WSN
# CREATE DATE: 17 July 2024
import sys
import pandas
import re
import matplotlib.pyplot as plt

###########################################################################
# MAIN
###########################################################################
def main():

    # Check the command line.
    if len(sys.argv) != 3:
        print("USAGE: make_bar_chart <TSV file> <plot> ", file=sys.stderr)
        sys.exit(1)

    # Read the statistics into a dataframe.
    stats_df = pandas.read_csv(sys.argv[1], sep="\t", index_col=0)
    species = list(stats_df.index)[:-1] # Remove total in last column.
    print(species, file=sys.stderr)

    # Remove hyphens.
    species = [re.sub("-", " ", one_species) for one_species in species]
    print(species, file=sys.stderr)

    num_psms = stats_df[["#PSMs"]].values.tolist()
    num_spectra = stats_df[["#spectra"]].values.tolist()
    percents = [100 * (my_tuple[0][0]/my_tuple[1][0])
                for my_tuple in zip(num_psms, num_spectra)][:-1] # Remove total
    print(percents, file=sys.stderr)

    # Make the barchart.
    fig, ax = plt.subplots()
    fig.set_size_inches(4,8)
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
              'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive']
    ax.barh(species, percents, color=colors[::-1])
    ax.set_xlabel("Percent identified")
    ax.set_yticklabels([])
    plt.tight_layout()
    plt.savefig(sys.argv[2])

if __name__ == "__main__":
    main()
