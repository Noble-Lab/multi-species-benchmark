#!/usr/bin/env python
import sys
import argparse
import matplotlib.pyplot as plt
import spectrum_utils.spectrum
import spectrum_utils.plot
import pyteomics.mgf

DESCRIPTION = """Given an annotated MGF, compute the proportion of intensity
matched to b- and y-ions.  Produces a tsv file (<root>.txt) and a plot
(<root>.pdf)."""

###########################################################################
# MAIN
###########################################################################
def main():

    # Parse the command line.
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--mgf', type=str, required=True,
                        help="Annotated MGF file")
    parser.add_argument('--root', type=str, required=True,
                        help="Output file root")
    args = parser.parse_args()

    title = f"controllerType=0 controllerNumber=1 scan={args.scan}"
    spectrum = pyteomics.mgf.get_spectrum(args.mgf, title)

    # Convert from pyteomics to spectrum_utils format.
    spectrum = spectrum_utils.spectrum.MsmsSpectrum(
        title,
        float(spectrum["params"]["pepmass"][0]), 
        int(spectrum["params"]["charge"][0]), 
        spectrum["m/z array"], 
       spectrum["intensity array"]
     )

    # Annotate the spectrum.
    spectrum = (
        spectrum.set_mz_range(min_mz=100, max_mz=1400)
        .remove_precursor_peak(args.tolerance, "ppm")
        .filter_intensity(min_intensity=0.05, max_num_peaks=50)
        .scale_intensity("root")
        .annotate_proforma(
            args.peptide, args.tolerance, "ppm", ion_types="aby"
        )
    )


if __name__ == "__main__":
    main()
