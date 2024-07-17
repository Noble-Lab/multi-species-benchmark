#!/usr/bin/env python
import sys
import argparse
import matplotlib.pyplot as plt
import spectrum_utils.spectrum
import spectrum_utils.plot
import pyteomics.mgf

DESCRIPTION = """Given an annotated MGF, compute the proportion of intensity
matched to b- and y-ions.  Prints a list of values to stdout."""

def get_percent_matched(spectrum):
    
    # Extract the peptide sequence.
    peptide = spectrum['params']['seq']

    # Convert from pyteomics to spectrum_utils format.
    spectrum = spectrum_utils.spectrum.MsmsSpectrum(
        spectrum['params']['title'],
        float(spectrum['params']['pepmass'][0]), 
        int(spectrum['params']['charge'][0]), 
        spectrum['m/z array'], 
        spectrum['intensity array']
    )
    
    # Annotate the spectrum.
    spectrum = (
        spectrum.set_mz_range(min_mz=100, max_mz=1400)
        .remove_precursor_peak(0.05, "Da")
        .filter_intensity(min_intensity=0.05, max_num_peaks=50)
        .scale_intensity("root")
        .annotate_proforma(peptide, 0.05, "Da", ion_types="by")
    )

    # Compute percent matched intensity.
    total_intensity = 0.0
    matched_intensity = 0.0
    for annotation, intensity in zip(
            spectrum.annotation, spectrum.intensity
    ):
        total_intensity += intensity
        if len(annotation.fragment_annotations) > 0:
#            print(f"Matching {annotation} with intensity {intensity}.")
            matched_intensity += intensity
#        else:
#            print(f"No matching for intensity {intensity}.")

    return 100 * matched_intensity / total_intensity

###########################################################################
# MAIN
###########################################################################
def main():

    # Parse the command line.
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--mgf', type=str, required=True,
                        help="Annotated MGF file")
    args = parser.parse_args()

    with pyteomics.mgf.read(args.mgf, use_index=False) as reader:
        for mgf_index, spectrum in enumerate(reader):
            print(f"{get_percent_matched(spectrum):.4f}")

if __name__ == "__main__":
    main()

#############################################################################
# TESTING
#############################################################################
import pytest
import tempfile

# Mass of +1 b1 ion is m(V) + m(proton)
# 99.06841391299 + 1.00727646677 = 100.0757

# Mass of +1 y1 ion is m(K) + m(H2O) + m(proton)
# 128.09496301399997 + (1.00782503207 + 17.00274) + 1.00727646677 = 147.1128
# https://github.com/crux-toolkit/crux-toolkit/blob/master/src/util/mass.h

# Just the y1 and b1 ions, plus one other peak.
bitty_mgf = """BEGIN IONS
TITLE=controllerType=0 controllerNumber=1 scan=2475
PEPMASS=561.798400878906
RTINSECONDS=832.4328
CHARGE=2+
SCANS=2475
SEQ=VVQEQGTHPK
100.0757 10
147.1124878 10
1033.7346191 10
END IONS
"""

my_mgf = """BEGIN IONS
TITLE=controllerType=0 controllerNumber=1 scan=2475
PEPMASS=561.798400878906
RTINSECONDS=832.4328
CHARGE=2+
SCANS=2475
SEQ=VVQEQGTHPK
84.0448303 3936.5886230469 
84.0814972 2155.9631347656 
97.1241684 675.0009765625 
101.0713272 8273.31640625 
106.5969009 551.2675170898 
110.071434 3480.9711914063 
112.0875931 655.4925537109 
117.1025467 787.0622558594 
124.0218506 573.0434570313 
129.1023407 4307.3549804688 
136.0753784 513.6068725586 
147.1124878 2242.1789550781 
171.1489716 17505.900390625 
175.118515 1430.5400390625 
199.1435089 5907.8701171875 
226.154892 1230.7514648438 
228.1339722 1189.3963623047 
230.1139221 730.3195800781 
238.2206421 578.6478271484 
239.6330261 573.1657104492 
241.0816956 1016.6735229492 
244.1650391 8036.1098632813 
327.2015991 2643.5959472656 
329.436676 615.9159545898 
369.1423035 615.8500976563 
381.2232971 1044.3002929688 
445.2155762 714.327331543 
445.7083435 1164.3920898438 
453.7232666 4702.521484375 
454.2202454 5765.0966796875 
454.717926 1017.1290893555 
462.7281799 1381.1580810547 
474.7731628 584.8037719727 
482.273468 749.9706420898 
512.2681274 999.9737548828 
539.2930298 3000.0100097656 
553.2305298 773.0104980469 
560.7866211 1212.0534667969 
561.2929688 10707.52734375 
561.829834 1013.0922241211 
562.2373657 1954.6754150391 
562.3233032 6883.9736328125 
562.7788696 1303.7845458984 
650.3289795 1264.9660644531 
664.2758179 1195.4071044922 
667.3475342 1121.3416748047 
778.3859253 1206.9233398438 
796.3920288 6032.0810546875 
797.3983765 848.6907958984 
802.911499 672.137878418 
906.4327393 937.2236938477 
907.4239502 4685.7861328125 
924.4491577 5191.5551757813 
925.4563599 2268.0749511719 
1033.7346191 725.4644775391 
END IONS
"""

def test_get_percent_matched():
    global my_peptide, bitty_mgf, my_mgf

    for sample_mgf in [bitty_mgf, my_mgf]:
    
        my_tempfile = tempfile.NamedTemporaryFile(delete=False, mode='w')
        my_tempfile.write(sample_mgf)
        my_tempfile.close()

        # There is just one spectrum in there.
        with pyteomics.mgf.read(my_tempfile.name, use_index=False) as reader:
            for mgf_index, spectrum in enumerate(reader):
                percent_matched = get_percent_matched(spectrum)
                assert((percent_matched == pytest.approx(2/3)) # bitty_mgf
                       or
                       # my_mgf
                       (percent_matched == pytest.approx(0.28868434441953056)))
            
