package main.java.util;

import com.compomics.util.experiment.mass_spectrometry.spectra.Precursor;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;

public class JSpectrum extends Spectrum {

    public String title;
    public String SEQ;
    public int SCANS;

    public JSpectrum(Precursor precursor, double[] mzArray, double[] intensityArray) {
        super(precursor, mzArray, intensityArray);
    }
}
