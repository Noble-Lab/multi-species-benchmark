package main.java.util;


import com.compomics.util.experiment.biology.ions.Ion;
import com.compomics.util.experiment.biology.ions.NeutralLoss;
import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.biology.modifications.Modification;
import com.compomics.util.experiment.biology.modifications.ModificationCategory;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.biology.modifications.ModificationType;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.matches.IonMatch;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.spectrum_annotation.AnnotationParameters;
import com.compomics.util.experiment.identification.spectrum_annotation.SpecificAnnotationParameters;
import com.compomics.util.experiment.identification.spectrum_annotation.SpectrumAnnotator;
import com.compomics.util.experiment.identification.spectrum_annotation.spectrum_annotators.PeptideSpectrumAnnotator;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.gui.waiting.waitinghandlers.WaitingHandlerCLIImpl;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.parameters.identification.search.ModificationParameters;
import com.compomics.util.waiting.WaitingHandler;
import org.apache.commons.lang3.StringUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MGFUtils {

    public double itol = 0.05;
    public String itolu = "Da";
    public boolean lossWaterNH3 = true;
    public double fragment_ion_intensity_cutoff = 0.01;
    public String out_dir = "./";
    public boolean export_input_for_pdv = true;

    public static void main(String[] args) throws IOException {
        MGFUtils mgfUtils = new MGFUtils();
        File F = new File(args[0]);
        List<File> files = new ArrayList<>();
        if(F.isDirectory()){
            mgfUtils.out_dir = args[1];
            files.addAll(getMgfFiles(args[0]));
            for (File file : files) {
                System.out.println("Processing " + file.getAbsolutePath());
                String out_mzlib = mgfUtils.out_dir+"/"+file.getName().replaceAll(".mgf$", ".mzlib.txt");
                mgfUtils.mgf2mzSpecLib(file.getAbsolutePath(), out_mzlib);
            }
        }else{
            System.out.println("Processing " + F.getAbsolutePath());
            String out_mzlib = args[1];
            System.out.println("Output library " + out_mzlib);
            File OF = new File(args[1]);
            OF = new File(OF.getAbsolutePath());
            mgfUtils.out_dir = OF.getParent();
            mgfUtils.mgf2mzSpecLib(F.getAbsolutePath(), out_mzlib);
        }
    }

    public static List<File> getMgfFiles(String inputFolder) {
        List<File> mgfFiles = new ArrayList<>();
        File folder = new File(inputFolder);
        if (folder.isDirectory()) {
            File[] files = folder.listFiles();
            if (files != null) {
                for (File file : files) {
                    if (file.isFile() && file.getName().endsWith(".mgf")) {
                        mgfFiles.add(file);
                    }
                }
            }
        }

        return mgfFiles;
    }

    public void mgf2mzSpecLib(String spectraFile, String mzSpecLibFile) throws IOException {
        File F = new File(spectraFile);
        BufferedWriter writer = new BufferedWriter(new FileWriter(mzSpecLibFile));
        BufferedWriter pdv_mgf_writer = null;
        BufferedWriter pdv_pep_writer = null;
        if(export_input_for_pdv){
            pdv_mgf_writer = new BufferedWriter(new FileWriter(out_dir +"/PDV_"+F.getName()));
            pdv_pep_writer = new BufferedWriter(new FileWriter(out_dir +"/PDV_"+F.getName().replaceAll(".mgf$", ".tsv")));
            pdv_pep_writer.write("peptide\tmodification\tspectrum_title\tcharge\n");
        }

        writer.write("<mzSpecLib 1.0>\n");

        writer.write("MS:1003188|library name="+F.getName().replaceAll(".mgf$", "")+"\n");
        File mgfFile = new File(spectraFile);
        WaitingHandler waitingHandler = new WaitingHandlerCLIImpl();
        waitingHandler.setDisplayProgress(false);

        JMgfFileIterator mgfFileIterator = new JMgfFileIterator(mgfFile, waitingHandler);
        String title;
        JSpectrum spectrum;
        int i = 0;
        String charge_str = "";
        while ((title = mgfFileIterator.next()) != null) {

            i++;
            spectrum = mgfFileIterator.getSpectrum();
            if(export_input_for_pdv && pdv_mgf_writer != null){
                pdv_mgf_writer.write(asMgf(spectrum, String.valueOf(i), spectrum.getPrecursor().possibleCharges[0], spectrum.SCANS)+"\n");
            }
            // copy the mz and intensity to double array
            // because the mz and intensity in JSpectrum may be changed in PSM annotation
            double[] mz = new double[spectrum.mz.length];
            double[] intensity = new double[spectrum.intensity.length];
            for(int j=0;j<spectrum.mz.length;j++){
                mz[j] = spectrum.mz[j];
                intensity[j] = spectrum.intensity[j];
            }
            Peptide modPeptide = getPeptide(spectrum);
            String header = get_spectrum_header4mzSpecLib(spectrum, modPeptide, i);
            writer.write(header);
            HashMap<Double, JPeak> fragment_mz2peak = psm_annotate(spectrum,modPeptide);

            writer.write("<Interpretation=1>\n");
            writer.write("MS:1001121|number of matched peaks="+fragment_mz2peak.size()+"\n");
            writer.write("MS:1001362|number of unmatched peaks="+(mz.length-fragment_mz2peak.size())+"\n");
            writer.write("<Peaks>\n");
            for(int k=0;k<mz.length;k++){
                if(fragment_mz2peak.containsKey(mz[k])){
                    if(fragment_mz2peak.get(mz[k]).charge > 1){
                        charge_str = "^"+fragment_mz2peak.get(mz[k]).charge;
                    }else{
                        charge_str = "";
                    }
                    writer.write(mz[k] + " " + intensity[k] + " "+fragment_mz2peak.get(mz[k]).ion_type+fragment_mz2peak.get(mz[k]).ion_number+fragment_mz2peak.get(mz[k]).neutral_loss_type+charge_str+"/"+String.format("%.4f",fragment_mz2peak.get(mz[k]).error)+"\n");
                }else{
                    writer.write(mz[k] + " " + intensity[k]+" ?\n");
                }
            }
            writer.write("\n");
            if(export_input_for_pdv && pdv_pep_writer != null){
                pdv_pep_writer.write(modPeptide.getSequence()+"\t"+getPDVModificationString(modPeptide)+"\t"+i+"\t"+spectrum.getPrecursor().possibleCharges[0]+"\n");
            }
        }
        writer.close();
        if(export_input_for_pdv && pdv_mgf_writer != null){
            pdv_mgf_writer.close();
            pdv_pep_writer.close();
        }
    }

    public Peptide getPeptide(JSpectrum spectrum){
        HashMap<Integer,Double> pos2mod = get_mod_from_peptide(spectrum.SEQ);
        String sequence = get_peptide(spectrum.SEQ);
        Peptide modPeptide = new Peptide(sequence);
        for (int pos: pos2mod.keySet()) {
            double monoMass = pos2mod.get(pos);
            if(pos == 0){
                add_modification("", monoMass, pos);
                Modification ptm = get_modification("", monoMass, pos);
                modPeptide.addVariableModification(new ModificationMatch(ptm.getName(), pos));
            }else{
                String aa = String.valueOf(sequence.charAt(pos - 1));
                add_modification(aa, monoMass, pos);
                Modification ptm = get_modification(aa, monoMass, pos);
                modPeptide.addVariableModification(new ModificationMatch(ptm.getName(), pos));
            }
        }
        return modPeptide;
    }

    public HashMap<Double, JPeak> psm_annotate(JSpectrum spectrum, Peptide modPeptide){
        HashMap<Double, JPeak> fragment_mz2peak = new HashMap<>();
        remove_precursor_ions(spectrum, 1.0);
        PeptideSpectrumAnnotator peptideSpectrumAnnotator = new PeptideSpectrumAnnotator();
        int precursor_charge = spectrum.getPrecursor().possibleCharges[0];
        PeptideAssumption peptideAssumption = new PeptideAssumption(modPeptide, precursor_charge);
        SpecificAnnotationParameters specificAnnotationPreferences = new SpecificAnnotationParameters();

        HashSet<Integer> charges = new HashSet<>(4);
        int precursorCharge = peptideAssumption.getIdentificationCharge();
        if (precursorCharge <= 1) {
            charges.add(precursorCharge);
        } else {
            int cur_max_fragment_ion_charge;
            if(precursorCharge==2){
                cur_max_fragment_ion_charge = 2;
            }else{
                cur_max_fragment_ion_charge = precursorCharge - 1;
            }
            for (int c = 1; c <= cur_max_fragment_ion_charge; c++) {
                charges.add(c);
            }
        }
        specificAnnotationPreferences.setSelectedCharges(charges);

        specificAnnotationPreferences.addIonType(Ion.IonType.PEPTIDE_FRAGMENT_ION, PeptideFragmentIon.B_ION);
        specificAnnotationPreferences.addIonType(Ion.IonType.PEPTIDE_FRAGMENT_ION, PeptideFragmentIon.Y_ION);
        specificAnnotationPreferences.setFragmentIonAccuracy(this.itol);
        specificAnnotationPreferences.setFragmentIonPpm(this.itolu.startsWith("ppm"));
        specificAnnotationPreferences.setNeutralLossesAuto(false);
        specificAnnotationPreferences.clearNeutralLosses();
        specificAnnotationPreferences.setPrecursorCharge(precursorCharge);

        if (lossWaterNH3) {
            specificAnnotationPreferences.addNeutralLoss(NeutralLoss.H2O);
            specificAnnotationPreferences.addNeutralLoss(NeutralLoss.NH3);
        }

        AnnotationParameters annotationSettings = new AnnotationParameters();
        annotationSettings.setTiesResolution(SpectrumAnnotator.TiesResolution.mostAccurateMz);
        annotationSettings.setFragmentIonAccuracy(this.itol);
        annotationSettings.setFragmentIonPpm(this.itolu.startsWith("p"));
        annotationSettings.setIntensityLimit(0.0001);
        annotationSettings.setNeutralLossesSequenceAuto(false);
        annotationSettings.setIntensityThresholdType(AnnotationParameters.IntensityThresholdType.percentile);


        ModificationParameters modificationParameters = new ModificationParameters();
        SequenceMatchingParameters sequenceMatchingParameters = new SequenceMatchingParameters();
        JSequenceProvider jSequenceProvider = new JSequenceProvider();
        IonMatch[] matched_ions = peptideSpectrumAnnotator.getSpectrumAnnotation(annotationSettings,
                specificAnnotationPreferences,
                "",
                "",
                spectrum,
                modPeptide,
                modificationParameters,
                jSequenceProvider,
                sequenceMatchingParameters);


        if (matched_ions == null || matched_ions.length == 0) {
            System.err.println("No ions matched!");
        }else{
            int ion_number;
            String ion_type = "-"; // b or y ion
            double max_fragment_ion_intensity = spectrum.getMaxIntensity();
            for (IonMatch ionMatch : matched_ions) {
                if (ionMatch.ion.getSubType() == PeptideFragmentIon.B_ION || ionMatch.ion.getSubType() == PeptideFragmentIon.Y_ION) {
                    PeptideFragmentIon fragmentIon = ((PeptideFragmentIon) ionMatch.ion);
                    ion_number = fragmentIon.getNumber();
                    if(ionMatch.ion.getSubType() == PeptideFragmentIon.Y_ION){
                        ion_type = "y";
                    }else if (ionMatch.ion.getSubType() == PeptideFragmentIon.B_ION){
                        ion_type = "b";
                    }else{
                        System.err.println("Unrecognized fragment ion type:"+ionMatch.ion.getSubType()+","+ionMatch.ion.getSubTypeAsString());
                        System.exit(1);
                    }

                    if(ionMatch.peakIntensity < this.fragment_ion_intensity_cutoff * max_fragment_ion_intensity){
                        continue;
                    }
                    JPeak peak = new JPeak();
                    peak.mz = ionMatch.peakMz;
                    peak.intensity = ionMatch.peakIntensity;
                    peak.ion_number = ion_number;
                    peak.ion_type = ion_type;
                    peak.charge = ionMatch.charge;
                    peak.error = ionMatch.getError(this.itolu.startsWith("ppm"));
                    fragment_mz2peak.put(ionMatch.peakMz, peak);
                    if(ionMatch.ion.hasNeutralLosses()){
                        peak.neutral_loss_type = ionMatch.ion.getNeutralLossesAsString();
                    }

                }
            }
        }
        return fragment_mz2peak;
    }

    private void remove_precursor_ions(JSpectrum spectrum, double mz_tol){
        ArrayList<Double> mz_list = new ArrayList<>();
        ArrayList<Double> intensity_list = new ArrayList<>();
        double precursor_mz = spectrum.getPrecursor().mz;
        for(int i=0;i<spectrum.mz.length;i++){
            if(Math.abs(spectrum.mz[i]-precursor_mz) > mz_tol){
                mz_list.add(spectrum.mz[i]);
                intensity_list.add(spectrum.intensity[i]);
            }
        }
        spectrum.mz = new double[mz_list.size()];
        spectrum.intensity = new double[intensity_list.size()];
        for(int i=0;i<mz_list.size();i++){
            spectrum.mz[i] = mz_list.get(i);
            spectrum.intensity[i] = intensity_list.get(i);
        }
    }

    static Pattern pattern = Pattern.compile("([A-Z]?)([-+]?\\d+\\.\\d+)");
    public static HashMap<Integer,Double> get_mod_from_peptide(String peptide){
        // -17.027QILLTQSPAIM+15.995SASPGQ+0.984K
        // +43.006-17.027IAHQTIAN+0.984MQAR: combine the two as one modification for annotation
        HashMap<Integer,Double> pos2mod = new HashMap<>();
        Matcher matcher = pattern.matcher(peptide);
        int sequencePosition = 0;
        int mod_length = 0;
        double mod_aa;
        while (matcher.find()) {
            String aminoAcid = matcher.group(1);
            String modificationName = matcher.group(2);
            if(aminoAcid.isEmpty()){
                sequencePosition = 0;
                mod_aa = Double.parseDouble(modificationName);
            }else{
                mod_aa = Double.parseDouble(modificationName);
                sequencePosition = peptide.indexOf(aminoAcid+modificationName, sequencePosition) + 1;
            }
            if(aminoAcid.isEmpty()){
                if(pos2mod.containsKey(sequencePosition)) {
                    pos2mod.put(sequencePosition, mod_aa+pos2mod.get(sequencePosition));
                }else{
                    pos2mod.put(sequencePosition, mod_aa);
                }
            }else{
                if(pos2mod.containsKey(sequencePosition-mod_length)){
                    pos2mod.put(sequencePosition - mod_length, mod_aa+pos2mod.get(sequencePosition - mod_length));
                }else {
                    pos2mod.put(sequencePosition - mod_length, mod_aa);
                }
            }
            mod_length = mod_length + modificationName.length();
        }
        return pos2mod;
    }

    private static String get_peptide(String mod_pep){
        return mod_pep.replaceAll("[-+]?\\d+\\.\\d+", "");
    }

    ModificationFactory ptmFactory = ModificationFactory.getInstance();
    private void add_modification(String aa, double monoMass, int pos){
        Modification ptm = null;
        String ptmName = get_modification_name(aa, monoMass, pos);
        if(pos==0){
            ptm = new Modification(ModificationType.modn_peptide,ptmName,monoMass,null, ModificationCategory.Other);
            ptm.setShortName(String.valueOf(monoMass));
        }else {
            ArrayList<String > residues = new ArrayList<>();
            residues.add(aa);
            ptm = new Modification(ModificationType.modaa, ptmName, monoMass, residues,ModificationCategory.Other);
            ptm.setShortName(String.valueOf(monoMass));
        }

        if(!ptmFactory.containsModification(ptm.getName())){
            ptmFactory.addUserModification(ptm);
        }
    }

    private Modification get_modification(String aa, double monoMass, int pos){
        return ptmFactory.getModification(get_modification_name(aa, monoMass, pos));
    }

    private String get_modification_name(String aa, double monoMass, int pos){
        String ptmName;
        if(pos==0){
            ptmName = String.valueOf(monoMass);
        }else {
            ptmName = monoMass + " of " + aa;
        }
        return ptmName;
    }

    public String getSkylineFormatPeptide(Peptide peptide){
        ModificationMatch []modificationMatchs = peptide.getVariableModifications();
        if(modificationMatchs !=null && modificationMatchs.length>=1){
            String [] aa = peptide.getSequence().split("");
            String d[] = new String[modificationMatchs.length];
            for(int i=0;i<d.length;i++){
                Modification ptm = ptmFactory.getModification(modificationMatchs[i].getModification());
                int pos = modificationMatchs[i].getSite();
                if(pos==0){
                    aa[0] = "["+String.format("%.4f",ptm.getMass())+"]" + aa[0];
                }else if(pos > aa.length){
                    aa[pos-1] = aa[pos-1] + "["+String.format("%.4f",ptm.getMass())+"]";
                }else{
                    aa[pos-1] = aa[pos-1] + "["+String.format("%.4f",ptm.getMass())+"]";
                }
            }
            return StringUtils.join(aa,"");
        }else {
            return peptide.getSequence();
        }
    }

    private String get_spectrum_header4mzSpecLib(JSpectrum spectrum, Peptide peptide, int spectrum_index){
        int other_attribute_index = 1;
        StringBuilder sb = new StringBuilder();
        sb.append("<Spectrum=").append(spectrum_index).append(">\n");
        sb.append("MS:1003061|library spectrum name=").append(spectrum.title).append("\n");
        sb.append("MS:1003208|experimental precursor monoisotopic m/z=").append(spectrum.getPrecursor().mz).append("\n");
        sb.append("MS:1000894|retention time=").append(String.format("%.4f",spectrum.getPrecursor().getRtInMinutes())).append("\n");
        sb.append("MS:1003059|number of peaks=").append(spectrum.mz.length).append("\n");

        sb.append("<Analyte=1>\n");
        sb.append("MS:1000888|stripped peptide sequence=").append(peptide.getSequence()).append("\n");
        sb.append("[").append(other_attribute_index).append("]MS:1003275|other attribute name=SEQ\n");
        sb.append("[").append(other_attribute_index).append("]MS:1003276|other attribute value=").append(spectrum.SEQ).append("\n");
        other_attribute_index++;
        sb.append("[").append(other_attribute_index).append("]MS:1003275|other attribute name=Modified Peptide\n");
        sb.append("[").append(other_attribute_index).append("]MS:1003276|other attribute value=").append(getSkylineFormatPeptide(peptide)).append("\n");
        sb.append("MS:1000041|charge state=").append(spectrum.getPrecursor().possibleCharges[0]).append("\n");
        return sb.toString();
    }

    public static String asMgf(Spectrum spectrum, String spectrumTitle, int charge, int scan_number){

        double intensity = spectrum.getPrecursor().intensity;
        double mz = spectrum.getPrecursor().mz;
        double[] mzArray = spectrum.mz;
        double[] intensityArray = spectrum.intensity;

        StringBuilder stringBuilder = new StringBuilder();

        stringBuilder.append("BEGIN IONS\n");
        stringBuilder.append("TITLE=").append(spectrumTitle).append("\n");
        stringBuilder.append("PEPMASS=").append(mz).append(" ").append(intensity).append("\n");
        if(spectrum.getPrecursor().rt >= 0.0) {
            double rt = spectrum.getPrecursor().rt;
            stringBuilder.append("RTINSECONDS=").append(rt).append("\n");
        }
        stringBuilder.append("CHARGE=").append(charge).append("+\n");
        if(scan_number >=0){
            stringBuilder.append("SCANS=").append(scan_number).append("\n");
        }
        for(int i=0;i<mzArray.length;i++){
            if(intensityArray[i] > 0.000000000001) {
                stringBuilder.append(mzArray[i]).append(" ").append(intensityArray[i]).append("\n");
            }
        }
        stringBuilder.append("END IONS\n");
        return(stringBuilder.toString());

    }

    public String getPDVModificationString(Peptide peptide){
        String mod = "-";
        ModificationMatch [] modificationMatchs = peptide.getVariableModifications();
        if(modificationMatchs !=null && modificationMatchs.length>=1){
            String d[] = new String[modificationMatchs.length];
            for(int i=0;i<d.length;i++){
                Modification ptm = ptmFactory.getModification(modificationMatchs[i].getModification());
                if(modificationMatchs[i].getSite()==0){
                    d[i] = ptm.getName() + " of N-term@" + modificationMatchs[i].getSite() + "[" + String.format("%.4f", ptm.getMass()) + "]";
                }else {
                    d[i] = ptm.getName() + "@" + modificationMatchs[i].getSite() + "[" + String.format("%.4f", ptm.getMass()) + "]";
                }
            }
            mod = StringUtils.join(d,';');
        }

        return(mod);
    }

}
