package no.uib.probe.optprot.util;

import com.compomics.util.experiment.biology.modifications.Modification;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.biology.proteins.Protein;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.identification.protein_inference.fm_index.FMIndex;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
import com.compomics.util.experiment.identification.spectrum_assumptions.TagAssumption;
import com.compomics.util.experiment.io.biology.protein.Header;
import com.compomics.util.experiment.io.biology.protein.iterators.FastaIterator;
import com.compomics.util.experiment.io.identification.IdfileReader;
import com.compomics.util.experiment.io.identification.IdfileReaderFactory;
import com.compomics.util.experiment.io.mass_spectrometry.MsFileHandler;
import com.compomics.util.experiment.io.mass_spectrometry.mgf.MgfFileWriter;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.io.IoUtil;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.parameters.identification.tool_specific.DirecTagParameters;
import com.compomics.util.parameters.identification.tool_specific.MyriMatchParameters;
import eu.isas.searchgui.SearchHandler;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;
import no.uib.probe.optprot.model.RawScoreModel;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.search.SearchExecuter;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author yfa041
 */
public class SpectraUtilities {

    private final SearchInputSetting searchOptimizerParameters;
    /**
     * The compomics PTM factory.
     */
    private static final ModificationFactory ptmFactory = ModificationFactory.getInstance();

    public SpectraUtilities() {
        searchOptimizerParameters = new SearchInputSetting();
    }

    private File resultsFolder;
    private File oreginalFastaFile;
    private File oreginalMsFile;
    private File oreginalIdentificationFile;
    private IdentificationParameters identificationParameters;
    private int startIndex = 0;
    private final Map<String, Spectrum> spectrumMap = new LinkedHashMap<>();

    ;

    public File[] initInputSubSetFiles(File msFile, File fastaFile, File identificationParametersFile, int fullSpectrumSize) {
        try {

            long start1 = System.currentTimeMillis();
            String spectraFileName = IoUtil.removeExtension(msFile.getName());
            oreginalFastaFile = fastaFile;
            oreginalMsFile = msFile;
            oreginalIdentificationFile = identificationParametersFile;
            identificationParameters = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            SearchParameters searchParameters = identificationParameters.getSearchParameters();
            searchParameters.getModificationParameters().clearVariableModifications();
            searchParameters.getModificationParameters().clearFixedModifications();
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + spectraFileName + Configurations.get_current_file_fingerprent() + "_" + fullSpectrumSize;
            File filteredSubMsFile = null;
            for (File f : oreginalMsFile.getParentFile().listFiles()) {
                if (f.getName().startsWith(Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_") && f.getName().endsWith("_" + oreginalMsFile.getName())) {
                    filteredSubMsFile = f;
                    break;
                }
            }
            if (filteredSubMsFile == null) {
                filteredSubMsFile = new File(msFile.getParent(), Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_" + fullSpectrumSize + "_" + msFile.getName());
            }
//            System.out.println("file name ms " + filteredSubMsFile.getName() + "   " + filteredSubMsFile.exists());
            if (!filteredSubMsFile.exists()) {

                /**
                 * ***
                 */
                MsFileHandler msFileHandler = new MsFileHandler();
                msFileHandler.register(msFile, new OptProtWaitingHandler());
                String fileNameWithoutExtension = IoUtil.removeExtension(msFile.getName());
                String[] spectrumTitles = msFileHandler.getSpectrumTitles(fileNameWithoutExtension);

                DirecTagParameters direcTagParameters = (DirecTagParameters) searchParameters.getIdentificationAlgorithmParameter(Advocate.direcTag.getIndex());
                direcTagParameters.setMaxTagCount(1);
                direcTagParameters.setTagLength(4);
                direcTagParameters.setNumChargeStates(4);
                searchOptimizerParameters.setRunDirecTag(true);
//               

                int step = spectrumTitles.length / (fullSpectrumSize);
                step = step / 4;
                while (spectrumMap.size() < fullSpectrumSize) {
                    int currentSpectrumSize = fullSpectrumSize - spectrumMap.size();
                    if (currentSpectrumSize <= 0) {
                        break;
                    }
                    currentSpectrumSize += (currentSpectrumSize * 0.1);
                    spectrumMap.putAll(substractSpectraWithConfidentTag(startIndex, currentSpectrumSize, msFileHandler, spectrumTitles, fileNameWithoutExtension));
                    startIndex += step;
//                    System.out.println("final size is " + spectrumMap.size() + "   " + currentSpectrumSize);

                }

                writeSpectraToFile(spectrumMap, filteredSubMsFile);

//                System.exit(0);
                /**
                 * **
                 */
//                File unfilteredSubMsFile = subsetSpectraFile(msFile, Configurations.EXTRACT_MS_TYPE.equals("TA"));
//                long end1st = System.currentTimeMillis();
//                double total = (end1st - start1) / 1000.0;
//                System.out.println("process I ( subsetSpectraFile) in seconds: " + total + "  generated file size ");
//                long start2 = System.currentTimeMillis();
//                DirecTagParameters direcTagParameters = (DirecTagParameters) searchParameters.getIdentificationAlgorithmParameter(Advocate.direcTag.getIndex());
//                direcTagParameters.setMaxTagCount(1);
//                direcTagParameters.setTagLength(4);
//                direcTagParameters.setNumChargeStates(4);
//                searchOptimizerParameters.setRunDirecTag(true);
//                resultsFolder = SearchExecuter.executeSearch(updatedName, searchOptimizerParameters, unfilteredSubMsFile, fastaFile, tempIdParam, identificationParametersFile);
//                File direcTagFile = new File(resultsFolder, IoUtil.removeExtension(unfilteredSubMsFile.getName()) + ".tags");
//                if (!direcTagFile.exists()) {
//                    //delete previos sub mgf and cms files 
//                    unfilteredSubMsFile.delete();
//                    File cms = new File(unfilteredSubMsFile.getParent(), unfilteredSubMsFile.getName().replace(".mgf", ".cms"));
//                    cms.delete();
//                    Configurations.EXTRACT_MAX_MS_SIZE += 1000;
//                    MainUtilities.cleanOutputFolder();
//                    return initInputSubSetFiles(msFile, fastaFile, identificationParametersFile);
//                }
//
//                IdfileReader idReader = IdfileReaderFactory.getInstance().getFileReader(direcTagFile);
//                MsFileHandler msFileHandler = new MsFileHandler();
//                msFileHandler.register(unfilteredSubMsFile, MainUtilities.OptProt_Waiting_Handler);
//                ArrayList<SpectrumMatch> matches = idReader.getAllSpectrumMatches(msFileHandler, MainUtilities.OptProt_Waiting_Handler, searchParameters);
//                Map<String, String> specTagMap = new LinkedHashMap<>();
//                for (SpectrumMatch sm : matches) {
//                    TagAssumption tag = sm.getAllTagAssumptions().toList().get(0);
//                    if (tag.getScore() < 0.01) {
//                        specTagMap.put(sm.getSpectrumTitle(), tag.getTag().getContent().get(1).asSequence());
//                    }
//                }
//                System.out.println("tag size " + specTagMap.size());
//                if (specTagMap.size() < Configurations.MIN_TAG_SIZE) {
//                    //delete previos sub mgf and cms files 
//                    unfilteredSubMsFile.delete();
//                    msFileHandler.close();
//                    File cms = new File(unfilteredSubMsFile.getParent(), unfilteredSubMsFile.getName().replace(".mgf", ".cms"));
//                    cms.delete();
////                    Configurations.EXTRACT_MAX_MS_SIZE += 500;
//                    Configurations.EXTRACT_MIN_MS_SIZE += 500;
//                    if (Configurations.EXTRACT_MAX_MS_SIZE <= Configurations.EXTRACT_MIN_MS_SIZE) {
//                        Configurations.EXTRACT_MAX_MS_SIZE = Configurations.EXTRACT_MIN_MS_SIZE + 500;
//                    }
//                    MainUtilities.cleanOutputFolder();
//                    return initInputSubSetFiles(msFile, fastaFile, identificationParametersFile);
//
//                }
//                filteredSubMsFile = subsetSpectraFile(unfilteredSubMsFile, specTagMap.keySet());
//                long end2ndst = System.currentTimeMillis();
//                total = (end2ndst - start2) / 1000.0;
//                System.out.println("process II ( DirecTag+Generated subMGF) in seconds: " + total);
//                File cms = new File(unfilteredSubMsFile.getParent(), unfilteredSubMsFile.getName().replace(".mgf", ".cms"));
//                cms.delete();
//                unfilteredSubMsFile.delete();
//                MainUtilities.deleteFolder(direcTagFile.getParentFile());
            }

            File subFasta = null;

            for (File f : oreginalFastaFile.getParentFile().listFiles()) {
                if (f.getName().startsWith(Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_") && f.getName().endsWith("_" + oreginalFastaFile.getName())) {
                    subFasta = f;
                    break;
                }
            }
            if (subFasta == null) {
                subFasta = new File(fastaFile.getParent(), Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_" + fullSpectrumSize + "_" + fastaFile.getName());
            }

            if (!subFasta.exists()) {
                long start3 = System.currentTimeMillis();
                searchOptimizerParameters.setRunNovor(true);
                resultsFolder = SearchExecuter.executeSearch(updatedName, searchOptimizerParameters, filteredSubMsFile, fastaFile, identificationParameters, identificationParametersFile);
                File NovorFile = new File(resultsFolder, IoUtil.removeExtension(filteredSubMsFile.getName()) + ".novor.csv");
                Set<String> sequences = SpectraUtilities.getSequences(NovorFile);
                System.out.println("sequence from nover " + sequences.size());
                long end3rd = System.currentTimeMillis();
                double total = (end3rd - start3) / 1000.0;
                System.out.println("process III ( Novor ) in seconds: " + total);
                long start4 = System.currentTimeMillis();
                subFasta = initSubFastaFile(fastaFile, sequences);
                long end4th = System.currentTimeMillis();
                total = (end4th - start4) / 1000.0;
                System.out.println("process IV ( Generated Fasta) in seconds: " + total);
                long end = System.currentTimeMillis();
                total = (end - start1) / 1000.0;
                System.out.println("Total Elapsed Time for initInputSubSetFiles in seconds: " + total);
                MainUtilities.deleteFolder(NovorFile.getParentFile());
            }

            return new File[]{filteredSubMsFile, subFasta};
        } catch (IOException ex) {  //ex) {//
            ex.printStackTrace();
            MainUtilities.cleanOutputFolder();
//            System.exit(0);
        }
        return null;

    }

    public static Map<String, Integer> initSpectraIndexMap(File msFile) {
        Map<String, Integer> spectraIndexMap = new HashMap<>();
        try {

            String fileNameWithoutExtension = IoUtil.removeExtension(msFile.getName());
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(msFile, new OptProtWaitingHandler());
            String[] spectrumTitles = msFileHandler.getSpectrumTitles(fileNameWithoutExtension);
            for (int i = 0; i < spectrumTitles.length; i++) {
                spectraIndexMap.put(spectrumTitles[i], i);
            }
        } catch (IOException ex) {
            Logger.getLogger(SpectraUtilities.class.getName()).log(Level.SEVERE, null, ex);
        }
        return spectraIndexMap;

    }

    private File initSubFastaFile(File fastaFile, Set<String> sequences) {

        File subFastaFile = new File(fastaFile.getParent(), Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_" + fastaFile.getName());
        if (subFastaFile.exists()) {
            subFastaFile.delete();
        }
        SpectraUtilities.createSubFastaFile(fastaFile, subFastaFile, sequences);
//        SpectraUtilities.createSubFastaFile(fastaFile, subFastaFile,true);
        return subFastaFile;

    }

    private File subsetSpectraFile(File oreginalMsFile, Set<String> spectraTags) {

        try {
            String name = oreginalMsFile.getName().replace(Configurations.DEFAULT_RESULT_NAME + "_sub" + Configurations.get_current_file_fingerprent() + "_", "");
            File subMsFile = new File(oreginalMsFile.getParent(), Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_" + name);
            if (subMsFile.exists()) {
                subMsFile.delete();
                File subSampleCMS = new File(oreginalMsFile.getParent(), subMsFile.getName().replace(".mgf", ".cms"));
                subSampleCMS.delete();
            } else {
                subMsFile.createNewFile();
            }
            SpectraUtilities.writeFilteredMgfFile(oreginalMsFile, subMsFile, spectraTags);
            return subMsFile;
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return oreginalMsFile;
    }

//    private File subsetSpectraFile(File oreginalMsFile, boolean targtedArea) {
//        try {
//            File subMsFile = new File(oreginalMsFile.getParent(), Configurations.DEFAULT_RESULT_NAME + "_sub" + Configurations.get_current_file_fingerprent() + "_" + oreginalMsFile.getName());
//            if (subMsFile.exists()) {
//                subMsFile.delete();
//                File subSampleCMS = new File(oreginalMsFile.getParent(), subMsFile.getName().replace(".mgf", ".cms"));
//                subSampleCMS.delete();
//            } else {
//                subMsFile.createNewFile();
//            }
//
//            if (targtedArea) {
//                SpectraUtilities.writeSubSetTargtedFileAreaMgfFile(oreginalMsFile, subMsFile, Configurations.EXTRACT_MIN_MS_SIZE, Configurations.EXTRACT_MAX_MS_SIZE);
//            } else {
//                SpectraUtilities.writeSubSetEveryNMgfFile(oreginalMsFile, subMsFile, Configurations.EXTRACT_MIN_MS_SIZE, Configurations.EXTRACT_MAX_MS_SIZE);
//            }
//            return subMsFile;
//        } catch (IOException ex) {
//            ex.printStackTrace();
//        }
//        return oreginalMsFile;
//    }
    /**
     * Writes the spectra of a file in the Mascot Generic File (mgf) format.
     *
     * @param mgfFile The spectrum file to use to get the spectra.
     * @param destinationFile The file where to write.
     * @param maxSpectNum number of spectra in the final file
     */
    public static void writeSubSetEveryNMgfFile(
            File mgfFile,
            File destinationFile,
            double minSpectNum,
            double maxSpectNum
    ) {
        try {
            String fileNameWithoutExtension = IoUtil.removeExtension(mgfFile.getName());
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(mgfFile, new OptProtWaitingHandler());
            String[] spectrumTitles = msFileHandler.getSpectrumTitles(fileNameWithoutExtension);

            int initiSize = Math.min((int) maxSpectNum, (int) Math.round(spectrumTitles.length * 20 / 100));
            initiSize = Math.max(initiSize, (int) minSpectNum);
//            System.out.println("initial zer " + initiSize + "  vs " + maxSpectNum);

            int n = (int) Math.round(spectrumTitles.length / initiSize);
            if (spectrumTitles == null) {
                throw new IllegalArgumentException(
                        fileNameWithoutExtension + " not loaded."
                );
            }
            MgfFileWriter writer = new MgfFileWriter(destinationFile);
            for (int i = 0; i < spectrumTitles.length;) {
                String spectrumTitle = spectrumTitles[i];
                Spectrum spectrum = msFileHandler.getSpectrum(fileNameWithoutExtension, spectrumTitle);
                writer.writeSpectrum(spectrumTitle, spectrum);
                i += n;
            }

            writer.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

    /**
     * Writes the spectra of a file in the Mascot Generic File (mgf) format.
     *
     * @param mgfFile The spectrum file to use to get the spectra.
     * @param destinationFile The file where to write.
     * @param ids
     */
    public static void writeFilteredMgfFile(
            File mgfFile,
            File destinationFile,
            Set<String> ids
    ) {
        try {
            String fileNameWithoutExtension = IoUtil.removeExtension(mgfFile.getName());
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(mgfFile, new OptProtWaitingHandler());
            int counter = 0;
            String[] spectrumTitles = msFileHandler.getSpectrumTitles(fileNameWithoutExtension);

            if (spectrumTitles == null) {
                throw new IllegalArgumentException(
                        fileNameWithoutExtension + " not loaded."
                );
            }
            MgfFileWriter writer = new MgfFileWriter(destinationFile);
            for (int i = 0; i < (spectrumTitles.length); i++) {
                String spectrumTitle = spectrumTitles[i];
                if (ids.contains(spectrumTitle)) {
                    Spectrum spectrum = msFileHandler.getSpectrum(fileNameWithoutExtension, spectrumTitle);
                    writer.writeSpectrum(spectrumTitle, spectrum);
                    counter++;
                }
            }
//            System.out.println("total filtered subset mgf " + counter);
            writer.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

    /**
     * Writes the spectra of a file in the Mascot Generic File (mgf) format.
     *
     * @param mgfFile The spectrum file to use to get the spectra.
     * @param destinationFile The file where to write.
     * @param minSpecNum
     *
     * @param maxSpectNum number of spectra in the final file
     */
    public static void writeSubSetTargtedFileAreaMgfFile(
            File mgfFile,
            File destinationFile,
            double minSpecNum,
            double maxSpectNum
    ) {
        try {
            String fileNameWithoutExtension = IoUtil.removeExtension(mgfFile.getName());
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(mgfFile, new OptProtWaitingHandler());

            String[] spectrumTitles = msFileHandler.getSpectrumTitles(fileNameWithoutExtension);
            int targetedAreaStep = spectrumTitles.length / 4;

            int initiSize = Math.min((int) maxSpectNum, (int) Math.round(targetedAreaStep * 20 / 100));
            initiSize = Math.max(initiSize, (int) minSpecNum);
//            System.out.println("initial zer " + initiSize + "  vs " + maxSpectNum);

            int n = (int) Math.round(targetedAreaStep * 2 / initiSize);

            if (spectrumTitles == null) {
                throw new IllegalArgumentException(
                        fileNameWithoutExtension + " not loaded."
                );
            }
            MgfFileWriter writer = new MgfFileWriter(destinationFile);
            for (int i = targetedAreaStep; i < (spectrumTitles.length - targetedAreaStep);) {
                String spectrumTitle = spectrumTitles[i];
                Spectrum spectrum = msFileHandler.getSpectrum(fileNameWithoutExtension, spectrumTitle);
                writer.writeSpectrum(spectrumTitle, spectrum);
                i += n;
            }

            writer.close();
//            System.out.println("at n: " + n + "  " + targetedAreaStep * 2 + "  " + spectrumTitles.length);
        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

    /**
     * Writes the spectra of a file in the Mascot Generic File (mgf) format.
     *
     * @param spectrumResource
     * @param destinationFile The file where to write.
     */
    public static void writeSpectraToFile(
            Map<String, Spectrum> spectrumResource,
            File destinationFile) {
        MgfFileWriter writer = new MgfFileWriter(destinationFile);
        for (String spectrumTitle : spectrumResource.keySet()) {
            Spectrum spectrum = spectrumResource.get(spectrumTitle);
            writer.writeSpectrum(spectrumTitle, spectrum);
        }
        writer.close();

    }

    /**
     * Appends decoy sequences to the provided FASTA file.
     *
     * @param fastaIn the FASTA file to read
     * @param fastaOut the FASTA file to write displaying progress
     * @param micro
     */
    public static void createSubFastaFile(
            File fastaIn,
            File fastaOut,
            boolean micro
    ) {

        try {
            int total = 0;
            int left = 0;
            int specIndex = 0;
            FastaIterator fastaIterator = new FastaIterator(fastaIn);

            try (BufferedWriter bw = new BufferedWriter(new FileWriter(fastaOut))) {
                Protein protein;
                while ((protein = fastaIterator.getNextProtein()) != null) {
                    Header header = fastaIterator.getLastHeader();
                    total++;
                    if (header.getProteinEvidence() != 1 || protein.getLength() <= 400) { //&& !header.asGenericHeader().contains("SV=2")
                        continue;
                    }
                    specIndex++;
                    if (micro) {
                        if (specIndex == 3) {
                            specIndex = 0;
                            continue;
                        }
                    }
                    left++;
                    String sequence = protein.getSequence();
                    String rawHeader = header.getRawHeader();
                    bw.write(rawHeader);
                    bw.newLine();
                    bw.write(sequence);

                    bw.newLine();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        }
    }

    public static void createSubFastaFile(
            File fastaIn,
            File fastaOut,
            Set<String> sequences
    ) {

        try {
            int total = 0;
            int included = 0;
            int specIndex = 0;
            FastaIterator fastaIterator = new FastaIterator(fastaIn);

            try (BufferedWriter bw = new BufferedWriter(new FileWriter(fastaOut))) {
                Protein protein;
                while ((protein = fastaIterator.getNextProtein()) != null) {
                    Header header = fastaIterator.getLastHeader();
                    total++;
                    if (header.getProteinEvidence() != 1 || protein.getAccession().endsWith("_REVERSED")) { //&& !header.asGenericHeader().contains("SV=2")|| protein.getLength() <= 400
                        continue;
                    }
                    boolean escape = true;
                    for (String seq : sequences) {
                        if (protein.getSequence().replace("I", "L").contains(seq.replace("I", "L"))) {
                            escape = false;
                            break;
                        }
                    }
                    if (escape) {
                        specIndex++;
                        continue;
                    }

                    included++;
                    String sequence = protein.getSequence();
                    String rawHeader = header.getRawHeader();
                    bw.write(rawHeader);
                    bw.newLine();
                    bw.write(sequence);
                    bw.newLine();

                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
            if (included < 50) {
                copySubFastaFile(fastaIn, fastaOut);
            }

        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        }

    }

    public static void copySubFastaFile(
            File fastaIn,
            File fastaOut
    ) {

        try {
            fastaOut.delete();
            fastaOut.createNewFile();
            int total = 0;
            int left = 0;
            FastaIterator fastaIterator = new FastaIterator(fastaIn);

            try (BufferedWriter bw = new BufferedWriter(new FileWriter(fastaOut))) {
                Protein protein;
                while ((protein = fastaIterator.getNextProtein()) != null) {
                    Header header = fastaIterator.getLastHeader();
                    total++;
                    if (header.getProteinEvidence() != 1 || protein.getAccession().endsWith("_REVERSED")) { //&& !header.asGenericHeader().contains("SV=2")|| protein.getLength() <= 400
                        continue;
                    }

                    left++;
                    String sequence = protein.getSequence();
                    String rawHeader = header.getRawHeader();
                    bw.write(rawHeader);
                    bw.newLine();
                    bw.write(sequence);
                    bw.newLine();

                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            Logger.getLogger(SpectraUtilities.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public static Set<String> getSequences(File novorOutputFile) {
        String line;
        String splitBy = ",";
        Set<String> sequences = new HashSet<>();
        try {
//parsing a CSV file into BufferedReader class constructor  
            BufferedReader br = new BufferedReader(new FileReader(novorOutputFile));
            while ((line = br.readLine()) != null) //returns a Boolean value 
            {
                if (line.startsWith("# id,")) {
                    break;
                }
            }

            while ((line = br.readLine()) != null) //returns a Boolean value 
            {
                String[] values = line.split(splitBy);    // use comma as separator  
                sequences.add(values[9].replaceAll("[()0123456]", "").trim());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return sequences;
    }

    public static Set<String> getTaggedSequences(File novorOutputFile) {
        String line;
        String splitBy = ",";
        Set<String> sequences = new HashSet<>();
        try {
//parsing a CSV file into BufferedReader class constructor  
            BufferedReader br = new BufferedReader(new FileReader(novorOutputFile));
            while ((line = br.readLine()) != null) //returns a Boolean value 
            {
                if (line.startsWith("# id,")) {
                    break;
                }
            }
            while ((line = br.readLine()) != null) //returns a Boolean value 
            {
                String[] values = line.split(splitBy);    // use comma as separator  
                sequences.add(values[9].trim());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return sequences;
    }

    /**
     * Lets the user select an output folder and starts the recalibration of
     * spectra.
     *
     * @param sourceMgfFile
     * @param fastaFile
     * @param spectrumMatches
     * @param identificationParameters
     * @return
     */
    public static File recalibrateSpectra(File sourceMgfFile, File fastaFile, ArrayList<SpectrumMatch> spectrumMatches, IdentificationParameters identificationParameters) {

        try {
            final File selectedFolder = new File(sourceMgfFile.getParent());
            MsFileHandler spectrumProvider = new MsFileHandler();
            spectrumProvider.register(sourceMgfFile, new OptProtWaitingHandler());
            final boolean precursors = true;
            final boolean fragments = true;
            FMIndex sequenceProvider = new FMIndex(fastaFile, null, new OptProtWaitingHandler(), false, identificationParameters);
            try {
                return SpectralRecalibrationUtilitiy.writeRecalibratedSpectra(
                        precursors,
                        fragments,
                        selectedFolder,
                        spectrumMatches,
                        sequenceProvider,
                        spectrumProvider,
                        identificationParameters
                ).get(0);

            } catch (IOException e) {
                e.printStackTrace();
            }

        } catch (IOException ex) {
            Logger.getLogger(SpectraUtilities.class.getName()).log(Level.SEVERE, null, ex);
        }
        return null;
    }

    public Map<String, Spectrum> substractSpectraWithConfidentTag(int startIndex, int maxSpectraNumber, MsFileHandler msFileHandler, String[] spectrumTitles, String msFileNameWithoutExtension) {
        Map<String, Spectrum> spectraMap = new LinkedHashMap<>();
        try {

            int targetedAreaStep = spectrumTitles.length / 4;
            int n = targetedAreaStep * 2 / (maxSpectraNumber);
            File destinationFile = new File(oreginalMsFile.getParentFile(), Configurations.DEFAULT_RESULT_NAME + "_" + msFileNameWithoutExtension + ".mgf");
            destinationFile.createNewFile();
            try (MgfFileWriter writer = new MgfFileWriter(destinationFile)) {
                for (int i = startIndex + targetedAreaStep; i < (spectrumTitles.length - targetedAreaStep);) {
                    Spectrum spectrum = msFileHandler.getSpectrum(msFileNameWithoutExtension, spectrumTitles[i]);
                    writer.writeSpectrum(spectrumTitles[i], spectrum);
                    spectraMap.put(spectrumTitles[i], spectrum);
                    i += n;
                }
            }
            String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + msFileNameWithoutExtension + Configurations.get_current_file_fingerprent();

            File tempResultsFolder = SearchExecuter.executeSearch(updatedName, searchOptimizerParameters, destinationFile, oreginalFastaFile, identificationParameters, oreginalIdentificationFile);

            File direcTagFile = new File(tempResultsFolder, IoUtil.removeExtension(destinationFile.getName()) + ".tags");

            if (!direcTagFile.exists()) {
                //delete previos sub mgf and cms files
                destinationFile.delete();
                File cms = new File(destinationFile.getParent(), destinationFile.getName().replace(".mgf", ".cms"));
                cms.delete();
                substractSpectraWithConfidentTag(startIndex, maxSpectraNumber + 500, msFileHandler, spectrumTitles, msFileNameWithoutExtension);
                return spectraMap;
            }

            IdfileReader idReader = IdfileReaderFactory.getInstance().getFileReader(direcTagFile);
            MsFileHandler subMsFileHandler = new MsFileHandler();
            subMsFileHandler.register(destinationFile, MainUtilities.OptProt_Waiting_Handler);
            ArrayList<SpectrumMatch> matches = idReader.getAllSpectrumMatches(subMsFileHandler, MainUtilities.OptProt_Waiting_Handler, identificationParameters.getSearchParameters());
            for (SpectrumMatch sm : matches) {
                TagAssumption tag = sm.getAllTagAssumptions().toList().get(0);
                if (tag.getScore() > 0.01) {
                    spectraMap.remove(sm.getSpectrumTitle());
                }
            }
            destinationFile.delete();
            File cms = new File(destinationFile.getParent(), destinationFile.getName().replace(".mgf", ".cms"));
            cms.delete();

        } catch (IOException | SQLException | ClassNotFoundException | InterruptedException | JAXBException | XmlPullParserException | XMLStreamException ex) {
            ex.printStackTrace();
        }
        return spectraMap;
    }

    public static double[] getValidatedIdentificationReferenceData(List<SpectrumMatch> matches, Advocate searchEngine) {
        double[] data = new double[matches.size()];
        int index = 0;
        if (searchEngine.getIndex() == Advocate.myriMatch.getIndex()) {

            for (SpectrumMatch sm : matches) {
                for (PeptideAssumption peptideAssumtion : sm.getAllPeptideAssumptions().toList()) {
                    if (peptideAssumtion.getRank() != 1) {
                        continue;
                    }
                    sm.setBestPeptideAssumption(peptideAssumtion);
                    double score = -Math.log(peptideAssumtion.getScore());
                    data[index++] = score;
                    break;
                }
            }
        } else if (searchEngine.getIndex() == Advocate.sage.getIndex()) {

            for (SpectrumMatch sm : matches) {
                for (PeptideAssumption peptideAssumtion : sm.getAllPeptideAssumptions().toList()) {
                    if (peptideAssumtion.getRank() != 1) {
                        continue;
                    }
                    double score = peptideAssumtion.getRawScore();
                    sm.setBestPeptideAssumption(peptideAssumtion);
                    data[index++] = score;
                    break;
                }
            }
        } else if (searchEngine.getIndex() == Advocate.xtandem.getIndex()) {

            for (SpectrumMatch sm : matches) {

                for (PeptideAssumption peptideAssumtion : sm.getAllPeptideAssumptions().toList()) {
                    if (peptideAssumtion.getRank() != 1) {
                        continue;
                    }
                    double score = peptideAssumtion.getRawScore();
                    data[index++] = score;
                    sm.setBestPeptideAssumption(peptideAssumtion);
                    break;
                }
            }
        }
        return data;
    }

    public static RawScoreModel getComparableRawScore(SearchingSubDataset optProtDataset, List<SpectrumMatch> matches, Advocate searchEngine, boolean addData, String scoreId) {

        List<Double> sharedToData = new ArrayList<>();
        List<Double> sharedFromData = new ArrayList<>();
        List<Double> onlyFromData = new ArrayList<>();
        List<Double> onlyToData = new ArrayList<>();
        Map<String, Double> matchScores = new HashMap<>();

        RawScoreModel rawScore = new RawScoreModel(scoreId);
        rawScore.setFullSpectraSize(optProtDataset.getTotalSpectraNumber());
        if (searchEngine.getIndex() == Advocate.xtandem.getIndex()) {
            for (SpectrumMatch sm : matches) {
                for (PeptideAssumption peptideAssumtion : sm.getAllPeptideAssumptions().toList()) {
                    if (peptideAssumtion.getRank() != 1) {
                        continue;
                    }
                    sm.setBestPeptideAssumption(peptideAssumtion);
                    matchScores.put(sm.getSpectrumTitle(), peptideAssumtion.getRawScore());
                    break;
                }
            }

        } else if (searchEngine.getIndex() == Advocate.myriMatch.getIndex()) {
            for (SpectrumMatch sm : matches) {
                for (PeptideAssumption peptideAssumtion : sm.getAllPeptideAssumptions().toList()) {
                    if (peptideAssumtion.getRank() != 1) {
                        continue;
                    }
                    sm.setBestPeptideAssumption(peptideAssumtion);
                    matchScores.put(sm.getSpectrumTitle(), peptideAssumtion.getRawScore());
                    break;
                }
            }
        } else if (searchEngine.getIndex() == Advocate.sage.getIndex()) {
            for (SpectrumMatch sm : matches) {
                for (PeptideAssumption peptideAssumtion : sm.getAllPeptideAssumptions().toList()) {
                    if (peptideAssumtion.getRank() != 1) {
                        continue;
                    }
                    sm.setBestPeptideAssumption(peptideAssumtion);
                    if (peptideAssumtion.getRawScore() > 0) {
                        matchScores.put(sm.getSpectrumTitle(), peptideAssumtion.getRawScore());
                    }
                    break;
                }
            }
        }

        for (String titile : optProtDataset.getFullSpectraScore().keySet()) {
            double fromScore = optProtDataset.getFullSpectraScore().get(titile);
            if (matchScores.containsKey(titile)) {
                double toScore = matchScores.get(titile);
                if (fromScore > 0) {
                    sharedToData.add(toScore);
                    sharedFromData.add(fromScore);
                } else {
                    onlyToData.add(toScore);
                }

            } else if (fromScore > 0) {
                onlyFromData.add(fromScore);
            }
        }
        if (sharedFromData.isEmpty() && matches.size() < ((double) onlyFromData.size() * 0.2)) {
            rawScore.setTotalNumber(matches.size());
            rawScore.setFinalScore(-10);
            return rawScore;
        }

        Collections.sort(sharedFromData);
        //devide the TODATA data into 3 category         

        double[] fromSharedData = sharedFromData.stream().mapToDouble(Double::doubleValue).toArray();
        double[] toSharedData = sharedToData.stream().mapToDouble(Double::doubleValue).toArray();
        double[] fromOnlyData = onlyFromData.stream().mapToDouble(Double::doubleValue).toArray();
        double[] toOnlyData = onlyToData.stream().mapToDouble(Double::doubleValue).toArray();

        double[] improved1 = SpectraUtilities.compareData(fromSharedData, toSharedData);
        double ratio = Math.max(onlyToData.size(), onlyFromData.size()) / (double) optProtDataset.getTotalSpectraNumber();
        double score2;
        score2 = SpectraUtilities.compareData(fromOnlyData, toOnlyData)[4];
        double fs1 = improved1[4];
        double fs2 = score2 * ratio;
        double fs = fs1 + fs2;
        rawScore.setS1(fs1);
        rawScore.setS2(fs2);
        rawScore.setTotalNumber(matches.size());
        rawScore.setFinalScore(fs);
        rawScore.setSameData(rawScore.getFinalScore() == 0.0 && matches.size() == ((onlyFromData.size() + toSharedData.length)));
        boolean accepted = rawScore.getFinalScore() >= optProtDataset.getComparisonsThreshold();
        rawScore.setSignificatChange(accepted);
        boolean senstive = false;
        if (!rawScore.isSameData()) {
            if (rawScore.getFinalScore() > 0) {
                senstive = true;
            } else if (optProtDataset.getCurrentScoreModel() != null && (rawScore.getFinalScore() >= -0.05 && matches.size() > optProtDataset.getCurrentScoreModel().getSpectrumMatchResult().size())) {
                senstive = true;
            } else if ((rawScore.getFinalScore() >= -0.05 && matches.size() > optProtDataset.getActiveIdentificationNum())) {
                senstive = true;
            } else if ((fs1 > 0 && matches.size() > optProtDataset.getActiveIdentificationNum())) {
                senstive = true;
            } else if ((fs1 - fs2 > 0)) {
                senstive = true;
            }
        }
        rawScore.setSensitiveChange(senstive);
        if (rawScore.isSignificatChange() || rawScore.isSensitiveChange() || rawScore.isSameData() || addData) {
            rawScore.setSpectrumMatchResult(matches);
        }
        return rawScore;

    }

    public static List<Double> getQuartileRatio(String[] titiles, List<SpectrumMatch> matches) {

        double partCount = 4;
        Map<String, Double> fullSpectraMap = new LinkedHashMap<>();
        for (String titel : titiles) {
            fullSpectraMap.put(titel, 100.0);
        }
        for (SpectrumMatch sm : matches) {
            TagAssumption tag = sm.getAllTagAssumptions().toList().get(0);
            fullSpectraMap.replace(sm.getSpectrumTitle(), tag.getScore());
        }
        double countId = 0;
        double countUnId = 0;
        double partI = Math.round((double) titiles.length / 4.0);
        double partII = 2 * partI;
        double partIII = 3 * partI;
        double[][] quartileData = new double[(int) partCount][2];
        int counter = 0;
        int quartileIndex = 0;
        int left = titiles.length;
        for (String titile : titiles) {
            double v = fullSpectraMap.get(titile);
            if (v <= 0.01) {
                countId++;
            } else {
                countUnId++;
            }
            counter++;
            left--;
            if ((left == partI) || (left == 0) || (left == partII) || (left == partIII)) {
                //initpart I 
                quartileData[quartileIndex++] = new double[]{countId, countUnId};
                counter = 0;
                countId = 0;
                countUnId = 0;

            }
        }

        List<Double> quartileRatios = new ArrayList<>();
        double total = 0;
        for (double[] qData : quartileData) {
            double qRatio = qData[0] / (qData[0] + qData[1]);
            total += qData[0];
            quartileRatios.add(qRatio);
        }
        total = total / (double) titiles.length;
        quartileRatios.add(total);

        return quartileRatios;
    }

    public static List<Double> getQuartileRatioT(String[] titiles, List<SpectrumMatch> matches) {

        double partCount = 4;
        Map<String, Double> fullSpectraMap = new LinkedHashMap<>();
        for (String titel : titiles) {
            fullSpectraMap.put(titel, 100.0);
        }
        for (SpectrumMatch sm : matches) {
            TagAssumption tag = sm.getAllTagAssumptions().toList().get(0);
            fullSpectraMap.replace(sm.getSpectrumTitle(), tag.getScore());
        }
        double countId = 0;
        double countUnId = 0;
        double partI = Math.round((double) titiles.length / 4.0);
        double partII = 2 * partI;
        double partIII = 3 * partI;
        double[][] quartileData = new double[(int) partCount][2];
        int counter = 0;
        int quartileIndex = 0;
        int left = titiles.length;
        for (String titile : titiles) {
            double v = fullSpectraMap.get(titile);
            if (v <= 0.1) {
                countId++;
            } else {
                countUnId++;
            }
            counter++;
            left--;
            if ((left == partI) || (left == 0) || (left == partII) || (left == partIII)) {
                //initpart I 
                quartileData[quartileIndex++] = new double[]{countId, countUnId};
                counter = 0;
                countId = 0;
                countUnId = 0;

            }
        }

        double totalCount = 0;
        List<Double> quartileRatios = new ArrayList<>();
        for (double[] qData : quartileData) {
            totalCount += qData[0];
            quartileRatios.add(qData[0]);
        }
        double quat_3_sum = 0;
        for (int j = 0; j < quartileRatios.size() - 1; j++) {
            double ratio = quartileRatios.get(j) / totalCount;
            ratio = Math.round(ratio * 100.0) / 100.0;
            quartileRatios.set(j, ratio);
            quat_3_sum += ratio;
        }
        double ratio = 1.0 - quat_3_sum;
        ratio = Math.round(ratio * 100.0) / 100.0;
        quartileRatios.set(quartileRatios.size() - 1, ratio);

        double totalRatio = totalCount / titiles.length;
        System.out.println("total ratio is " + totalRatio);
        quartileRatios.add(totalRatio);

        return quartileRatios;
    }

    public static List<SpectrumMatch> getValidatedIdentificationResults(File resultOutput, File msFile, Advocate searchEngine, IdentificationParameters identificationParameters) {
        List<SpectrumMatch> validatedMaches = Collections.synchronizedList(new ArrayList<>());
        if (resultOutput == null) {
            System.out.println("output result file was null");
            return validatedMaches;
        }
        try {
            File idResultsFile = null;
            if (searchEngine.getIndex() == Advocate.myriMatch.getIndex()) {
                MyriMatchParameters myriMatchParameters = (MyriMatchParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.myriMatch.getIndex());
//                myriMatchParameters.setTicCutoffPercentage();
                idResultsFile = new File(resultOutput, SearchHandler.getMyriMatchFileName(IoUtil.removeExtension(msFile.getName()), myriMatchParameters));
            } else if (searchEngine.getIndex() == Advocate.xtandem.getIndex()) {
                idResultsFile = new File(resultOutput, SearchHandler.getXTandemFileName(IoUtil.removeExtension(msFile.getName())));
            } else if (searchEngine.getIndex() == Advocate.sage.getIndex()) {
                idResultsFile = new File(resultOutput, (IoUtil.removeExtension(msFile.getName()) + ".sage.tsv"));
            }
            if (idResultsFile == null || !idResultsFile.exists()) {
                System.out.println("id result file was null " + resultOutput.getAbsolutePath() + "   ---  " + searchEngine.getName() + "  " + idResultsFile);
                System.exit(0);
                return validatedMaches;
            }
            IdfileReader idReader = IdfileReaderFactory.getInstance().getFileReader(idResultsFile);
            MsFileHandler subMsFileHandler = new MsFileHandler();
            subMsFileHandler.register(msFile, MainUtilities.OptProt_Waiting_Handler);
            ArrayList<SpectrumMatch> matches = idReader.getAllSpectrumMatches(subMsFileHandler, MainUtilities.OptProt_Waiting_Handler, identificationParameters.getSearchParameters());
            for (SpectrumMatch sm : matches) {
                for (PeptideAssumption peptideAssumtion : sm.getAllPeptideAssumptions().toList()) {
                    if (peptideAssumtion.getRank() != 1) {
                        continue;
                    }
                    sm.setBestPeptideAssumption(peptideAssumtion);
                    if (peptideAssumtion.getRawScore() > 0) {
                        validatedMaches.add(sm);
                    }
                    break;
                }
            }

//            validatedMaches.addAll(matches);
        } catch (IOException | SQLException | ClassNotFoundException | InterruptedException | JAXBException | XmlPullParserException | XMLStreamException ex) {
            System.out.println("no.uib.probe.optprot.util.SpectraFileUtilities.getValidatedIdentificationResults() " + ex.getLocalizedMessage());
            ex.printStackTrace();
        }
//         MainUtilities.cleanOutputFolder();
        return validatedMaches;
    }

    public static Set<String> getModifiedSpectrumSet(List<SpectrumMatch> spectraResults) {
//        System.out.println("mod id is ---------- "+mod);
        Set<String> modSpectrumSet = new HashSet<>();
        for (SpectrumMatch sm : spectraResults) {
            if (sm == null) {
                continue;
            }
            PeptideAssumption peptideAssm = sm.getBestPeptideAssumption();
            if (peptideAssm.getPeptide().getVariableModifications().length > 0) {
//                for (ModificationMatch mm : peptideAssm.getPeptide().getVariableModifications()) {
//                    if (mm.getModification().endsWith(mod.getPattern().toString()) && Precision.equals(mod.getMass(), Double.parseDouble(mm.getModification().split("@")[0]), 0.01)) {
                modSpectrumSet.add(sm.getSpectrumTitle());

//                    }
//                }
            }
        }

        return modSpectrumSet;

    }

    public static Set<String> getIntersectionSet(Set<String> setI, Set<String> setII) {
        Set<String> intersection = new HashSet<>();
        Set<String> full = new HashSet<>(setI);
        full.addAll(setII);
        for (String title : full) {
            if (setI.contains(title) && setII.contains(title)) {
                intersection.add(title);
            }
        }

        return intersection;

    }

    public static double[] isBetterScore(double[] fromData, double[] toData, boolean pairData) {
        ScoreComparison sc = new ScoreComparison();
        if (pairData) {
            if (toData.length > fromData.length) {
                double[] updatedData = new double[toData.length];
                System.arraycopy(fromData, 0, updatedData, 0, fromData.length);
                fromData = updatedData;
            } else if (toData.length < fromData.length) {
                double[] updatedData = new double[fromData.length];
                System.arraycopy(toData, 0, updatedData, 0, toData.length);
                toData = updatedData;
            }
        }
        double[] finalScores = sc.calculateFinalScore(fromData, toData);
        return finalScores;

    }

    public static String compareScoresSet(Map<String, RawScoreModel> resultsMap, int totalSpecNumber) {
//        TreeMap<RawScoreModel, String> sorter = new TreeMap<>(Collections.reverseOrder());
        List<RawScoreModel> scoreModelSorter = new ArrayList<>();

        for (String rs1Key : resultsMap.keySet()) {
//            sorter.put(resultsMap.get(rs1Key), rs1Key);
            scoreModelSorter.add(resultsMap.get(rs1Key));
        }
        Collections.sort(scoreModelSorter);
        Collections.reverse(scoreModelSorter);
//        return sorter.lastEntry().getValue();
        int index1 = -1;
        List<String> topScoreSet = new ArrayList<>(resultsMap.keySet());
        for (RawScoreModel score1 : scoreModelSorter) {
            String rs1Key = score1.getComparisonId();
            index1++;
            if (!topScoreSet.contains(rs1Key)) {
                continue;
            }
            int index2 = -1;
            for (RawScoreModel score2 : scoreModelSorter) {
                String rs2Key = score2.getComparisonId();
                index2++;
                if (index2 <= index1 || !topScoreSet.contains(rs2Key)) {
                    continue;
                }
                double key1Better = isBetterScore(resultsMap.get(rs2Key).getSpectrumMatchResult(), resultsMap.get(rs1Key).getSpectrumMatchResult(), totalSpecNumber);
//                System.out.println("resultsMap.get(rs1Key) " + resultsMap.get(rs1Key).getFinalScore() + "   " + resultsMap.get(rs2Key).getFinalScore());
//                System.out.println(rs1Key + " better " + rs2Key + "  " + key1Better);
                if ((key1Better > 0) || Double.isNaN(key1Better)) {// && (resultsMap.get(rs1Key).getFinalScore() > resultsMap.get(rs2Key).getFinalScore())
                    topScoreSet.remove(rs2Key);
                } else {
                    System.out.println("-----------------------------------remove " + rs1Key + "-------------------------------");
                    topScoreSet.remove(rs1Key);
                    break;
                }
            }
        }
        if (topScoreSet.size() > 1) {
            System.out.println("error should be only oine value " + topScoreSet);
        } else if (topScoreSet.isEmpty()) {
            System.out.println("error should be only oine value " + topScoreSet);
            return "";
        }
        String topSelection = topScoreSet.get(0);
        if (!scoreModelSorter.get(0).getComparisonId().equalsIgnoreCase(topSelection)) {
            System.out.println("-------------------------------->> disagree with both methods " + scoreModelSorter.get(0).getComparisonId() + "  " + topSelection);
        }
        return topSelection;
    }

    public static String getTopScoresSet(Map<String, RawScoreModel> resultsMap, Set<String> refMatches) {
//        TreeMap<RawScoreModel, String> sorter = new TreeMap<>(Collections.reverseOrder());
        List<RawScoreModel> scoreModelSorter = new ArrayList<>();

        for (String rs1Key : resultsMap.keySet()) {
            scoreModelSorter.add(resultsMap.get(rs1Key));
            List<SpectrumMatch> updatedAddOnly= SpectraUtilities.getModificationFrequentScore(rs1Key, resultsMap.get(rs1Key).getSpecTitles(), refMatches, resultsMap.get(rs1Key).getSpectrumMatchResult());
            resultsMap.get(rs1Key).setSpectrumMatchResult(updatedAddOnly);
        }
        Collections.sort(scoreModelSorter);
        Collections.reverse(scoreModelSorter);
        int index1 = -1;
        List<String> topScoreSet = new ArrayList<>(resultsMap.keySet());
        for (RawScoreModel score1 : scoreModelSorter) {
            String rs1Key = score1.getComparisonId();
            index1++;
            if (!topScoreSet.contains(rs1Key)) {
                continue;
            }
            int index2 = -1;
            for (RawScoreModel score2 : scoreModelSorter) {
                String rs2Key = score2.getComparisonId();
                index2++;
                if (index2 <= index1 || !topScoreSet.contains(rs2Key)) {
                    continue;
                }
                Set<String> total = new HashSet<>(resultsMap.get(rs1Key).getSpecTitles());
                total.addAll(resultsMap.get(rs2Key).getSpecTitles());
                double key1Better = isBetterScore(resultsMap.get(rs2Key).getSpectrumMatchResult(), resultsMap.get(rs1Key).getSpectrumMatchResult(), total.size());
                System.out.println(rs1Key + " better " + rs2Key + "  " + key1Better);
                System.out.println("resultsMap.get(rs1Key) " + resultsMap.get(rs1Key).getTotalNumber() + "(" + resultsMap.get(rs1Key).getFinalScore() + ")   vs " + resultsMap.get(rs2Key).getTotalNumber() + "(" + resultsMap.get(rs2Key).getFinalScore() + ")" + "   " + key1Better);

                if ((key1Better > 0) || Double.isNaN(key1Better)) {// && (resultsMap.get(rs1Key).getFinalScore() > resultsMap.get(rs2Key).getFinalScore())
                    System.out.println("-----------------------------------remove " + rs2Key + "-------------------------------");
                    topScoreSet.remove(rs2Key);
                } else {
                    System.out.println("-----------------------------------remove " + rs1Key + "-------------------------------");
                    topScoreSet.remove(rs1Key);
                    break;
                }
            }
        }
        if (topScoreSet.size() > 1) {
            System.out.println("error should be only oine value " + topScoreSet);
        } else if (topScoreSet.isEmpty()) {
            System.out.println("error should be only oine value " + topScoreSet);
            return "";
        }
        String topSelection = topScoreSet.get(0);
        if (!scoreModelSorter.get(0).getComparisonId().equalsIgnoreCase(topSelection)) {
            topSelection = topSelection + "_-_" + scoreModelSorter.get(0).getComparisonId();
            System.out.println("-------------------------------->> disagree with both methods " + scoreModelSorter.get(0).getComparisonId() + "  " + topSelection);
        }
        return topSelection;
    }

    public static double isBetterScore(List<SpectrumMatch> fromData, List<SpectrumMatch> toData, int totalSpeNumber) {
        List<Double> sharedToData = new ArrayList<>();
        List<Double> sharedFromData = new ArrayList<>();
        List<Double> onlyToData = new ArrayList<>();
        List<Double> onlyFromData = new ArrayList<>();
        Map<String, Double> fullMatchScores = new HashMap<>();
        Map<String, Double> toMatchScores = new HashMap<>();
        if (fromData == null && toData != null) {
            return 1;
        }
        if (fromData != null && toData == null) {
            return -1;
        }
        if (fromData == null && toData == null) {
            return -1;
        }
        toData.stream().map(sm -> {
            fullMatchScores.put(sm.getSpectrumTitle(), 0.0);
            return sm;
        }).forEachOrdered(sm -> {
            toMatchScores.put(sm.getSpectrumTitle(), sm.getBestPeptideAssumption().getRawScore());
        });
        fromData.forEach(sm -> {
            fullMatchScores.put(sm.getSpectrumTitle(), sm.getBestPeptideAssumption().getRawScore());
        });

        for (String titile : fullMatchScores.keySet()) {
            double fromScore = fullMatchScores.get(titile);
            if (toMatchScores.containsKey(titile)) {
                double toScore = toMatchScores.get(titile);
                if (fromScore > 0) {
                    sharedToData.add(toScore);
                    sharedFromData.add(fromScore);
                } else {
                    onlyToData.add(toScore);
                }
            } else if (fromScore > 0) {
                onlyFromData.add(fromScore);
            }
        }
        Collections.sort(sharedFromData);
//        System.out.println("at total score is better is " + " from: " + onlyFromData.size() + " - share: " + sharedToData.size() + " - to: " + onlyToData.size());
        //devide the TODATA data into 3 category         
        double[] fromSharedData = sharedFromData.stream().mapToDouble(Double::doubleValue).toArray();
        double[] toSharedData = sharedToData.stream().mapToDouble(Double::doubleValue).toArray();
        double[] toOnlyData = onlyToData.stream().mapToDouble(Double::doubleValue).toArray();
        double[] fromOnlyData = onlyFromData.stream().mapToDouble(Double::doubleValue).toArray();
        double[] improved1 = SpectraUtilities.compareData(fromSharedData, toSharedData);
        double score2;
//        if (onlyFromData.size() == onlyToData.size()) {
        score2 = SpectraUtilities.compareData(fromOnlyData, toOnlyData)[4];
//        } else {
//            score2 = SpectraUtilities.compareUnPairedData(fromOnlyData, toOnlyData);
//        }

//       double
        double ratio = ((double) onlyToData.size() / totalSpeNumber);
//         double total =sharedFromData.size()+onlyFromData.size()+onlyToData.size();  
//         ratio = (double)(total-onlyFromData.size()-onlyToData.size())/total;

        double fs1 = improved1[4];//* (1.0 - ratio);
        double fs2 = score2 * ratio;

        double fs = fs1 + fs2;
//        if (Double.isNaN(fs)) {
        System.out.println("score " + fs1 + "   " + fs2 + "  " + score2 + "*" + ratio);
//            System.exit(0);
//        }
        return fs;

    }

    public static double[] compareData(double[] fromData, double[] toData) {
        ScoreComparison sc = new ScoreComparison();
        double[] finalScores = sc.calculateFinalScore(fromData, toData);
        return finalScores;

    }

    public static double calculateDatasetScoreThreshold(double oreginalDataSize, double subDataSize, double identificationRate, double idNumber) {
        double comparisonsThreshold = oreginalDataSize / subDataSize;
        comparisonsThreshold = comparisonsThreshold / 10.0;

        double idToFullFactor = oreginalDataSize / (idNumber * 100.0);

//        System.out.println("at level one " + comparisonsThreshold + "   identificationRate " + identificationRate + "  " + idToFullFactor);
        comparisonsThreshold = ScoreComparison.logScaleNormalize(comparisonsThreshold, Math.E);
        idToFullFactor = ScoreComparison.logScaleNormalize(idToFullFactor, Math.E);
        identificationRate = ScoreComparison.logScaleNormalize(identificationRate, Math.E);
//        System.out.println("at level one normalised " + comparisonsThreshold + "   identificationRate " + identificationRate + "  " + idToFullFactor);
        comparisonsThreshold = (comparisonsThreshold + identificationRate + idToFullFactor) / 3.0;
//        System.out.println("at level two normalised " + comparisonsThreshold);
        comparisonsThreshold = Math.floor(comparisonsThreshold * 100.0) / 100.0;
//        System.out.println("at finalnormalised " + comparisonsThreshold);
        return comparisonsThreshold;

    }

    public static double compareUnPairedData(double[] fromData, double[] toData) {//, List<Double> sharedFromData, List<Double> sharedToData
//
//        double[] updatedFrom;
//        double[] updatedTo;
//        double[] updatedleftFrom;
//        double[] updatedleftTo;
        double finalScore;
//
//        if (fromData.length == toData.length) {
//            System.out.println("case I same length ");
        finalScore = SpectraUtilities.compareData(fromData, toData)[4];
//        }

//        else if (fromData.length < toData.length) {
//            int medianIndex = SpectraUtilities.getMedianIndex(sharedFromData);
////            System.out.println("case II fromData.length < toData.length");
//            updatedFrom = fromData;
//            updatedTo = new double[fromData.length];
//
//            updatedleftTo = new double[toData.length - fromData.length];
//            updatedleftFrom = new double[updatedleftTo.length];
//            // new double[fromData.length + toData.length];
//            Arrays.sort(toData);
//            Arrays.sort(updatedFrom);
//            for (int i = 0; i < updatedTo.length; i++) {
//                updatedTo[i] = toData[toData.length - 1 - i];
//            }
////            System.out.println("at median index " + medianIndex + "  " + sharedFromData.size() + "  " + updatedleftTo.length);
//            for (int i = 0; i < updatedleftTo.length; i++) {
//                updatedleftTo[i] = toData[i];
//                if (medianIndex + i < sharedFromData.size() && medianIndex >= 0) {
//                    updatedleftFrom[i] = sharedFromData.get(medianIndex + i);
//                }
//            }
//
//            double score1 = SpectraUtilities.compareData(updatedFrom, updatedTo)[4];
//            double score2 = SpectraUtilities.compareData(updatedleftFrom, updatedleftTo)[4];
//            double ratio1 = (double) updatedFrom.length / (double) (toData.length);
//            double ratio2 = 1.0 - ratio1;
//            finalScore = (score1 * 1) + (score2 * ratio2);
////            System.out.println("--------------------------------->>>  " + score1 + " * " + 1 + " + " + score2 + " * " + ratio2 + "  " + (updatedFrom.length) + "   " + (toData.length));
////              System.out.println("--------------------------------->>>  " + score1 + " * " + 1 + " + " + score2 + " * " + ratio2 + "  " + (updatedleftFrom.length) + "   " + (updatedleftTo.length));
//
//        } else {
//            System.out.println(" here we go time to compare ");
//            
//            if (fromData.length * 0.9 > toData.length + sharedToData.size()) {
//                return -1.0;
//            }
//            int medianIndex = SpectraUtilities.getMedianIndex(sharedFromData);
////            System.out.println("case III fromData.length > toData.length");
//            updatedFrom = new double[toData.length];
//            updatedleftFrom = new double[fromData.length - toData.length];
//            updatedleftTo = new double[updatedleftFrom.length];
//            updatedTo = toData;// new double[fromData.length + toData.length];
//            Arrays.sort(fromData);
//            Arrays.sort(updatedTo);
//            for (int i = 0; i < updatedFrom.length; i++) {
//                updatedFrom[i] = fromData[updatedFrom.length - i];
//            }
//
////            System.out.println("at median index " + medianIndex + "  " + sharedToData.size() + "  " + updatedleftFrom.length);
//            for (int i = 0; i < updatedleftFrom.length; i++) {
//                updatedleftFrom[i] = fromData[i];
//                if (medianIndex + i < sharedFromData.size()) {
//                    updatedleftTo[i] = sharedFromData.get(medianIndex + i);
//                }
//            }
//
//            double score1 = SpectraUtilities.compareData(updatedFrom, updatedTo)[4];
//            //for score 2 the higher the worest is the score should be ...mean we get rid of bas spectra
//            double score2 = SpectraUtilities.compareData(updatedleftFrom, updatedleftTo)[4] * -1.0;
//
//            double ratio1 = (double) updatedFrom.length / (double) (fromData.length);
//            double ratio2 = 1.0 - ratio1;
//            finalScore = (score1 * 1) + (score2 * ratio2);
////            System.out.println("--------------------------------->>> other case from data is larger " + updatedFrom.length + "  " + updatedTo.length + " --- >> " + updatedleftFrom.length + "  " + updatedleftTo.length + "  " + score1 + "  " + score2);
//
//        }
        if (Double.isNaN(finalScore)) {
            if (fromData.length > toData.length) {
                finalScore = -1;
            } else if (fromData.length < toData.length) {
                finalScore = 1;
            } else {
                finalScore = 0;
            }
        }
        return finalScore;

    }

    public static int scaleSubsetSize(int oreginalFileSize, Advocate se) {
        double minVal = 2000;
        double maxVal = 60000;
        double a = 1000;
        double b = 1500;
//        if (se.getIndex() == Advocate.sage.getIndex()) {
//            a = 2000;
//            b = 8000;
//        }

        if (maxVal <= oreginalFileSize) {
            return (int) 2000;
        }
        if (minVal >= oreginalFileSize) {
            return (int) oreginalFileSize;
        }
        // Normalize to [0, 1]
        double normalized = (oreginalFileSize - minVal) / (maxVal - minVal);
        // Scale to [a, b]
        return (int) Math.round((a + (b - a) * normalized));
    }

    public static Set<String> getModifiedSpectraSubset(String mod, List<SpectrumMatch> matches) {

        Set<String> subset = new LinkedHashSet<>();
        Modification modification = ptmFactory.getModification(mod);

        double mass = Math.round(modification.getMass());
        String updatedModName = mass + "@" + modification.getPattern().toString();
        boolean terminal = (modification.getModificationType().isCTerm() || modification.getModificationType().isNTerm()) && modification.getPattern().toString().trim().equals("");
//        System.out.println("mod " + mod + "  " + modification.getMass() + "  " + updatedModName + "  ");   
//        System.out.println("at mod "+mod+"  "+modification.getPattern().toString());
        for (SpectrumMatch sm : matches) {
            if (sm.getBestPeptideAssumption().getPeptide().getNVariableModifications() > 0) {
                ModificationMatch[] mms = sm.getBestPeptideAssumption().getPeptide().getVariableModifications();
                for (ModificationMatch mm : mms) {
                    double rMass = Math.round(Double.parseDouble(mm.getModification().split("@")[0]));//    mm.getModification().split("@")[0].split("\\.")[1].length();
                    String modAsString;
                    if (terminal) {
                        modAsString = rMass + "@";
                    } else {
                        modAsString = rMass + "@" + mm.getModification().split("@")[1];
                    }
//                    System.out.println("mod maching is " + modAsString + "  oreginal--> " + mm.getModification());
                    if (modAsString.equalsIgnoreCase(updatedModName)) {
                        subset.add(sm.getSpectrumTitle());
                        break;
                    }
                }

            }
        }
        return subset;
    }

    public static List<SpectrumMatch> getModificationFrequentScore(String mod, Set<String> modMaches, Set<String> refMatches, List<SpectrumMatch> specList) {

        double newAdded = 0;
        double lost = 0;
        double shared = 0;
        Set<String> fullSpecTitiles = new HashSet<>(refMatches);
        fullSpecTitiles.addAll(modMaches);
        Set<String> sharedSpec = new HashSet<>();
        Set<String> newAddedSpec = new HashSet<>();
        List<SpectrumMatch> newAddedspecList = new ArrayList<>();
        for (String sm : fullSpecTitiles) {
            if (!refMatches.contains(sm)) {
                newAdded++;
                newAddedSpec.add(sm);
            } else if (!modMaches.contains(sm)) {
                lost++;
            } else {
                shared++;
                sharedSpec.add(sm);
            }
        }
        double effect = ((double) (newAdded) / (double) (shared + lost));
//        if (mod.contains("of K")) {
        System.out.println(mod + "---->>> \t" + effect);
//        }
        for (SpectrumMatch sm : specList) {
            if (newAddedSpec.contains(sm.getSpectrumTitle())) {
                newAddedspecList.add(sm);
            }

        }

        return newAddedspecList;//((lost) / (double) shared)*(newAdded);
    }

    public static double[] isPotintialVariableModification(String modPattern, Map<String, SpectrumMatch> modSubSet, Map<String, SpectrumMatch> nonModRefsubSet) {

        double frequency = 0;
        char m = modPattern.charAt(0);
        for (String key : nonModRefsubSet.keySet()) {
            SpectrumMatch refMach = nonModRefsubSet.get(key);
            String refSequence = refMach.getBestPeptideAssumption().getPeptide().getSequence();
            if (refSequence.contains(modPattern)) {
//                frequency++;
            }

        }
        double potintialPSM = 0;
        for (String key : modSubSet.keySet()) {
            SpectrumMatch vmMach = modSubSet.get(key);
            String vmodSequence = vmMach.getBestPeptideAssumption().getPeptide().getSequence();
            if (vmodSequence.contains(modPattern) && !nonModRefsubSet.containsKey(key)) {
                int freq = SpectraUtilities.getFraquent(modPattern.charAt(0), vmodSequence);
                if (freq == 2 || freq == 1) {
                    potintialPSM++;
                }
                frequency += freq;
            }
        }
        return new double[]{potintialPSM, frequency};
    }

    private static int getFraquent(char pattern, String seq) {
        int freq = 0;
        for (char c : seq.toCharArray()) {
            if (c == pattern) {
                freq++;
            }
        }

        return freq;
    }

    public static double getQuartileValue(List<Double> sortedList, double percentile) {
        // Sort the list to find the median
        Collections.sort(sortedList);
        int n = sortedList.size();
        double position = (n + 1) * percentile;
        int integerPart = (int) position;
        double fractionalPart = position - integerPart;

        if (integerPart == 0) {
            return sortedList.get(0);
        } else if (integerPart >= n) {
            return sortedList.get(n - 1);
        } else {
            return sortedList.get(integerPart - 1) + fractionalPart
                    * (sortedList.get(integerPart) - sortedList.get(integerPart - 1));
        }
    }

    public static int getMedianIndex(List<Double> list) {
        // Sort the list to find the median
        Collections.sort(list);
        int size = list.size();
        int medianIndex;
        if (size % 2 == 1) {
            // Odd number of elements, median is the middle one
            medianIndex = size / 2;
        } else {
            // Even number of elements, median is the average of the two middle ones
            // Here, we'll return the lower median index (size / 2 - 1)
            medianIndex = (size / 2) - 1;
        }

        return medianIndex;
    }

    public static void main(String[] args) {
        double[] from = new double[]{5, 10, 15, 20, 25, 30, 40, 50};
        double[] to = new double[]{10, 15, 20, 25, 30, 40, 50, 60};
        compareData(to, to);
    }

    public static int findMainDrop(int[] numbers) {
        if (numbers.length < 2) {
            System.out.println("Not enough data points to detect a drop.");
            return 0;
        }

        double maxDrop = 0;  // Largest drop value
        double previous = 0; // Previous number in the sequence
        double next = 0;     // Next number in the sequence
        int selectedIndex = 0;
        for (int i = 1; i < numbers.length; i++) {
            double drop = numbers[i - 1] - numbers[i];  // Calculate the drop between consecutive elements

            if (drop > maxDrop) {
                maxDrop = drop;
                previous = numbers[i - 1];
                selectedIndex = (i - 1);
                next = numbers[i];
            }
        }

        if (maxDrop > 0) {
            System.out.println("The main drop is from " + previous + " to " + next + " with a drop of " + maxDrop + "  " + selectedIndex);
        } else {
            System.out.println("No drop detected.");
        }
//        System.exit(0);
        return selectedIndex;
    }

}
