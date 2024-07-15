package no.uib.probe.optprot.util;

import com.compomics.util.experiment.biology.modifications.Modification;
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
import java.util.HashSet;
import java.util.LinkedHashMap;
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
import org.apache.commons.math3.util.Precision;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author yfa041
 */
public class SpectraFileUtilities {

    private final SearchInputSetting searchOptimizerParameters;

    public SpectraFileUtilities() {
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
            System.out.println("file name ms " + filteredSubMsFile.getName() + "   " + filteredSubMsFile.exists());
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
                    System.out.println("final size is " + spectrumMap.size() + "   " + currentSpectrumSize);

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
                Set<String> sequences = SpectraFileUtilities.getSequences(NovorFile);
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

    private File initSubFastaFile(File fastaFile, Set<String> sequences) {

        File subFastaFile = new File(fastaFile.getParent(), Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_" + fastaFile.getName());
        if (subFastaFile.exists()) {
            subFastaFile.delete();
        }
        SpectraFileUtilities.createSubFastaFile(fastaFile, subFastaFile, sequences);
//        SpectraFileUtilities.createSubFastaFile(fastaFile, subFastaFile,true);
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
            SpectraFileUtilities.writeFilteredMgfFile(oreginalMsFile, subMsFile, spectraTags);
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
//                SpectraFileUtilities.writeSubSetTargtedFileAreaMgfFile(oreginalMsFile, subMsFile, Configurations.EXTRACT_MIN_MS_SIZE, Configurations.EXTRACT_MAX_MS_SIZE);
//            } else {
//                SpectraFileUtilities.writeSubSetEveryNMgfFile(oreginalMsFile, subMsFile, Configurations.EXTRACT_MIN_MS_SIZE, Configurations.EXTRACT_MAX_MS_SIZE);
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
            System.out.println("initial zer " + initiSize + "  vs " + maxSpectNum);

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
            System.out.println("total filtered subset mgf " + counter);
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
            System.out.println("initial zer " + initiSize + "  vs " + maxSpectNum);

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
            System.out.println("at n: " + n + "  " + targetedAreaStep * 2 + "  " + spectrumTitles.length);
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
                    String accession = protein.getAccession();
                    String sequence = protein.getSequence();

                    String rawHeader = header.getRawHeader();

                    bw.write(rawHeader);
                    bw.newLine();
                    bw.write(sequence);

                    bw.newLine();
//                    bw.newLine();
//
//                    int accessionEndIndex = rawHeader.indexOf(accession) + accession.length();
//
//                    String part0 = rawHeader.substring(0, rawHeader.indexOf(accession));
//                    String part1 = rawHeader.substring(rawHeader.indexOf(accession), accessionEndIndex);
//                    String part2 = rawHeader.substring(accessionEndIndex);
//
//                    bw.write(part0);
//
//                    bw.write(part1);
//
//                    bw.write(part2);
//
//                    bw.newLine();
//
//                    char[] sequenceAsArray = protein.getSequence().toCharArray();
//
//                    for (int i = sequenceAsArray.length - 1; i >= 0; i--) {
//
//                        char aa = sequenceAsArray[i];
//
//                        bw.write(aa);
//
//                    }

//                    bw.newLine();
//                    bw.newLine();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
            System.out.println("total: " + total + "  left " + left);
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
                System.out.println("Fasta copy all total: " + total + "  included " + included + "  seq not exist ");
                copySubFastaFile(fastaIn, fastaOut);
            }
            System.out.println("Fasta total: " + total + "  included " + included + "  seq not exist " + specIndex);

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

            System.out.println("Fasta total: " + total + "  left " + left + "  seq not exist " + specIndex);

        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            Logger.getLogger(SpectraFileUtilities.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

//
//    public static Set<String> getTagsSequences(File directTag) {
//        String line;
//        String splitBy = ",";
//        Set<String> sequences = new HashSet<>();
//        try {
////parsing a CSV file into BufferedReader class constructor  
//            BufferedReader br = new BufferedReader(new FileReader(directTag));
////            while ((line = br.readLine()) != null) //returns a Boolean value 
////            {
////                if (line.startsWith("# id,")) {
////                    break;
////                }
////            }
//
//            while ((line = br.readLine()) != null) //returns a Boolean value 
//            {
//                System.out.println("line: " + line);
//            }
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//        return sequences;
//    }
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

    /**
     * Lets the user select an output folder and starts the recalibration of
     * spectra.
     *
     * @param sourceMgfFile
     * @param fastaFile
     * @param spectrumMatches
     * @param identificationParameters
     */
    public static File recalibrateSpectra(File sourceMgfFile, File fastaFile, ArrayList<SpectrumMatch> spectrumMatches, IdentificationParameters identificationParameters) {

        try {
            final File selectedFolder = new File(sourceMgfFile.getParent());
            MsFileHandler spectrumProvider = new MsFileHandler();
            spectrumProvider.register(sourceMgfFile, new OptProtWaitingHandler());
            final boolean precursors = true;
            final boolean fragments = true;
//            ModificationParameters mp =  new ModificationParameters();//identificationParameters.getSearchParameters().getModificationParameters();
//            identificationParameters.getSearchParameters().setModificationParameters(null);
            FMIndex sequenceProvider = new FMIndex(fastaFile, null, new OptProtWaitingHandler(), false, identificationParameters);
//            identificationParameters.getSearchParameters().setModificationParameters(mp);
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

            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                System.out.println("clean after");
            }

        } catch (IOException ex) {
            Logger.getLogger(SpectraFileUtilities.class.getName()).log(Level.SEVERE, null, ex);
        }
        return null;
    }

    public Map<String, Spectrum> substractSpectraWithConfidentTag(int startIndex, int maxSpectraNumber, MsFileHandler msFileHandler, String[] spectrumTitles, String msFileNameWithoutExtension) {
        Map<String, Spectrum> spectraMap = new LinkedHashMap<>();
        try {

            int targetedAreaStep = spectrumTitles.length / 4;
//            targetedAreaStep=0;
//            int n = spectrumTitles.length / (maxSpectraNumber);
            int n = targetedAreaStep * 2 / (maxSpectraNumber);
            File destinationFile = new File(oreginalMsFile.getParentFile(), Configurations.DEFAULT_RESULT_NAME + "_" + msFileNameWithoutExtension + ".mgf");
            destinationFile.createNewFile();
            MgfFileWriter writer = new MgfFileWriter(destinationFile);

            for (int i = startIndex + targetedAreaStep; i < (spectrumTitles.length - targetedAreaStep);) {
                Spectrum spectrum = msFileHandler.getSpectrum(msFileNameWithoutExtension, spectrumTitles[i]);
                writer.writeSpectrum(spectrumTitles[i], spectrum);
                spectraMap.put(spectrumTitles[i], spectrum);
                i += n;
            }
            writer.close();
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
            System.out.println("tag size before" + spectraMap.size());
            for (SpectrumMatch sm : matches) {
                TagAssumption tag = sm.getAllTagAssumptions().toList().get(0);
                if (tag.getScore() > 0.01) {
                    spectraMap.remove(sm.getSpectrumTitle());
                }
            }
            System.out.println("tag size after" + spectraMap.size());
            destinationFile.delete();
            File cms = new File(destinationFile.getParent(), destinationFile.getName().replace(".mgf", ".cms"));
            cms.delete();

        } catch (IOException | SQLException | ClassNotFoundException | InterruptedException | JAXBException | XmlPullParserException | XMLStreamException ex) {
            ex.printStackTrace();
        }
        return spectraMap;
    }

//    public static double[] getValidatedIdentificationReferenceData(File resultOutput, File msFile, Advocate searchEngine, IdentificationParameters identificationParameters, int validatedIdNumber) {
//
//        try {
//            File idResultsFile = null;
//            if (searchEngine.getIndex() == Advocate.myriMatch.getIndex()) {
//                MyriMatchParameters myriMatchParameters = (MyriMatchParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.myriMatch.getIndex());
//                idResultsFile = new File(resultOutput, SearchHandler.getMyriMatchFileName(IoUtil.removeExtension(msFile.getName()), myriMatchParameters));
//            } else if (searchEngine.getIndex() == Advocate.xtandem.getIndex()) {
//                idResultsFile = new File(resultOutput, SearchHandler.getXTandemFileName(IoUtil.removeExtension(msFile.getName())));
//            }
//            if (idResultsFile == null || !idResultsFile.exists()) {
//                System.out.println("id result file was null " + resultOutput.getAbsolutePath());
//                return new double[]{};
//            }
//            IdfileReader idReader = IdfileReaderFactory.getInstance().getFileReader(idResultsFile);
//            MsFileHandler subMsFileHandler = new MsFileHandler();
//            subMsFileHandler.register(msFile, MainUtilities.OptProt_Waiting_Handler);
//            ArrayList<SpectrumMatch> matches = idReader.getAllSpectrumMatches(subMsFileHandler, MainUtilities.OptProt_Waiting_Handler, identificationParameters.getSearchParameters());
//            if (searchEngine.getIndex() == Advocate.myriMatch.getIndex()) {
//                double[] data = new double[matches.size()];
//                int index = 0;
//                for (SpectrumMatch sm : matches) {
//                    for (PeptideAssumption peptideAssumtion : sm.getAllPeptideAssumptions().toList()) {
//                        if (peptideAssumtion.getRank() != 1) {
//                            continue;
//                        }
//                        double score = -Math.log(peptideAssumtion.getScore());
//                        data[index++] = score;
//                        break;
//                    }
//                }
//                return data;
//            }
////            else if (searchEngine.getIndex() == Advocate.xtandem.getIndex()) {
////                double[] data = new double[matches.size()];
////                int index = 0;
////                for (SpectrumMatch sm : matches) {
////                    for (PeptideAssumption peptideAssumtion : sm.getAllPeptideAssumptions().toList()) {
////                        if (peptideAssumtion.getRank() != 1) {
////                            continue;
////                        }
////                        double score = -Math.log(peptideAssumtion.getScore());
////                        data[index++] = score;
////                        break;
////                    }
////                }
////                return data;
////            }
//
//        } catch (IOException | SQLException | ClassNotFoundException | InterruptedException | JAXBException | XmlPullParserException | XMLStreamException ex) {
//            System.out.println("no.uib.probe.optprot.util.SpectraFileUtilities.getValidatedIdentificationResults() " + ex.getLocalizedMessage());
//            ex.printStackTrace();
//        }
//        return new double[]{};
//    }
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
//                    if(score<21)
//                        score=0;
//                    else 
//                        score=1;
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

    public static RawScoreModel getComparableRawScore(SearchingSubDataset optProtDataset, List<SpectrumMatch> matches, Advocate searchEngine, boolean pairData) {

        double[] data = new double[matches.size()];
        int index = 0;
        RawScoreModel rawScore = new RawScoreModel();
        if (searchEngine.getIndex() == Advocate.xtandem.getIndex()) {
            for (SpectrumMatch sm : matches) {
                for (PeptideAssumption peptideAssumtion : sm.getAllPeptideAssumptions().toList()) {
                    if (peptideAssumtion.getRank() != 1) {
                        continue;
                    }
//                    double score = -Math.log(peptideAssumtion.getScore());
                    sm.setBestPeptideAssumption(peptideAssumtion);
                    data[index++] = peptideAssumtion.getRawScore();
                    break;
                }
            }

        } else if (searchEngine.getIndex() == Advocate.myriMatch.getIndex()) {
            for (SpectrumMatch sm : matches) {
                for (PeptideAssumption peptideAssumtion : sm.getAllPeptideAssumptions().toList()) {
                    if (peptideAssumtion.getRank() != 1) {
                        continue;
                    }
                    double score = -Math.log(peptideAssumtion.getScore());
                    sm.setBestPeptideAssumption(peptideAssumtion);
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
        }
        double[] fromData;
        double[] toData;
//        if (data.length == optProtDataset.getValidatedIdRefrenceData().length) {
        fromData = optProtDataset.getValidatedIdRefrenceData();
        toData = data;
//        } else if (data.length > optProtDataset.getValidatedIdRefrenceData().length) {
//            double[] updatedData = new double[data.length];
//            System.arraycopy(optProtDataset.getValidatedIdRefrenceData(), 0, updatedData, 0, optProtDataset.getValidatedIdRefrenceData().length);
//            fromData = updatedData;
//            toData = data;
//            System.out.println("updated from data to " + data.length + "  from become " + fromData.length + " instead of " + optProtDataset.getValidatedIdRefrenceData().length);
//
//        } else {
//            double[] updatedData = new double[optProtDataset.getValidatedIdRefrenceData().length];
//            System.arraycopy(data, 0, updatedData, 0, data.length);
//            fromData = optProtDataset.getValidatedIdRefrenceData();
//            toData = updatedData;
//
//        }

        double[] improved = SpectraFileUtilities.compareData(fromData, toData, pairData);
        rawScore.setTotalNumber(data.length);
        rawScore.setImprovmentScore(improved[0]);
        rawScore.settTestStat(improved[1]);
        rawScore.setpValue(improved[2]);
        rawScore.setFinalScore(improved[3]);
        boolean accepted = true;
        if (improved[0] < 1) {
            accepted = false;
        } else if (improved[1] < 0 && improved[2] < optProtDataset.getpValueThresholds()) {
            accepted = false;
        } else if (Double.isNaN(improved[3])) {
            accepted = false;
        }
//        else {
//            System.out.println("pass all check time to compare score " + improved[3]);
//        }
//        System.out.println("rawScore.getTotalNumber() > optProtDataset.getValidatedIdRefrenceData().length " + rawScore.getTotalNumber() + " > " + optProtDataset.getValidatedIdRefrenceData().length + "  " + optProtDataset.getpValueThresholds());
        rawScore.setSensitiveChange(rawScore.getImprovmentScore() > 0 && rawScore.gettTestStat() > 0 && rawScore.getTotalNumber() >= optProtDataset.getValidatedIdRefrenceData().length);
        rawScore.setSameData(rawScore.getImprovmentScore() == 0 && Double.isNaN(rawScore.gettTestStat()));
        if (optProtDataset.isFullDataSpectaInput()) {
            accepted = rawScore.isSensitiveChange();
        }
        rawScore.setSignificatChange(accepted);

//        System.out.println("is it accepted score " + accepted + "  " + rawScore);
        if (rawScore.isSignificatChange() || rawScore.isSensitiveChange() || rawScore.isSameData()) {
            rawScore.setData(data);
        }
        return rawScore;

    }

//    public static double[] compare2RawScoresData(RawScoreModel rs1, RawScoreModel rs2) {
//
////        double ttest = StatisticsTests.tTest(rs1.getData(), rs2.getData());
////        double pValue = StatisticsTests.pValueForTTest(rs1.getData(), rs2.getData());
////        boolean significatChange = StatisticsTests.isSignificatChange(ttest, pValue, rs1.getData().length, rs2.getData().length);
////        ScoreComparison sc = new ScoreComparison();
////        double improvmentScore;
////        double finalScore;
////        if (rs1.getData().length == rs2.getData().length) {
////            improvmentScore = sc.percentageImprovement(rs2.getData(), rs1.getData());
////            finalScore = sc.calculateFinalScore(rs2.getData(), rs1.getData());
////        } else if (rs1.getData().length > rs2.getData().length) {
////            double[] updatedData = new double[rs1.getData().length];
////            System.arraycopy(rs2.getData(), 0, updatedData, 0, rs2.getData().length);
////            improvmentScore = sc.percentageImprovement(updatedData, rs1.getData());//  
////            finalScore = sc.calculateFinalScore(updatedData, rs1.getData());// 
////        } else {
////            double[] updatedData = new double[rs2.getData().length];
////            System.arraycopy(rs1.getData(), 0, updatedData, 0, rs1.getData().length);
////            improvmentScore = sc.percentageImprovement(rs2.getData(), updatedData);// 
////            finalScore = sc.calculateFinalScore(rs2.getData(), updatedData);
////        }
////        System.out.println("final score 2 " + improvmentScore + "  " + "  " + finalScore);
//        double[] improved = SpectraFileUtilities.compareData(rs1.getData(), rs2.getData());
//
//        double r = (rs2.getImprovmentScore() - rs1.getImprovmentScore()) / rs1.getImprovmentScore();
//
//        System.out.println("result of score comparisons " + improved + " vs " + r);
//        return improved;
//
//    }
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
//                System.exit(0);
                return validatedMaches;
            }
            IdfileReader idReader = IdfileReaderFactory.getInstance().getFileReader(idResultsFile);
            MsFileHandler subMsFileHandler = new MsFileHandler();
            subMsFileHandler.register(msFile, MainUtilities.OptProt_Waiting_Handler);
            ArrayList<SpectrumMatch> matches = idReader.getAllSpectrumMatches(subMsFileHandler, MainUtilities.OptProt_Waiting_Handler, identificationParameters.getSearchParameters());
            validatedMaches.addAll(matches);

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

    public static boolean isBetterScore(double[] data1, double[] data2, boolean pairData) {
        ScoreComparison sc = new ScoreComparison();
        if (pairData) {
            if (data2.length > data1.length) {
                double[] updatedData = new double[data2.length];
                System.arraycopy(data1, 0, updatedData, 0, data1.length);
                data1 = updatedData;
            } else if (data2.length < data1.length) {
                double[] updatedData = new double[data1.length];
                System.arraycopy(data2, 0, updatedData, 0, data2.length);
                data2 = updatedData;
            }
        }
//        double p = sc.(fromData, toData);
//        double pvalue = sc.performTTest(data1, data2);
//        double tvalue = sc.performTStatTest(data1, data2);
        double improvmentScore = sc.percentageImprovement(data1, data2);
//        double finalScore = sc.calculateFinalScore(data1, data2);
//        System.out.println("final improvment score " + improvmentScore + "  t " + tvalue + " Pvalue " + pvalue + "  Final:   " + finalScore);
////        if (improvmentScore <= 0) {
//            return new double[]{-1, -1, -1, 1};
//        }
//         if (tvalue < 0 && pvalue <= 0.1) {
//            return false;
//        } 
//        if (tvalue<0 && pvalue < 0.05) {
//             System.out.println("----------------------supposed to be false but -------------------------finalScore "+finalScore+", improvmentScore "+improvmentScore+"----");
//             return new double[]{-1, -1,-1,1};
//            
//        }
//        else if (finalScore < 10) {
//            return false;
//        } else {
//            System.out.println("----------------------return true-----------------------------");
//
//        }

        return improvmentScore > 0;

    }

    public static double[] compareData(double[] fromData, double[] toData, boolean pairData) {
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
//        double p = sc.(fromData, toData);
        double pvalue = sc.performTTest(fromData, toData);
        double tvalue = sc.performTStatTest(fromData, toData);
        double improvmentScore = sc.percentageImprovement(fromData, toData);
        double finalScore = sc.calculateFinalScore(fromData, toData);
//        if (improvmentScore == 0) {
//            System.out.println("final improvment score " + improvmentScore + "  t " + tvalue + " Pvalue " + pvalue + "  Final:   " + finalScore);
//            return new double[]{-1, -1, 1, -1};
//        }
//        if (improvmentScore <= 0) {
//            return new double[]{-1, -1, -1, 1};
//        }
//         if (tvalue < 0 && pvalue <= 0.1) {
//            return false;
//        } 
//        if (tvalue<0 && pvalue < 0.05) {
//             System.out.println("----------------------supposed to be false but -------------------------finalScore "+finalScore+", improvmentScore "+improvmentScore+"----");
//             return new double[]{-1, -1,-1,1};
//            
//        }
//        else if (finalScore < 10) {
//            return false;
//        } else {
//            System.out.println("----------------------return true-----------------------------");
//
//        }

        return new double[]{improvmentScore, tvalue, pvalue, finalScore};

    }

    public static void main(String[] args) {
        double[] from = new double[]{5, 10, 15, 20, 25, 30, 40, 50};
        double[] to = new double[]{10, 15, 20, 25, 30, 40, 50, 60};
        compareData(to, to, true);
    }

}
