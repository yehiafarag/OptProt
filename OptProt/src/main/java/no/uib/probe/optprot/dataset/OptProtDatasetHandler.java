/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.dataset;

import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.identification.spectrum_assumptions.TagAssumption;
import com.compomics.util.experiment.io.identification.IdfileReader;
import com.compomics.util.experiment.io.identification.IdfileReaderFactory;
import com.compomics.util.experiment.io.mass_spectrometry.MsFileHandler;
import com.compomics.util.experiment.io.mass_spectrometry.mgf.MgfFileWriter;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.io.IoUtil;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.parameters.identification.tool_specific.DirecTagParameters;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.search.SearchExecuter;
import no.uib.probe.optprot.util.MainUtilities;
import no.uib.probe.optprot.util.OptProtWaitingHandler;
import no.uib.probe.optprot.util.SpectraFileUtilities;
import static no.uib.probe.optprot.util.SpectraFileUtilities.writeSpectraToFile;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author yfa041
 */
public class OptProtDatasetHandler {

    private final SearchInputSetting searchInputSetting = new SearchInputSetting();
    private final Map<String, Spectrum> spectrumMap = new LinkedHashMap<>();
    private int startIndex = 0;
    private File subFastaFile = null;
    private File filteredSubMsFile = null;
    private final Random randomGenerator = new Random(200000);

    public SearchingSubDataset generateOptProtDataset(File msFile, File fastaFile, Advocate searchEngineToOptimise, File identificationParametersFile) {
        SearchingSubDataset subDataset = new SearchingSubDataset();
        subDataset.setSubMsFile(msFile);
        subDataset.setSubFastaFile(fastaFile);
        subFastaFile = null;
        filteredSubMsFile = null;
        try {
            String fileNameWithoutExtension = IoUtil.removeExtension(msFile.getName());
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(msFile, new OptProtWaitingHandler());
            double startRatio = 0.04;
            String[] spectrumTitles = msFileHandler.getSpectrumTitles(fileNameWithoutExtension);
            int finalSpectraNumberInGeneratedFile = Math.max((int) (Configurations.REFINED_MS_SIZE * 1.5), (int) (spectrumTitles.length * startRatio));
            subDataset.setOreginalDatasize(spectrumTitles.length);
//            System.out.println("----------------------------------------------------------------------------------------------------------------->>>first round ");
            while (true) {
                int[] idNums = runOptProtSubDataset(msFile, fastaFile, Configurations.ACTIVE_SEARCH_SETTINGS_FILE, searchEngineToOptimise, fileNameWithoutExtension, finalSpectraNumberInGeneratedFile, spectrumTitles.length, msFileHandler);
                subDataset.setIdentificationNum(idNums[0]);
                subDataset.setTotalSpectraNumber(idNums[1]);
//                System.out.println("id number " + idNums[0] + " / " + (idNums[1] * 0.085) + "   " + finalSpectraNumberInGeneratedFile);
                if (idNums[0] >= (idNums[1] * 0.085)) {
                    break;
                }

                startRatio += 0.02;
//                System.out.println("----------------------------------------------------------------------------------------------------------------->>>second round " + finalSpectraNumberInGeneratedFile + "   ");
                finalSpectraNumberInGeneratedFile = (int) (finalSpectraNumberInGeneratedFile * (1.2));//Math.max((int) (Configurations.REFINED_MS_SIZE * 1.5), (int) (spectrumTitles.length * startRatio));
                MainUtilities.cleanOutputFolder();
                File cms = new File(filteredSubMsFile.getParent(), filteredSubMsFile.getName().replace(".mgf", ".cms"));
                cms.delete();
                filteredSubMsFile.delete();
                subFastaFile.delete();
//                System.out.println("filteredSubMsFile " + filteredSubMsFile.exists() + "  subFastaFile " + subFastaFile.exists() + "   updated finalSpectraNumberInGeneratedFile  " + finalSpectraNumberInGeneratedFile);
            }
//            System.out.println("-------------------------------------------head to user reference-----------------------------------------------------------------");
            int[] userReference = runOptProtSubDataset(msFile, fastaFile, identificationParametersFile, searchEngineToOptimise, fileNameWithoutExtension, startIndex, startIndex, msFileHandler);

            subDataset.setIdentificationNum(userReference[0]);
            subDataset.setTotalSpectraNumber(userReference[1]);
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        if (subFastaFile != null) {
            subDataset.setSubFastaFile(subFastaFile);
        }
        if (filteredSubMsFile != null) {
            subDataset.setSubMsFile(filteredSubMsFile);
        }
        //re run with user input search to get the reference number
//        System.out.println("final ds handleing " + subDataset.getIdentificationNum() + "/" + subDataset.getTotalSpectraNumber());
        return subDataset;
    }

    private int[] runOptProtSubDataset(File msFile, File fastaFile, File identificationParametersFile, Advocate searchEngineToOptimise, String fileNameWithoutExtension, int finalSpectraNumberInGeneratedFile, int spectrumTitlesLength, MsFileHandler msFileHandler) {

        try {
            long start1 = System.currentTimeMillis();

            IdentificationParameters identificationParameters = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            SearchParameters searchParameters = identificationParameters.getSearchParameters();
            searchParameters.getModificationParameters().clearVariableModifications();
            searchParameters.getModificationParameters().clearFixedModifications();

            for (File f : msFile.getParentFile().listFiles()) {
                if (f.getName().startsWith(Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_") && f.getName().endsWith("_" + msFile.getName())) {
                    filteredSubMsFile = f;
//                    System.out.println("sub ms file exist " + f.getAbsolutePath());
                    break;
                }
            }
            if (filteredSubMsFile == null) {
                filteredSubMsFile = new File(msFile.getParent(), Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_" + finalSpectraNumberInGeneratedFile + "_" + msFile.getName());
            }
            if (!filteredSubMsFile.exists()) {
                filteredSubMsFile = new File(msFile.getParent(), Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_" + finalSpectraNumberInGeneratedFile + "_" + msFile.getName());
                DirecTagParameters direcTagParameters = (DirecTagParameters) searchParameters.getIdentificationAlgorithmParameter(Advocate.direcTag.getIndex());
                direcTagParameters.setMaxTagCount(1);
                direcTagParameters.setTagLength(3);
                direcTagParameters.setNumChargeStates(4);
                direcTagParameters.setDuplicateSpectra(false);
                direcTagParameters.setUseChargeStateFromMS(false);
                 
                searchInputSetting.setRunDirecTag(true);

                double devidValue = 0.5;
                startIndex = (int) Math.round((double) spectrumTitlesLength * 0.25);
                if (spectrumTitlesLength < 20000) {
                    devidValue = 0.75;
                    startIndex = (int) Math.round((double) spectrumTitlesLength * 0.125);
                }
                int coverageSize = (int) Math.round((double) spectrumTitlesLength * devidValue);
                int step = coverageSize / (finalSpectraNumberInGeneratedFile);
                int lastIndex = startIndex + coverageSize;

                while (spectrumMap.size() < finalSpectraNumberInGeneratedFile) {
                    int currentSpectrumSize = finalSpectraNumberInGeneratedFile - spectrumMap.size();
                    if (currentSpectrumSize <= 0) {
                        break;
                    }
//                    currentSpectrumSize += (currentSpectrumSize * 0.1);
                    spectrumMap.putAll(substractSpectraWithConfidentTag(msFile, fastaFile, startIndex, step, lastIndex, msFileHandler, identificationParameters, identificationParametersFile));
                    if (spectrumMap.isEmpty()) {
                        return runOptProtSubDataset(msFile, fastaFile, identificationParametersFile, searchEngineToOptimise, fileNameWithoutExtension, finalSpectraNumberInGeneratedFile + 500, spectrumTitlesLength, msFileHandler);
                    } else {
                        startIndex += (step / 2);
                    }

                }
                writeSpectraToFile(refineSpectrumMap(spectrumMap, Configurations.REFINED_MS_SIZE), filteredSubMsFile);
            }

            for (File f : fastaFile.getParentFile().listFiles()) {
                if (f.getName().startsWith(Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_") && f.getName().endsWith("_" + fastaFile.getName())) {
                    subFastaFile = f;
                    break;
                }
            }
            if (subFastaFile == null) {
                subFastaFile = new File(fastaFile.getParent(), Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_" + finalSpectraNumberInGeneratedFile + "_" + fastaFile.getName());
            }

            if (!subFastaFile.exists()) {
                long start3 = System.currentTimeMillis();
                searchInputSetting.setRunNovor(true);
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + fileNameWithoutExtension + Configurations.get_current_file_fingerprent() + "_" + finalSpectraNumberInGeneratedFile;
                File resultsFolder = SearchExecuter.executeSearch(updatedName, searchInputSetting, filteredSubMsFile, fastaFile, identificationParameters, identificationParametersFile);
                File NovorFile = new File(resultsFolder, IoUtil.removeExtension(filteredSubMsFile.getName()) + ".novor.csv");
                Set<String> sequences = SpectraFileUtilities.getSequences(NovorFile);
                System.out.println("sequence from nover " + sequences.size());
                long end3rd = System.currentTimeMillis();
                double total = (end3rd - start3) / 1000.0;
                System.out.println("process III ( Novor ) in seconds: " + total);
                long start4 = System.currentTimeMillis();
                subFastaFile = initSubFastaFile(fastaFile, sequences);
                long end4th = System.currentTimeMillis();
                total = (end4th - start4) / 1000.0;
                System.out.println("process IV ( Generated Fasta) in seconds: " + total);
                long end = System.currentTimeMillis();
                total = (end - start1) / 1000.0;
                System.out.println("Total Elapsed Time for initInputSubSetFiles in seconds: " + total);
                MainUtilities.deleteFolder(NovorFile.getParentFile());
            }
            //run reference search
            identificationParameters = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            searchInputSetting.setSelectedSearchEngine(searchEngineToOptimise);
            final String option = "reference_run_" + searchEngineToOptimise;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + fileNameWithoutExtension;
            File resultsFolder = SearchExecuter.executeSearch(updatedName, searchInputSetting, filteredSubMsFile, subFastaFile, identificationParameters, identificationParametersFile);

            ArrayList<SpectrumMatch> validatedMaches = SpectraFileUtilities.readIdentificationResults(resultsFolder, filteredSubMsFile, searchEngineToOptimise, identificationParameters);

            fileNameWithoutExtension = IoUtil.removeExtension(filteredSubMsFile.getName());
            msFileHandler = new MsFileHandler();
            msFileHandler.register(filteredSubMsFile, new OptProtWaitingHandler());
            MainUtilities.deleteFolder(resultsFolder);
            return new int[]{validatedMaches.size(), msFileHandler.getSpectrumTitles(fileNameWithoutExtension).length};

//           
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        return new int[2];

    }

    private Map<String, Spectrum> substractSpectraWithConfidentTag(File msFile, File fastaFile, int startIndex, int stepSize, int lastIndex, MsFileHandler msFileHandler, IdentificationParameters identificationParameters, File identificationParametersFile) {
        Map<String, Spectrum> spectraMap = new LinkedHashMap<>();
        try {
//            System.out.println("substractSpectraWithConfidentTag() " + startIndex + "  " + stepSize + "  " + lastIndex);
            String msFileNameWithoutExtension = IoUtil.removeExtension(msFile.getName());
            String[] spectrumTitles = msFileHandler.getSpectrumTitles(msFileNameWithoutExtension);

            File destinationFile = new File(msFile.getParentFile(), Configurations.DEFAULT_RESULT_NAME + "_" + msFileNameWithoutExtension + ".mgf");
            if (destinationFile.exists()) {
                destinationFile.delete();
            }
            destinationFile.createNewFile();
            try (MgfFileWriter writer = new MgfFileWriter(destinationFile)) {
                for (int i = startIndex; i < lastIndex;) {
                    Spectrum spectrum = msFileHandler.getSpectrum(msFileNameWithoutExtension, spectrumTitles[i]);
                    writer.writeSpectrum(spectrumTitles[i], spectrum);
                    spectraMap.put(spectrumTitles[i], spectrum);
                    i += stepSize;
                }
                writer.close();
            }

            String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + msFileNameWithoutExtension + Configurations.get_current_file_fingerprent();

            File tempResultsFolder = SearchExecuter.executeSearch(updatedName, searchInputSetting, destinationFile, fastaFile, identificationParameters, identificationParametersFile);
            File direcTagFile = new File(tempResultsFolder, IoUtil.removeExtension(destinationFile.getName()) + ".tags");

            if (!direcTagFile.exists()) {
//                //delete previos sub mgf and cms files
//                destinationFile.delete();
//                File cms = new File(destinationFile.getParent(), destinationFile.getName().replace(".mgf", ".cms"));
//                cms.delete();
//                //not enough confident tag , increase the number 
//                substractSpectraWithConfidentTag(msFile, fastaFile, startIndex, maxSpectraNumber + 500, msFileHandler, identificationParameters, identificationParametersFile);
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

    private File initSubFastaFile(File fastaFile, Set<String> sequences) {

        File tempSubFastaFile = new File(fastaFile.getParent(), Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_" + fastaFile.getName());
        if (tempSubFastaFile.exists()) {
            tempSubFastaFile.delete();
        }
        SpectraFileUtilities.createSubFastaFile(fastaFile, tempSubFastaFile, sequences);
        return tempSubFastaFile;

    }

    private Map<String, Spectrum> refineSpectrumMap(Map<String, Spectrum> spectrumMap, int refineSize) {
        if (spectrumMap.size() <= refineSize) {
            return spectrumMap;
        }
        Map<String, Spectrum> refineSpectraMap = new LinkedHashMap<>();
        Set<Integer> filter = new HashSet<>();
        for (int i = 0; i < refineSize; i++) {
            int keyIndex = this.randomGenerator.nextInt(0, spectrumMap.size() - 1);
            while (filter.contains(keyIndex)) {
                keyIndex = this.randomGenerator.nextInt(0, spectrumMap.size() - 1);
            }
            String key = spectrumMap.keySet().toArray()[keyIndex].toString();
            refineSpectraMap.put(key, spectrumMap.get(key));
            filter.add(keyIndex);

        }
        spectrumMap.clear();
        spectrumMap.putAll(refineSpectraMap);
        return refineSpectraMap;
    }
}
