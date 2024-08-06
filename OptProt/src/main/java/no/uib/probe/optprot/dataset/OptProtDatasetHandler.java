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
import com.compomics.util.parameters.identification.tool_specific.SageParameters;
import com.compomics.util.parameters.identification.tool_specific.XtandemParameters;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.dataset.model.ConfidentTagSorter;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.search.SearchExecuter;
import no.uib.probe.optprot.util.MainUtilities;
import no.uib.probe.optprot.util.OptProtWaitingHandler;
import no.uib.probe.optprot.util.SpectraUtilities;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author yfa041
 */
public class OptProtDatasetHandler {

    private final SearchInputSetting searchInputSetting = new SearchInputSetting();
    private int startIndex = 0;
    private File subFastaFile = null;

    private File subMsFile = null;
//    private int counter = 0;
    private double acceptedTagEvalue;
//    private boolean smallDataset;

    public SearchingSubDataset generateOptProtDataset(File msFile, File fastaFile, Advocate searchEngineToOptimise, File subDataFolder, File identificationParametersFile, boolean wholeDataTest) {
        acceptedTagEvalue = Configurations.ACCEPTED_TAG_EVALUE;
        long start1 = System.currentTimeMillis();
        Advocate standeredReferenceSearchEngine = Advocate.xtandem;
        SearchingSubDataset subDataset = new SearchingSubDataset();
        subDataset.setSubMsFile(msFile);
        subDataset.setSubFastaFile(fastaFile);
        TreeMap<String, File> subFilesMap = new TreeMap<>(Collections.reverseOrder());

        if (!wholeDataTest) {
            subFastaFile = new File(subDataFolder, Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_" + fastaFile.getName());
            subMsFile = new File(subDataFolder, Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_" + msFile.getName());

        } else {
            subFastaFile = new File(subDataFolder, Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_Full_" + fastaFile.getName());
            subMsFile = new File(subDataFolder, Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_Full_" + msFile.getName());
        }
        for (File f : subDataFolder.listFiles()) {
            if (f.getName().toLowerCase().endsWith(".par") || f.getName().toLowerCase().endsWith(".mgf") || f.getName().toLowerCase().endsWith(".fasta") || f.getName().toLowerCase().endsWith(".txt")) {
                subFilesMap.put(f.getName().toLowerCase(), f);
            } else {
                f.delete();
            }
        }
        final String fileNameWithoutExtension = IoUtil.removeExtension(msFile.getName());
        MsFileHandler msFileHandler = new MsFileHandler();
        try {
            msFileHandler.register(msFile, new OptProtWaitingHandler());

        } catch (IOException ex) {
            if (subMsFile != null) {
                subMsFile.delete();
            }
            if (subFastaFile != null) {
                subFastaFile.delete();
            }
            ex.printStackTrace();
        }

        String[] spectrumTitles = msFileHandler.getSpectrumTitles(fileNameWithoutExtension);
        subDataset.setOreginalDatasize(spectrumTitles.length);
        int step;
        int lastIndex;
        boolean update = false;
        try {
            final IdentificationParameters identificationParameters = IdentificationParameters.getIdentificationParameters(new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE));
            if (!subMsFile.exists()) {
                update = true;
                System.out.println("sub ms not exist");
                subMsFile.createNewFile();
                //initial param to handel only one time
                SearchParameters searchParameters = identificationParameters.getSearchParameters();
                searchParameters.getModificationParameters().clearVariableModifications();
                searchParameters.getModificationParameters().clearFixedModifications();
                //initial direcTag param
                DirecTagParameters direcTagParameters = (DirecTagParameters) searchParameters.getIdentificationAlgorithmParameter(Advocate.direcTag.getIndex());
                direcTagParameters.setMaxTagCount(1);
                direcTagParameters.setTagLength(3);
                direcTagParameters.setNumChargeStates(4);
                direcTagParameters.setDuplicateSpectra(false);
                direcTagParameters.setUseChargeStateFromMS(false);
                //generate subset of spectra 
                Map<String, Spectrum> spectraMap = new LinkedHashMap<>();
                //init Q1
                int subSize;
                int left;
                List<Double> quartileRatios;
                if (wholeDataTest) {

                    quartileRatios = this.getQuartileRatios(msFile, msFileHandler, 0, 1, spectrumTitles.length, fastaFile, identificationParameters);
                    subSize = (int) Math.round(quartileRatios.get(4) * spectrumTitles.length);
                    left = subSize;
//                    spectraMap.putAll(initSubDatasetPart(msFile, msFileHandler, 0, 1, spectrumTitles.length, spectrumTitles.length, fastaFile, identificationParameters, 0.8));
//                    MainUtilities.cleanOutputFolder();
//                    System.exit(0);
                } else {
                    left = Configurations.EXTRACT_MAX_MS_SIZE;
                    subSize = SpectraUtilities.scaleSubsetSize(spectrumTitles.length);
                    System.out.println("at subsetSize " + subSize);
                    step = (int) Math.round((double) spectrumTitles.length / (double) subSize);
                    quartileRatios = this.getQuartileRatios(msFile, msFileHandler, 0, step, spectrumTitles.length, fastaFile, identificationParameters);
                }
                startIndex = 0;

                int Q_1_4_coverageSize = (int) Math.round((double) spectrumTitles.length / 4.0);
                lastIndex = startIndex + Q_1_4_coverageSize;
                if (wholeDataTest) {
                    step = 1;
                } else {
                    step = Math.max(Q_1_4_coverageSize / 3000, 1);
                }

                int qSize = (int) Math.round(Configurations.EXTRACT_MAX_MS_SIZE * quartileRatios.get(0));
                Map<String, Spectrum> Q1_spectraMap = initSubDatasetPart(msFile, msFileHandler, startIndex, step, lastIndex, qSize, fastaFile, identificationParameters, 1.0 - quartileRatios.get(0));
                spectraMap.putAll(Q1_spectraMap);
                left = left - Q1_spectraMap.size();

                startIndex = lastIndex + 1;
                lastIndex = startIndex + +Q_1_4_coverageSize;
                if (wholeDataTest) {
                    step = 1;
                } else {
                    step = Math.max(Q_1_4_coverageSize / 3000, 1);
                }
                qSize = (int) Math.round(Configurations.EXTRACT_MAX_MS_SIZE * quartileRatios.get(1));
                Map<String, Spectrum> Q2_spectraMap = initSubDatasetPart(msFile, msFileHandler, startIndex, step, lastIndex, qSize, fastaFile, identificationParameters, 1.0 - quartileRatios.get(1));
                spectraMap.putAll(Q2_spectraMap);
                left = left - Q2_spectraMap.size();
                startIndex = lastIndex + 1;
                lastIndex = startIndex + +Q_1_4_coverageSize;
                if (wholeDataTest) {
                    step = 1;
                } else {
                    step = Math.max(Q_1_4_coverageSize / 3000, 1);
                }
                qSize = (int) Math.round(Configurations.EXTRACT_MAX_MS_SIZE * quartileRatios.get(2));
                Map<String, Spectrum> Q3_spectraMap = initSubDatasetPart(msFile, msFileHandler, startIndex, step, lastIndex, qSize, fastaFile, identificationParameters, 1.0 - quartileRatios.get(2));
                spectraMap.putAll(Q3_spectraMap);
                left = left - Q3_spectraMap.size();
                startIndex = lastIndex + 1;
                lastIndex = startIndex + +Q_1_4_coverageSize;
                if (wholeDataTest) {
                    step = 1;
                } else {
                    step = Math.max(Q_1_4_coverageSize / 3000, 1);
                }
                qSize = left;
                Map<String, Spectrum> Q4_spectraMap = initSubDatasetPart(msFile, msFileHandler, startIndex, step, lastIndex, qSize, fastaFile, identificationParameters, 1.0 - quartileRatios.get(3));
                spectraMap.putAll(Q4_spectraMap);
//                }
                MainUtilities.cleanOutputFolder();
                //create stabkle subMs file
                subMsFile = generateMsSubFile(spectraMap, subMsFile);
                subFastaFile.delete();
            }
            //create stabkle subfasta file
            if (!subFastaFile.exists() || update) {
                subFastaFile.createNewFile();
                long start3 = System.currentTimeMillis();
                searchInputSetting.setRunNovor(true);
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + fileNameWithoutExtension + Configurations.get_current_file_fingerprent();
                File resultsFolder = SearchExecuter.executeSearch(updatedName, searchInputSetting, subMsFile, fastaFile, identificationParameters, new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE));
                File NovorFile = new File(resultsFolder, IoUtil.removeExtension(subMsFile.getName()) + ".novor.csv");
                Set<String> sequences = SpectraUtilities.getSequences(NovorFile);
                System.out.println("sequence from nover " + sequences.size());
                long end3rd = System.currentTimeMillis();
                double total = (end3rd - start3) / 1000.0;
                System.out.println("process III ( Novor ) in seconds: " + total);
                long start4 = System.currentTimeMillis();
                subFastaFile = initSubFastaFile(subFastaFile, fastaFile, sequences);
                long end4th = System.currentTimeMillis();
                total = (end4th - start4) / 1000.0;
                System.out.println("process IV ( Generated Fasta) in seconds: " + total);
                long end = System.currentTimeMillis();
                total = (end - start1) / 1000.0;
                System.out.println("Total Elapsed Time for initInputSubSetFiles in seconds: " + total);
            }
        } catch (IOException ex) {
            if (subMsFile != null) {
                subMsFile.delete();
            }
            if (subFastaFile != null) {
                subFastaFile.delete();
            }
            ex.printStackTrace();
        }

        if (subFastaFile != null) {
            subDataset.setSubFastaFile(subFastaFile);
        }
        if (subMsFile != null) {
            subDataset.setSubMsFile(subMsFile);
        }
        try {

            //run initial identification with user selected SE
            final IdentificationParameters identificationParameters = IdentificationParameters.getIdentificationParameters(new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE));
            searchInputSetting.setSelectedSearchEngine(searchEngineToOptimise);
            final String option = "reference_run_default_" + searchEngineToOptimise;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + fileNameWithoutExtension;

            if (standeredReferenceSearchEngine.getIndex() == Advocate.sage.getIndex()) {
                SageParameters sageParameters = (SageParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
                sageParameters.setMaxVariableMods(0);
                sageParameters.setNumPsmsPerSpectrum(1);
                sageParameters.setGenerateDecoys(false);

            } else if (standeredReferenceSearchEngine.getIndex() == Advocate.xtandem.getIndex()) {
                XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
                xtandemParameters.setProteinQuickAcetyl(false);
                xtandemParameters.setQuickPyrolidone(false);
                xtandemParameters.setStpBias(false);
                xtandemParameters.setRefine(false);
                xtandemParameters.setOutputResults("all");

            }
            String subfileNameWithoutExtension = IoUtil.removeExtension(subMsFile.getName());
            MsFileHandler subMsFileHandler = new MsFileHandler();
            subMsFileHandler.register(subMsFile, new OptProtWaitingHandler());
            //init spectra map
            subDataset.setSpectraTitiles(subMsFileHandler.getSpectrumTitles(subfileNameWithoutExtension));

            File resultsFolder = SearchExecuter.executeSearch(updatedName, searchInputSetting, subMsFile, subFastaFile, identificationParameters, new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE));
          System.out.println("system excuted "+resultsFolder.getName()+"  and done "+MainUtilities.OptProt_Waiting_Handler.isRunFinished()+"  "+Arrays.asList(resultsFolder.list())); 
//          if(MainUtilities.OptProt_Waiting_Handler.isRunFinished())
          
          List<SpectrumMatch> validatedMaches = SpectraUtilities.getValidatedIdentificationResults(resultsFolder, subMsFile, searchEngineToOptimise, identificationParameters);
            
            if (validatedMaches==null || validatedMaches.isEmpty()) {
                System.out.println("Error in the system please restart!");
                System.exit(0);
            }
            subDataset.setDefaultSettingIdentificationNum(validatedMaches.size());
            subDataset.updateValidatedIdRefrenceData(validatedMaches);

            MainUtilities.deleteFolder(resultsFolder);
            int total = subMsFileHandler.getSpectrumTitles(subfileNameWithoutExtension).length;
            subDataset.setTotalSpectraNumber(total);
            MainUtilities.cleanOutputFolder();
        } catch (IOException ex) {
            if (subMsFile != null) {
                subMsFile.delete();
            }
            if (subFastaFile != null) {
                subFastaFile.delete();
            }
            ex.printStackTrace();
        }
        return subDataset;
    }

    private Map<String, Spectrum> substractSpectraFirstLevelDataFiltering(File msFile, MsFileHandler msFileHandler, int startIndex, int stepSize, int lastIndex) {
        Map<String, Spectrum> spectraMap = new LinkedHashMap<>();
        String msFileNameWithoutExtension = IoUtil.removeExtension(msFile.getName());
        String[] spectrumTitles = msFileHandler.getSpectrumTitles(msFileNameWithoutExtension);
        for (int i = startIndex; i < lastIndex && i < spectrumTitles.length;) {
            Spectrum spectrum = msFileHandler.getSpectrum(msFileNameWithoutExtension, spectrumTitles[i]);
            spectraMap.put(spectrumTitles[i], spectrum);
            i += stepSize;
        }
        return spectraMap;
    }

    private File generateMsSubFile(Map<String, Spectrum> spectraMap, File destinationFile) {

        try {

            if (destinationFile.exists()) {
                destinationFile.delete();
            }
            destinationFile.createNewFile();
            try (MgfFileWriter writer = new MgfFileWriter(destinationFile)) {
                for (String spectrumTitle : spectraMap.keySet()) {
                    Spectrum spectrum = spectraMap.get(spectrumTitle);
                    writer.writeSpectrum(spectrumTitle, spectrum);
                }
                writer.close();
            }

        } catch (IOException ex) {
            if (subMsFile != null) {
                subMsFile.delete();
            }
            if (subFastaFile != null) {
                subFastaFile.delete();
            }
            ex.printStackTrace();
        }
        return destinationFile;

    }

    private List<Double> getQuartileRatios(File msFile, MsFileHandler msFileHandler, int startIndex, int step, int lastIndex, File fastaFile, IdentificationParameters identificationParameters) {
        ArrayList<Double> arrayList = new ArrayList<>();
        try {
            arrayList.add(0.25);
            arrayList.add(0.25);
            arrayList.add(0.25);
            arrayList.add(0.25);
            arrayList.add(1.0);

            Map<String, Spectrum> spectraMap = substractSpectraFirstLevelDataFiltering(msFile, msFileHandler, startIndex, step, lastIndex);
            System.out.println("sub file size " + spectraMap.size());
            final String fileNameWithoutExtension = IoUtil.removeExtension(msFile.getName());
            File destinationFile = new File(Configurations.GET_OUTPUT_FOLDER_PATH(), Configurations.DEFAULT_RESULT_NAME + "_temp_" + startIndex + "_" + lastIndex + "_" + spectraMap.size() + "_-_" + fileNameWithoutExtension + ".mgf");
            if (destinationFile.exists()) {
                destinationFile.delete();
            }
            destinationFile.createNewFile();
            destinationFile = generateMsSubFile(spectraMap, destinationFile);
            final String subfileNameWithoutExtension = IoUtil.removeExtension(destinationFile.getName());

            ArrayList<SpectrumMatch> matches = getTagMaches(destinationFile, fastaFile, identificationParameters, new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE), subfileNameWithoutExtension);
            if (matches.isEmpty()) {
                System.out.println("there is no tags in the file ...very poor data " + subfileNameWithoutExtension);
                return arrayList;

            }
            MsFileHandler subMsFileHandler = new MsFileHandler();
            subMsFileHandler.register(destinationFile, MainUtilities.OptProt_Waiting_Handler);

//            DataAnalysisHandler.showTagDistribution(subfileNameWithoutExtension.split("_-_")[1], subMsFileHandler.getSpectrumTitles(IoUtil.removeExtension(destinationFile.getName())), matches);
            arrayList.clear();
            arrayList.addAll(SpectraUtilities.getQuartileRatio(subMsFileHandler.getSpectrumTitles(IoUtil.removeExtension(destinationFile.getName())), matches));
            arrayList.add((double) matches.size() / (double) spectraMap.size());
        } catch (IOException ex) {
            Logger.getLogger(OptProtDatasetHandler.class.getName()).log(Level.SEVERE, null, ex);
        }
        return arrayList;
    }

    private Map<String, Spectrum> initSubDatasetPart(File msFile, MsFileHandler msFileHandler, int startIndex, int step, int lastIndex, int refinSize, File fastaFile, IdentificationParameters identificationParameters, double highQualityRatio) {
        Map<String, Spectrum> spectraMap = substractSpectraFirstLevelDataFiltering(msFile, msFileHandler, startIndex, step, lastIndex);
        System.out.println("sub file size " + spectraMap.size());
        try {

            //generate submsFile
            final String fileNameWithoutExtension = IoUtil.removeExtension(msFile.getName());
            File destinationFile = new File(Configurations.GET_OUTPUT_FOLDER_PATH(), Configurations.DEFAULT_RESULT_NAME + "_temp_" + startIndex + "_" + lastIndex + "_" + spectraMap.size() + "_-_" + fileNameWithoutExtension + ".mgf");
            if (destinationFile.exists()) {
                destinationFile.delete();
            }
            destinationFile.createNewFile();
            destinationFile = generateMsSubFile(spectraMap, destinationFile);
            final String subfileNameWithoutExtension = IoUtil.removeExtension(destinationFile.getName());
            //run direct tag and get confident tags    File destinationFile, File fastaFile, IdentificationParameters identificationParameters, File identificationParametersFile, String msFileNameWithoutExtension, int spectraSizeLimit
            Map<String, Spectrum> confidentSpectraSet = getSubSpectraWithConfidentTag(destinationFile, fastaFile, identificationParameters, new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE), subfileNameWithoutExtension, refinSize, highQualityRatio);

            return confidentSpectraSet;

        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return new HashMap<>();
    }

    private ArrayList<SpectrumMatch> getTagMaches(File destinationFile, File fastaFile, IdentificationParameters identificationParameters, File identificationParametersFile, String msFileNameWithoutExtension) {
        try {
            searchInputSetting.setRunDirecTag(true);
            String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + msFileNameWithoutExtension + Configurations.get_current_file_fingerprent();
            File tempResultsFolder = SearchExecuter.executeSearch(updatedName, searchInputSetting, destinationFile, fastaFile, identificationParameters, identificationParametersFile);
            File direcTagFile = new File(tempResultsFolder, IoUtil.removeExtension(destinationFile.getName()) + ".tags");

            if (!direcTagFile.exists()) {
                System.out.println("there is no tags in the file ...very poor data " + msFileNameWithoutExtension);
//                //delete previos sub mgf and cms files
//                destinationFile.delete();
//                File cms = new File(destinationFile.getParent(), destinationFile.getName().replace(".mgf", ".cms"));
//                cms.delete();
//                //not enough confident tag , increase the number 
//                substractSpectraWithConfidentTag(msFile, fastaFile, startIndex, maxSpectraNumber + 500, msFileHandler, identificationParameters, identificationParametersFile);
                return new ArrayList<>();
            }

            IdfileReader idReader = IdfileReaderFactory.getInstance().getFileReader(direcTagFile);
            System.out.println("file name " + msFileNameWithoutExtension + "   " + IoUtil.removeExtension(destinationFile.getName()));
            MsFileHandler subMsFileHandler = new MsFileHandler();
            subMsFileHandler.register(destinationFile, MainUtilities.OptProt_Waiting_Handler);
            ArrayList<SpectrumMatch> matches = idReader.getAllSpectrumMatches(subMsFileHandler, MainUtilities.OptProt_Waiting_Handler, identificationParameters.getSearchParameters());
            return matches;
        } catch (IOException | SQLException | ClassNotFoundException | InterruptedException | JAXBException | XmlPullParserException | XMLStreamException ex) {
            Logger.getLogger(OptProtDatasetHandler.class.getName()).log(Level.SEVERE, null, ex);
        }
        return new ArrayList<>();

    }

    private Map<String, Spectrum> getSubSpectraWithConfidentTag(File destinationFile, File fastaFile, IdentificationParameters identificationParameters, File identificationParametersFile, String msFileNameWithoutExtension, int spectraSizeLimit, double highQualityRatio) {
        Set<ConfidentTagSorter> confidentSpectraSet = new LinkedHashSet<>();
        Map<String, Spectrum> subSpectraMap = new LinkedHashMap<>();
        try {
            ArrayList<SpectrumMatch> matches = getTagMaches(destinationFile, fastaFile, identificationParameters, identificationParametersFile, msFileNameWithoutExtension);
            if (matches.isEmpty()) {
                System.out.println("there is no tags in the file ...very poor data " + msFileNameWithoutExtension);
//                //delete previos sub mgf and cms files
//                destinationFile.delete();
//                File cms = new File(destinationFile.getParent(), destinationFile.getName().replace(".mgf", ".cms"));
//                cms.delete();
//                //not enough confident tag , increase the number 
//                substractSpectraWithConfidentTag(msFile, fastaFile, startIndex, maxSpectraNumber + 500, msFileHandler, identificationParameters, identificationParametersFile);
                return subSpectraMap;
            }

            MsFileHandler subMsFileHandler = new MsFileHandler();
            subMsFileHandler.register(destinationFile, MainUtilities.OptProt_Waiting_Handler);

//            DataAnalysisHandler.showTagDistribution(msFileNameWithoutExtension.split("_-_")[1], subMsFileHandler.getSpectrumTitles(IoUtil.removeExtension(destinationFile.getName())), matches);
            for (SpectrumMatch sm : matches) {
                TagAssumption tag = sm.getAllTagAssumptions().toList().get(0);
//                if (tag.getScore() < acceptedTagEvalue) {
                confidentSpectraSet.add(new ConfidentTagSorter(tag.getScore(), sm.getSpectrumTitle(), subMsFileHandler.getSpectrum(msFileNameWithoutExtension, sm.getSpectrumTitle())));
//                }
            }
            File cms = new File(destinationFile.getParent(), destinationFile.getName().replace(".mgf", ".cms"));
            cms.delete();
            int counter = 0;
            double highQualitylimit = (double) spectraSizeLimit * highQualityRatio;
            double avgQualityLimit = spectraSizeLimit - highQualitylimit;
            int n = (int) (confidentSpectraSet.size() / highQualitylimit);
            int step = 0;
            for (ConfidentTagSorter tag : confidentSpectraSet) {
                step++;
                if (tag.getValue() > 0.01 || step < n) {
                    continue;
                }
                subSpectraMap.put(tag.getTitle(), tag.getSpectrum());
                counter++;
                step = 0;
                if (counter >= highQualitylimit) {
                    break;
                }
            }
//            highQualitylimit += avgQualityLimit - counter;
//            counter = 0;
            step = 0;
            n = (int) (confidentSpectraSet.size() / avgQualityLimit);
            for (ConfidentTagSorter tag : confidentSpectraSet) {
                step++;
                if (tag.getValue() > 0.1 || tag.getValue() <= 0.01 || step < n) {
                    continue;
                }
                subSpectraMap.put(tag.getTitle(), tag.getSpectrum());
                step = 0;
                counter++;
                if (counter >= avgQualityLimit) {
                    break;
                }
            }
            if (subSpectraMap.size() < spectraSizeLimit) {
                for (ConfidentTagSorter tag : confidentSpectraSet) {
                    if (tag.getValue() > 0.01 || subSpectraMap.containsKey(tag.getTitle())) {
                        continue;
                    }
                    subSpectraMap.put(tag.getTitle(), tag.getSpectrum());
                    if (subSpectraMap.size() >= spectraSizeLimit) {
                        break;
                    }

                }
            }

            if (subSpectraMap.size() < spectraSizeLimit) {
                for (ConfidentTagSorter tag : confidentSpectraSet) {
                    if (tag.getValue() < 0.1 || tag.getValue() >= 1.0) {
                        continue;
                    }
                    subSpectraMap.put(tag.getTitle(), tag.getSpectrum());
                    if (subSpectraMap.size() >= spectraSizeLimit) {
                        break;
                    }
                }
            }
            if (subSpectraMap.size() < spectraSizeLimit) {
                for (ConfidentTagSorter tag : confidentSpectraSet) {
                    if (tag.getValue() < 1) {
                        continue;
                    }
                    subSpectraMap.put(tag.getTitle(), tag.getSpectrum());
                    if (subSpectraMap.size() >= spectraSizeLimit) {
                        break;
                    }
                }
            }

            if (subSpectraMap.size() < spectraSizeLimit) {
                String[] titiles = subMsFileHandler.getSpectrumTitles(msFileNameWithoutExtension);
                for (String str : titiles) {
                    if (!subSpectraMap.containsKey(str)) {
                        subSpectraMap.put(str, subMsFileHandler.getSpectrum(msFileNameWithoutExtension, str));
                    }
                    if (subSpectraMap.size() >= spectraSizeLimit) {
                        break;
                    }
                }
//                  subMsFileHandler.getSpectrum(msFileNameWithoutExtension, updatedName)
            }

        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return subSpectraMap;
    }

    private File initSubFastaFile(File tempSubFastaFile, File fastaFile, Set<String> sequences) {

        if (tempSubFastaFile.exists()) {
            tempSubFastaFile.delete();
        }
        SpectraUtilities.createSubFastaFile(fastaFile, tempSubFastaFile, sequences);
        return tempSubFastaFile;

    }

}
