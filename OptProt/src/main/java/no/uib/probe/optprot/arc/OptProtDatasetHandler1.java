/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.arc;

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
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
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
public class OptProtDatasetHandler1 {

    private final SearchInputSetting searchInputSetting = new SearchInputSetting();
    private int startIndex = 0;
    private File subFastaFile = null;

    private File subMsFile = null;
    private int counter = 0;
    private double acceptedTagEvalue;
    private boolean smallDataset;

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
        smallDataset = spectrumTitles.length < 15000;
        final Map<String, Spectrum> fullSpectrumMap = new LinkedHashMap<>();
        final Set<ConfidentTagSorter> fullConfidentSpectraSet = new TreeSet<>();

        final double coverageRatio;
        final int initialIndex;
        if (spectrumTitles.length < 30000) {
            initialIndex = 0;
            coverageRatio = 1;
        } else {
//            initialIndex = (int) Math.round((double) spectrumTitles.length * 0.25);
//            coverageRatio = 0.5;
            initialIndex = 0;
            coverageRatio = 1;
        }
        startIndex = initialIndex;

        int coverageSize;
        int step;
        int lastIndex;
        int extendor;
        boolean update;

        if (wholeDataTest || smallDataset) {
            coverageSize = spectrumTitles.length;
            step = 1;
            lastIndex = startIndex + coverageSize;
            extendor = 0;
            update = false;
        } else {
            double startRatio = 0.2;
            int finalSpectraNumberInGeneratedFile = Math.max((int) (Configurations.REFINED_MS_SIZE * 1.5), (int) (spectrumTitles.length * startRatio));
            coverageSize = (int) Math.round((double) spectrumTitles.length * coverageRatio);
            step = coverageSize / (finalSpectraNumberInGeneratedFile);
            lastIndex = startIndex + coverageSize;
            extendor = 0;
            update = false;
        }

//        if (!wholeDataTest) {
        while (true) {
            try {
                final IdentificationParameters identificationParameters = IdentificationParameters.getIdentificationParameters(new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE));
                if (!subMsFile.exists()) {
                    update = true;
                    System.out.println("sub ms not exist");
//                    subMsFile = new File(subDataFolder, Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_" + counter + "_" + msFile.getName());
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

                    while (true) {
//                            System.out.println("----------------------------------------------------------------------------------------------------------------->>>" + counter + " round " + subDataset.getOreginalDatasize() + "  final size " + finalSpectraNumberInGeneratedFile + "  start index " + startIndex);
                        //generate subset of spectra 
                        System.out.println("start index updated " + startIndex + "   " + lastIndex + "  " + step);
                        Map<String, Spectrum> spectraMap = substractSpectraFirstLevelDataFiltering(msFile, msFileHandler, startIndex, step, lastIndex);
                        System.out.println("spectraMap " + spectraMap.size());
                        fullSpectrumMap.putAll(spectraMap);

                        //generate submsFile
                        File destinationFile = new File(Configurations.GET_OUTPUT_FOLDER_PATH(), Configurations.DEFAULT_RESULT_NAME + "_" + counter + fileNameWithoutExtension + ".mgf");
                        if (destinationFile.exists()) {
                            destinationFile.delete();
                        }
                        destinationFile.createNewFile();
                        destinationFile = generateMsSubFile(spectraMap, destinationFile);
                        final String subfileNameWithoutExtension = IoUtil.removeExtension(destinationFile.getName());
                        //run direct tag and get confident tags                
                        Set<ConfidentTagSorter> confidentSpectraSet = getSpectraWithConfidentTagSecondStageFiltering(destinationFile, fastaFile, identificationParameters, new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE), subfileNameWithoutExtension);
                        fullConfidentSpectraSet.addAll(confidentSpectraSet);
                        System.out.println("full size " + fullConfidentSpectraSet.size() + "  " + confidentSpectraSet.size());
                        if (counter == 2) {
                            counter = 0;
                            subFastaFile.delete();
                            subMsFile.delete();
                            subFastaFile = null;
                            subMsFile = null;
                            startIndex = 0;
                            return generateOptProtDataset(msFile, fastaFile, searchEngineToOptimise, subDataFolder, identificationParametersFile, true);
                        }
                        if (!wholeDataTest && fullConfidentSpectraSet.size() < Configurations.REFINED_MS_SIZE + extendor) {
                            startIndex++;
                            if (startIndex == initialIndex + step || counter == 8) {
                                System.out.println("the file is totally covered and data is so bad " + confidentSpectraSet.size() + "/" + Configurations.REFINED_MS_SIZE);
                                break;
                            }
                            counter++;
                        } else {
                            break;
                        }
                    }
                    int refinedIndex = 0;
                    final Map<String, Spectrum> refinedSpectrumMap = new LinkedHashMap<>();
                    //first refine to run 
                    for (ConfidentTagSorter tag : fullConfidentSpectraSet) {
                        refinedSpectrumMap.put(tag.getTitle(), fullSpectrumMap.get(tag.getTitle()));
                        if (!wholeDataTest && refinedIndex > Configurations.EXTRACT_MAX_MS_SIZE + (extendor * 0.5)) {//
                            break;
                        }
                        refinedIndex++;
                    }
                    System.out.println("refined spectrum map size " + refinedSpectrumMap.size());
                    //create stabkle subMs file
                    subMsFile = generateMsSubFile(refinedSpectrumMap, subMsFile);
                    subFastaFile.delete();
                }
                //create stabkle subfasta file
                if (!subFastaFile.exists() || update) {
                    update = true;
                    System.out.println("sub fasta not exist");
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
                    subFastaFile = initSubFastaFile(subDataFolder, subFastaFile, fastaFile, sequences);
                    long end4th = System.currentTimeMillis();
                    total = (end4th - start4) / 1000.0;
                    System.out.println("process IV ( Generated Fasta) in seconds: " + total);
                    long end = System.currentTimeMillis();
                    total = (end - start1) / 1000.0;
                    System.out.println("Total Elapsed Time for initInputSubSetFiles in seconds: " + total);
//                    MainUtilities.deleteFolder(NovorFile.getParentFile());
                }
                if (update && !wholeDataTest) {
                    //run initial identification with user selected SE
                    long s1 = System.currentTimeMillis();
                    searchInputSetting.setSelectedSearchEngine(standeredReferenceSearchEngine);
                    final String option = "init_input_files" + standeredReferenceSearchEngine;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + fileNameWithoutExtension;
//                    XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//                    xtandemParameters.setProteinQuickAcetyl(false);
//                    xtandemParameters.setQuickPyrolidone(false);
//                    xtandemParameters.setStpBias(false);
//                    xtandemParameters.setRefine(false);
//                    xtandemParameters.setOutputResults("valid");//"valid"
//                    xtandemParameters.setMaxEValue(0.01);

                    File resultsFolder = SearchExecuter.executeSearch(updatedName, searchInputSetting, subMsFile, subFastaFile, identificationParameters, new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE));
                    List<SpectrumMatch> validatedMaches = SpectraUtilities.getValidatedIdentificationResults(resultsFolder, subMsFile, standeredReferenceSearchEngine, identificationParameters);
                    long e = System.currentTimeMillis();
                    double t = (e - s1) / 1000.0;
                    subDataset.setProcessDelay(t);
                    String subfileNameWithoutExtension = IoUtil.removeExtension(subMsFile.getName());
                    MsFileHandler subMsFileHandler = new MsFileHandler();
                    subMsFileHandler.register(subMsFile, new OptProtWaitingHandler());
                    MainUtilities.deleteFolder(resultsFolder);
                    int total = subMsFileHandler.getSpectrumTitles(subfileNameWithoutExtension).length;
//                       subDataset.setTotalSpectraNumber(total);
                    MainUtilities.cleanOutputFolder();
                    if (validatedMaches.size() >= (subMsFileHandler.getSpectrumTitles(subfileNameWithoutExtension).length * Configurations.ACCEPTED_REFERENCE_ID_RATIO) || total >= Configurations.EXTRACT_MAX_MS_SIZE) {
                        subDataset.setDefaultSettingIdentificationNum(validatedMaches.size());
                        subMsFileHandler.close();
                        subMsFileHandler.getCmsFilePaths().clear();
                        break;
                    } else {
                        subMsFileHandler.close();
                        subMsFile.delete();
                        subFastaFile.delete();
                        File cms = new File(subDataFolder, subMsFile.getName().replace(".mgf", ".cms"));
                        cms.delete();
                        counter++;
                        startIndex++;
                        extendor += 500;
                    }
                } else {
                    break;
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
            File resultsFolder = SearchExecuter.executeSearch(updatedName, searchInputSetting, subMsFile, subFastaFile, identificationParameters, new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE));

            List<SpectrumMatch> validatedMaches = SpectraUtilities.getValidatedIdentificationResults(resultsFolder, subMsFile, searchEngineToOptimise, identificationParameters);
            System.out.println("refrence default " + validatedMaches.size());
            if (validatedMaches.isEmpty()) {
                System.out.println("Error in the system please restart!");
                System.exit(0);
            }
            subDataset.setDefaultSettingIdentificationNum(validatedMaches.size());//                }
//            subDataset.setActiveScoreModel(SpectraUtilities.getValidatedIdentificationReferenceData(resultsFolder, subMsFile, searchEngineToOptimise, identificationParameters, subDataset.getDefaultSettingIdentificationNum()));
 //@replace later           subDataset.setValidatedIdRefrenceData(SpectraUtilities.getValidatedIdentificationReferenceData(validatedMaches, searchEngineToOptimise));
            String subfileNameWithoutExtension = IoUtil.removeExtension(subMsFile.getName());
            MsFileHandler subMsFileHandler = new MsFileHandler();
            subMsFileHandler.register(subMsFile, new OptProtWaitingHandler());
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

//        }
        //run nonEvalue score search engine
//        if (searchEngineToOptimise.getIndex() == Advocate.myriMatch.getIndex()) {
//            try {
//
//                //run initial identification with user selected SE
//                final IdentificationParameters identificationParameters = IdentificationParameters.getIdentificationParameters(new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE));
//                searchInputSetting.setSelectedSearchEngine(searchEngineToOptimise);
//                final String option = "reference_run_default_" + searchEngineToOptimise;
//                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + fileNameWithoutExtension;
//                File resultsFolder = SearchExecuter.executeSearch(updatedName, searchInputSetting, subMsFile, subFastaFile, identificationParameters, new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE));
//                Configurations.VALIDATED_ID_REF_DATA = SpectraUtilities.getValidatedIdentificationReferenceData(resultsFolder, subMsFile, searchEngineToOptimise, identificationParameters, subDataset.getDefaultSettingIdentificationNum());
//                System.out.println("refrence default threshold " + Configurations.VALIDATED_ID_REF_DATA.length + "  " + resultsFolder.getAbsolutePath());
//                MainUtilities.cleanOutputFolder();
//            } catch (IOException ex) {
//                ex.printStackTrace();
//            }
//
//        }
        //re run with user input search to get the reference number
//        try {
//            final String option = "reference_run_user_" + searchEngineToOptimise;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + fileNameWithoutExtension;
//            final IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
//            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
//            tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
//            tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
//            System.out.println("at searchEngineToOptimise " + searchEngineToOptimise);
//            File resultsFolder = SearchExecuter.executeSearch(updatedName, searchInputSetting, subMsFile, subFastaFile, tempIdParam, identificationParametersFile);
//            List<SpectrumMatch> validatedMaches = SpectraUtilities.getValidatedIdentificationResults(resultsFolder, subMsFile, searchEngineToOptimise, tempIdParam);
//            double ratio = (validatedMaches.size() * 100.0 / subDataset.getTotalSpectraNumber());
//            subDataset.setHighResolutionMassSpectrometers(ratio >= 20);
//            String subfileNameWithoutExtension = IoUtil.removeExtension(subMsFile.getName());
//            MsFileHandler subMsFileHandler = new MsFileHandler();
//            subMsFileHandler.register(subMsFile, new OptProtWaitingHandler());
//            subDataset.setUserReferenceIdentificationNum(validatedMaches.size());
//            subDataset.setTotalSpectraNumber(subMsFileHandler.getSpectrumTitles(subfileNameWithoutExtension).length);
//            System.out.println("refrence user " + validatedMaches.size());
//            MainUtilities.cleanOutputFolder();
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
        //re run quality with user input search to get the reference number
//        double ratio = (subDataset.getDefaultSettingIdentificationNum() * 100.0 / subDataset.getTotalSpectraNumber());
//        subDataset.setHighResolutionMassSpectrometers(ratio >= 20);
//        try {
//            final String option = "reference_run_user_Quality_" + searchEngineToOptimise;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + fileNameWithoutExtension;
//            final IdentificationParameters identificationParameters;
////            if (subDataset.getUserReferenceIdentificationNum() > subDataset.getDefaultSettingIdentificationNum()) {
////                identificationParameters = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
////            } else {
//            identificationParameters = IdentificationParameters.getIdentificationParameters(new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE));
////            }
//            identificationParameters.getSearchParameters().getModificationParameters().clearFixedModifications();
//            identificationParameters.getSearchParameters().getModificationParameters().clearVariableModifications();
//            identificationParameters.getSearchParameters().getModificationParameters().clearRefinementModifications();
//            identificationParameters.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
//            File resultsFolder = SearchExecuter.executeSearch(updatedName, searchInputSetting, subMsFile, fastaFile, identificationParameters, identificationParametersFile);
//            ArrayList<SpectrumMatch> validatedMaches = SpectraUtilities.getValidatedIdentificationResults(resultsFolder, subMsFile, searchEngineToOptimise, identificationParameters);
//            double ratio = (validatedMaches.size() * 100.0 / subDataset.getTotalSpectraNumber());
//            subDataset.setHighResolutionMassSpectrometers(ratio >= 20);
//            MainUtilities.cleanOutputFolder();
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//        System.out.println("final ds handleing " + subDataset.getActiveIdentificationNum() + "/" + subDataset.getTotalSpectraNumber());
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

    private Set<ConfidentTagSorter> getSpectraWithConfidentTagSecondStageFiltering(File destinationFile, File fastaFile, IdentificationParameters identificationParameters, File identificationParametersFile, String msFileNameWithoutExtension) {
        Set<ConfidentTagSorter> confidentSpectraSet = new TreeSet<>();
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
                return confidentSpectraSet;
            }

            IdfileReader idReader = IdfileReaderFactory.getInstance().getFileReader(direcTagFile);

            MsFileHandler subMsFileHandler = new MsFileHandler();
            subMsFileHandler.register(destinationFile, MainUtilities.OptProt_Waiting_Handler);
            ArrayList<SpectrumMatch> matches = idReader.getAllSpectrumMatches(subMsFileHandler, MainUtilities.OptProt_Waiting_Handler, identificationParameters.getSearchParameters());

            for (SpectrumMatch sm : matches) {
                TagAssumption tag = sm.getAllTagAssumptions().toList().get(0);
                if (tag.getScore() < acceptedTagEvalue) {
//                    confidentSpectraSet.add(new ConfidentTagSorter(tag.getScore(), sm.getSpectrumTitle()));
                }
            }
            File cms = new File(destinationFile.getParent(), destinationFile.getName().replace(".mgf", ".cms"));
            cms.delete();

        } catch (IOException | SQLException | ClassNotFoundException | InterruptedException | JAXBException | XmlPullParserException | XMLStreamException ex) {
            ex.printStackTrace();
        }
        return confidentSpectraSet;
    }

    private File initSubFastaFile(File subDataFolder, File tempSubFastaFile, File fastaFile, Set<String> sequences) {

        if (tempSubFastaFile.exists()) {
            tempSubFastaFile.delete();
        }
        SpectraUtilities.createSubFastaFile(fastaFile, tempSubFastaFile, sequences);
        return tempSubFastaFile;

    }

}
