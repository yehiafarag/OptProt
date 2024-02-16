package no.uib.probe.optprot.search;

import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.io.identification.IdfileReader;
import com.compomics.util.experiment.io.identification.IdfileReaderFactory;
import com.compomics.util.experiment.io.mass_spectrometry.MsFileHandler;
import com.compomics.util.experiment.personalization.ExperimentObject;
import com.compomics.util.io.IoUtil;
import com.compomics.util.parameters.UtilitiesUserParameters;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.parameters.searchgui.OutputParameters;
import com.compomics.util.parameters.tools.ProcessingParameters;
import eu.isas.searchgui.SearchHandler;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import no.uib.probe.optprot.model.Configurations;
import no.uib.probe.optprot.model.OptimisedSearchParameters;
import no.uib.probe.optprot.model.SearchOptimizerParameters;
import no.uib.probe.optprot.util.OptProtWaitingHandler;
import no.uib.probe.optprot.util.SpectraFileUtilities;
import org.xmlpull.v1.XmlPullParserException;

/**
 * This class responsible for performing sub-searches in order to optimize
 * different search parameters
 *
 * @author Yehia Mokhtar Farag
 */
public class SearchOptimizerHandler {

    /**
     * The identification settings file.
     */
    private File identificationParametersFile;
    private File optimizedIdentificationParametersFile;
    private File subMsFile;
    private ArrayList<File> msFileInList;
    private File subFastaFile;

    public void setSubMsFile(File msFile) {
        this.subMsFile = msFile;
        this.msFileInList = new ArrayList<>();
        this.msFileInList.add(msFile);
    }

    public void setSubFastaFile(File subFastaFile) {
        this.subFastaFile = subFastaFile;
    }
    private File resultsOutput;
    /**
     * The processing preferences.
     */
    private ProcessingParameters processingParameters;
    private UtilitiesUserParameters userParameters;

    private ExecutorService executor;

    /**
     * The identification file reader factory of compomics utilities.
     */
    private final IdfileReaderFactory readerFactory = IdfileReaderFactory.getInstance();
    private SearchOptimizerParameters searchEngineParameters;
    private final OptimisedSearchParameters optimisedSearchParameter;

    public SearchOptimizerHandler() {
        this.initOutputFolder();
        this.initParameters();
        this.optimisedSearchParameter = new OptimisedSearchParameters();

    }

    public void executeParametersOptimization(SearchOptimizerParameters searchEngineParameters) {
        try {
            executor = Executors.newFixedThreadPool(5);
            this.searchEngineParameters = searchEngineParameters;
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(optimizedIdentificationParametersFile);
            if (searchEngineParameters.isOptimizeDigestionParameter()) {
                this.optimizeDigestionParameter();

                tempIdParam.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf(optimisedSearchParameter.getDigestionParameter()));
                IdentificationParameters.saveIdentificationParameters(tempIdParam, optimizedIdentificationParametersFile);

            }
            if (optimisedSearchParameter.getDigestionParameter().equalsIgnoreCase(DigestionParameters.CleavageParameter.enzyme.name()) && (searchEngineParameters.isOptimizeDigestionParameter() || searchEngineParameters.isOptimizeEnzymeParameter())) {
                this.optimizeEnzymeParameter();
                tempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().clear();
                tempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().add(EnzymeFactory.getInstance().getEnzyme(optimisedSearchParameter.getEnzymeName()));
                tempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(optimisedSearchParameter.getEnzymeName(), optimisedSearchParameter.getMaxMissedCleavage());
                IdentificationParameters.saveIdentificationParameters(tempIdParam, optimizedIdentificationParametersFile);
            }
            if (optimisedSearchParameter.getDigestionParameter().equalsIgnoreCase(DigestionParameters.CleavageParameter.enzyme.name()) && (searchEngineParameters.isOptimizeDigestionParameter() || searchEngineParameters.isOptimizeSpecificityParameter())) {
                this.optimizeSpecificityParameter();
                tempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(optimisedSearchParameter.getEnzymeName(), DigestionParameters.Specificity.valueOf(optimisedSearchParameter.getEnzymeSpecificity()));
                IdentificationParameters.saveIdentificationParameters(tempIdParam, optimizedIdentificationParametersFile);

            }
            if (optimisedSearchParameter.getDigestionParameter().equalsIgnoreCase(DigestionParameters.CleavageParameter.enzyme.name()) && (searchEngineParameters.isOptimizeDigestionParameter() || searchEngineParameters.isOptimizeMaxMissCleavagesParameter())) {
                this.optimizeMaxMissCleavagesParameter();
                tempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(optimisedSearchParameter.getEnzymeName(), optimisedSearchParameter.getMaxMissedCleavage());
                IdentificationParameters.saveIdentificationParameters(tempIdParam, optimizedIdentificationParametersFile);
            }
            if (searchEngineParameters.isOptimizeFragmentIonTypesParameter()) {
                this.optimizeFragmentIonTypesParameter();
                tempIdParam.getSearchParameters().setForwardIons(optimisedSearchParameter.getSelectedForwardIons());
                tempIdParam.getSearchParameters().setRewindIons(optimisedSearchParameter.getSelectedRewindIons());
                IdentificationParameters.saveIdentificationParameters(tempIdParam, optimizedIdentificationParametersFile);
            }
            if (this.searchEngineParameters.isOptimizePrecursorToleranceParameter()) {
                this.optimizePrecursorToleranceParameter();
                tempIdParam.getSearchParameters().setPrecursorAccuracy(optimisedSearchParameter.getPrecursorTolerance());
                tempIdParam.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
                IdentificationParameters.saveIdentificationParameters(tempIdParam, optimizedIdentificationParametersFile);
            }
             if (this.searchEngineParameters.isOptimizeFragmentToleranceParameter()) {
                this.optimizeFragmentToleranceParameter();
                tempIdParam.getSearchParameters().setFragmentIonAccuracy(optimisedSearchParameter.getFragmentTolerance());
                tempIdParam.getSearchParameters().setFragmentAccuracyType(SearchParameters.MassAccuracyType.DA);
                IdentificationParameters.saveIdentificationParameters(tempIdParam, optimizedIdentificationParametersFile);
            }

        } catch (IOException ex) {
            ex.printStackTrace();
        } finally {
            executor.shutdown();
            while (!executor.isTerminated()) {
            }
            cleanFolder(resultsOutput);

            System.out.println("-------------------------------------------------------------------------------------------");
            System.out.println("Digestion           :\t" + optimisedSearchParameter.getDigestionParameter());
            System.out.println("Enzyme              :\t" + optimisedSearchParameter.getEnzymeName());
            System.out.println("Specificity         :\t" + optimisedSearchParameter.getEnzymeSpecificity());
            System.out.println("Max Missed Cleavages:\t" + optimisedSearchParameter.getMaxMissedCleavage());
            System.out.println("Fragment Ion Types  :\t" + optimisedSearchParameter.getSelectedForwardIons() + "-" + optimisedSearchParameter.getSelectedRewindIons());
            System.out.println("Precursor Accuracy  :\t" + optimisedSearchParameter.getPrecursorTolerance()+" ppm");
            System.out.println("Fragment Accuracy   :\t" + optimisedSearchParameter.getFragmentTolerance()+" Da");
            System.out.println("-------------------------------------------------------------------------------------------");

        }

    }

    private void optimizeDigestionParameter() throws IOException {
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        // Creates a pool of 3 threads  
        for (int i = 0; i < DigestionParameters.CleavageParameter.values().length; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(optimizedIdentificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.getCleavageParameters(i));
            final String option = DigestionParameters.CleavageParameter.getCleavageParameters(i).name();
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            File resourceFolder = new File(resultsOutput, updatedName + "_temp");
            resourceFolder.mkdir();
            File resultOutput = new File(resultsOutput, updatedName);
            resultOutput.mkdir();
            Future f = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, "Digestion_Parameters", option, tempIdParam));
            });
            while (!f.isDone()) {
            }

        }
        System.out.println("------------->> I work is done " + resultsMap);
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;
            }
        }
        optimisedSearchParameter.setDigestionParameter(selectedOption);

    }

    private void optimizeEnzymeParameter() throws IOException {
        // Processing        
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        //optimise enzyme   
        for (Enzyme enzyme : EnzymeFactory.getInstance().getEnzymes()) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(optimizedIdentificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
            tempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(enzyme);
            final String option = enzyme.getName();
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            File resourceFolder = new File(resultsOutput, updatedName + "_temp");
            resourceFolder.mkdir();
            File resultOutput = new File(resultsOutput, updatedName);
            resultOutput.mkdir();
            Future f = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, "Digestion_Parameters", option, tempIdParam));
            });
            while (!f.isDone()) {
            }

        }
        System.out.println("------------->> II work is done " + resultsMap);
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;
            }
        }
        optimisedSearchParameter.setEnzymeName(selectedOption);

    }

    private void optimizeSpecificityParameter() throws IOException {
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        for (int i = 0; i < DigestionParameters.Specificity.values().length; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(optimizedIdentificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(optimisedSearchParameter.getEnzymeName(), DigestionParameters.Specificity.getSpecificity(i));
            final String option = DigestionParameters.Specificity.getSpecificity(i).name();
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            File resourceFolder = new File(resultsOutput, updatedName + "_temp");
            resourceFolder.mkdir();
            File resultOutput = new File(resultsOutput, updatedName);
            resultOutput.mkdir();
            Future f = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, "Digestion_Parameters", option, tempIdParam));
            });
            while (!f.isDone()) {
            }
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;
            }
        }
        System.out.println("------------->> III work is done " + resultsMap);
        optimisedSearchParameter.setEnzymeSpecificity(selectedOption);

    }

    private void optimizeMaxMissCleavagesParameter() throws IOException {
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        for (int i = 0; i < 7; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(optimizedIdentificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(optimisedSearchParameter.getEnzymeName(), i);
            final String option = i + "_missed_cleavages";
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            File resourceFolder = new File(resultsOutput, updatedName + "_temp");
            resourceFolder.mkdir();
            File resultOutput = new File(resultsOutput, updatedName);
            resultOutput.mkdir();
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, "Digestion_Parameters", option, tempIdParam));
            });
            while (!future.isDone()) {
            }
//            if (resultsMap.get(option) > idRate) {
//
//            } else if (resultsMap.get(option) == idRate) {
//                if (cont) {
//                    cont = false;
//                } else {
//                    break;
//                }
//            } else if (resultsMap.get(option) < idRate) {
//                break;
//            }
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;

            }
        }
        int selectedNum = Integer.parseInt(selectedOption.split("_")[0]);
        optimisedSearchParameter.setMaxMissedCleavage(selectedNum);
        System.out.println("------------->> IV work is done " + resultsMap + "  " + selectedNum);

    }

    private void optimizeFragmentIonTypesParameter() throws IOException {
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(optimizedIdentificationParametersFile);
        SearchParameters tempSearchParameters = tempIdParam.getSearchParameters();
        ArrayList<Integer> selectedForwardIons = tempSearchParameters.getForwardIons();
        String[] forwardIons = new String[]{"a", "b", "c"};
        String[] rewindIons = new String[]{"x", "y", "z"};

        for (String forwardIon : forwardIons) {
            selectedForwardIons.clear();
            Integer forwardIonType = PeptideFragmentIon.getIonType(forwardIon);
            selectedForwardIons.add(forwardIonType);
            for (String rewindIon : rewindIons) {
                Integer rewindIonType = PeptideFragmentIon.getIonType(rewindIon);
                ArrayList<Integer> selectedRewindIons = new ArrayList<>();
                selectedRewindIons.add(rewindIonType);
                tempSearchParameters.setRewindIons(selectedRewindIons);
                String option = selectedForwardIons + "-" + selectedRewindIons;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
                File resourceFolder = new File(resultsOutput, updatedName + "_temp");
                resourceFolder.mkdir();
                File resultOutput = new File(resultsOutput, updatedName);
                resultOutput.mkdir();
                Future f = executor.submit(() -> {
                    resultsMap.put(option, excuteXTandomSearches(updatedName, "Digestion_Parameters", option, tempIdParam));
                });
                while (!f.isDone()) {
                }
            }
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;
            }
        }
        System.out.println("------------->> V work is done " + resultsMap);

        optimisedSearchParameter.setSelectedForwardIons(selectedForwardIons);
        optimisedSearchParameter.setSelectedRewindIons(selectedForwardIons);

    }

    private void optimizePrecursorToleranceParameter() throws IOException {
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        for (double i = 5; i < 30;) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(optimizedIdentificationParametersFile);
            tempIdParam.getSearchParameters().setPrecursorAccuracy(i);
            tempIdParam.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
            final String option = i + "_precursor_accuracy";
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            File resourceFolder = new File(resultsOutput, updatedName + "_temp");
            resourceFolder.mkdir();
            File resultOutput = new File(resultsOutput, updatedName);
            resultOutput.mkdir();
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, "precursor_accuracy", option, tempIdParam));
            });
            while (!future.isDone()) {
            }
            i += 5;
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;

            }
        }
        double selectedNum = Double.parseDouble(selectedOption.split("_")[0]);
        optimisedSearchParameter.setPrecursorTolerance(selectedNum);
        System.out.println("------------->> IV work is done " + resultsMap + "  " + selectedNum);

    }
      private void optimizeFragmentToleranceParameter() throws IOException {
        Map<String, Double> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        double[]values=new double[]{0.02,0.05,0.2,0.5};
        for (double i :values) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(optimizedIdentificationParametersFile);
            tempIdParam.getSearchParameters().setFragmentIonAccuracy(i);
            tempIdParam.getSearchParameters().setFragmentAccuracyType(SearchParameters.MassAccuracyType.DA);
            final String option = i + "_fragment_accuracy";
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            File resourceFolder = new File(resultsOutput, updatedName + "_temp");
            resourceFolder.mkdir();
            File resultOutput = new File(resultsOutput, updatedName);
            resultOutput.mkdir();
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, "fragment_accuracy", option, tempIdParam));
            });
            while (!future.isDone()) {
            }
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;

            }
        }
        double selectedNum = Double.parseDouble(selectedOption.split("_")[0]);
        optimisedSearchParameter.setFragmentTolerance(selectedNum);
        System.out.println("------------->> VII work is done " + resultsMap + "  " + selectedNum);

    }

    private synchronized Double excuteXTandomSearches(String defaultOutputFileName, String paramName, String paramOption, IdentificationParameters tempIdParam) {
        try {
            SearchParameters searchParameters = tempIdParam.getSearchParameters();
            OptProtWaitingHandler waitingHandlerCLIImpl = new OptProtWaitingHandler();
            File resultOutput = new File(resultsOutput, defaultOutputFileName);
            File configFolder = new File(resultOutput, defaultOutputFileName + "_temp");
            SearchHandler searchHandler = new SearchHandler(tempIdParam, resultOutput, configFolder, defaultOutputFileName, msFileInList,
                    subFastaFile, new ArrayList<>(),
                    optimizedIdentificationParametersFile,
                    searchEngineParameters.isRunOmssa(),
                    searchEngineParameters.isRunXTandem(),
                    searchEngineParameters.isRunMsgf(),
                    searchEngineParameters.isRunMsAmanda(),
                    this.searchEngineParameters.isRunMyriMatch(),
                    this.searchEngineParameters.isRunComet(),
                    this.searchEngineParameters.isRunTide(),
                    this.searchEngineParameters.isRunAndromeda(),
                    this.searchEngineParameters.isRunMetaMorpheus(),
                    this.searchEngineParameters.isRunSage(),
                    this.searchEngineParameters.isRunNovor(),
                    this.searchEngineParameters.isRunDirecTag(),
                    this.searchEngineParameters.getOmssaFolder(),
                    this.searchEngineParameters.getxTandemFolder(),
                    this.searchEngineParameters.getMsgfFolder(),
                    this.searchEngineParameters.getMsAmandaFolder(),
                    this.searchEngineParameters.getMyriMatchFolder(),
                    this.searchEngineParameters.getCometFolder(),
                    this.searchEngineParameters.getTideFolder(),
                    this.searchEngineParameters.getTideIndexLocation(),
                    this.searchEngineParameters.getAndromedaFolder(),
                    this.searchEngineParameters.getMetaMorpheusFolder(),
                    this.searchEngineParameters.getSageFolder(),
                    this.searchEngineParameters.getNovorFolder(),
                    this.searchEngineParameters.getDirecTagFolder(),
                    this.searchEngineParameters.getMakeblastdbFolder(),
                    processingParameters
            );
            searchHandler.startSearch(waitingHandlerCLIImpl);
            ArrayList<File> xTandFiles = searchHandler.getXTandemFiles(resultOutput, IoUtil.removeExtension(subMsFile.getName()));
            IdfileReader idReader = readerFactory.getFileReader(xTandFiles.get(0));

//            String cometFileName = searchHandler.getCometFileName(IoUtil.removeExtension(sampleMgf.getName()));
//            IdfileReader idReader = readerFactory.getFileReader(new File(searchHandler.getResultsFolder(), cometFileName));
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(subMsFile, waitingHandlerCLIImpl);
            ArrayList<SpectrumMatch> matches = idReader.getAllSpectrumMatches(msFileHandler, waitingHandlerCLIImpl, searchParameters);
            double leng = msFileHandler.getSpectrumTitles(IoUtil.removeExtension(subMsFile.getName())).length;
            System.out.println("total Param : " + paramName + " - option " + paramOption + " total (%) " + matches.size() + " (" + ((matches.size() * 100.0) / leng) + "%)");
            return ((matches.size() * 100.0) / leng);
        } catch (IOException | InterruptedException | SQLException | ClassNotFoundException | JAXBException | XMLStreamException | XmlPullParserException ex) {
            ex.printStackTrace();
        }
        return -1.0;
    }

    public synchronized Set<String> excuteNovorSearches() {

        try {
            String spectraFileName = IoUtil.removeExtension(subMsFile.getName());
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + spectraFileName;
            File configFolder = new File(resultsOutput, updatedName + "_temp");
            configFolder.mkdir();
            File resultOutput = new File(resultsOutput, updatedName);
            resultOutput.mkdir();
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            SearchParameters searchParameters = tempIdParam.getSearchParameters();
            String error = SearchHandler.loadModifications(searchParameters);
            if (error != null) {
                System.out.println(error);
            }
            OptProtWaitingHandler waitingHandlerCLIImpl = new OptProtWaitingHandler();
            searchEngineParameters.setRunXTandem(false);
            searchEngineParameters.setRunNovor(true);
            SearchHandler searchHandler = new SearchHandler(tempIdParam, resultOutput, configFolder, spectraFileName, msFileInList,
                    subFastaFile, new ArrayList<>(),
                    identificationParametersFile,
                    searchEngineParameters.isRunOmssa(),
                    searchEngineParameters.isRunXTandem(),
                    searchEngineParameters.isRunMsgf(),
                    searchEngineParameters.isRunMsAmanda(),
                    this.searchEngineParameters.isRunMyriMatch(),
                    this.searchEngineParameters.isRunComet(),
                    this.searchEngineParameters.isRunTide(),
                    this.searchEngineParameters.isRunAndromeda(),
                    this.searchEngineParameters.isRunMetaMorpheus(),
                    this.searchEngineParameters.isRunSage(),
                    this.searchEngineParameters.isRunNovor(),
                    this.searchEngineParameters.isRunDirecTag(),
                    this.searchEngineParameters.getOmssaFolder(),
                    this.searchEngineParameters.getxTandemFolder(),
                    this.searchEngineParameters.getMsgfFolder(),
                    this.searchEngineParameters.getMsAmandaFolder(),
                    this.searchEngineParameters.getMyriMatchFolder(),
                    this.searchEngineParameters.getCometFolder(),
                    this.searchEngineParameters.getTideFolder(),
                    this.searchEngineParameters.getTideIndexLocation(),
                    this.searchEngineParameters.getAndromedaFolder(),
                    this.searchEngineParameters.getMetaMorpheusFolder(),
                    this.searchEngineParameters.getSageFolder(),
                    this.searchEngineParameters.getNovorFolder(),
                    this.searchEngineParameters.getDirecTagFolder(),
                    this.searchEngineParameters.getMakeblastdbFolder(),
                    processingParameters
            );

            searchHandler.startSearch(waitingHandlerCLIImpl);
            File resultsFile = searchHandler.getResultsFolder();
            File NovorFile = new File(resultsFile, spectraFileName + ".novor.csv");
            Set<String> sequences = SpectraFileUtilities.getSequences(NovorFile);
            return sequences;
        } catch (IOException | InterruptedException ex) {
            ex.printStackTrace();
        }
        return new HashSet<>();
    }

    private void initOutputFolder() {
        this.resultsOutput = new File(Configurations.OUTPUT_FOLDER_PATH);
        if (resultsOutput.exists()) {
            cleanFolder(resultsOutput);
        }
        resultsOutput.mkdir();
    }

    private void cleanFolder(File folder) {
        for (File file : folder.listFiles()) {
            if (file.isDirectory()) {
                for (File f : file.listFiles()) {
                    f.delete();
                }
                file.delete();
            } else {
                file.delete();
            }

        }
        folder.delete();

    }

    private void initParameters() {

        SearchHandler.setCloseProcessWhenDone(false);
        // set the processing preferences
        this.processingParameters = new ProcessingParameters();
        processingParameters.setnThreads(Runtime.getRuntime().availableProcessors());
        // Processing
        processingParameters.setnThreads(15);

        this.userParameters = UtilitiesUserParameters.loadUserParameters();
        userParameters.setGzip(false);
        userParameters.setSearchGuiOutputParameters(OutputParameters.no_zip);
        userParameters.setRenameXTandemFile(true);
        UtilitiesUserParameters.saveUserParameters(userParameters);

    }

    public void setIdentificationParametersFile(File identificationParametersFile) {
        try {
            this.identificationParametersFile = identificationParametersFile;
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);

            optimizedIdentificationParametersFile = new File(identificationParametersFile.getParent(), "optimized_" + identificationParametersFile.getName());
            optimizedIdentificationParametersFile.createNewFile();
            IdentificationParameters.saveIdentificationParameters(tempIdParam, optimizedIdentificationParametersFile);

            optimisedSearchParameter.setDigestionParameter(tempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name());
            if (optimisedSearchParameter.getDigestionParameter().equalsIgnoreCase(DigestionParameters.CleavageParameter.enzyme.name())) {
                optimisedSearchParameter.setEnzymeName(tempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName());
                optimisedSearchParameter.setEnzymeSpecificity(tempIdParam.getSearchParameters().getDigestionParameters().getSpecificity(optimisedSearchParameter.getEnzymeName()).name());
                optimisedSearchParameter.setMaxMissedCleavage(tempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(optimisedSearchParameter.getEnzymeName()));
            }
            optimisedSearchParameter.setSelectedForwardIons(tempIdParam.getSearchParameters().getForwardIons());
            optimisedSearchParameter.setSelectedRewindIons(tempIdParam.getSearchParameters().getRewindIons());
            optimisedSearchParameter.setPrecursorTolerance(tempIdParam.getSearchParameters().getPrecursorAccuracy());
            optimisedSearchParameter.setFragmentTolerance(tempIdParam.getSearchParameters().getFragmentIonAccuracyInDaltons());
        } catch (IOException ex) {
            Logger.getLogger(SearchOptimizerHandler.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private IdentificationParameters initSearchParameter() {

        return null;
    }

}
