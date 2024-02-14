package no.uib.probe.optprot.search;

import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
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
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import no.uib.probe.optprot.model.Configurations;
import no.uib.probe.optprot.model.OptimisedSearchParameters;
import no.uib.probe.optprot.model.SearchEngineParameters;
import no.uib.probe.optprot.util.OptProtWaitingHandler;
import org.xmlpull.v1.XmlPullParserException;

/**
 * This class responsible for performing sub-searches in order to optimize
 * different search parameters
 *
 * @author Yehia Mokhtar Farag
 */
public class SubSearchHandler {

    /**
     * The identification settings file.
     */
    private final File identificationParametersFile;
    private final ArrayList<File> msFiles;
    private final File fastaFile;
    private File resultsOutput;
    /**
     * The processing preferences.
     */
    private final ProcessingParameters processingParameters;
    private final UtilitiesUserParameters userParameters;

    private ExecutorService executor;

    /**
     * The identification file reader factory of compomics utilities.
     */
    private final IdfileReaderFactory readerFactory = IdfileReaderFactory.getInstance();
    private File sampleMgf;
    private SearchEngineParameters searchEngineParameters;
    private final OptimisedSearchParameters optimisedSearchParameter;

    public SubSearchHandler(File identificationParametersFile, ArrayList<File> msFiles, File fastaFile) {
        this.identificationParametersFile = identificationParametersFile;
        this.msFiles = msFiles;
        this.fastaFile = fastaFile;
        // set the processing preferences
        processingParameters = new ProcessingParameters();
        processingParameters.setnThreads(Runtime.getRuntime().availableProcessors());
        // Processing
        processingParameters.setnThreads(15);
        this.sampleMgf = msFiles.get(0);
        userParameters = UtilitiesUserParameters.loadUserParameters();
        userParameters.setGzip(false);
        userParameters.setSearchGuiOutputParameters(OutputParameters.no_zip);
        userParameters.setRenameXTandemFile(true);
        UtilitiesUserParameters.saveUserParameters(userParameters);
        initOutputFolder();
        this.initParameters();
        this.optimisedSearchParameter = new OptimisedSearchParameters();

    }

    public void execute(ExperimentObject parameter) {
        if (parameter instanceof DigestionParameters) {
            this.optimiseDigestionParameters();
        }

    }

    private void optimiseDigestionParameters() {
        try {
            //first optimise digestion 

            // Processing        
            Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
            int idRate = -1;
            String selectedOption = "";
            //optimise enzyme   
            if (optimisedSearchParameter.getCleavageParameters() == null) {
                executor = Executors.newFixedThreadPool(DigestionParameters.CleavageParameter.values().length + 1); // Creates a pool of 3 threads  
                for (int i = 0; i < DigestionParameters.CleavageParameter.values().length; i++) {
                    IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
                    tempIdParam.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.getCleavageParameters(i));
                    final String option = DigestionParameters.CleavageParameter.getCleavageParameters(i).name();
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
                    File resourceFolder = new File(resultsOutput, updatedName + "_temp");
                    resourceFolder.mkdir();
                    File resultOutput = new File(resultsOutput, updatedName);
                    resultOutput.mkdir();
                    executor.submit(() -> {
                        resultsMap.put(option, excuteSearches(updatedName, "Digestion_Parameters", option, tempIdParam));
                    });
                }
                executor.shutdown();
                while (!executor.isTerminated()) {
                }
                System.out.println("------------->> I work is done " + resultsMap);

                for (String option : resultsMap.keySet()) {
                    if (resultsMap.get(option) > idRate) {
                        idRate = resultsMap.get(option);
                        selectedOption = option;
                    }
                }
                optimisedSearchParameter.setCleavageParameters(selectedOption);
                resultsMap.clear();
            }

            //next stage of optimisation is to choose enzyme
            if (optimisedSearchParameter.getCleavageParameters() != null && optimisedSearchParameter.getCleavageParameters().equalsIgnoreCase(DigestionParameters.CleavageParameter.enzyme.name()) && optimisedSearchParameter.getEnzymeName() == null) {
                executor = Executors.newFixedThreadPool(EnzymeFactory.getInstance().getEnzymes().size() + 1); // Creates a pool of 3 threads  
                for (Enzyme enzyme : EnzymeFactory.getInstance().getEnzymes()) {
                    IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
                    tempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
                    tempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(enzyme);
                    final String option = enzyme.getName();
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
                    File resourceFolder = new File(resultsOutput, updatedName + "_temp");
                    resourceFolder.mkdir();
                    File resultOutput = new File(resultsOutput, updatedName);
                    resultOutput.mkdir();
                    executor.submit(() -> {
                        resultsMap.put(option, excuteSearches(updatedName, "Digestion_Parameters", option, tempIdParam));
                    });

                }
                executor.shutdown();
                while (!executor.isTerminated()) {
                }
                System.out.println("------------->> II work is done " + resultsMap);
                idRate = -1;
                selectedOption = "";
                for (String option : resultsMap.keySet()) {
                    if (resultsMap.get(option) > idRate) {
                        idRate = resultsMap.get(option);
                        selectedOption = option;
                    }
                }
                resultsMap.clear();
                optimisedSearchParameter.setEnzymeName(selectedOption);
                System.out.println("Output results is " + selectedOption + " with idRate : " + idRate);
            }
            if (optimisedSearchParameter.getEnzymeName() != null && optimisedSearchParameter.getEnzymeSpecificity() == null) {
                System.out.println("Search for specifty type start ----------------------------->");
                executor = Executors.newFixedThreadPool(DigestionParameters.Specificity.values().length + 1);
                for (int i = 0; i < DigestionParameters.Specificity.values().length; i++) {
                    IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
                    tempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(optimisedSearchParameter.getEnzymeName(), DigestionParameters.Specificity.getSpecificity(i));
                    final String option = DigestionParameters.Specificity.getSpecificity(i).name();
                    System.out.println("at specifty is " + option);
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
                    File resourceFolder = new File(resultsOutput, updatedName + "_temp");
                    resourceFolder.mkdir();
                    File resultOutput = new File(resultsOutput, updatedName);
                    resultOutput.mkdir();
                    executor.submit(() -> {
                        resultsMap.put(option, excuteSearches(updatedName, "Digestion_Parameters", option, tempIdParam));
                    });
                }

                executor.shutdown();
                while (!executor.isTerminated()) {
                }
                System.out.println("------------->> III work is done " + resultsMap);
                idRate = -1;
                selectedOption = "";

                for (String option : resultsMap.keySet()) {
                    if (resultsMap.get(option) > idRate) {
                        idRate = resultsMap.get(option);
                        selectedOption = option;
                    }
                }
                resultsMap.clear();
                System.out.println("Output results is " + selectedOption + " with idRate : " + idRate);
                optimisedSearchParameter.setEnzymeSpecificity(selectedOption);
            }
            if (optimisedSearchParameter.getMaxMissedCleavage() == -1) {
                executor = Executors.newFixedThreadPool(5 + 1);
                idRate = -1;
                selectedOption = "";
                int selectedNum = -1;
                boolean cont = true;
                for (int i = 0; i < 6; i++) {
                    IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
                    tempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(optimisedSearchParameter.getEnzymeName(), i);
                    final String option = i + "_missed_cleavages";
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
                    File resourceFolder = new File(resultsOutput, updatedName + "_temp");
                    resourceFolder.mkdir();
                    File resultOutput = new File(resultsOutput, updatedName);
                    resultOutput.mkdir();
                    Future future = executor.submit(() -> {
                        resultsMap.put(option, excuteSearches(updatedName, "Digestion_Parameters", option, tempIdParam));
                    });
                    while (!future.isDone()) {
                    }
                    if (resultsMap.get(option) > idRate) {
                        idRate = resultsMap.get(option);
                        selectedOption = option;
                        selectedNum = i;
                    } else if (resultsMap.get(option) == idRate) {
                        if (cont) {
                            cont = false;
                        } else {
                            break;
                        }
                    } else if (resultsMap.get(option) < idRate) {
                        break;
                    }
                }
                executor.shutdown();
                while (!executor.isTerminated()) {
                }
//                System.out.println("result map "+resultsMap);
//                for (String option : resultsMap.keySet()) {
//                    if (resultsMap.get(option) > idRate) {
//                        idRate = resultsMap.get(option);
//                        selectedOption = option;
//                    }
//                }
                System.out.println("------------->> IV work is done " + resultsMap);
                resultsMap.clear();
                optimisedSearchParameter.setMaxMissedCleavage(selectedNum);
                System.out.println("Output results is " + selectedOption + " with idRate : " + idRate);

            }
        } catch (Exception ex) {
        } finally {
            cleanFolder(resultsOutput);
        }

    }

    private synchronized int excuteSearches(String defaultOutputFileName, String paramName, String paramOption, IdentificationParameters tempIdParam) {
        try {
            SearchParameters searchParameters = tempIdParam.getSearchParameters();
            String error = SearchHandler.loadModifications(searchParameters);
            if (error != null) {
                System.out.println(error);
            }

            OptProtWaitingHandler waitingHandlerCLIImpl = new OptProtWaitingHandler();
            File resultOutput = new File(resultsOutput, defaultOutputFileName);
            File configFolder = new File(resultOutput, defaultOutputFileName + "_temp");
            SearchHandler searchHandler = new SearchHandler(tempIdParam, resultOutput, configFolder, defaultOutputFileName, msFiles,
                    fastaFile, new ArrayList<>(),
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
            ArrayList<File> xTandFiles = searchHandler.getXTandemFiles(resultOutput, IoUtil.removeExtension(sampleMgf.getName()));
            IdfileReader idReader = readerFactory.getFileReader(xTandFiles.get(0));
//            System.out.println("xtandom file " + xTandFiles.get(0).getAbsolutePath());
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(sampleMgf, waitingHandlerCLIImpl);
            ArrayList<SpectrumMatch> matches = idReader.getAllSpectrumMatches(msFileHandler, waitingHandlerCLIImpl, searchParameters);
            int leng=msFileHandler.getSpectrumTitles(IoUtil.removeExtension(sampleMgf.getName())).length;
            System.out.println("------------------------------------------------------------------------------------------------------------------------------ total Param : " + paramName + " - option " + paramOption + " total (%) " +matches.size()+" ("+ ((matches.size()*100)/leng)+"%)");
            return matches.size();
        } catch (IOException | InterruptedException | SQLException | ClassNotFoundException | JAXBException | XMLStreamException | XmlPullParserException ex) {
            ex.printStackTrace();
        }
        return -1;
    }

    private void initOutputFolder() {
        this.resultsOutput = new File(Configurations.OUTPUT_FOLDER_PATH);
        if (!resultsOutput.exists()) {
            resultsOutput.mkdir();
        } else {
            cleanFolder(resultsOutput);
        }

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

    }

    private void initParameters() {
        this.searchEngineParameters = new SearchEngineParameters();
        searchEngineParameters.setxTandemFolder(new File("D:\\Apps\\searchgui\\resources\\XTandem\\windows\\windows_64bit"));
        searchEngineParameters.setRunXTandem(true);

    }

}
