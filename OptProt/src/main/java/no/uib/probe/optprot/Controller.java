package no.uib.probe.optprot;

import java.io.File;
import java.util.List;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.search.OptProtSearchHandler;
import no.uib.probe.optprot.util.MainUtilities;
import no.uib.probe.optprot.dataset.OptProtDatasetHandler;
import no.uib.probe.optprot.util.ReportExporter;
import no.uib.probe.optprot.util.SpectraUtilities;

/**
 *
 * @author yfa041
 */
public class Controller {
    
    private final OptProtDatasetHandler optProtDatasetHandler;
    
    public Controller() {
        this.optProtDatasetHandler = new OptProtDatasetHandler();
        
    }
    
    public void processDataset(String datasetId, File oreginalMsFile, File oreginalFastaFile, File identificationParametersFile, SearchInputSetting optProtSearchSettings, boolean wholeDataTest, List<String> paramOrder) {
        File subDataFolder = new File(Configurations.GET_DATA_FOLDER() + datasetId, optProtSearchSettings.getSelectedSearchEngine().getName());
        if (subDataFolder.exists()) {
//            for (File f : subDataFolder.listFiles()) {
//                System.out.println("File found  " + f.getName());
//                if (f.getName().endsWith("_optAll.par") && optProtSearchSettings.isOptimizeAllParameters()) {
//                    try {
//                        IdentificationParameters identificationParameters = IdentificationParameters.getIdentificationParameters(f);
//                        ReportExporter.printFullReport(identificationParameters, null, optProtSearchSettings.getSelectedSearchEngine());
////                        System.exit(0);
//                    } catch (IOException ex) {
//                        Logger.getLogger(Controller.class.getName()).log(Level.SEVERE, null, ex);
//                    }
//
//                }
//
//            }
        } else {
            subDataFolder.mkdir();
        }
        MainUtilities.cleanOutputFolder();
        long startDsInit = System.currentTimeMillis();
        SearchingSubDataset optProtDataset = optProtDatasetHandler.generateOptProtDataset(oreginalMsFile, oreginalFastaFile, optProtSearchSettings.getSelectedSearchEngine(), subDataFolder, identificationParametersFile, wholeDataTest);
        long endDsInit = System.currentTimeMillis();
        String totalDsTime = MainUtilities.msToTime(endDsInit - startDsInit);
        optProtDataset.setSubDataFolder(subDataFolder);
        optProtDataset.setFullDataSpectaInput(wholeDataTest);
        File selectedSearchSettingsFile;
        
        if (optProtSearchSettings.isOptimizeAllParameters()) {
            selectedSearchSettingsFile = new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE);
        } else {
            selectedSearchSettingsFile = identificationParametersFile;            
        }
        optProtDataset.setSearchSettingsFile(selectedSearchSettingsFile);
        
        if (!optProtDataset.isFullDataSpectaInput() || true) {
            
            double comparisonsThreshold = SpectraUtilities.calculateDatasetScoreThreshold((double) optProtDataset.getOreginalDatasize(), (double) optProtDataset.getTotalSpectraNumber(), (optProtDataset.getIdentificationRate() / 100.0), (double) optProtDataset.getActiveIdentificationNum());
            optProtDataset.setComparisonsThreshold(comparisonsThreshold);
        }
        
        MainUtilities.cleanOutputFolder();
        
        OptProtSearchHandler optProtSearchHandler = new OptProtSearchHandler();
        long start = System.currentTimeMillis();
        File generatedFile = optProtSearchHandler.optimizeSearchEngine(optProtDataset, optProtSearchSettings, paramOrder);
        long end = System.currentTimeMillis();
        String totalTime = MainUtilities.msToTime(end - start);
        if (generatedFile != null) {
            ReportExporter.exportFullReport(generatedFile, optProtDataset, optProtSearchSettings.getSelectedSearchEngine(), datasetId, totalTime, totalDsTime,optProtDataset.getParameterScoreMap());
            ReportExporter.printFullReport(generatedFile, optProtDataset, optProtSearchSettings.getSelectedSearchEngine(), datasetId);
        }
        System.out.println("Total Elapsed Time for Init dataset : " + totalDsTime);
        System.out.println("Total Elapsed Time for optimizing the data in : " + totalTime);
    }
}
