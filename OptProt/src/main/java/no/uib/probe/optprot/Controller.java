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

/**
 *
 * @author yfa041
 */
public class Controller {

    private final OptProtDatasetHandler optProtDatasetHandler;

    public Controller() {
        this.optProtDatasetHandler = new OptProtDatasetHandler();

    }

    String v = "";

    public void processDataset(String datasetId, File oreginalMsFile, File oreginalFastaFile, File identificationParametersFile, SearchInputSetting optProtSearchSettings, boolean fullDataTest, List<String> paramOrder) {

        File subDataFolder = new File(Configurations.DATA_FOLDER + datasetId, optProtSearchSettings.getSelectedSearchEngine().getName());
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
//        System.exit(0);
        MainUtilities.cleanOutputFolder();
        SearchingSubDataset optProtDataset = optProtDatasetHandler.generateOptProtDataset(oreginalMsFile, oreginalFastaFile, optProtSearchSettings.getSelectedSearchEngine(), subDataFolder, identificationParametersFile, fullDataTest);
        optProtDataset.setSubDataFolder(subDataFolder);
        File selectedSearchSettingsFile;
//        optProtDataset.setSubFastaFile(oreginalFastaFile);
//        optProtDataset.setSubMsFile(oreginalMsFile);
        if (optProtSearchSettings.isOptimizeAllParameters()) {//&& (optProtDataset.getDefaultSettingIdentificationNum() >= optProtDataset.getUserReferenceIdentificationNum())) {
            selectedSearchSettingsFile = new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE);
            optProtDataset.setActiveIdentificationNum(optProtDataset.getDefaultSettingIdentificationNum());
        } else {
            selectedSearchSettingsFile = identificationParametersFile;
//             Configurations.ACTIVE_SEARCH_SETTINGS_FILE = new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE);            
            optProtDataset.setActiveIdentificationNum(optProtDataset.getUserReferenceIdentificationNum());
        }
        optProtDataset.setSearchSettingsFile(selectedSearchSettingsFile);

        OptProtSearchHandler optProtSearchHandler = new OptProtSearchHandler();
        long start = System.currentTimeMillis();
        File generatedFile = optProtSearchHandler.optimizeSearchEngine(optProtDataset, optProtSearchSettings, paramOrder);
        long end = System.currentTimeMillis();
        double total = (end - start) / 60000.0;

//        System.exit(0);
        if (generatedFile != null) {
            ReportExporter.printFullReport(generatedFile, optProtDataset, optProtSearchSettings.getSelectedSearchEngine(), datasetId);
        }
        System.out.println("Total Elapsed Time for optimizing the data in min: " + total);
    }
}
