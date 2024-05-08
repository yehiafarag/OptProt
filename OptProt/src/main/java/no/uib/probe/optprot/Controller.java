package no.uib.probe.optprot;

import java.io.File;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.dataset.OptProtDatasetHandler;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.search.OptProtSearchHandler;
import no.uib.probe.optprot.util.MainUtilities;

/**
 *
 * @author yfa041
 */
public class Controller {

    private final OptProtDatasetHandler optProtDatasetHandler;

    public Controller() {
        this.optProtDatasetHandler = new OptProtDatasetHandler();
    }

    public void processDataset(String datasetId, File oreginalMsFile, File oreginalFastaFile, File identificationParametersFile, SearchInputSetting optProtSearchSettings) {
        MainUtilities.cleanOutputFolder();
        Configurations.ACTIVE_SEARCH_SETTINGS_FILE = new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE);
        SearchingSubDataset optProtDataset = optProtDatasetHandler.generateOptProtDataset(oreginalMsFile, oreginalFastaFile, optProtSearchSettings.getSelectedSearchEngine(), identificationParametersFile);
//        System.out.println("----------------------------------------------------------------<<<<subMs : " + optProtDataset.getSubMsFile().getName() + "  " + optProtDataset.getSubFastaFile().getName());
 
//            optProtDataset.setSubFastaFile(oreginalFastaFile);        
//            optProtDataset.setSubMsFile(oreginalMsFile);
        
        
        OptProtSearchHandler optProtSearchHandler = new OptProtSearchHandler();
        long start = System.currentTimeMillis();
        optProtSearchHandler.optimizeSearchEngine(optProtDataset, optProtSearchSettings);
        long end = System.currentTimeMillis();
        double total = (end - start) / 60000.0;
        System.out.println("Total Elapsed Time for optimizing the data in min: " + total);
        System.exit(0);

    }
}
