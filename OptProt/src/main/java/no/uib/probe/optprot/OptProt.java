package no.uib.probe.optprot;

import com.compomics.util.experiment.identification.Advocate;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import javax.swing.SwingUtilities;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.search.sage.OptProtSageParameterSettings;
import no.uib.probe.optprot.search.xtandam.OptProtXtandemParameterSettings;
import no.uib.probe.optprot.util.MainUtilities;

/**
 * This app is search settings optimization workflow that aim to optimize search
 * settings for different proteomics search engines
 *
 * @author Yehia Mokhtar Farag
 */
public class OptProt {

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            try {
                Map<Advocate, List<String>> paramOrderMap = new HashMap<>();

//                List<String> paramOrder = new ArrayList<>();
                Set<Advocate> supportedSearchEngine = new LinkedHashSet<>();

                supportedSearchEngine.add(Advocate.sage);
                paramOrderMap.put(Advocate.sage, OptProtSageParameterSettings.Get_Sage_Parameters_List());

//                supportedSearchEngine.add(Advocate.xtandem);
                paramOrderMap.put(Advocate.xtandem, OptProtXtandemParameterSettings.Get_Xtandem_Parameters_List());
//                supportedSearchEngine.add(Advocate.myriMatch);

//////   
//                paramOrder.add("MyriMatchAdvancedParameter");
//               
//                
                Set<String> datasettoTestSet = new LinkedHashSet<>();
                if (args == null || args.length == 0) {
                    datasettoTestSet.add("PXD028427");    //1
                    datasettoTestSet.add("PXD000561");    //2
                    datasettoTestSet.add("PXD001468");          //3
                    datasettoTestSet.add("PXD047036");        //4
                    datasettoTestSet.add("PXD009340");        //5
                    datasettoTestSet.add("PXD001250");        //6
////////////////            datasettoTestSet.add("PXD000815");        
////////////////            datasettoTestSet.add("PXD054727");    
                } else {
                    datasettoTestSet.addAll(Arrays.asList(args));
                    System.exit(0);
                }

//             
                boolean cleanAll = false;
                SearchInputSetting searchOpParameter = new SearchInputSetting();
                boolean all = true;
                boolean useFullFasta = false;
                boolean useOreginalInputs = true;
                searchOpParameter.setOptimizeAllParameters(all);
                searchOpParameter.setOptimizeDigestionParameter(true || all);
                searchOpParameter.setOptimizeCleavageParameter(false);
                searchOpParameter.setOptimizeEnzymeParameter(true);
                searchOpParameter.setOptimizeMaxMissCleavagesParameter(false || all);
                searchOpParameter.setOptimizeSpecificityParameter(false);
                searchOpParameter.setOptimizeFragmentIonTypesParameter(false || all);
                searchOpParameter.setOptimizePrecursorToleranceParameter(false || all);
                searchOpParameter.setOptimizeFragmentToleranceParameter(false || all);
                searchOpParameter.setOptimizePrecursorChargeParameter(false || all);
                searchOpParameter.setOptimizeIsotopsParameter(false || all);
                searchOpParameter.setOptimizeModificationParameter(false || all);
                searchOpParameter.setOptimizeSageAdvancedParameter(false||all);
                 searchOpParameter.setOptimizeXtandemAdvancedParameter(false||all);
//            searchOpParameter.setRecalibrateSpectraParameter(false);

                for (Advocate se : supportedSearchEngine) {
//                    if (se.getIndex() == Advocate.xtandem.getIndex()) {
//                        System.out.println("---------------------------------------------------------full-" + "PXD001250" + "----------------------------------------------");
//                        System.gc();
//                        MainUtilities.cleanOutputFolder();
//                        runDataset("PXD001250", cleanAll, paramOrderMap.get(se), searchOpParameter, true);
//                        continue;
//                    }

                    searchOpParameter.setSelectedSearchEngine(se);
                    for (String datasetId : datasettoTestSet) {
                        System.out.println("--------------------------------------------------------- ds " + datasetId + "----------------------------------------------");
//                        cleanAll = true;
                        MainUtilities.cleanOutputFolder();
                        runDataset(datasetId, cleanAll, paramOrderMap.get(se), searchOpParameter, false, useFullFasta, false);
//                        cleanAll = false;
//                        System.out.println("---------------------------------------------------------full-" + datasetId + "----------------------------------------------");
//                        System.gc();
//                        MainUtilities.cleanOutputFolder();
//                        runDataset(datasetId, cleanAll, paramOrderMap.get(se), searchOpParameter, true, useFullFasta, useOreginalInputs);
                    }
                }
                MainUtilities.cleanOutputFolder();
                System.exit(0);

            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                MainUtilities.cleanOutputFolder();
                System.exit(0);
            }
        }
        );
    }

    private static void runDataset(String datasetId, boolean cleanAll, List<String> paramOrder, SearchInputSetting searchOpParameter, boolean wholeDataTest, boolean fullFasta, boolean useOreginalInputs) {
        ArrayList<File> msFiles = new ArrayList<>();
        File datasetFolder = new File(Configurations.GET_DATA_FOLDER() + datasetId);//  
        File searchParamFile = null;
        File fastaFile = null;
        for (File f : datasetFolder.listFiles()) {
            if (cleanAll) {
                if (f.isDirectory() && f.getName().equals(searchOpParameter.getSelectedSearchEngine().getName())) {
                    for (File ff : f.listFiles()) {
                        if (ff.getName().startsWith(Configurations.DEFAULT_RESULT_NAME) && !ff.getName().endsWith(".txt")) {
                            System.out.println("to be deleted files " + ff.getAbsolutePath());
                            ff.delete();

                        }
                    }
                }
            }

            if (f.getName().endsWith(".mgf")) {
                msFiles.add(f);
            } else if (f.getName().endsWith(".fasta")) {
                fastaFile = f;
            } else if (f.getName().endsWith(".par")) {
                searchParamFile = f;
            }
        }

        Controller controller = new Controller();
        controller.processDataset(datasetId, msFiles.get(0), fastaFile, searchParamFile, searchOpParameter, wholeDataTest, fullFasta, paramOrder, useOreginalInputs);
        MainUtilities.cleanOutputFolder();

    }
}
