package no.uib.probe.optprot;

import com.compomics.util.experiment.identification.Advocate;
import static com.sun.tools.xjc.reader.Ring.add;
import java.io.File;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import javax.swing.SwingUtilities;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.model.SearchInputSetting;

/**
 * This app is search settings optimization workflow that aim to optimize search
 * settings for different proteomics search engines
 *
 * @author Yehia Mokhtar Farag
 */
public class OptProt {

//    public OptProt(String datasetId, File identificationParametersFile, File oreginalFastaFile, ArrayList<File> msFiles, Advocate searchEngineToOptimise) {
//        try {
//            long start = System.currentTimeMillis();
//            Configurations.Dataset_Id = datasetId;
//            File seFolder = new File(Configurations.DATA_FOLDER + datasetId, searchEngineToOptimise.getName());
//            if (seFolder.exists()) {
//                for (File f : seFolder.listFiles()) {
//                    System.out.println("File found  " + f.getName());
//                }
//            } else {
//                seFolder.mkdir();
//            }
//
//            System.exit(0);
//            System.out.println("opt prot still on :( ");
//            File resultsOutput = new File(Configurations.OUTPUT_FOLDER_PATH);
//            if (resultsOutput.exists()) {
//                MainUtilities.cleanOutputFolder();
//            }
//            SearchInputSetting searchOpParameter = new SearchInputSetting();
////        MainUtilities.cleanOptProtData(oreginalFastaFile.getParentFile());
//            SearchOptimizerHandler searchOptimizerHandler = new SearchOptimizerHandler();
//            searchOptimizerHandler.setIdentificationParametersFile(identificationParametersFile);
//            SpectraFileUtilities spectraFileUtilities = new SpectraFileUtilities();
//            MsFileHandler msFileHandler = new MsFileHandler();
//            msFileHandler.register(msFiles.get(0), new OptProtWaitingHandler());
//            String fileNameWithoutExtension = IoUtil.removeExtension(msFiles.get(0).getName());
//            String[] spectrumTitles = msFileHandler.getSpectrumTitles(fileNameWithoutExtension);
//            double startRatio = 0.04;
//            int totalTagsNumb = Math.max(1000, (int) (spectrumTitles.length * startRatio));
//            while (true) {
//                File[] subSetInputs = spectraFileUtilities.initInputSubSetFiles(msFiles.get(0), oreginalFastaFile, identificationParametersFile, totalTagsNumb);
//                searchOptimizerHandler.setSubMsFile(subSetInputs[0]);
//                searchOptimizerHandler.setSubFastaFile(subSetInputs[1]);
//                searchOpParameter.setSelectedSearchEngine(searchEngineToOptimise);
//                searchOptimizerHandler.setOptProtSearchParameters(searchOpParameter);
//                searchOptimizerHandler.runReferenceSearch(false);
//                System.out.println("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<id rate for reference run is " + " " + searchOptimizerHandler.getReferenceIdRate() + " >=  " + (totalTagsNumb * 0.085) + ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
//                if (searchOptimizerHandler.getReferenceIdRate() >= (totalTagsNumb * 0.085)) {
//                    break;
//                }
//                File cms = new File(subSetInputs[0].getParent(), subSetInputs[0].getName().replace(".mgf", ".cms"));
//                cms.deleteOnExit();
//                cms.delete();
//                subSetInputs[0].delete();
//                subSetInputs[1].delete();
//                startRatio += 0.02;
//                totalTagsNumb = Math.max(1000, (int) (spectrumTitles.length * startRatio));
//
//            }
//            boolean all = true;
//            searchOpParameter.setOptimizeDigestionParameter(true || all);
////            searchOpParameter.setOptimizeEnzymeParameter(false || all);
////            searchOpParameter.setOptimizeSpecificityParameter(false || all);
////            searchOpParameter.setOptimizeMaxMissCleavagesParameter(false || all);
//            searchOpParameter.setOptimizeFragmentIonTypesParameter(false || all);
//            searchOpParameter.setOptimizePrecursorToleranceParameter(false || all);
//            searchOpParameter.setOptimizeFragmentToleranceParameter(false || all);
//            searchOpParameter.setOptimizePrecursorChargeParameter(false || all);
//            searchOpParameter.setOptimizeIsotopsParameter(false || all);
//            searchOpParameter.setOptimizeModificationParameter(false || all);
//            searchOpParameter.setRecalibrateSpectraParameter(false);
////            searchOpParameter.setRunXTandem(true);
//            searchOpParameter.setSelectedSearchEngine(searchEngineToOptimise);
////            searchOpParameter.setOptimizeXtandemAdvancedParameter(true || all);
////            searchOpParameter.getXtandemOptProtAdvancedSearchParameters().setOptAll(true || all);
//            long start1 = System.currentTimeMillis();
//            searchOptimizerHandler.executeParametersOptimization();//
//            long end1 = System.currentTimeMillis();
//
//            double full = (end1 - start) / 1000.0;
//            double total = (end1 - start1) / 1000.0;
//            System.out.println("-------------------------------------------------------------------->>>>>>>>>>>>>>>>>> Elapsed Time for the optimization of " + datasetId + "  in seconds: " + total + "  and in mins " + (total / 60.0) + "  full process in seconds: " + full + "  and in mins " + (full / 60.0));
//            if (resultsOutput.exists()) {
//                MainUtilities.cleanOutputFolder();
//            }
//
////        MainUtilities.cleanOptProtData(oreginalFastaFile.getParentFile());
//        } catch (IOException ex) {
//            Logger.getLogger(OptProt.class.getName()).log(Level.SEVERE, null, ex);
//        } finally {
//            MainUtilities.cleanOutputFolder();
//            System.exit(0);
//        }
//
//    }
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            try {
                List<String> paramOrder = new ArrayList<>();
                Set<Advocate> supportedSearchEngine = new LinkedHashSet<>();
//                supportedSearchEngine.add(Advocate.xtandem);
                supportedSearchEngine.add(Advocate.myriMatch);
//                supportedSearchEngine.add(Advocate.sage);
//                paramOrder.add("ModificationParameter");
//                paramOrder.add("DigestionParameter_1");
//                paramOrder.add("FragmentIonTypesParameter");
//                paramOrder.add("DigestionParameter_2");
//                paramOrder.add("FragmentToleranceParameter");
//                paramOrder.add("PrecursorChargeParameter");
                paramOrder.add("IsotopParameter");
//
//                paramOrder.add("XtandemAdvancedParameter");
//                paramOrder.add("MyriMatchAdvancedParameter");
//                paramOrder.add("DigestionParameter_3");
//                paramOrder.add("PrecursorToleranceParameter");
                String datasetId = "PXD028427";//PXD028427 PXD047036   PXD009340 PXD000561   PXD000815  PXD001250 PXD001468
                boolean cleanAll = false;
                SearchInputSetting searchOpParameter = new SearchInputSetting();
                boolean all = true;
                searchOpParameter.setOptimizeAllParameters(all);
                searchOpParameter.setOptimizeDigestionParameter(false || all);

                searchOpParameter.setOptimizeFragmentIonTypesParameter(false || all);
                searchOpParameter.setOptimizePrecursorToleranceParameter(false || all);
                searchOpParameter.setOptimizeFragmentToleranceParameter(false || all);
                searchOpParameter.setOptimizePrecursorChargeParameter(true || all);
                searchOpParameter.setOptimizeIsotopsParameter(false || all);
                searchOpParameter.setOptimizeModificationParameter(false || all);
//            searchOpParameter.setRecalibrateSpectraParameter(false);
//             
//               //|| all
//            searchOpParameter.setOptimizeMyriMatchAdvancedParameter(false || all);
//            searchOpParameter.getXtandemOptProtAdvancedSearchParameters().setOptAll(false );//|| all
//            runDataset(datasetId, cleanAll, paramOrder, searchOpParameter);
//            File oreginalMGFFile = new File("D:\\Apps\\OptProt\\data\\sample.mgf");
//            msFiles.add(oreginalMGFFile);
//            searchOpParameter.setSelectedSearchEngine(Advocate.xtandem);
//            searchOpParameter.setOptimizeXtandemAdvancedParameter(false || all);//|| all
//            runDataset(datasetId, cleanAll, paramOrder, searchOpParameter);

                for (Advocate se : supportedSearchEngine) {
                    searchOpParameter.setSelectedSearchEngine(se);
//                searchOpParameter.setSelectedSearchEngine(Advocate.sage);
//                searchOpParameter.setSelectedSearchEngine(Advocate.myriMatch);
                    runDataset(datasetId, cleanAll, paramOrder, searchOpParameter, false);
                }
                System.exit(0);

            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                System.exit(0);
            }
        }
        );
    }

    private static void runDataset(String datasetId, boolean cleanAll, List<String> paramOrder, SearchInputSetting searchOpParameter, boolean wholeDataTest) {
        ArrayList<File> msFiles = new ArrayList<>();
        File datasetFolder = new File(Configurations.DATA_FOLDER + datasetId);//  
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
        controller.processDataset(datasetId, msFiles.get(0), fastaFile, searchParamFile, searchOpParameter, wholeDataTest, paramOrder);

    }
}
