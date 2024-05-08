package no.uib.probe.optprot;

import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.io.mass_spectrometry.MsFileHandler;
import com.compomics.util.io.IoUtil;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.SwingUtilities;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.search.SearchOptimizerHandler;
import no.uib.probe.optprot.util.MainUtilities;
import no.uib.probe.optprot.util.OptProtWaitingHandler;
import no.uib.probe.optprot.util.SpectraFileUtilities;

/**
 * This app is search settings optimization workflow that aim to optimize search
 * settings for different proteomics search engines
 *
 * @author Yehia Mokhtar Farag
 */
public class OptProt {

    public OptProt(String datasetId, File identificationParametersFile, File oreginalFastaFile, ArrayList<File> msFiles, Advocate searchEngineToOptimise) {
        try {
            long start = System.currentTimeMillis();
            Configurations.Dataset_Id = datasetId;
            File resultsOutput = new File(Configurations.OUTPUT_FOLDER_PATH);
            if (resultsOutput.exists()) {
                MainUtilities.cleanOutputFolder();
            }
            SearchInputSetting searchOpParameter = new SearchInputSetting();
//        MainUtilities.cleanOptProtData(oreginalFastaFile.getParentFile());
            SearchOptimizerHandler searchOptimizerHandler = new SearchOptimizerHandler();
            searchOptimizerHandler.setIdentificationParametersFile(identificationParametersFile);
            SpectraFileUtilities spectraFileUtilities = new SpectraFileUtilities();
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(msFiles.get(0), new OptProtWaitingHandler());
            String fileNameWithoutExtension = IoUtil.removeExtension(msFiles.get(0).getName());
            String[] spectrumTitles = msFileHandler.getSpectrumTitles(fileNameWithoutExtension);
            double startRatio = 0.04;
            int totalTagsNumb = Math.max(1000, (int) (spectrumTitles.length * startRatio));
            while (true) {
                File[] subSetInputs = spectraFileUtilities.initInputSubSetFiles(msFiles.get(0), oreginalFastaFile, identificationParametersFile, totalTagsNumb);
                searchOptimizerHandler.setSubMsFile(subSetInputs[0]);
                searchOptimizerHandler.setSubFastaFile(subSetInputs[1]);
                searchOpParameter.setSelectedSearchEngine(searchEngineToOptimise);
                searchOptimizerHandler.setOptProtSearchParameters(searchOpParameter);
                searchOptimizerHandler.runReferenceSearch(false);
                System.out.println("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<id rate for reference run is " + " " + searchOptimizerHandler.getReferenceIdRate() + " >=  " + (totalTagsNumb * 0.085) + ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
                if (searchOptimizerHandler.getReferenceIdRate() >= (totalTagsNumb * 0.085)) {
                    break;
                }
                File cms = new File(subSetInputs[0].getParent(), subSetInputs[0].getName().replace(".mgf", ".cms"));
                cms.deleteOnExit();
                cms.delete();
                subSetInputs[0].delete();
                subSetInputs[1].delete();
                startRatio += 0.02;
                totalTagsNumb = Math.max(1000, (int) (spectrumTitles.length * startRatio));

            }
            boolean all = true;
            searchOpParameter.setOptimizeDigestionParameter(true || all);
            searchOpParameter.setOptimizeEnzymeParameter(false || all);
            searchOpParameter.setOptimizeSpecificityParameter(false || all);
            searchOpParameter.setOptimizeMaxMissCleavagesParameter(false || all);
            searchOpParameter.setOptimizeFragmentIonTypesParameter(false || all);
            searchOpParameter.setOptimizePrecursorToleranceParameter(false || all);
            searchOpParameter.setOptimizeFragmentToleranceParameter(false || all);
            searchOpParameter.setOptimizePrecursorChargeParameter(false || all);
            searchOpParameter.setOptimizeIsotopsParameter(false || all);
            searchOpParameter.setOptimizeModificationParameter(false || all);
            searchOpParameter.setRecalibrateSpectraParameter(false);
//            searchOpParameter.setRunXTandem(true);
            searchOpParameter.setSelectedSearchEngine(searchEngineToOptimise);
//            searchOpParameter.setOptimizeXtandemAdvancedParameter(true || all);
//            searchOpParameter.getXtandemOptProtAdvancedSearchParameters().setOptAll(true || all);
            long start1 = System.currentTimeMillis();
            searchOptimizerHandler.executeParametersOptimization();//
            long end1 = System.currentTimeMillis();
            
            
             double full = (end1 - start) / 1000.0;
            double total = (end1 - start1) / 1000.0;
            System.out.println("-------------------------------------------------------------------->>>>>>>>>>>>>>>>>> Elapsed Time for the optimization of " + datasetId + "  in seconds: " + total + "  and in mins " + (total / 60.0)+ "  full process in seconds: " + full + "  and in mins " + (full / 60.0));
            if (resultsOutput.exists()) {
                MainUtilities.cleanOutputFolder();
            }

//        MainUtilities.cleanOptProtData(oreginalFastaFile.getParentFile());
        } catch (IOException ex) {
            Logger.getLogger(OptProt.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            MainUtilities.cleanOutputFolder();
            System.exit(0);
        }

    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            ArrayList<File> msFiles = new ArrayList<>();
//            File oreginalMGFFile = new File("D:\\Apps\\OptProt\\data\\sample.mgf");
//            msFiles.add(oreginalMGFFile);
            String datasetId = "PXD028427";//PXD047036  PXD028427  PXD009340
            File datasetFolder = new File("D:\\Apps\\OptProt\\data\\" + datasetId);//  
            File searchParamFile = null;
            File fastaFile = null;
            boolean cleanAll = false;
            for (File f : datasetFolder.listFiles()) {
                if (f.getName().startsWith(Configurations.DEFAULT_RESULT_NAME)) {
                    if (cleanAll || (!f.getName().startsWith(Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_"))) {
                        f.delete();
                    } else if (f.getName().startsWith(Configurations.DEFAULT_RESULT_NAME + Configurations.get_current_file_fingerprent() + "_") && f.getName().endsWith(".cms")) {
                        f.delete();

                    }
                    continue;
                }
                if (f.getName().endsWith(".mgf")) {
                    msFiles.add(f);
                } else if (f.getName().endsWith(".fasta")) {
                    fastaFile = f;
                } else if (f.getName().endsWith(".par")) {
                    searchParamFile = f;
                }
            }

            SearchInputSetting searchOpParameter = new SearchInputSetting();
            boolean all = true;
            searchOpParameter.setOptimizeDigestionParameter(false || all);
            searchOpParameter.setOptimizeEnzymeParameter(false || all);
            searchOpParameter.setOptimizeSpecificityParameter(false || all);
            searchOpParameter.setOptimizeMaxMissCleavagesParameter(false || all);
            searchOpParameter.setOptimizeFragmentIonTypesParameter(false || all);
            searchOpParameter.setOptimizePrecursorToleranceParameter(false || all);
            searchOpParameter.setOptimizeFragmentToleranceParameter(false || all);
            searchOpParameter.setOptimizePrecursorChargeParameter(false || all);
            searchOpParameter.setOptimizeIsotopsParameter(false || all);
            searchOpParameter.setOptimizeModificationParameter(false || all);
//            searchOpParameter.setRecalibrateSpectraParameter(false);
//            searchOpParameter.setRunXTandem(true);
            searchOpParameter.setSelectedSearchEngine(Advocate.xtandem);
            searchOpParameter.setOptimizeXtandemAdvancedParameter(true || all);
            searchOpParameter.getXtandemOptProtAdvancedSearchParameters().setOptAll(true || all);

            Controller controller = new Controller();
//            controller.processDataset(datasetId, msFiles.get(0), fastaFile, searchParamFile, Advocate.myriMatch);
//            searchParamFile = new File(Configurations.DEFAULT_OPTPROT_SEARCH_SETTINGS_FILE);
            controller.processDataset(datasetId, msFiles.get(0), fastaFile, searchParamFile, searchOpParameter);
//            OptProt optProt = new OptProt(datasetId, searchParamFile, fastaFile, msFiles,Advocate.myriMatch);

        });
    }
}
