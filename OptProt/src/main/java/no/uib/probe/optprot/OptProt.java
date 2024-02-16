package no.uib.probe.optprot;

import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.parameters.tools.ProcessingParameters;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Set;
import javax.swing.SwingUtilities;
import no.uib.probe.optprot.model.SearchOptimizerParameters;
import no.uib.probe.optprot.search.SearchOptimizerHandler;
import no.uib.probe.optprot.util.SpectraFileUtilities;

/**
 * This app is search settings optimization workflow that aim to optimize search
 * settings for different proteomics search engines
 *
 * @author Yehia Mokhtar Farag
 */
public class OptProt {

    private final SearchOptimizerHandler subSearchHandler;

    public OptProt(File identificationParametersFile, File oreginalFastaFile, ArrayList<File> msFiles) {

        this.subSearchHandler = new SearchOptimizerHandler();
        this.subSearchHandler.setIdentificationParametersFile(identificationParametersFile);

        long start1 = System.currentTimeMillis();
         File subMGFFile = new File("D:\\Apps\\OptProt\\data\\sub_mgf_file.mgf");
        subSearchHandler.setSubMsFile(subMGFFile);//initSpectraSubSet(msFiles));
        subSearchHandler.setSubFastaFile(oreginalFastaFile);
        File subFasta = new File("D:\\Apps\\OptProt\\data\\sub_fasta_file.fasta");// initSubFastaFile(oreginalFastaFile, subSearchHandler.excuteNovorSearches());
        subSearchHandler.setSubFastaFile(subFasta);
//     
        SearchOptimizerParameters searchOpParameter=this.initSearchOptimizerParameters();
        searchOpParameter.setRunXTandem(true);
        searchOpParameter.setOptimizeDigestionParameter(true);
        searchOpParameter.setOptimizeEnzymeParameter(false);
        searchOpParameter.setOptimizeSpecificityParameter(false);
        searchOpParameter.setOptimizeMaxMissCleavagesParameter(false);
        searchOpParameter.setOptimizeFragmentIonTypesParameter(false);
        searchOpParameter.setOptimizePrecursorToleranceParameter(true);
         searchOpParameter.setOptimizeFragmentToleranceParameter(true);
        this.subSearchHandler.executeParametersOptimization(searchOpParameter);
        long end1 = System.currentTimeMillis();
        double total = (end1 - start1) / 1000.0;
        System.out.println("-------------------------------------------------------------------->>>>>>>>>>>>>>>>>> Elapsed Time in nano seconds: " + total + "   " + (total / 3));
    }

    private File initSpectraSubSet(ArrayList<File> msFiles) {

        File oreginalMGFFile = msFiles.get(0);
        try {
            File subMGFFile = new File("D:\\Apps\\OptProt\\data\\sub_mgf_file.mgf");
            if (subMGFFile.exists()) {
                subMGFFile.delete();
                File subSampleCMS = new File("D:\\Apps\\OptProt\\data\\sub_mgf_file.cms");
                subSampleCMS.delete();
            } else {
                subMGFFile.createNewFile();
            }
//            SpectraFileUtilities.writeSubSetEveryNMgfFile(sampleMgf, subSampleMgf,500.0);
            SpectraFileUtilities.writeSubSetTargtedFileAreaMgfFile(oreginalMGFFile, subMGFFile, 250);
            return subMGFFile;
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return oreginalMGFFile;
    }

    private File initSubFastaFile(File fastaFile, Set<String> sequences) {

        File subFastaFile = new File("D:\\Apps\\OptProt\\data\\sub_fasta_file.fasta");
        if (subFastaFile.exists()) {
            subFastaFile.delete();
        }
        SpectraFileUtilities.createSubFastaFile(fastaFile, subFastaFile, sequences);
//        SpectraFileUtilities.createSubFastaFile(fastaFile, subFastaFile,true);
        return subFastaFile;

    }
    
    private SearchOptimizerParameters initSearchOptimizerParameters (){
      SearchOptimizerParameters searchOptimizerParameters = new SearchOptimizerParameters();
        searchOptimizerParameters.setxTandemFolder(new File("D:\\Apps\\searchgui\\resources\\XTandem\\windows\\windows_64bit"));
        searchOptimizerParameters.setNovorFolder(new File("D:\\Apps\\searchgui\\resources\\Novor"));
        searchOptimizerParameters.setCometFolder(new File("D:\\Apps\\searchgui\\resources\\Comet\\windows"));
        return searchOptimizerParameters;
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            ArrayList<File> msFiles = new ArrayList<>();
            File oreginalMGFFile = new File("D:\\Apps\\OptProt\\data\\sample.mgf");
            msFiles.add(oreginalMGFFile);
            OptProt optProt = new OptProt(new File("D:\\Apps\\OptProt\\data\\searchparam.par"), new File("D:\\Apps\\OptProt\\data\\sample.fasta"), msFiles);

        });
    }
}
