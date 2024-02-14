package no.uib.probe.optprot;

import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.io.identification.IdfileReader;
import com.compomics.util.experiment.io.identification.IdfileReaderFactory;
import com.compomics.util.experiment.io.mass_spectrometry.MsFileHandler;
import com.compomics.util.io.IoUtil;
import com.compomics.util.io.file.filefilters.FastaFileFilter;
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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.SwingUtilities;
import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import no.uib.probe.optprot.model.Configurations;
import no.uib.probe.optprot.model.SearchEngineParameters;
import no.uib.probe.optprot.search.SubSearchHandler;
import no.uib.probe.optprot.util.OptProtWaitingHandler;
import no.uib.probe.optprot.util.SpectraFileUtilities;
import org.xmlpull.v1.XmlPullParserException;

/**
 * This app is search settings optimization workflow that aim to optimize search
 * settings for different proteomics search engines
 *
 * @author Yehia Mokhtar Farag
 */
public class OptProt {

    /**
     * The identification settings file.
     */
    private File identificationParametersFile;

    /**
     * The search parameters.
     */
    private IdentificationParameters identificationParameters = null;

    /**
     * The processing preferences.
     */
    private ProcessingParameters processingParameters;

    private final SubSearchHandler subSearchHandler;

    public OptProt() {
        this.identificationParametersFile = new File("D:\\Apps\\OptProt\\data\\searchparam.par");
        if (this.identificationParametersFile == null) {
            if (identificationParametersFile == null) {
                try {
                    String name = identificationParameters.getName();
                    if (name == null) {
                        name = "SearchCLI.par";
                    } else {
                        name += ".par";

                    }
//                    identificationParametersFile = new File(Configurations.CONFIG_FOLDER, name);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
            }

        }
        ArrayList<File> msList = initSubSet();
        File fasta = initFastaFile();
        this.subSearchHandler = new SubSearchHandler(identificationParametersFile, msList, fasta);
        long start1 = System.currentTimeMillis();
        this.subSearchHandler.execute(DigestionParameters.getDefaultParameters());
        long end1 = System.currentTimeMillis();
        double total = (end1 - start1) / 1000.0;
        System.out.println("-------------------------------------------------------------------->>>>>>>>>>>>>>>>>> Elapsed Time in nano seconds: " + total + "   " + (total / 3));
    }

    private ArrayList<File> initSubSet() {
        ArrayList<File> msFiles = new ArrayList<>();
        File sampleMgf = new File("D:\\Apps\\OptProt\\data\\sample.mgf");
        try {
            File subSampleMgf = new File("D:\\Apps\\OptProt\\data\\sub_sample.mgf");
            if (subSampleMgf.exists()) {
                subSampleMgf.delete();
                File subSampleCMS = new File("D:\\Apps\\OptProt\\data\\sub_sample.cms");
                subSampleCMS.delete();
            } else {
                subSampleMgf.createNewFile();
            }
//            SpectraFileUtilities.writeSubSetEveryNMgfFile(sampleMgf, subSampleMgf,500.0);
            SpectraFileUtilities.writeSubSetTargtedFileAreaMgfFile(sampleMgf, subSampleMgf, 250.0);
            msFiles.add(subSampleMgf);
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return msFiles;
    }

    private File initFastaFile() {
        File fastaFile = new File("D:\\Apps\\OptProt\\data\\sample.fasta");
        File subFastaFile = new File("D:\\Apps\\OptProt\\data\\sub_sample.fasta");
        if (subFastaFile.exists()) {
            subFastaFile.delete();
        }
        SpectraFileUtilities.createSubFastaFile(fastaFile, subFastaFile,false);
        return subFastaFile;

    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            OptProt optProt = new OptProt();

        });
    }
}
