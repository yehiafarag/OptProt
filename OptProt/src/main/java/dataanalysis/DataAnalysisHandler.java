/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package dataanalysis;

import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
import com.compomics.util.experiment.identification.spectrum_assumptions.TagAssumption;
import com.compomics.util.experiment.io.identification.IdfileReader;
import com.compomics.util.experiment.io.identification.IdfileReaderFactory;
import com.compomics.util.experiment.io.mass_spectrometry.MsFileHandler;
import com.compomics.util.experiment.mass_spectrometry.spectra.Spectrum;
import com.compomics.util.io.IoUtil;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import no.uib.probe.optprot.util.MainUtilities;
import no.uib.probe.optprot.util.ScoreComparison;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author yfa041
 */
public class DataAnalysisHandler {

    public DataAnalysisHandler(String dataFolderURL, String datasetName) {
        try {
            double partCount = 4;
            Map<String, File> files = new LinkedHashMap<>();
            File dataFolder = new File(dataFolderURL);
            for (File f : dataFolder.listFiles()) {
                if (f.getName().endsWith(".mgf")) {
                    files.put("MGF", f);
                } else if (f.getName().endsWith(".fasta")) {
                    files.put("FASTA", f);
                } else if (f.getName().endsWith(".par")) {
                    files.put("SS", f);
                } else if (f.getName().endsWith(".t.xml")) {
                    files.put("XTan", f);
                } else if (f.getName().endsWith(".tags")) {
                    files.put("Tag", f);
                }
            }

            File msFile = files.get("MGF");
            IdentificationParameters identificationParam = IdentificationParameters.getIdentificationParameters(files.get("SS"));
            SearchParameters searchParameters = identificationParam.getSearchParameters();
            File xTandemFile = files.get("XTan");
            IdfileReader idReader = readerFactory.getFileReader(xTandemFile);
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(msFile, MainUtilities.OptProt_Waiting_Handler);

            Map<String, Double> fullSpectraMap = new LinkedHashMap<>();
            String[] titiles = msFileHandler.getSpectrumTitles(IoUtil.removeExtension(msFile.getName()));
            for (String titel : titiles) {
                fullSpectraMap.put(titel, 0.0);
            }

            final List<SpectrumMatch> matches = Collections.synchronizedList(idReader.getAllSpectrumMatches(msFileHandler, MainUtilities.OptProt_Waiting_Handler, searchParameters));

            for (SpectrumMatch sm : matches) {
                PeptideAssumption peptide = sm.getAllPeptideAssumptions().toList().get(0);
                fullSpectraMap.replace(sm.getSpectrumTitle(), peptide.getRawScore());
            }
//            double[] fullSpectraData = new double[titiles.length];
            double countId = 0;
            double countUnId = 0;
            double partI = Math.round((double) titiles.length / 4.0);
            double partII_I = partI;// Math.round((titiles.length - partI - partI) / 2.0);
            double partVI = titiles.length - (3.0 * partI);// - partII_I - partII_I;
//            System.out.println("total " + titiles.length + "   " + QuartileSize + "  " + partCount);
            double[][] quartileData = new double[(int) partCount][2];
            int counter = 0;
            int quartileIndex = 0;
            int left = titiles.length;
            for (int i = 0; i < titiles.length; i++) {
                double v = fullSpectraMap.get(titiles[i]);
//                fullSpectraData[i] = v;
                if (v > 0) {
                    countId++;
                } else {
                    countUnId++;
                }
                counter++;
                left--;
//                System.out.println("left " + left + " counter " + counter);
                if ((left == titiles.length - partI) || (left == 0) || (left == titiles.length - partI - partII_I) || (left == titiles.length - partI - partII_I - partII_I)) {
                    //initpart I 
                    System.out.println("A part" + quartileIndex + "total " + counter + "   " + partI);
                    System.out.println("A add q " + quartileIndex + "  " + counter);
                    quartileData[quartileIndex++] = new double[]{countId, countUnId};
                    counter = 0;
                    countId = 0;
                    countUnId = 0;

                }

//                if (counter >= QuartileSize || left == 0 && quartileIndex < quartileData.length) {
//                    System.out.println("add q " + quartileIndex + "  " + counter);
//                    quartileData[quartileIndex++] = new double[]{countId, countUnId};
//                    counter = 0;
//                    countId = 0;
//                    countUnId = 0;
//                }
            }

            BarChartView bcv = new BarChartView("PSM-Distribution " + datasetName, quartileData);
            System.out.println("A End " + "  " + counter);
            File tagFile = files.get("Tag");
            idReader = readerFactory.getFileReader(tagFile);
            ArrayList<SpectrumMatch> tagMatches = idReader.getAllSpectrumMatches(msFileHandler, MainUtilities.OptProt_Waiting_Handler, identificationParam.getSearchParameters());

            for (String titile : fullSpectraMap.keySet()) {
                fullSpectraMap.replace(titile, 100.0);
            }

            for (SpectrumMatch sm : tagMatches) {
                TagAssumption tag = sm.getAllTagAssumptions().toList().get(0);
                fullSpectraMap.replace(sm.getSpectrumTitle(), tag.getScore());
            }

            

            double[][] quartileEntData = new double[(int) partCount][];
            double quartIntens = 0;
            counter = 0;
            quartileIndex = 0;
            double totalFileInt = 0;
            String msFileName = IoUtil.removeExtension(msFile.getName());
            left = titiles.length;

            for (int i = 0; i < titiles.length; i++) {
                Spectrum spec = msFileHandler.getSpectrum(msFileName, titiles[i]);
                double score1 = ((spec.getMaxMz() - spec.getMinMz()));
                quartIntens += (score1);
//                quartIntens += (spec.getMaxMz() - spec.getMinMz())+(spec.getTotalIntensity()/(double)spec.getNPeaks());
                counter++;
                left--;

                if (left == titiles.length - partI) {
                    //initpart I 
                    System.out.println("part 1 total " + counter + "   " + partI);
                    quartileEntData[quartileIndex++] = new double[]{quartIntens / (double) counter, 0};
                    totalFileInt += quartileEntData[quartileIndex - 1][0];

                    counter = 0;
                    quartIntens = 0;

                }
                if (left == titiles.length - partI - partII_I) {
                    //initpart II 
                    System.out.println("part 2 total " + counter + "   " + partII_I);
                    quartileEntData[quartileIndex++] = new double[]{quartIntens / (double) counter, 0};
                    totalFileInt += quartileEntData[quartileIndex - 1][0];

                    counter = 0;
                    quartIntens = 0;

                }
                if (left == titiles.length - partI - partII_I - partII_I) {
                    //initpart II 
                    System.out.println("part 3 total " + counter + "   " + partII_I);
                    quartileEntData[quartileIndex++] = new double[]{quartIntens / (double) counter, 0};
                    totalFileInt += quartileEntData[quartileIndex - 1][0];

                    counter = 0;
                    quartIntens = 0;

                }
                if (left == 0) {
                    //initpart 4 
                    System.out.println("part 4 total " + counter + "   " + partII_I);
                    quartileEntData[quartileIndex++] = new double[]{quartIntens / (double) counter, 0};
                    totalFileInt += quartileEntData[quartileIndex - 1][0];

                    counter = 0;
                    quartIntens = 0;

                }

//                if (counter >= QuartileSize || left == 0 && quartileIndex < quartileData.length) {
//                    quartileEntData[quartileIndex++] = new double[]{quartIntens / (double) counter, 0};
//                    totalFileInt += quartileEntData[quartileIndex - 1][0];
//
//                    counter = 0;
//                    quartIntens = 0;
//                }
            }
//            quartileIndex = 1;
//            for (double[] ratio : quartileEntData) {
////                System.out.println("at Q " + quartileIndex + "  ratio is " + Math.round(ratio[0] / totalFileInt));
//                quartileIndex++;
//            }
            BarChartView bcv3 = new BarChartView("Int-Distribution " + datasetName, quartileEntData);

            System.out.println("at matches size " + matches.size());
        } catch (IOException | SQLException | ClassNotFoundException | InterruptedException | JAXBException | XmlPullParserException | XMLStreamException ex) {
            Logger.getLogger(DataAnalysisHandler.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public static void showTagDistribution(String datasetName, String[] titiles, List<SpectrumMatch> matches) {
        double partCount = 4;

        Map<String, Double> fullSpectraMap = new LinkedHashMap<>();
        for (String titel : titiles) {
            fullSpectraMap.put(titel, 100.0);
        }
        for (SpectrumMatch sm : matches) {
            TagAssumption tag = sm.getAllTagAssumptions().toList().get(0);
            fullSpectraMap.replace(sm.getSpectrumTitle(), tag.getScore());
        }
        double countId = 0;
        double countUnId = 0;
        double partI = Math.round((double) titiles.length / 4.0);
        double partII = 2 * partI;
        double partIII = 3 * partI;
        double[][] quartileData = new double[(int) partCount][2];
        int counter = 0;
        int quartileIndex = 0;
        int left = titiles.length;
        for (String titile : titiles) {
            double v = fullSpectraMap.get(titile);
            if (v <=0.1) {
                countId++;
            } else {
                countUnId++;
            }
            counter++;
            left--;
            if ((left == partI) || (left == 0) || (left == partII) || (left == partIII)) {
                //initpart I 
                quartileData[quartileIndex++] = new double[]{countId, countUnId};
                counter = 0;
                countId = 0;
                countUnId = 0;

            }
        }

        BarChartView bcv = new BarChartView("TAG-Distribution " + datasetName, quartileData);

    }
    /**
     * The identification file reader factory of compomics utilities.
     */
    private final IdfileReaderFactory readerFactory = IdfileReaderFactory.getInstance();

    public static void main(String[] args) {
//        String dataFolderURL = "D:\\Apps\\OptProt\\pepshaker\\data_ps\\PXD028427";
//        new DataAnalysisHandler(dataFolderURL, "PXD028427");
        String dataFolderURL2 = "D:\\Apps\\OptProt\\pepshaker\\data_ps\\PXD000561";
        new DataAnalysisHandler(dataFolderURL2, "PXD000561");
//        String dataFolderURL3 = "D:\\Apps\\OptProt\\pepshaker\\data_ps\\PXD001468";
//        new DataAnalysisHandler(dataFolderURL3, "PXD001468");
    }
}
