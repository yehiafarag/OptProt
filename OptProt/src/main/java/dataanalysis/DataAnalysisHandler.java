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
import com.compomics.util.io.IoUtil;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import java.awt.Rectangle;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JLabel;
import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import no.uib.probe.optprot.util.MainUtilities;
import no.uib.probe.optprot.util.StatisticsTests;
import org.apache.commons.math.stat.StatUtils;
import org.jfree.svg.SVGGraphics2D;
import org.jfree.svg.SVGUtils;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author yfa041
 */
public class DataAnalysisHandler {

    public DataAnalysisHandler(String dataFolderURL, String datasetName) {
        double quartiles = 4;
        Map<String, File> files = new LinkedHashMap<>();
        File dataFolder = new File(dataFolderURL);
        for (File f : dataFolder.listFiles()) {

            if (f.isDirectory() && f.getName().equalsIgnoreCase("optprot")) {
                for (File ff : f.listFiles()) {
                    if (ff.getName().startsWith("optsearch_results") && ff.getName().endsWith(".mgf")) {
                        files.put("OPT_MGF", ff);
                    } else if (ff.getName().startsWith("optsearch_results") && ff.getName().endsWith(".tags")) {
                        files.put("OPT_Tag", ff);
                    } else if (ff.getName().endsWith("_fullFasta.t.xml")) {
                        files.put("OPT_FullFasta_XTan", ff);
                    } else if (ff.getName().endsWith(".t.xml")) {
                        files.put("OPT_XTan", ff);
                    } else if (ff.getName().endsWith(".fasta")) {
                        files.put("OPT_FASTA", ff);
                    } else if (ff.getName().endsWith("optprot.par")) {
                        files.put("OPT_SS", ff);
                    }
                }

            } else if (f.getName().endsWith(".mgf")) {
                files.put("MGF", f);
            } else if (f.getName().endsWith(".fasta")) {
                files.put("FASTA", f);
            } else if (f.getName().endsWith("default.par")) {
                files.put("SS", f);
            } else if (f.getName().endsWith(".par")) {
                files.put("DS", f);
            } else if (f.getName().endsWith(".t.xml")) {
                files.put("XTan", f);
            } else if (f.getName().endsWith(".tags")) {
                files.put("Tag", f);
            }
        }
        File msFile = files.get("MGF");
        File xTandemFile = files.get("XTan");
        String title = "Identified spectra in full dataset quartiles (" + datasetName + ")";
        BarChartView bcv = DataAnalysisHandler.getIdentifiedSpectrumChart(msFile, files.get("SS"), xTandemFile, title);
        //
        File tagFile = files.get("Tag");
        title = "Spectra with confident tag in full dataset quartiles (" + datasetName + ")";
        BarChartView bcv2 = DataAnalysisHandler.getSpectrumTagChart(msFile, files.get("DS"), tagFile, title);
        System.out.println("tag file " + files.get("DS").getAbsolutePath());
        double agreeId = 0;
        double disagreeId = 0;
        int i = 0;
        for (String spectitle : bcv2.getConfidentScores().keySet()) {
            if (bcv.getConfidentScores().get(spectitle) > 0) {
                agreeId++;
            } else {
                disagreeId++;
            }
            i++;
        }
        double idAgreeRatio = agreeId / (agreeId + disagreeId);
        idAgreeRatio = Math.round(idAgreeRatio * 100.0);
        System.out.println("agree " + agreeId + "  disagree " + disagreeId + "  " + idAgreeRatio);
        Map<String, Double[]> dataMap = new LinkedHashMap<>();
        dataMap.put("Identified Spectra", bcv.getValues());
        dataMap.put("Confident Tag Spectra", bcv2.getValues());
        LineChartView lcView = new LineChartView(title, new String[]{"Q1", "Q2", "Q3", "Q4"}, dataMap);

//        double correlation = StatisticsTests.calculatePearsonCorrelationTest(confidentTag, comparableIdSpectra);
//        double correlation = StatisticsTests.calculatePearsonCorrelationTest(bcv.getScores(), bcv2.getScores());
        SVGGraphics2D g2 = new SVGGraphics2D((30 + 1400), 480);
        Rectangle r1 = new Rectangle(0, 20, 500, 400);
        bcv.getChart().draw(g2, r1);
        Rectangle r2 = new Rectangle(510, 20, 500, 400);
        bcv2.getChart().draw(g2, r2);
        Rectangle r3 = new Rectangle(1020, 20, 400, 400);
        lcView.getChart().draw(g2, r3);
        g2.drawString("Spectra with confidentTag and Identified ratio " + (int) idAgreeRatio + "  %", 40, 440);
        File f = new File(dataFolder, "full dataset(" + datasetName + ").svg");
        try {
            SVGUtils.writeToSVG(f, g2.getSVGElement());
        } catch (IOException ex) {
            Logger.getLogger(BarChartView.class.getName()).log(Level.SEVERE, null, ex);
        }
        //opt prot analysys

        msFile = files.get("OPT_MGF");
        xTandemFile = files.get("OPT_FullFasta_XTan");
        title = "Identified spectra in sub dataset with full fasta (" + datasetName + ")";
        BarChartView sub_bcv = DataAnalysisHandler.getIdentifiedSpectrumChart(msFile, files.get("OPT_SS"), xTandemFile, title);

        xTandemFile = files.get("OPT_XTan");
        title = "Identified spectra in sub dataset with sub fasta (" + datasetName + ")";
        BarChartView sub_bcv2 = DataAnalysisHandler.getIdentifiedSpectrumChart(msFile, files.get("OPT_SS"), xTandemFile, title);

        agreeId = 0;
        disagreeId = 0;
        for (String spectitle : sub_bcv2.getConfidentScores().keySet()) {
            if (sub_bcv.getConfidentScores().get(spectitle) > 0) {
                agreeId++;
            } else {
                disagreeId++;
            }
            i++;
        }

        idAgreeRatio = agreeId / (agreeId + disagreeId);
        idAgreeRatio = Math.round(idAgreeRatio * 100.0);
        System.out.println("agree " + agreeId + "  disagree " + disagreeId + "  " + idAgreeRatio);
        dataMap = new LinkedHashMap<>();
        dataMap.put("Identified Spectra(Full FASTA)", sub_bcv.getValues());
        dataMap.put("Identified Spectra(SUB FASTA)", sub_bcv2.getValues());
        LineChartView sub_lcView = new LineChartView(title, new String[]{"Q1", "Q2", "Q3", "Q4"}, dataMap);

//        double correlation = StatisticsTests.calculatePearsonCorrelationTest(confidentTag, comparableIdSpectra);
//        double correlation = StatisticsTests.calculatePearsonCorrelationTest(bcv.getScores(), bcv2.getScores());
        SVGGraphics2D sub_g2 = new SVGGraphics2D((30 + 1400), 480);
         r1 = new Rectangle(0, 20, 500, 400);
        sub_bcv.getChart().draw(sub_g2, r1);
         r2 = new Rectangle(510, 20, 500, 400);
        sub_bcv2.getChart().draw(sub_g2, r2);
        r3 = new Rectangle(1020, 20, 400, 400);
        sub_lcView.getChart().draw(sub_g2, r3);
        sub_g2.drawString("Identified with Full and SUB FASTA file ratio " + (int) idAgreeRatio + "  %", 40, 440);
        File sub_f = new File(dataFolder, "Sub_dataset(" + datasetName + ").svg");
        try {
            SVGUtils.writeToSVG(sub_f, sub_g2.getSVGElement());
        } catch (IOException ex) {
            Logger.getLogger(BarChartView.class.getName()).log(Level.SEVERE, null, ex);
        }

//        msFile = files.get("OPT_MGF");
//        File opttagFile = files.get("OPT_Tag");
//        title = "Spectra with confident tag in sub-dataset quartiles (" + datasetName + ")";
//        BarChartView bcv3 = DataAnalysisHandler.getSpectrumTagChart(msFile, files.get("DS"), opttagFile, title);
//
//        correlation = StatisticsTests.calculatePearsonCorrelationTest(bcv2.getValues(), bcv3.getValues());
//        g2 = new SVGGraphics2D(1040, 600);
//        r1 = new Rectangle(20, 20, 500, 500);
//        bcv2.getChart().draw(g2, r1);
//        r2 = new Rectangle(540, 20, 500, 500);
//        bcv3.getChart().draw(g2, r2);
//        g2.drawString("Correlation between Confident tag in full and sub dataset is " + correlation, 20, 560);
//        f = new File(dataFolder, "sub dataset(" + datasetName + ").svg");
//        try {
//            SVGUtils.writeToSVG(f, g2.getSVGElement());
//        } catch (IOException ex) {
//            Logger.getLogger(BarChartView.class.getName()).log(Level.SEVERE, null, ex);
//        }
//
//
//            double[][] quartileEntData = new double[(int) partCount][];
//            double quartIntens = 0;
//            counter = 0;
//            quartileIndex = 0;
//            double totalFileInt = 0;
//            String msFileName = IoUtil.removeExtension(msFile.getName());
//            left = titiles.length;
//
//            for (int i = 0; i < titiles.length; i++) {
//                Spectrum spec = msFileHandler.getSpectrum(msFileName, titiles[i]);
//                double score1 = ((spec.getMaxMz() - spec.getMinMz()));
//                quartIntens += (score1);
////                quartIntens += (spec.getMaxMz() - spec.getMinMz())+(spec.getTotalIntensity()/(double)spec.getNPeaks());
//                counter++;
//                left--;
//
//                if (left == titiles.length - partI) {
//                    //initpart I 
//                    System.out.println("part 1 total " + counter + "   " + partI);
//                    quartileEntData[quartileIndex++] = new double[]{quartIntens / (double) counter, 0};
//                    totalFileInt += quartileEntData[quartileIndex - 1][0];
//
//                    counter = 0;
//                    quartIntens = 0;
//
//                }
//                if (left == titiles.length - partI - partII_I) {
//                    //initpart II 
//                    System.out.println("part 2 total " + counter + "   " + partII_I);
//                    quartileEntData[quartileIndex++] = new double[]{quartIntens / (double) counter, 0};
//                    totalFileInt += quartileEntData[quartileIndex - 1][0];
//
//                    counter = 0;
//                    quartIntens = 0;
//
//                }
//                if (left == titiles.length - partI - partII_I - partII_I) {
//                    //initpart II 
//                    System.out.println("part 3 total " + counter + "   " + partII_I);
//                    quartileEntData[quartileIndex++] = new double[]{quartIntens / (double) counter, 0};
//                    totalFileInt += quartileEntData[quartileIndex - 1][0];
//
//                    counter = 0;
//                    quartIntens = 0;
//
//                }
//                if (left == 0) {
//                    //initpart 4 
//                    System.out.println("part 4 total " + counter + "   " + partII_I);
//                    quartileEntData[quartileIndex++] = new double[]{quartIntens / (double) counter, 0};
//                    totalFileInt += quartileEntData[quartileIndex - 1][0];
//
//                    counter = 0;
//                    quartIntens = 0;
//
//                }
//
////                if (counter >= QuartileSize || left == 0 && quartileIndex < quartileData.length) {
////                    quartileEntData[quartileIndex++] = new double[]{quartIntens / (double) counter, 0};
////                    totalFileInt += quartileEntData[quartileIndex - 1][0];
////
////                    counter = 0;
////                    quartIntens = 0;
////                }
//            }
////            quartileIndex = 1;
////            for (double[] ratio : quartileEntData) {
//////                System.out.println("at Q " + quartileIndex + "  ratio is " + Math.round(ratio[0] / totalFileInt));
////                quartileIndex++;
////            }
//            BarChartView bcv3 = new BarChartView("Int-Distribution " + datasetName, quartileEntData);
//
//            System.out.println("at matches size " + matches.size());
    }

    public static BarChartView getIdentifiedSpectrumChart(File msFile, File searchParam, File idFile, String title) {
        try {
            IdentificationParameters identificationParam = IdentificationParameters.getIdentificationParameters(searchParam);
            SearchParameters searchParameters = identificationParam.getSearchParameters();
            IdfileReader idReader = readerFactory.getFileReader(idFile);
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(msFile, MainUtilities.OptProt_Waiting_Handler);

            Map<String, Double> fullSpectraMap = new LinkedHashMap<>();
            String[] titiles = msFileHandler.getSpectrumTitles(IoUtil.removeExtension(msFile.getName()));
//            double[] scores = new double[titiles.length];
            for (String titel : titiles) {
                fullSpectraMap.put(titel, 0.0);
            }

            final List<SpectrumMatch> matches = Collections.synchronizedList(idReader.getAllSpectrumMatches(msFileHandler, MainUtilities.OptProt_Waiting_Handler, searchParameters));
            int index = 0;
            for (SpectrumMatch sm : matches) {
                PeptideAssumption peptide = sm.getAllPeptideAssumptions().toList().get(0);
//                if (peptide.getScore() <= 0.01) {
                fullSpectraMap.replace(sm.getSpectrumTitle(), peptide.getRawScore());
//                }
            }
            double countId = 0;
            double countUnId = 0;
            double quartileSize = Math.round((double) titiles.length / 4.0);
            Double[][] quartileData = new Double[4][2];
            int quartileCounter = 0;
            int quartileIndex = 0;
            int left = titiles.length;
            for (String titile : titiles) {
                double v = fullSpectraMap.get(titile);
                if (v > 0) {
                    countId++;

                } else {
                    countUnId++;
                }
                quartileCounter++;
                left--;
                if (quartileCounter >= quartileSize) {
                    quartileData[quartileIndex] = new Double[]{countId, countUnId};
                    System.out.println("Q" + quartileIndex + "  " + countId + "  " + countUnId + "  " + left);
                    quartileIndex++;
                    quartileCounter = 0;
                    countId = 0;
                    countUnId = 0;
                    if (left < quartileSize) {
                        quartileSize = left;
                    }

                }

            }
            BarChartView bcv = new BarChartView(title, quartileData);
            bcv.setConfidentScores(fullSpectraMap);
            System.out.println("Full spec map " + fullSpectraMap.size());
            return bcv;

        } catch (IOException | SQLException | ClassNotFoundException | InterruptedException | JAXBException | XmlPullParserException | XMLStreamException ex) {
            Logger.getLogger(DataAnalysisHandler.class.getName()).log(Level.SEVERE, null, ex);
        }
        return null;

    }

    public static BarChartView getSpectrumTagChart(File msFile, File searchParam, File idFile, String title) {
        try {
            IdentificationParameters identificationParam = IdentificationParameters.getIdentificationParameters(searchParam);
            SearchParameters searchParameters = identificationParam.getSearchParameters();
            System.out.println("id file exist " + idFile.exists() + "  " + idFile.getAbsolutePath());
            IdfileReader idReader = readerFactory.getFileReader(idFile);
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(msFile, MainUtilities.OptProt_Waiting_Handler);

            Map<String, Double> fullSpectraMap = new LinkedHashMap<>();
            String[] titiles = msFileHandler.getSpectrumTitles(IoUtil.removeExtension(msFile.getName()));
            Map<String, Double> confidentTagMap = new HashMap<>();
            for (String titel : titiles) {
                fullSpectraMap.put(titel, 100.0);
            }

            ArrayList<SpectrumMatch> tagMatches = idReader.getAllSpectrumMatches(msFileHandler, MainUtilities.OptProt_Waiting_Handler, searchParameters);

            for (SpectrumMatch sm : tagMatches) {
                TagAssumption tag = sm.getAllTagAssumptions().toList().get(0);
                fullSpectraMap.replace(sm.getSpectrumTitle(), tag.getScore());
            }
            double countId = 0;
            double countUnId = 0;
            double quartileSize = Math.round((double) titiles.length / 4.0);
            Double[][] quartileData = new Double[4][2];
            int quartileCounter = 0;
            int quartileIndex = 0;
            int left = titiles.length;
            for (String specTitle : titiles) {
                double v = fullSpectraMap.get(specTitle);
                if (v <= 0.01) {
                    countId++;
                    confidentTagMap.put(specTitle, v);
                } else {
                    countUnId++;
                }

//                  if (v < 1) {
//                    scores[index++] = (1.0 - v);
//                } else {
//                    scores[index++] = (0);
//                }
//                
                quartileCounter++;
                left--;
                if (quartileCounter >= quartileSize) {
                    quartileData[quartileIndex] = new Double[]{countId, countUnId};
                    System.out.println("Q" + quartileIndex + "  " + countId + "  " + countUnId + "  " + left);
                    quartileIndex++;
                    quartileCounter = 0;
                    countId = 0;
                    countUnId = 0;
                    if (left < quartileSize) {
                        quartileSize = left;
                    }

                }

            }
            BarChartView bcv = new BarChartView(title, quartileData);
            bcv.setConfidentScores(confidentTagMap);
            System.out.println("confidentTagMap " + confidentTagMap.size() + "  ");
            return bcv;

        } catch (IOException | SQLException | ClassNotFoundException | InterruptedException | JAXBException | XmlPullParserException | XMLStreamException ex) {
            Logger.getLogger(DataAnalysisHandler.class.getName()).log(Level.SEVERE, null, ex);
        }
        return null;

    }

    /**
     * The identification file reader factory of compomics utilities.
     */
    private static final IdfileReaderFactory readerFactory = IdfileReaderFactory.getInstance();

    public static void main(String[] args) {
        String dataFolderURL = "D:\\Apps\\OptProt\\pepshaker\\data_ps\\PXD028427";
        new DataAnalysisHandler(dataFolderURL, "PXD028427");
        String dataFolderURL2 = "D:\\Apps\\OptProt\\pepshaker\\data_ps\\PXD047036";
        new DataAnalysisHandler(dataFolderURL2, "PXD047036");
//        String dataFolderURL3 = "D:\\Apps\\OptProt\\pepshaker\\data_ps\\PXD001468";
//        new DataAnalysisHandler(dataFolderURL3, "PXD001468");
// System.out.println(2 + 2 + "Hello" + 2 + 2);
    }
}
