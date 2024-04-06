/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.search;

import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.biology.proteins.Peptide;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.identification.modification.search_engine_mapping.ModificationLocalizationMapper;
import com.compomics.util.experiment.identification.protein_inference.fm_index.FMIndex;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
import com.compomics.util.experiment.io.identification.IdfileReader;
import com.compomics.util.experiment.io.identification.IdfileReaderFactory;
import com.compomics.util.experiment.io.mass_spectrometry.MsFileHandler;
import com.compomics.util.io.IoUtil;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.advanced.SequenceMatchingParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.parameters.identification.tool_specific.XtandemParameters;
import eu.isas.searchgui.SearchHandler;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.model.OptProtSearchParameters;
import no.uib.probe.optprot.util.MainUtilities;
import no.uib.probe.optprot.util.OptProtWaitingHandler;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author yfa041
 */
public class XtandemAdvancedSearchOptimizerHandler {

    private ExecutorService executor;
    private final File identificationParametersFile;
    private final IdentificationParameters identificationParam;
//    private double referenceIdRate;
    private final XtandemParameters xtandemParameters;
    private final File subMsFile, subFastaFile;
    private final OptProtSearchParameters searchEngineParameters;
    private final SearchParameters searchParameters;
    /**
     * The identification file reader factory of compomics utilities.
     */
    private final IdfileReaderFactory readerFactory = IdfileReaderFactory.getInstance();
    /**
     * The compomics PTM factory.
     */
    private final ModificationFactory ptmFactory = ModificationFactory.getInstance();

    public XtandemAdvancedSearchOptimizerHandler(File subMsFile, File subFastaFile, File optimizedIdentificationParametersFile, OptProtSearchParameters searchEngineParameters, double referenceIdRate) throws IOException {
        executor = Executors.newFixedThreadPool(2);
        this.identificationParametersFile = optimizedIdentificationParametersFile;
//        this.referenceIdRate = referenceIdRate;
        identificationParam = IdentificationParameters.getIdentificationParameters(optimizedIdentificationParametersFile);
        searchParameters = identificationParam.getSearchParameters();
        xtandemParameters = (XtandemParameters) searchParameters.getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        xtandemParameters.setQuickPyrolidone(true);
        xtandemParameters.setProteinQuickAcetyl(true);
        xtandemParameters.setStpBias(true);
        xtandemParameters.setRefine(true);

        this.subFastaFile = subFastaFile;
        this.subMsFile = subMsFile;
        this.searchEngineParameters = searchEngineParameters;
    }

    public double optimizeSpectrumDynamicRange() throws IOException {
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        for (double i = 100.0; i > 0;) {
            final String option = "spectrumDR_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setDynamicRange(i);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());
            });
            while (!future.isDone()) {
            }
            i -= 20;
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;

            }
        }
        xtandemParameters.setDynamicRange(100.0);
        return Double.parseDouble(selectedOption.split("_")[1]);
    }

    public int optimizeSpectrumPeaksNumber() throws IOException {
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        for (int i = 30; i <= 70;) {
            final String option = "peaksNum_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setnPeaks(i);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());
            });
            while (!future.isDone()) {
            }
            i += 10;
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;

            }
        }
        xtandemParameters.setnPeaks(50);
        return Integer.parseInt(selectedOption.split("_")[1]);
    }

    public double optimizeMinimumFragmentMz() throws IOException {
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        for (int i = 100; i <= 300;) {
            final String option = "minimumFragmentMz_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setMinFragmentMz(i);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());
            });
            while (!future.isDone()) {
            }
            i += 50;
        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;

            }
        }
        xtandemParameters.setMinFragmentMz(Integer.parseInt(selectedOption.split("_")[1]));
        return Integer.parseInt(selectedOption.split("_")[1]);
    }

    public int optimizeMinimumPeaks() throws IOException {
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        int lastValue = 0;
        for (int i = 5; i <= 100;) {
            final String option = "minpeaksNum_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setMinPeaksPerSpectrum(i);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());
            });
            while (!future.isDone()) {
            }
            if (lastValue > resultsMap.get(option)) {
                break;
            }
            lastValue = resultsMap.get(option);
            i += 10;

        }
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) >= idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;

            }
        }
        xtandemParameters.setMinPeaksPerSpectrum(Integer.parseInt(selectedOption.split("_")[1]));
        return Integer.parseInt(selectedOption.split("_")[1]);
    }

    public double optimizeNoiseSuppression() throws IOException {
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        final String option = "noiseSupression_" + false;
        double minPrecursorMass = -1;
        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
        xtandemParameters.setUseNoiseSuppression(false);
        Future future = executor.submit(() -> {
            resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());

        });
        while (!future.isDone()) {
        }
        xtandemParameters.setUseNoiseSuppression(true);
        for (double j = 500; j < 1600;) {
            final String suboption = "noiseSupression_" + true + "_" + j;
            final String subupdatedName = Configurations.DEFAULT_RESULT_NAME + "_" + suboption;
            xtandemParameters.setMinPrecursorMass(j);
            future = executor.submit(() -> {
                resultsMap.put(suboption, excuteXTandomSearches(subupdatedName, false).size());

            });
            while (!future.isDone()) {
            }
            j += 350;

        }

        for (String availableOption : resultsMap.keySet()) {
            if (resultsMap.get(availableOption) > idRate && availableOption.split("_").length == 3) {
                minPrecursorMass = Double.parseDouble(availableOption.split("_")[2]);
                idRate = resultsMap.get(availableOption);
            }
        }
        xtandemParameters.setMinPrecursorMass(minPrecursorMass);
        return minPrecursorMass;// Boolean.parseBoolean(selectedOption.split("_")[1]);
    }

    public boolean optimizeParentIsotopExpansion() throws IOException {
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        boolean useParentIsotopExpansion = false;
        for (int i = 0; i < 2; i++) {
            final String option = useParentIsotopExpansion + "";

            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setParentMonoisotopicMassIsotopeError(useParentIsotopExpansion);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());

            });
            while (!future.isDone()) {
            }
            useParentIsotopExpansion = true;
        }

        for (String availableOption : resultsMap.keySet()) {
            if (resultsMap.get(availableOption) > idRate) {
                idRate = resultsMap.get(availableOption);
                selectedOption = availableOption;

            }
        }
        xtandemParameters.setParentMonoisotopicMassIsotopeError(Boolean.parseBoolean(selectedOption));
        return Boolean.parseBoolean(selectedOption);
    }

    public boolean optimizeQuickAcetyl() throws IOException {
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        boolean useQuickAcetyl = false;
        for (int i = 0; i < 2; i++) {
            final String option = useQuickAcetyl + "";

            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setProteinQuickAcetyl(useQuickAcetyl);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());

            });
            while (!future.isDone()) {
            }
            useQuickAcetyl = true;
        }

        for (String availableOption : resultsMap.keySet()) {
            if (resultsMap.get(availableOption) > idRate) {
                idRate = resultsMap.get(availableOption);
                selectedOption = availableOption;

            }
        }
        xtandemParameters.setProteinQuickAcetyl(Boolean.parseBoolean(selectedOption));
        return Boolean.parseBoolean(selectedOption);
    }

    public boolean optimizeQuickPyrolidone() throws IOException {
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        boolean useQuickPyrolidone = false;
        for (int i = 0; i < 2; i++) {
            final String option = useQuickPyrolidone + "";

            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setQuickPyrolidone(useQuickPyrolidone);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());

            });
            while (!future.isDone()) {
            }
            useQuickPyrolidone = true;
        }

        for (String availableOption : resultsMap.keySet()) {
            if (resultsMap.get(availableOption) > idRate) {
                idRate = resultsMap.get(availableOption);
                selectedOption = availableOption;

            }
        }
        xtandemParameters.setQuickPyrolidone(Boolean.parseBoolean(selectedOption));
        return Boolean.parseBoolean(selectedOption);
    }

    public boolean optimizeStPBias() throws IOException {
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        boolean useStPBias = false;
        for (int i = 0; i < 2; i++) {
            final String option = useStPBias + "";

            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setStpBias(useStPBias);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());

            });
            while (!future.isDone()) {
            }
            useStPBias = true;
        }

        for (String availableOption : resultsMap.keySet()) {
            if (resultsMap.get(availableOption) > idRate) {
                idRate = resultsMap.get(availableOption);
                selectedOption = availableOption;

            }
        }
        xtandemParameters.setStpBias(Boolean.parseBoolean(selectedOption));
        return Boolean.parseBoolean(selectedOption);
    }

    public boolean optimizeUseRefine() throws IOException {
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        boolean useRefine = false;
        for (int i = 0; i < 2; i++) {
            final String option = useRefine + "";

            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setRefine(useRefine);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());

            });
            while (!future.isDone()) {
            }
            useRefine = true;
        }

        for (String availableOption : resultsMap.keySet()) {
            if (resultsMap.get(availableOption) > idRate) {
                idRate = resultsMap.get(availableOption);
                selectedOption = availableOption;

            }
        }
        xtandemParameters.setStpBias(Boolean.parseBoolean(selectedOption));
        return Boolean.parseBoolean(selectedOption);
    }

    public boolean optimizeRefineUnanticipatedCleavage() throws IOException {
        if (!xtandemParameters.isRefine()) {
            return false;
        }
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        boolean useUnanticipatedCleavage = false;
        for (int i = 0; i < 2; i++) {
            final String option = useUnanticipatedCleavage + "";

            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setRefineUnanticipatedCleavages(useUnanticipatedCleavage);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());

            });
            while (!future.isDone()) {
            }
            useUnanticipatedCleavage = true;
        }

        for (String availableOption : resultsMap.keySet()) {
            if (resultsMap.get(availableOption) > idRate) {
                idRate = resultsMap.get(availableOption);
                selectedOption = availableOption;

            }
        }
        xtandemParameters.setRefineUnanticipatedCleavages(Boolean.parseBoolean(selectedOption));
        return Boolean.parseBoolean(selectedOption);
    }

    public boolean optimizeRefineSimiEnzymaticCleavage() throws IOException {
        if (!xtandemParameters.isRefine()) {
            return false;
        }
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        boolean simiEnzymaticCleavage = false;
        for (int i = 0; i < 2; i++) {
            final String option = simiEnzymaticCleavage + "";

            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setRefineSemi(simiEnzymaticCleavage);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());

            });
            while (!future.isDone()) {
            }
            simiEnzymaticCleavage = true;
        }

        for (String availableOption : resultsMap.keySet()) {
            if (resultsMap.get(availableOption) > idRate) {
                idRate = resultsMap.get(availableOption);
                selectedOption = availableOption;

            }
        }
        xtandemParameters.setRefineSemi( Boolean.parseBoolean(selectedOption));
        return Boolean.parseBoolean(selectedOption);
    }

    public boolean optimizeRefinePotintialModification() throws IOException {
        if (!xtandemParameters.isRefine()) {
            return false;
        }
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        boolean useRefinePotintialModification = false;
        for (int i = 0; i < 2; i++) {
            final String option = useRefinePotintialModification + "";

            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setPotentialModificationsForFullRefinment(useRefinePotintialModification);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());

            });
            while (!future.isDone()) {
            }
            useRefinePotintialModification = true;
        }

        for (String availableOption : resultsMap.keySet()) {
            if (resultsMap.get(availableOption) > idRate) {
                idRate = resultsMap.get(availableOption);
                selectedOption = availableOption;

            }
        }
        xtandemParameters.setPotentialModificationsForFullRefinment( Boolean.parseBoolean(selectedOption));
        return Boolean.parseBoolean(selectedOption);
    }

    public boolean optimizeRefinePointMutations() throws IOException {
        if (!xtandemParameters.isRefine()) {
            return false;
        }
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        boolean useUnanticipatedCleavage = false;
        for (int i = 0; i < 2; i++) {
            final String option = useUnanticipatedCleavage + "";

            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setRefinePointMutations(useUnanticipatedCleavage);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());

            });
            while (!future.isDone()) {
            }
            useUnanticipatedCleavage = true;
        }

        for (String availableOption : resultsMap.keySet()) {
            if (resultsMap.get(availableOption) > idRate) {
                idRate = resultsMap.get(availableOption);
                selectedOption = availableOption;

            }
        }
        xtandemParameters.setRefinePointMutations( Boolean.parseBoolean(selectedOption));
        return Boolean.parseBoolean(selectedOption);
    }

    public boolean optimizeRefineSnAPs() throws IOException {
        if (!xtandemParameters.isRefine()) {
            return false;
        }
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        boolean useRefineSnaps = false;
        for (int i = 0; i < 2; i++) {
            final String option = useRefineSnaps + "";

            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setRefineSnaps(useRefineSnaps);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());

            });
            while (!future.isDone()) {
            }
            useRefineSnaps = true;
        }

        for (String availableOption : resultsMap.keySet()) {
            if (resultsMap.get(availableOption) > idRate) {
                idRate = resultsMap.get(availableOption);
                selectedOption = availableOption;

            }
        }
        xtandemParameters.setRefineSnaps( Boolean.parseBoolean(selectedOption));
        return Boolean.parseBoolean(selectedOption);
    }

    public boolean optimizeRefineSpectrumSynthesis() throws IOException {
        if (!xtandemParameters.isRefine()) {
            return false;
        }
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double idRate = -1;
        String selectedOption = "";
        boolean useSpectrumSynthesis = false;
        for (int i = 0; i < 2; i++) {
            final String option = useSpectrumSynthesis + "";

            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option;
            xtandemParameters.setRefineSpectrumSynthesis(useSpectrumSynthesis);
            Future future = executor.submit(() -> {
                resultsMap.put(option, excuteXTandomSearches(updatedName, false).size());

            });
            while (!future.isDone()) {
            }
            useSpectrumSynthesis = true;
        }

        for (String availableOption : resultsMap.keySet()) {
            if (resultsMap.get(availableOption) > idRate) {
                idRate = resultsMap.get(availableOption);
                selectedOption = availableOption;

            }
        }
        xtandemParameters.setRefineSpectrumSynthesis( Boolean.parseBoolean(selectedOption));
        return Boolean.parseBoolean(selectedOption);
    }

    private synchronized ArrayList<SpectrumMatch> excuteXTandomSearches(String defaultOutputFileName, boolean addPeptideMasses) {

        try {
            File resultOutput = SearchExecuter.executeSearch(defaultOutputFileName, searchEngineParameters, subMsFile, subFastaFile, identificationParam, identificationParametersFile);
            File xTandemFile = new File(resultOutput, SearchHandler.getXTandemFileName(IoUtil.removeExtension(subMsFile.getName())));
            IdfileReader idReader = readerFactory.getFileReader(xTandemFile);
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(subMsFile, MainUtilities.OptProt_Waiting_Handler);
            final List<SpectrumMatch> matches = Collections.synchronizedList(idReader.getAllSpectrumMatches(msFileHandler, MainUtilities.OptProt_Waiting_Handler, searchParameters));

            if (addPeptideMasses) {
                SequenceMatchingParameters modificationSequenceMatchingParameters = identificationParam.getModificationLocalizationParameters().getSequenceMatchingParameters();
                FMIndex sequenceProvider = new FMIndex(subFastaFile, null, new OptProtWaitingHandler(), false, identificationParam);
                executor = Executors.newFixedThreadPool(2);
                Future future = executor.submit(() -> {
                    for (SpectrumMatch sm : matches) {
                        for (PeptideAssumption pepAss : sm.getAllPeptideAssumptions().toList()) {
                            Peptide pep = pepAss.getPeptide();
                            ModificationLocalizationMapper.modificationLocalization(
                                    pep,
                                    identificationParam,
                                    idReader,
                                    ptmFactory,
                                    sequenceProvider
                            );
                            pepAss.getPeptide().setMass(pep.getMass(searchParameters.getModificationParameters(), sequenceProvider, modificationSequenceMatchingParameters));
                            sm.setBestPeptideAssumption(pepAss);
                        }
                    }
                });
                while (!future.isDone()) {
                }
//                executor.shutdown();
            }
            return new ArrayList<>(matches);
        } catch (IOException | InterruptedException | SQLException | ClassNotFoundException | JAXBException | XMLStreamException | XmlPullParserException ex) {
            ex.printStackTrace();
        }
        return new ArrayList<>();
    }

}
