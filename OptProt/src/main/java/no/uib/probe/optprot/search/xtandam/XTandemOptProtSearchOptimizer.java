/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.search.xtandam;

import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
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
import com.compomics.util.parameters.identification.search.DigestionParameters;
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
import java.util.Set;
import java.util.concurrent.Future;
import javax.xml.bind.JAXBException;
import javax.xml.stream.XMLStreamException;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;
import no.uib.probe.optprot.model.OptimisedSearchResults;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.search.DefaultOptProtSearchOptimizer;
import no.uib.probe.optprot.search.SearchExecuter;
import no.uib.probe.optprot.util.MainUtilities;
import static no.uib.probe.optprot.util.MainUtilities.executor;
import no.uib.probe.optprot.util.OptProtWaitingHandler;
import no.uib.probe.optprot.util.ReportExporter;
import org.xmlpull.v1.XmlPullParserException;

/**
 *
 * @author yfa041
 */
public class XTandemOptProtSearchOptimizer extends DefaultOptProtSearchOptimizer {

    /**
     * The identification file reader factory of compomics utilities.
     */
    private final IdfileReaderFactory readerFactory = IdfileReaderFactory.getInstance();
    /**
     * The compomics PTM factory.
     */
    private final ModificationFactory ptmFactory = ModificationFactory.getInstance();
    
    private final SearchingSubDataset searchingSubDataset;
    private final SearchInputSetting searchInputSetting;
    private final File identificationParametersFile;
    private final OptimisedSearchResults optimisedSearchResults;
    private final IdentificationParameters identificationParameters;
    
    public XTandemOptProtSearchOptimizer(SearchingSubDataset searchingSubDataset, SearchInputSetting searchInputSetting, File identificationParametersFile) throws IOException {
        
        this.searchingSubDataset = searchingSubDataset;
        this.searchInputSetting = searchInputSetting;
        this.identificationParametersFile = identificationParametersFile;
        this.identificationParameters = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        this.optimisedSearchResults = new OptimisedSearchResults();
        XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        xtandemParameters.setQuickPyrolidone(false);
        xtandemParameters.setProteinQuickAcetyl(false);
        xtandemParameters.setStpBias(false);
        IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
    }
    
    public void startProcess() throws IOException {
//        System.out.println("at------------------------------------------------------------------------->>> started id rate " + searchingSubDataset.getIdentificationNum());
        advancedParam = false;
        if (searchInputSetting.isOptimizeDigestionParameter()) {
            String value = this.optimizeDigestionParameter(searchingSubDataset, identificationParametersFile, searchInputSetting);
//            System.out.println("optimizeDigestionParameter " + value + "------------------------------------------------------------------------->>> 1 id rate " + searchingSubDataset.getIdentificationNum());
            optimisedSearchResults.setDigestionParameter(value);
            identificationParameters.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf(value));
            if (value.equalsIgnoreCase("enzyme")) {
                identificationParameters.getSearchParameters().getDigestionParameters().clearEnzymes();
                identificationParameters.getSearchParameters().getDigestionParameters().addEnzyme(EnzymeFactory.getInstance().getEnzyme("Trypsin"));
            }
            IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
        }
        
        if (optimisedSearchResults.getDigestionParameter().equalsIgnoreCase(DigestionParameters.CleavageParameter.enzyme.name()) && (searchInputSetting.isOptimizeDigestionParameter() || searchInputSetting.isOptimizeEnzymeParameter())) {
            String value = this.optimizeEnzymeParameter(searchingSubDataset, identificationParametersFile, searchInputSetting);

//            System.out.println("optimizeEnzymeParameter " + value + "------------------------------------------------------------------------->>> 2 id rate " + searchingSubDataset.getIdentificationNum());
            optimisedSearchResults.setEnzymeName(value);
            identificationParameters.getSearchParameters().getDigestionParameters().clearEnzymes();
            identificationParameters.getSearchParameters().getDigestionParameters().addEnzyme(EnzymeFactory.getInstance().getEnzyme(value));
            IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
        }
        if (optimisedSearchResults.getDigestionParameter().equalsIgnoreCase(DigestionParameters.CleavageParameter.enzyme.name()) && (searchInputSetting.isOptimizeDigestionParameter() || searchInputSetting.isOptimizeSpecificityParameter())) {
            String value = this.optimizeSpecificityParameter(searchingSubDataset, identificationParametersFile, searchInputSetting);
            identificationParameters.getSearchParameters().getDigestionParameters().setSpecificity(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName(), DigestionParameters.Specificity.valueOf(value));
            IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//            System.out.println("optimizeSpecificityParameter " + value + "------------------------------------------------------------------------->>> 3 id rate " + searchingSubDataset.getIdentificationNum());
        }
        
        if (searchInputSetting.isOptimizeFragmentIonTypesParameter()) {
            String value = this.optimizeFragmentIonTypesParameter(searchingSubDataset, identificationParametersFile, searchInputSetting);
            int forward = Integer.parseInt(value.split("-")[0]);
            int rewind = Integer.parseInt(value.split("-")[1]);
            if (!identificationParameters.getSearchParameters().getForwardIons().contains(forward)) {
                ArrayList<Integer> forwardIonsList = new ArrayList<>();
                forwardIonsList.add(forward);
                identificationParameters.getSearchParameters().setForwardIons(forwardIonsList);
            }
            if (!identificationParameters.getSearchParameters().getRewindIons().contains(rewind)) {
                ArrayList<Integer> rewindIonsList = new ArrayList<>();
                rewindIonsList.add(rewind);
                identificationParameters.getSearchParameters().setRewindIons(rewindIonsList);
            }
            IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//            System.out.println("OptimizeFragmentIonTypesParameter" + value + "------------------------------------------------------------------------->>> 5 id rate " + searchingSubDataset.getIdentificationNum());
        }
        
        if (searchInputSetting.isOptimizeIsotopsParameter()) {
            //not supported in xtandem search engine
//            int[] values = this.optimizeIsotopParameter(searchingSubDataset, identificationParametersFile, searchInputSetting);
//            System.out.println("optimizeIsotopParameter : " + values[0]+"--"+values[1]);
//            identificationParameters.getSearchParameters().setMinIsotopicCorrection(values[0]);
//             identificationParameters.getSearchParameters().setMaxIsotopicCorrection(values[1]);
//            IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
        }
        
        if (searchInputSetting.isOptimizePrecursorToleranceParameter()) {
            double value = this.optimizePrecursorToleranceParameter(searchingSubDataset, identificationParametersFile, searchInputSetting);
            identificationParameters.getSearchParameters().setPrecursorAccuracy(value);
            IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//            System.out.println("OptimizePrecursorToleranceParameter" + value + "------------------------------------------------------------------------->>> 6 id rate " + searchingSubDataset.getIdentificationNum());
        }
//confusing param
        if ((optimisedSearchResults.getDigestionParameter().equalsIgnoreCase(DigestionParameters.CleavageParameter.enzyme.name()) && (searchInputSetting.isOptimizeDigestionParameter()) || searchInputSetting.isOptimizeMaxMissCleavagesParameter())) {
            int value = this.optimizeMaxMissCleavagesParameter(searchingSubDataset, identificationParametersFile, searchInputSetting);
            identificationParameters.getSearchParameters().getDigestionParameters().setnMissedCleavages(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName(), value);
            IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//            System.out.println("optimizeMaxMissCleavagesParameter " + value + "------------------------------------------------------------------------->>> 4 id rate " + searchingSubDataset.getIdentificationNum());
        }
        if (searchInputSetting.isOptimizeFragmentToleranceParameter()) {
            double value = this.optimizeFragmentToleranceParameter(searchingSubDataset, identificationParametersFile, searchInputSetting);
            identificationParameters.getSearchParameters().setFragmentIonAccuracy(value);
            IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//            System.out.println("optimizeFragmentToleranceParameter" + value + "------------------------------------------------------------------------->>> 7 id rate " + searchingSubDataset.getIdentificationNum());

        }
        
        if (searchInputSetting.isOptimizePrecursorChargeParameter()) {
            int[] values = this.optimizePrecursorChargeParameter(searchingSubDataset, identificationParametersFile, searchInputSetting);
            identificationParameters.getSearchParameters().setMinChargeSearched(values[0]);
            identificationParameters.getSearchParameters().setMaxChargeSearched(values[1]);
            IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//            System.out.println("optimizePrecursorChargeParameter" + values[0] + "--" + values[1] + "------------------------------------------------------------------------->>> 8 id rate " + searchingSubDataset.getIdentificationNum());

        }
        
        if (searchInputSetting.isOptimizeModificationParameter()) {
            //not supported in xtandem search engine
//            int[] values = 
//            this.enableQuik = false;
            Map<String, Set<String>> modificationsResults = this.optimizeModificationsParameter(searchingSubDataset, identificationParametersFile, searchInputSetting);
            
            identificationParameters.getSearchParameters().getModificationParameters().clearFixedModifications();
            identificationParameters.getSearchParameters().getModificationParameters().clearVariableModifications();
            identificationParameters.getSearchParameters().getModificationParameters().clearRefinementModifications();
            identificationParameters.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
            for (String fixedMod : modificationsResults.get("fixedModifications")) {
                identificationParameters.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(fixedMod));
                identificationParameters.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(fixedMod));
            }
            for (String variableMod : modificationsResults.get("variableModifications")) {
                identificationParameters.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableMod));
            }
            for (String refinmentVariableMod : modificationsResults.get("refinmentVariableModifications")) {
                identificationParameters.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(refinmentVariableMod));
            }
            IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//            System.out.println("OptimizeModificationParameter " + modificationsResults + "------------------------------------------------------------------------->>> 7 id rate " + searchingSubDataset.getIdentificationNum());

        }
        
        if (searchInputSetting.isOptimizeXtandemAdvancedParameter()) {
            advancedParam = true;
            XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            xtandemParameters.setProteinQuickAcetyl(false);
            xtandemParameters.setQuickPyrolidone(false);
            xtandemParameters.setStpBias(false);
//                XtandemAdvancedSearchOptimizerHandler xtendemHandler = new XtandemAdvancedSearchOptimizerHandler(subMsFile, subFastaFile, optimizedIdentificationParametersFile, optProtSearchParameters, referenceIdRate);
            if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptSpectrumDynamicRange()) {
                double value = optimizeSpectrumDynamicRange(searchingSubDataset, identificationParametersFile, searchInputSetting);
                xtandemParameters.setDynamicRange(value);
                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                
            }
            if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptSpectrumNumbPeaks()) {
                int value = optimizeSpectrumPeaksNumber(searchingSubDataset, identificationParametersFile, searchInputSetting);
                xtandemParameters.setnPeaks(value);
                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                
            }
            if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptSpectrumMinimumFragment()) {
                double value = optimizeMinimumFragmentMz(searchingSubDataset, identificationParametersFile, searchInputSetting);
                xtandemParameters.setMinFragmentMz(value);
                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                
            }
            if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptSpectrumMinimumPeaks()) {
                int value = optimizeMinimumPeaks(searchingSubDataset, identificationParametersFile, searchInputSetting);
                xtandemParameters.setMinPeaksPerSpectrum(value);
                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                
            }
            if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptNoiseSuppression()) {
                double minPrecursorMass = optimizeNoiseSuppression(searchingSubDataset, identificationParametersFile, searchInputSetting);
                xtandemParameters.setUseNoiseSuppression(minPrecursorMass == 0.0);
                if (minPrecursorMass != 0.0) {
                    xtandemParameters.setMinPrecursorMass(minPrecursorMass);
                }
                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                
            }
            if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptParentIsotopeExpansion()) {
                boolean value = optimizeParentIsotopExpansion(searchingSubDataset, identificationParametersFile, searchInputSetting);
                xtandemParameters.setParentMonoisotopicMassIsotopeError(value);
                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
            }
            if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptQuickAcetyl()) {
                boolean value = optimizeQuickAcetyl(searchingSubDataset, identificationParametersFile, searchInputSetting);
                xtandemParameters.setProteinQuickAcetyl(value);
                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
            }
            if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptQuickPyrolidone()) {
                boolean value = optimizeQuickPyrolidone(searchingSubDataset, identificationParametersFile, searchInputSetting);
                xtandemParameters.setQuickPyrolidone(value);
                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
            }
            if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptstPBias()) {
                boolean value = optimizeStPBias(searchingSubDataset, identificationParametersFile, searchInputSetting);
                xtandemParameters.setStpBias(value);
                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
            }
            
            if (xtandemParameters.isRefine() || searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptUseRefine()) {
                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptRefineUnanticipatedCleavage()) {
                    boolean value = optimizeRefineUnanticipatedCleavage(searchingSubDataset, identificationParametersFile, searchInputSetting);
                    xtandemParameters.setRefineUnanticipatedCleavages(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptRefineSimiEnzymaticCleavage()) {
                    boolean value = optimizeRefineSimiEnzymaticCleavage(searchingSubDataset, identificationParametersFile, searchInputSetting);
                    xtandemParameters.setRefineSemi(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptRefinePotintialModification()) {
                    boolean value = optimizePotintialModification(searchingSubDataset, identificationParametersFile, searchInputSetting);
                    xtandemParameters.setPotentialModificationsForFullRefinment(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptRefinePointMutations()) {
                    boolean value = optimizeRefinePointMutations(searchingSubDataset, identificationParametersFile, searchInputSetting);
                    xtandemParameters.setRefinePointMutations(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptRefineSnAPs()) {
                    boolean value = optimizeRefineSnAPs(searchingSubDataset, identificationParametersFile, searchInputSetting);
                    xtandemParameters.setRefineSnaps(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptRefineSpectrumSynthesis()) {
                    boolean value = optimizeRefineSpectrumSynthesis(searchingSubDataset, identificationParametersFile, searchInputSetting);
                    xtandemParameters.setRefineSpectrumSynthesis(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                
            }
            
            if (searchInputSetting.getXtandemOptProtAdvancedSearchParameters().isOptUseRefine()) {
                boolean value = optimizeUseRefine(searchingSubDataset, identificationParametersFile, searchInputSetting);
                xtandemParameters.setRefine(value);
                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
            }
            
        }
        ReportExporter.printFullReport(identificationParameters, searchingSubDataset, Advocate.xtandem);
        
    }
    private boolean advancedParam = true;
    
    @Override
    public synchronized ArrayList<SpectrumMatch> excuteSearch(SearchingSubDataset optProtDataset, String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addPeptideMasses, SearchInputSetting optProtSearchSettings, File identificationParametersFile) {
        try {
            if (!optProtSearchSettings.getXTandemEnabledParameters().getParamsToOptimize().isEnabledParam(paramOption.split("_")[0])) {
                System.out.println("param " + paramOption + " is not supported " + paramOption);
                return new ArrayList<>();
            }
            if (!advancedParam && tempIdParam.getSearchParameters().getModificationParameters().getFixedModifications().size() == 1) {
                boolean terminalMod = ptmFactory.getModification(tempIdParam.getSearchParameters().getModificationParameters().getFixedModifications().get(0)).getModificationType().isCTerm() || ptmFactory.getModification(tempIdParam.getSearchParameters().getModificationParameters().getFixedModifications().get(0)).getModificationType().isNTerm();
                if (terminalMod) {
                    System.out.println("terminal mod only ");
                    return new ArrayList<>();
                }
                
            }
            SearchParameters searchParameters = tempIdParam.getSearchParameters();
            if (!advancedParam) {
                
                XtandemParameters xtandemParameters = (XtandemParameters) searchParameters.getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
                xtandemParameters.setQuickPyrolidone(false);
                xtandemParameters.setProteinQuickAcetyl(false);
                xtandemParameters.setStpBias(false);
                xtandemParameters.setRefine(false);
            }
            
            File resultOutput = SearchExecuter.executeSearch(defaultOutputFileName, optProtSearchSettings, optProtDataset.getSubMsFile(), optProtDataset.getSubFastaFile(), tempIdParam, identificationParametersFile);
            
            File xTandemFile = new File(resultOutput, SearchHandler.getXTandemFileName(IoUtil.removeExtension(optProtDataset.getSubMsFile().getName())));
            IdfileReader idReader = readerFactory.getFileReader(xTandemFile);
            MsFileHandler msFileHandler = new MsFileHandler();
            msFileHandler.register(optProtDataset.getSubMsFile(), MainUtilities.OptProt_Waiting_Handler);
            final List<SpectrumMatch> matches = Collections.synchronizedList(idReader.getAllSpectrumMatches(msFileHandler, MainUtilities.OptProt_Waiting_Handler, searchParameters));
            
            if (addPeptideMasses) {
                SequenceMatchingParameters modificationSequenceMatchingParameters = tempIdParam.getModificationLocalizationParameters().getSequenceMatchingParameters();
                FMIndex sequenceProvider = new FMIndex(optProtDataset.getSubFastaFile(), null, new OptProtWaitingHandler(), false, tempIdParam);
                for (SpectrumMatch sm : matches) {
                    Thread t = new Thread(() -> {
                        for (PeptideAssumption pepAss : sm.getAllPeptideAssumptions().toList()) {
                            try {
                                Peptide pep = pepAss.getPeptide();
                                ModificationLocalizationMapper.modificationLocalization(
                                        pep,
                                        tempIdParam,
                                        idReader,
                                        ptmFactory,
                                        sequenceProvider
                                );
                                pepAss.getPeptide().setMass(pep.getMass(searchParameters.getModificationParameters(), sequenceProvider, modificationSequenceMatchingParameters));
                                sm.setBestPeptideAssumption(pepAss);
                            } catch (Exception e) {
                            }
                        }
                    });
                    t.start();
                    while (t.isAlive()) {
                        Thread.currentThread().sleep(10);
                    }
                    
                }
            }
            
            MainUtilities.deleteFolder(resultOutput);
            
            return new ArrayList<>(matches);
        } catch (IOException | InterruptedException | SQLException | ClassNotFoundException | JAXBException | XMLStreamException | XmlPullParserException ex) {
            ex.printStackTrace();
        }
        return new ArrayList<>();
    }
    
    public double optimizeSpectrumDynamicRange(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Double, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedOption = xtandemParameters.getDynamicRange();
        for (double i = 100.0; i > 0;) {
            final String option = "spectrumDR_" + i;
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setDynamicRange(i);
            final double j = i;
            
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
            i -= 20;
        }
        
        int localId = -1;
        double localSelection = 0;
        for (double dRangeScore : resultsMap.keySet()) {
//            System.out.println(" option 1 " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
        }
        return selectedOption;
    }
    
    public int optimizeSpectrumPeaksNumber(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        Integer selectedOption = xtandemParameters.getnPeaks();
        int idRate = optProtDataset.getIdentificationNum();
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        
        for (int i = 30; i <= 70;) {
            final String option = "peaksNum_" + i;
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setnPeaks(i);
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
            i += 10;
        }
        int localId = -1;
        int localSelection = 0;
        for (int option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
        }
        return selectedOption;
    }
    
    public double optimizeMinimumFragmentMz(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Double, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedOption = xtandemParameters.getMinFragmentMz();
        for (int i = 100; i <= 300;) {
            final String option = "minimumFragmentMz_" + i;
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setMinFragmentMz(i);
            final double j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
            i += 50;
        }
        int localId = -1;
        double localSelection = 0;
        for (double dRangeScore : resultsMap.keySet()) {
//            System.out.println(" option 1 " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
    }
//

    public int optimizeMinimumPeaks(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        Integer selectedOption = xtandemParameters.getMinPeaksPerSpectrum();
        int idRate = optProtDataset.getIdentificationNum();
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int lastValue = 0;
        for (int i = 5; i <= 100;) {
            final String option = "minpeaksNum_" + i;
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setMinPeaksPerSpectrum(i);
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
            if (lastValue > resultsMap.get(i)) {
                break;
            }
            lastValue = resultsMap.get(i);
            i += 10;
            
        }
        int localId = -1;
        int localSelection = 0;
        for (int option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
        }
        return selectedOption;
    }
//

    public double optimizeNoiseSuppression(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Double, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedOption = xtandemParameters.getMinPrecursorMass();
        final String option = "noiseSupression_" + false;
        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
        xtandemParameters.setUseNoiseSuppression(false);
        Future future = executor.submit(() -> {
            resultsMap.put(0.0, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
        });
        while (!future.isDone()) {
        }
        xtandemParameters.setUseNoiseSuppression(true);
        for (double j = 500; j < 1600;) {
            final String suboption = "noiseSupression_" + true + "_" + j;
            final String subupdatedName = Configurations.DEFAULT_RESULT_NAME + "_" + suboption + "_" + msFileName;
            final double i = j;
            xtandemParameters.setMinPrecursorMass(j);
            future = executor.submit(() -> {
                resultsMap.put(i, excuteSearch(optProtDataset, subupdatedName, suboption, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
            j += 350;
            
        }
        
        int localId = -1;
        double localSelection = 0;
        for (double dRangeScore : resultsMap.keySet()) {
//            System.out.println(" option 1 " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
    }
//

    public boolean optimizeParentIsotopExpansion(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.getParentMonoisotopicMassIsotopeError();
        
        for (int i = 0; i < 2; i++) {
            final String option = "parentMonoisotopicMassIsotopeError_" + (i == 1);
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setParentMonoisotopicMassIsotopeError(i == 1);
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                
            });
            while (!future.isDone()) {
            }
        }
        
        int localId = -1;
        int localSelection = 0;
        for (int dRangeScore : resultsMap.keySet()) {
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection == 1;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
    }
//

    public boolean optimizeQuickAcetyl(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isProteinQuickAcetyl();
        
        for (int i = 0; i < 2; i++) {
            boolean useQuickAcetyl = (i == 1);
            final String option = "useQuickAcetyl_" + useQuickAcetyl;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setProteinQuickAcetyl(useQuickAcetyl);
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                
            });
            while (!future.isDone()) {
            }
            System.out.println("at quick acetyle " + useQuickAcetyl + " " + resultsMap.get(j));
        }
        
        int localId = -1;
        int localSelection = 0;
        for (int dRangeScore : resultsMap.keySet()) {
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection == 1;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
    }
    
    public boolean optimizeQuickPyrolidone(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isQuickPyrolidone();
        
        for (int i = 0; i < 2; i++) {
            boolean useQuickPyrolidone = (i == 1);
            final String option = "useQuickPyrolidone_" + useQuickPyrolidone;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setQuickPyrolidone(useQuickPyrolidone);
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                
            });
            while (!future.isDone()) {
            }
            System.out.println("at quick prolyien  " + useQuickPyrolidone + " " + resultsMap.get(j));
        }
        
        int localId = -1;
        int localSelection = 0;
        for (int dRangeScore : resultsMap.keySet()) {
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection == 1;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
    }
//

    public boolean optimizeStPBias(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isStpBias();
        
        for (int i = 0; i < 2; i++) {
            boolean useStpBias = (i == 1);
            final String option = "useStpBias_" + useStpBias;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setStpBias(useStpBias);
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                
            });
            while (!future.isDone()) {
            }
            
            System.out.println("at useStpBias  " + useStpBias + " " + resultsMap.get(j));
        }
        
        int localId = -1;
        int localSelection = 0;
        for (int dRangeScore : resultsMap.keySet()) {
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection == 1;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
    }
//

    public boolean optimizeUseRefine(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefine();
        
        for (int i = 0; i < 2; i++) {
            boolean useRefine = (i == 1);
            final String option = "useRefine_" + useRefine;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setRefine(useRefine);
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                
            });
            while (!future.isDone()) {
            }
        }
        
        int localId = -1;
        int localSelection = 0;
        for (int dRangeScore : resultsMap.keySet()) {
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection == 1;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
        
    }
//

    public boolean optimizeRefineUnanticipatedCleavage(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefineUnanticipatedCleavages();
        
        for (int i = 0; i < 2; i++) {
            boolean useRefineUnanticipatedCleavages = (i == 1);
            final String option = "useRefineUnanticipatedCleavages_" + useRefineUnanticipatedCleavages;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setRefineUnanticipatedCleavages(useRefineUnanticipatedCleavages);
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                
            });
            while (!future.isDone()) {
            }
        }
        
        int localId = -1;
        int localSelection = 0;
        for (int dRangeScore : resultsMap.keySet()) {
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection == 1;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
        
    }
//

    public boolean optimizeRefineSimiEnzymaticCleavage(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefineSemi();
        
        for (int i = 0; i < 2; i++) {
            boolean useRefineSimiEnzymaticCleavage = (i == 1);
            final String option = "useRefineSimiEnzymaticCleavage_" + useRefineSimiEnzymaticCleavage;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setRefineSemi(useRefineSimiEnzymaticCleavage);
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                
            });
            while (!future.isDone()) {
            }
        }
        
        int localId = -1;
        int localSelection = 0;
        for (int dRangeScore : resultsMap.keySet()) {
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection == 1;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
        
    }
//

    public boolean optimizePotintialModification(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isPotentialModificationsForFullRefinment();
        
        for (int i = 0; i < 2; i++) {
            boolean usePotintialModification = (i == 1);
            final String option = "usePotintialModification_" + usePotintialModification;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setPotentialModificationsForFullRefinment(usePotintialModification);
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                
            });
            while (!future.isDone()) {
            }
        }
        
        int localId = -1;
        int localSelection = 0;
        for (int dRangeScore : resultsMap.keySet()) {
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection == 1;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
        
    }
//

    public boolean optimizeRefinePointMutations(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefinePointMutations();
        
        for (int i = 0; i < 2; i++) {
            boolean useRefinePointMutations = (i == 1);
            final String option = "useRefinePointMutations_" + useRefinePointMutations;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setRefinePointMutations(useRefinePointMutations);
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                
            });
            while (!future.isDone()) {
            }
        }
        
        int localId = -1;
        int localSelection = 0;
        for (int dRangeScore : resultsMap.keySet()) {
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection == 1;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
        
    }
//
//

    public boolean optimizeRefineSnAPs(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefineSnaps();
        
        for (int i = 0; i < 2; i++) {
            boolean useRefineSnAPs = (i == 1);
            final String option = "useRefineSnAPs_" + useRefineSnAPs;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setRefineSnaps(useRefineSnAPs);
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                
            });
            while (!future.isDone()) {
            }
        }
        
        int localId = -1;
        int localSelection = 0;
        for (int dRangeScore : resultsMap.keySet()) {
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection == 1;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
        
    }
    
    public boolean optimizeRefineSpectrumSynthesis(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefineSpectrumSynthesis();
        
        for (int i = 0; i < 2; i++) {
            boolean useRefineSpectrumSynthesis = (i == 1);
            final String option = "useRefineSpectrumSynthesis_" + useRefineSpectrumSynthesis;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setRefineSpectrumSynthesis(useRefineSpectrumSynthesis);
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                
            });
            while (!future.isDone()) {
            }
        }
        
        int localId = -1;
        int localSelection = 0;
        for (int dRangeScore : resultsMap.keySet()) {
            if (resultsMap.get(dRangeScore) > localId) {
                localId = resultsMap.get(dRangeScore);
                localSelection = dRangeScore;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection == 1;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
        
    }
    
}
