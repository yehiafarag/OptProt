package no.uib.probe.optprot.search.sage;

import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.io.identification.IdfileReaderFactory;
import com.compomics.util.io.IoUtil;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
//import com.compomics.util.parameters.identification.tool_specific.MyriMatchParameters;
import com.compomics.util.parameters.identification.tool_specific.SageParameters;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;
import no.uib.probe.optprot.model.OptimisedSearchResults;
import no.uib.probe.optprot.model.ParameterScoreModel;
import no.uib.probe.optprot.model.RawScoreModel;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.search.DefaultOptProtSearchOptimizer;
import no.uib.probe.optprot.search.SearchExecuter;
import no.uib.probe.optprot.util.MainUtilities;
import no.uib.probe.optprot.util.SpectraUtilities;

/**
 *
 * @author yfa041
 */
public class SageOptProtSearchOptimizer extends DefaultOptProtSearchOptimizer {

    /**
     * The identification file reader factory of compomics utilities.
     */
    private final IdfileReaderFactory readerFactory = IdfileReaderFactory.getInstance();
    /**
     * The compomics PTM factory.
     */
    private final ModificationFactory ptmFactory = ModificationFactory.getInstance();

    private final SearchingSubDataset optProtDataset;
    private final SearchInputSetting searchInputSetting;
    private final File identificationParametersFile;
    private final OptimisedSearchResults optimisedSearchResults;
    private final IdentificationParameters identificationParameters;
    private final Map<String, TreeSet<ParameterScoreModel>> parameterScoreMap;

    public SageOptProtSearchOptimizer(SearchingSubDataset optProtDataset, SearchInputSetting searchInputSetting, File identificationParametersFile) throws IOException {

        this.optProtDataset = optProtDataset;
        this.searchInputSetting = searchInputSetting;
        this.identificationParametersFile = identificationParametersFile;
        this.identificationParameters = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameters = (SageParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        sageParameters.setMaxVariableMods(4);
        sageParameters.setNumPsmsPerSpectrum(1);
        sageParameters.setGenerateDecoys(false);
        this.optimisedSearchResults = new OptimisedSearchResults();
        this.parameterScoreMap = new LinkedHashMap<>();
        MainUtilities.cleanOutputFolder();
        parameterScoreMap.put("DigestionParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("EnzymeParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SpecificityParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("MaxMissCleavagesParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("FragmentIonTypesParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("PrecursorToleranceParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("FragmentToleranceParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("PrecursorChargeParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("IsotopParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("ModificationsParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SagePeptideLengthParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SagePeptideMassParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageBucketSizeParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageFragmentMzParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageIonMinIndexParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageMaxVariableModificationParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageGenerateDecoyParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageDeisotopParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageChimericSpectraParameter", new TreeSet<>(Collections.reverseOrder()));
//
        parameterScoreMap.put("SageWideWindowParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SagePredectRetentionTimeParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageNumberOfPeakParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageMinMatchedPeaksPeakParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageMaxFragmentChargeParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageNumberOfPeakParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageBatchSizeParameter", new TreeSet<>(Collections.reverseOrder()));

    }

    private String digestionParameterOpt;

    public void startProcess(List<String> paramOrder) throws IOException {

        digestionParameterOpt = identificationParameters.getSearchParameters().getDigestionParameters().getCleavageParameter().name();
        for (String param : paramOrder) {
            System.out.println("-------------------------------------------param " + param + "-------------------------------------------");
            if (param.equalsIgnoreCase("DigestionParameter_1") && searchInputSetting.isOptimizeDigestionParameter()) {
                optimisedSearchResults.setDigestionParameter("enzyme");
                String[] values = this.optimizeEnzymeParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("EnzymeParameter"));
                optimisedSearchResults.setEnzymeName(values[0]);
                if (!values[0].equalsIgnoreCase(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName())) {
                    int nMissesCleavages = identificationParameters.getSearchParameters().getDigestionParameters().getnMissedCleavages(values[0]);
                    identificationParameters.getSearchParameters().getDigestionParameters().clearEnzymes();
                    identificationParameters.getSearchParameters().getDigestionParameters().addEnzyme(EnzymeFactory.getInstance().getEnzyme(values[0]));
                    identificationParameters.getSearchParameters().getDigestionParameters().setnMissedCleavages(values[0], nMissesCleavages);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                values[0] = this.optimizeSpecificityParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SpecificityParameter"));
                if (!values[0].equalsIgnoreCase(identificationParameters.getSearchParameters().getDigestionParameters().getSpecificity(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName()).name())) {
                    identificationParameters.getSearchParameters().getDigestionParameters().setSpecificity(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName(), DigestionParameters.Specificity.valueOf(values[0]));
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                if (optProtDataset.getIdentificationRate() < 10) {
                    digestionParameterOpt = this.optimizeDigestionCleavageParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("DigestionParameter"));
                }

                continue;

            }

            if (param.equalsIgnoreCase("FragmentIonTypesParameter") && searchInputSetting.isOptimizeFragmentIonTypesParameter()) {
                String value = this.optimizeFragmentIonTypesParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("FragmentIonTypesParameter"));
                int forward = Integer.parseInt(value.split("-")[0]);
                int rewind = Integer.parseInt(value.split("-")[1]);
                boolean update = false;
                if (!identificationParameters.getSearchParameters().getForwardIons().contains(forward)) {
                    ArrayList<Integer> forwardIonsList = new ArrayList<>();
                    forwardIonsList.add(forward);
                    identificationParameters.getSearchParameters().setForwardIons(forwardIonsList);
                    update = true;
                }
                if (!identificationParameters.getSearchParameters().getRewindIons().contains(rewind)) {
                    ArrayList<Integer> rewindIonsList = new ArrayList<>();
                    rewindIonsList.add(rewind);
                    identificationParameters.getSearchParameters().setRewindIons(rewindIonsList);
                    update = true;
                }
                if (update) {
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                continue;

            }

//confusing param
            if (param.equalsIgnoreCase("DigestionParameter_2") && searchInputSetting.isOptimizeDigestionParameter()) {
                int value = this.optimizeMaxMissCleavagesParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("MaxMissCleavagesParameter"));
                if (value != identificationParameters.getSearchParameters().getDigestionParameters().getnMissedCleavages(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName())) {
                    identificationParameters.getSearchParameters().getDigestionParameters().setnMissedCleavages(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName(), value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                continue;

            }

            if (param.equalsIgnoreCase("FragmentToleranceParameter") && searchInputSetting.isOptimizeFragmentToleranceParameter()) {
                double value = this.optimizeFragmentToleranceParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("FragmentToleranceParameter"));
                if (value != identificationParameters.getSearchParameters().getFragmentIonAccuracy()) {
                    identificationParameters.getSearchParameters().setFragmentIonAccuracy(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                continue;
            }

            if (param.equalsIgnoreCase("PrecursorChargeParameter") && searchInputSetting.isOptimizePrecursorChargeParameter()) {

                int[] values = this.optimizePrecursorChargeParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("PrecursorChargeParameter"));
                if (values[1] != identificationParameters.getSearchParameters().getMaxChargeSearched()) {
                    identificationParameters.getSearchParameters().setMinChargeSearched(values[0]);
                    identificationParameters.getSearchParameters().setMaxChargeSearched(values[1]);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                continue;
            }

            if (param.equalsIgnoreCase("ModificationParameter") && searchInputSetting.isOptimizeModificationParameter()) {
                Map<String, Set<String>> modificationsResults = this.optimizeModificationsParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("ModificationsParameter"));
                identificationParameters.getSearchParameters().getModificationParameters().clearFixedModifications();
                identificationParameters.getSearchParameters().getModificationParameters().clearVariableModifications();
                identificationParameters.getSearchParameters().getModificationParameters().clearRefinementModifications();
                identificationParameters.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
                for (String fixedMod : modificationsResults.get("fixedModifications")) {
                    if (ptmFactory.getModification(fixedMod) != null) {
                        identificationParameters.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(fixedMod));
                        identificationParameters.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(fixedMod));
                    }
                }
                for (String variableMod : modificationsResults.get("variableModifications")) {
                    if (ptmFactory.getModification(variableMod) != null) {
                        identificationParameters.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableMod));
                    }
                }
                SageParameters sageParameters = (SageParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
                sageParameters.setMaxVariableMods(modificationsResults.get("variableModifications").size());
                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                continue;
            }

            if (param.equalsIgnoreCase("SageAdvancedParameter") && searchInputSetting.isOptimizeSageAdvancedParameter()) {
                SageParameters sageParameters = (SageParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());

                int[] values = optimizePeptideLengthParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SagePeptideLengthParameter"));
                if (values[1] != sageParameters.getMaxPeptideLength() || values[0] != sageParameters.getMinPeptideLength()) {
                    sageParameters.setMinPeptideLength(values[0]);
                    sageParameters.setMaxPeptideLength(values[1]);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                    System.out.println("peptide length " + values[0] + "--" + values[1] + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
                }

                double[] dvalues = optimizeFragmentMzParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SageFragmentMzParameter"));
                if (dvalues[1] != sageParameters.getMaxFragmentMz() || dvalues[0] != sageParameters.getMinFragmentMz()) {
                    sageParameters.setMinFragmentMz(dvalues[0]);
                    sageParameters.setMaxFragmentMz(dvalues[1]);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                    System.out.println("optimizeFragmentMzParameter " + dvalues[0] + "--" + dvalues[1] + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
                }
                dvalues = optimizePeptideMassParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SagePeptideMassParameter"));
                if (dvalues[1] != sageParameters.getMaxPeptideMass() || dvalues[0] != sageParameters.getMinPeptideMass()) {
                    sageParameters.setMinPeptideMass(dvalues[0]);
                    sageParameters.setMaxPeptideMass(dvalues[1]);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                    System.out.println("optimizePeptideMassParameter " + dvalues[0] + "--" + dvalues[1] + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
                }

                int value = optimizeIonMinIndexParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SageIonMinIndexParameter"));
                if (value != sageParameters.getMinIonIndex()) {
                    sageParameters.setMinIonIndex(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);

                }
                System.out.println("optimizeIonMinIndexParameter " + value + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());

                value = optimizeMaxVariableModificationParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SageMaxVariableModificationParameter"));
                if (value != sageParameters.getMaxVariableMods()) {
                    sageParameters.setMaxVariableMods(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);

                }
                System.out.println("optimizeMaxVariableModificationParameter " + value + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());

                value = optimizeMinMatchedPeaksParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SageMinMatchedPeaksPeakParameter"));
                if (value != sageParameters.getMinMatchedPeaks()) {
                    sageParameters.setMinMatchedPeaks(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);

                }
                System.out.println("optimizeMinMatchedPeaksParameter " + value + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());

                value = optimizeMaxFragmentChargeParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SageMaxFragmentChargeParameter"));
                if (value != 0) {
                    sageParameters.setMaxFragmentCharge(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                System.out.println("optimizeMaxFragmentChargeParameter " + value + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());

//                String valueStr = optimizeFragmentationMethod(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("MyriMatchFragmentationMethiodParameter"));
//                sageParameters.setFragmentationRule(valueStr);
//                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                value = optimizeMaxVarPTMsNumber(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("MyriMatchMaxVarPTMsParameter"));
//                sageParameters.setMaxVariableMods(value);
//                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//                String valueStr = optimizeFragmentationMethod(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("MyriMatchFragmentationMethiodParameter"));
//                sageParameters.setFragmentationRule(valueStr);
//                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//
//                value = optimizeEnzymaticTerminals(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("MyriMatchEnzymatricTerminalsParameter"));
//                sageParameters.setMinTerminiCleavages(value);
//                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//
                boolean valueBoolean = optimizeGenerateDecoyParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SageGenerateDecoyParameter"));
                if (valueBoolean != sageParameters.getGenerateDecoys()) {
                    sageParameters.setGenerateDecoys(valueBoolean);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                System.out.println("optimizeGenerateDecoyParameter " + valueBoolean + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());

                valueBoolean = optimizeDeisotopParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SageDeisotopParameter"));
                if (valueBoolean != sageParameters.getDeisotope()) {
                    sageParameters.setDeisotope(valueBoolean);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                System.out.println("optimizeDeisotopParameter " + valueBoolean + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());

                valueBoolean = optimizeChimericSpectraParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SageChimericSpectraParameter"));
                if (valueBoolean != sageParameters.getChimera()) {
                    sageParameters.setChimera(valueBoolean);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                System.out.println("SageChimericSpectraParameter " + valueBoolean + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());

                valueBoolean = optimizeWideWindowParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SageWideWindowParameter"));
                if (valueBoolean != sageParameters.getWideWindow()) {
                    sageParameters.setWideWindow(valueBoolean);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                System.out.println("SageWideWindowParameter " + valueBoolean + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());

                valueBoolean = optimizePredectRetentionTimeParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SagePredectRetentionTimeParameter"));
                if (valueBoolean != sageParameters.getPredictRt()) {
                    sageParameters.setPredictRt(valueBoolean);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
                System.out.println("SagePredectRetentionTimeParameter " + valueBoolean + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());

//
//                valueBoalen = optimizeComputeXCorr(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("MyriMatchComputeXCorrParameter"));
//                sageParameters.setComputeXCorr(valueBoalen);
//                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//
//                double dvalue = optimizeoptimizeTICCutoff(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("MyriMatchTICCutoffParameter"));
//                sageParameters.setTicCutoffPercentage(dvalue);
//                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//
//                value = optimizeNumberOfIntensityClasses(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("MyriMatchNumberOfIntensityClassesParameter"));
//                sageParameters.setNumIntensityClasses(value);
//                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//
//                value = optimizeClassSizeMultiplier(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("MyriMatchClassSizeMultiplierParameter"));
//                sageParameters.setClassSizeMultiplier(value);
//                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//
//                value = optimizeNumberOfBatches(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("MyriMatchNumberOfBatchesParameter"));
//                sageParameters.setNumberOfBatches(value);
//                IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
//
                values = optimizeNumberOfPeakParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("SageNumberOfPeakParameter"));
                if (values[1] != sageParameters.getMaxPeaks() || values[0] != sageParameters.getMinPeaks()) {
                    sageParameters.setMinPeaks(values[0]);
                    sageParameters.setMaxPeaks(values[1]);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                    System.out.println("peaks length " + values[0] + "--" + values[1] + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
                }
                continue;
            }

            if (param.equalsIgnoreCase("PrecursorToleranceParameter") && searchInputSetting.isOptimizePrecursorToleranceParameter()) {

                double value = this.optimizePrecursorToleranceParameter(optProtDataset, identificationParametersFile, searchInputSetting, parameterScoreMap.get("PrecursorToleranceParameter"));
                if (value != identificationParameters.getSearchParameters().getPrecursorAccuracy()) {
                    identificationParameters.getSearchParameters().setPrecursorAccuracy(value);
                    if (value > 1) {
                        identificationParameters.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
                    } else {
                        identificationParameters.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.DA);
                    }
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
                }
            }
        }
        if (!digestionParameterOpt.equalsIgnoreCase(identificationParameters.getSearchParameters().getDigestionParameters().getCleavageParameter().name())) {
            optimisedSearchResults.setDigestionParameter(digestionParameterOpt);
            identificationParameters.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf(digestionParameterOpt));
            IdentificationParameters.saveIdentificationParameters(identificationParameters, identificationParametersFile);
        }

    }

    @Override
    public synchronized RawScoreModel excuteSearch(SearchingSubDataset optProtDataset, String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addSpectraList, SearchInputSetting optProtSearchSettings, File identificationParametersFile, boolean pairData) {
        if (!optProtSearchSettings.getSageEnabledParameters().getParamsToOptimize().isEnabledParam(paramOption.split("_")[0])) {
            System.out.println("param " + paramOption + " is not supported " + paramOption);
            return new RawScoreModel();
        }
        if (defaultOutputFileName.contains("_resultsf_Carbamilation of protein N-term") || defaultOutputFileName.contains("_resultsf_Acetylation of protein N-term") || defaultOutputFileName.contains("_resultsf_Pyrolidone from carbamidomethylated C")) {
            System.out.println("param " + paramOption + " is not supported " + paramOption);
            return new RawScoreModel();
        }
        //            SearchParameters searchParameters = tempIdParam.getSearchParameters();
//        if (addSpectraList) {

//        }
        File resultOutput = SearchExecuter.executeSearch(defaultOutputFileName, optProtSearchSettings, optProtDataset.getSubMsFile(), optProtDataset.getSubFastaFile(), tempIdParam, identificationParametersFile);
        List<SpectrumMatch> validatedMaches = SpectraUtilities.getValidatedIdentificationResults(resultOutput, optProtDataset.getSubMsFile(), Advocate.sage, tempIdParam);
        System.out.println("at param name is " + paramOption + "  " + validatedMaches.size());

//            List<SpectrumMatch> updatedList = new ArrayList<>();
//        if (addSpectraList) {
//            System.out.println("----------------------------------------------------------olyyyy add use peptide mass----------------------------------------------------------");
//                File sageFile = new File(resultOutput, SearchHandler.getSageFileName(IoUtil.removeExtension(optProtDataset.getSubMsFile().getName())));
//                IdfileReader idReader = readerFactory.getFileReader(sageFile);
//                SequenceMatchingParameters modificationSequenceMatchingParameters = tempIdParam.getModificationLocalizationParameters().getSequenceMatchingParameters();
//                FMIndex sequenceProvider = new FMIndex(optProtDataset.getSubFastaFile(), null, new OptProtWaitingHandler(), false, tempIdParam);
//                for (SpectrumMatch sm : validatedMaches) {
//                    Thread t = new Thread(() -> {
//                        for (PeptideAssumption pepAss : sm.getAllPeptideAssumptions().toList()) {
//                            try {
//                                Peptide pep = pepAss.getPeptide();
//                                ModificationLocalizationMapper.modificationLocalization(
//                                        pep,
//                                        tempIdParam,
//                                        idReader,
//                                        ptmFactory,
//                                        sequenceProvider
//                                );
//                                pepAss.getPeptide().setMass(pep.getMass(searchParameters.getModificationParameters(), sequenceProvider, modificationSequenceMatchingParameters));
//                                sm.setBestPeptideAssumption(pepAss);
//                            } catch (Exception e) {
//                            }
//                        }
//                    });
//                    t.start();
//                    while (t.isAlive()) {
//                        Thread.currentThread().sleep(10);
//                    }
//
//                }
//        }
        RawScoreModel rawScore = SpectraUtilities.getComparableRawScore(optProtDataset, validatedMaches, Advocate.sage, pairData,addSpectraList);//(optProtDataset, resultOutput, optProtDataset.getSubMsFile(), Advocate.sage, tempIdParam, updateDataReference);
        //
//            if (rawScore.isSignificatChange()) {
//                System.out.println("at test ACCEPT CHANGE  " + paramOption + "  " + optProtDataset.getValidatedIdRefrenceData().length + "  " + validatedMaches.size());
//                System.out.println("maches--> " + validatedMaches.size() + "  the active size " + optProtDataset.getActiveIdentificationNum() + "   the reference :" + paramOption + "   " + optProtDataset.getValidatedIdRefrenceData().length);
//                if (validatedMaches.size() < optProtDataset.getActiveIdentificationNum()) {
//                    int size = validatedMaches.size();
//                    for (int i = size; i < size * 1.02; i++) {
//                        System.out.println("add compåonent");
//                        updatedList.add(null);
//                    }
//                }
//                updatedList.addAll(validatedMaches);
//                validatedMaches.clear();
//            }
//            if (!rawScore.isSignificatChange()) {
//                System.out.println("length before removing the maches--> " + validatedMaches.size() + "  the active size " + optProtDataset.getActiveIdentificationNum() + "   the reference :" + paramOption + "   " + optProtDataset.getValidatedIdRefrenceData().length);
//                validatedMaches.clear();
//            }
        MainUtilities.deleteFolder(resultOutput);
        if (addSpectraList && rawScore.isSignificatChange()) {
            rawScore.setSpectrumMatchResult(validatedMaches);
        }

//            validatedMaches.addAll(updatedList);
//            System.out.println("validated size after " + validatedMaches.size());
//
        return (rawScore);
    }

    public int[] optimizePeptideLengthParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameter = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());

        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxPeptideLengthOption = sageParameter.getMaxPeptideLength();
        int selectedMinPeptideLengthOption = sageParameter.getMinPeptideLength();
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters tempSageParameters = (SageParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        tempSageParameters.setMaxPeptideLength(selectedMaxPeptideLengthOption);

        for (int i = Math.min(10, selectedMinPeptideLengthOption + 2); i >= Math.max(selectedMinPeptideLengthOption - 2, 5); i--) {
            if (i == selectedMinPeptideLengthOption) {
                continue;
            }
            tempSageParameters.setMinPeptideLength(i);
            final String option = "minPeptideLength_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange() || scoreModel.isSameData()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        int indexI = 0;
        Set<Integer> filterSet = new LinkedHashSet<>();
        Set<Integer> finalSet = new LinkedHashSet<>();
        if (!resultsMap.isEmpty()) {

            for (int option : resultsMap.keySet()) {
                int indexII = 0;
                filterSet.clear();
                filterSet.add(option);
                for (int option2 : resultsMap.keySet()) {
                    if (indexII <= indexI) {
                        indexII++;
                        continue;
                    }
                    if (SpectraUtilities.isBetterScore(resultsMap.get(option).getSpectrumMatchResult(), resultsMap.get(option2).getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber())==0) {

                        filterSet.add(option2);
                    }

                    indexII++;
                }
                if (filterSet.size() > 1 && filterSet.contains(option)) {
                    filterSet.remove(option);
                    finalSet.remove(option);
                }
                finalSet.addAll(filterSet);
                indexI++;
            }
            selectedMinPeptideLengthOption = finalSet.iterator().next();
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedMinPeptideLengthOption));
        }

        resultsMap.clear();
        tempSageParameters.setMinPeptideLength(selectedMinPeptideLengthOption);
        for (int i = Math.max(25, selectedMaxPeptideLengthOption - 5); i < Math.min(35, selectedMaxPeptideLengthOption + 5); i++) {
            if (i == selectedMaxPeptideLengthOption) {
                continue;
            }
            tempSageParameters.setMaxPeptideLength(i);
            final String option = "maxPeptideLength_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("maxPeptideLength");

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        indexI = 0;
        filterSet.clear();
        finalSet.clear();
        if (!resultsMap.isEmpty()) {

            for (int option : resultsMap.keySet()) {
                int indexII = 0;
                filterSet.clear();
                filterSet.add(option);
                for (int option2 : resultsMap.keySet()) {
                    if (indexII <= indexI) {
                        indexII++;
                        continue;
                    }
                    if (SpectraUtilities.isBetterScore(resultsMap.get(option).getSpectrumMatchResult(), resultsMap.get(option2).getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber())==0) {
                        filterSet.add(option2);
                    }
                    indexII++;
                }
                if (filterSet.size() > 1 && filterSet.contains(option)) {
                    filterSet.remove(option);
                    finalSet.remove(option);
                }
                finalSet.addAll(filterSet);
                indexI++;
            }
            selectedMaxPeptideLengthOption = finalSet.iterator().next();
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedMaxPeptideLengthOption));
        }

//        if (!resultsMap.isEmpty()) {
//            for (int option : resultsMap.keySet()) {
//
//                sortedResultsMap.put(resultsMap.get(option), option);
//            }
//            selectedMaxPeptideLengthOption = sortedResultsMap.firstEntry().getValue();
//            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey().getSpectrumMatchResult());
//        }
        return new int[]{selectedMinPeptideLengthOption, selectedMaxPeptideLengthOption};
    }

    public double[] optimizeFragmentMzParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameter = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());

        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedMaxFragmentMzOption = sageParameter.getMaxFragmentMz();
        double selectedMinFragmentMzOption = sageParameter.getMinFragmentMz();
        Map<Double, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters tempSageParameters = (SageParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        tempSageParameters.setMaxFragmentMz(selectedMaxFragmentMzOption);
        for (double i = 300; i >= 100;) {
            if (i == selectedMinFragmentMzOption) {
                i -= 20.0;
                continue;
            }
            tempSageParameters.setMinFragmentMz(i);
            final String option = "minFragmentMz_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            double j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("minFragmentMz");

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            i -= 20.0;
        }

        int indexI = 0;
        Set<Double> filterSet = new LinkedHashSet<>();
        Set<Double> finalSet = new LinkedHashSet<>();
        if (!resultsMap.isEmpty()) {

            for (double option : resultsMap.keySet()) {
                int indexII = 0;
                filterSet.clear();
                filterSet.add(option);
                for (double option2 : resultsMap.keySet()) {
                    if (indexII <= indexI) {
                        indexII++;
                        continue;
                    }
                    if (SpectraUtilities.isBetterScore(resultsMap.get(option).getSpectrumMatchResult(), resultsMap.get(option2).getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber())==0) {
                        filterSet.add(option2);
                    }
                    indexII++;
                }
                if (filterSet.size() > 1 && filterSet.contains(option)) {
                    filterSet.remove(option);
                    finalSet.remove(option);
                }
                finalSet.addAll(filterSet);
                indexI++;
            }
            selectedMinFragmentMzOption = finalSet.iterator().next();
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedMinFragmentMzOption));
        }

//        TreeMap<RawScoreModel, Double> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
//        if (!resultsMap.isEmpty()) {
//            for (Double option : resultsMap.keySet()) {
//                sortedResultsMap.put(resultsMap.get(option), option);
//            }
//            selectedMinFragmentMzOption = sortedResultsMap.firstEntry().getValue();
//            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey().getSpectrumMatchResult());
//        }
        tempSageParameters.setMinFragmentMz(selectedMinFragmentMzOption);
        resultsMap.clear();
        for (double i = 1500; i <= 2500;) {
            if (i == selectedMaxFragmentMzOption) {
                i += 250;
                continue;
            }
            tempSageParameters.setMaxFragmentMz(i);
            final String option = "maxFragmentMz_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            double j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("maxFragmentMz");

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

            i += 250;
        }
        indexI = 0;
        filterSet.clear();
        finalSet.clear();
        if (!resultsMap.isEmpty()) {
            for (double option : resultsMap.keySet()) {
                int indexII = 0;
                filterSet.clear();
                filterSet.add(option);
                for (double option2 : resultsMap.keySet()) {
                    if (indexII <= indexI) {
                        indexII++;
                        continue;
                    }
                    if (SpectraUtilities.isBetterScore(resultsMap.get(option).getSpectrumMatchResult(), resultsMap.get(option2).getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber())==0) {
                        filterSet.add(option2);
                    }
                    indexII++;
                }
                if (filterSet.size() > 1 && filterSet.contains(option)) {
                    filterSet.remove(option);
                    finalSet.remove(option);
                }
                finalSet.addAll(filterSet);
                indexI++;
            }
            selectedMaxFragmentMzOption = finalSet.iterator().next();
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedMaxFragmentMzOption));
        }
        System.out.println("final FragmentMzOption " + selectedMinFragmentMzOption + "," + selectedMaxFragmentMzOption);
        return new double[]{selectedMinFragmentMzOption, selectedMaxFragmentMzOption};
    }

    public double[] optimizePeptideMassParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameter = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());

        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedMaxPeptideMassOption = sageParameter.getMaxPeptideMass();
        double selectedMinPeptideMassOption = sageParameter.getMinPeptideMass();
        Map<Double, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters tempSageParameters = (SageParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        tempSageParameters.setMaxPeptideMass(selectedMaxPeptideMassOption);
        for (double i = Math.min(1000, selectedMinPeptideMassOption + 200); i >= Math.max(200, selectedMinPeptideMassOption - 200);) {
            if (i == selectedMinPeptideMassOption) {
                i -= 100.0;
                continue;
            }
            tempSageParameters.setMinPeptideMass(i);
//            ((SageParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex())).setMinPeptideMass(i);
            final String option = "minPeptideMass_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            double j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("minPeptideMass");

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            i -= 100.0;
        }

        int indexI = 0;
        Set<Double> filterSet = new LinkedHashSet<>();
        Set<Double> finalSet = new LinkedHashSet<>();
        if (!resultsMap.isEmpty()) {

            for (double option : resultsMap.keySet()) {
                int indexII = 0;
                filterSet.clear();
                filterSet.add(option);
                for (double option2 : resultsMap.keySet()) {
                    if (indexII <= indexI) {
                        indexII++;
                        continue;
                    }
                    if (SpectraUtilities.isBetterScore(resultsMap.get(option).getSpectrumMatchResult(), resultsMap.get(option2).getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber())==0) {
                        filterSet.add(option2);
                    }
                    indexII++;
                }
                if (filterSet.size() > 1 && filterSet.contains(option)) {
                    filterSet.remove(option);
                    finalSet.remove(option);
                }
                finalSet.addAll(filterSet);
                indexI++;
            }
            selectedMinPeptideMassOption = finalSet.iterator().next();
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedMinPeptideMassOption));
        }
        tempSageParameters.setMinPeptideMass(selectedMinPeptideMassOption);
        resultsMap.clear();
        for (double i = Math.max(2000, selectedMaxPeptideMassOption - 2000); i <= Math.min(8000, selectedMaxPeptideMassOption + 2000);) {
            if (i == selectedMaxPeptideMassOption) {
                i += 500;
                continue;
            }
            tempSageParameters.setMaxPeptideMass(i);
            final String option = "maxPeptideMass_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            double j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("maxPeptideMass");

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

            i += 500;
        }
        indexI = 0;
        filterSet.clear();
        finalSet.clear();
        if (!resultsMap.isEmpty()) {
            for (double option : resultsMap.keySet()) {
                int indexII = 0;
                filterSet.clear();
                filterSet.add(option);
                for (double option2 : resultsMap.keySet()) {
                    if (indexII <= indexI) {
                        indexII++;
                        continue;
                    }
                    if (SpectraUtilities.isBetterScore(resultsMap.get(option).getSpectrumMatchResult(), resultsMap.get(option2).getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber())==0) {
                        filterSet.add(option2);
                    }
                    indexII++;
                }
                if (filterSet.size() > 1 && filterSet.contains(option)) {
                    filterSet.remove(option);
                    finalSet.remove(option);
                }
                finalSet.addAll(filterSet);
                indexI++;
            }
            selectedMaxPeptideMassOption = finalSet.iterator().next();
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedMaxPeptideMassOption));
        }

        return new double[]{selectedMinPeptideMassOption, selectedMaxPeptideMassOption};
    }

    public int[] optimizeNumberOfPeakParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameter = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());

        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxPeaksNumberOption = sageParameter.getMaxPeaks();
        int selectedMinPeaksNumberOption = sageParameter.getMinPeaks();
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters tempSageParameters = (SageParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        tempSageParameters.setMaxPeaks(selectedMaxPeaksNumberOption);

        for (int i = Math.min(20, selectedMinPeaksNumberOption + 5); i >= Math.max(selectedMinPeaksNumberOption - 5, 5); i--) {
            if (i == selectedMinPeaksNumberOption) {
                continue;
            }
            tempSageParameters.setMinPeaks(i);
            final String option = "minPeaks_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange() || scoreModel.isSameData()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        int indexI = 0;
        Set<Integer> filterSet = new LinkedHashSet<>();
        Set<Integer> finalSet = new LinkedHashSet<>();
        if (!resultsMap.isEmpty()) {

            for (int option : resultsMap.keySet()) {
                int indexII = 0;
                filterSet.clear();
                filterSet.add(option);
                for (int option2 : resultsMap.keySet()) {
                    if (indexII <= indexI) {
                        indexII++;
                        continue;
                    }

                    if (SpectraUtilities.isBetterScore(resultsMap.get(option).getSpectrumMatchResult(), resultsMap.get(option2).getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber())==0) {

                        filterSet.add(option2);
                    }
                    System.out.println("option 1 " + option + " in compare with " + option2 + "  " + filterSet);
                    indexII++;
                }
                if (filterSet.size() > 1 && filterSet.contains(option)) {
                    filterSet.remove(option);
                    finalSet.remove(option);
                }
                finalSet.addAll(filterSet);
                indexI++;
            }
            selectedMinPeaksNumberOption = finalSet.iterator().next();
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedMinPeaksNumberOption));
        }

        resultsMap.clear();
        tempSageParameters.setMinPeaks(selectedMinPeaksNumberOption);
        for (int i = Math.max(100, selectedMaxPeaksNumberOption - 50); i < Math.min(200, selectedMaxPeaksNumberOption + 50);) {
            if (i == selectedMaxPeaksNumberOption) {
                i += 10;
                continue;
            }
            tempSageParameters.setMaxPeaks(i);
            final String option = "maxPeaks_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("maxPeaks");

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            i += 10;
        }

        indexI = 0;
        filterSet.clear();
        finalSet.clear();
        if (!resultsMap.isEmpty()) {

            for (int option : resultsMap.keySet()) {
                int indexII = 0;
                filterSet.clear();
                filterSet.add(option);
                for (int option2 : resultsMap.keySet()) {
                    if (indexII <= indexI) {
                        indexII++;
                        continue;
                    }
                    if (SpectraUtilities.isBetterScore(resultsMap.get(option).getSpectrumMatchResult(), resultsMap.get(option2).getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber())
                            ==0) {
                        filterSet.add(option2);
                    }
                    indexII++;
                }
                if (filterSet.size() > 1 && filterSet.contains(option)) {
                    filterSet.remove(option);
                    finalSet.remove(option);
                }
                finalSet.addAll(filterSet);
                indexI++;
            }
            selectedMaxPeaksNumberOption = finalSet.iterator().next();
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedMaxPeaksNumberOption));
        }
//        System.out.println("final peptide mass " + selectedMinPeaksNumberOption + "," + selectedMaxPeaksNumberOption);

        return new int[]{selectedMinPeaksNumberOption, selectedMaxPeaksNumberOption};
    }

    public int optimizeIonMinIndexParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {

        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        int selectedOption = sageParameters.getMinIonIndex();
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        sageParameters = (SageParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        for (int i = 0; i <= 6; i++) {
            final String option = "ionMinIndex_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setMinIonIndex(i);
            final int j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("ionMinIndex");

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange() || scoreModel.isSameData()) {
                    System.out.println("ion was seg " + option);
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            if (selectedOption != (sortedResultsMap.firstEntry().getValue())) {
                selectedOption = sortedResultsMap.firstEntry().getValue();
                optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
            }
        }
        return selectedOption;
    }

    public int optimizeMaxVariableModificationParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {

        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//      sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        int selectedOption = oreginaltempIdParam.getSearchParameters().getModificationParameters().getVariableModifications().size();
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameters = (SageParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        for (int i = 0; i <= selectedOption; i++) {
            final String option = "maxVarPTMs_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setMaxVariableMods(i);
            final int j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("maxVarPTMs");

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange() || scoreModel.isSameData()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            if (selectedOption != (sortedResultsMap.firstEntry().getValue())) {
                selectedOption = sortedResultsMap.firstEntry().getValue();
                optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
            }
        }
        return selectedOption;
    }

    public int optimizeMinMatchedPeaksParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {

        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        int selectedOption = sageParameters.getMinMatchedPeaks();
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        sageParameters = (SageParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        for (int i = Math.max(2, selectedOption - 2); i <= Math.min(6, selectedOption + 2); i++) {
            final String option = "minMatchedPeaks_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setMinMatchedPeaks(i);
            final int j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("minMatchedPeaks");

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange() || scoreModel.isSameData()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            if (selectedOption != (sortedResultsMap.firstEntry().getValue())) {
                selectedOption = sortedResultsMap.firstEntry().getValue();
                optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
            }
        }
        return selectedOption;
    }

    public int optimizeMaxFragmentChargeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {

        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        int selectedOption = 0;
        if (sageParameters.getMaxFragmentCharge() != null) {
            selectedOption = sageParameters.getMaxFragmentCharge();
        }
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        sageParameters = (SageParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        for (int i = 1; i <= 5; i++) {
            final String option = "maxFragmentCharge_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setMaxFragmentCharge(i);
            final int j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("maxFragmentCharge");

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange() || scoreModel.isSameData()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            if (selectedOption != (sortedResultsMap.firstEntry().getValue())) {
                selectedOption = sortedResultsMap.firstEntry().getValue();
                optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
            }
        }
        return selectedOption;
    }
//
//    public int optimizeEnzymaticTerminals(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        MyriMatchParameters myriMatchParameters = (MyriMatchParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.myriMatch.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        int selectedOption = myriMatchParameters.getMinTerminiCleavages();
//
//        for (int i = 0; i <= 2; i++) {
//            final String option = "enzymatricTerminals_" + i;
//
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            myriMatchParameters.setMinTerminiCleavages(i);
//
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("enzymatricTerminals");
//            int j = i;
//            Future future = MainUtilities.getExecutorService().submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = -1;
//        for (int enzymatricTerminals : resultsMap.keySet()) {
//            System.out.println(" option enzymatricTerminals" + enzymatricTerminals + "  " + resultsMap.get(enzymatricTerminals) + " > " + localId + "  " + idRate);
//            if (resultsMap.get(enzymatricTerminals) > localId) {
//                localId = resultsMap.get(enzymatricTerminals);
//                localSelection = enzymatricTerminals;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//        return selectedOption;
//    }
//
//    public double optimizeoptimizeTICCutoff(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Double, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        MyriMatchParameters myriMatchParameters = (MyriMatchParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.myriMatch.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        double selectedOption = myriMatchParameters.getTicCutoffPercentage();
//        for (double i = 0.90; i <= 1;) {
//            final String option = "TICCutoff_" + i;
//
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            myriMatchParameters.setTicCutoffPercentage(i);
//            final double j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("TICCutoff");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//            i += 0.01;
//        }
//        int localId = -1;
//        double localSelection = 0;
//        for (double dRangeScore : resultsMap.keySet()) {
//            System.out.println(" option 1 " + dRangeScore + "  " + resultsMap.get(dRangeScore) + " > " + localId + "  " + idRate);
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//    }
//
//    public int optimizeNumberOfIntensityClasses(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        MyriMatchParameters myriMatchParameters = (MyriMatchParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.myriMatch.getIndex());
//        Integer selectedOption = myriMatchParameters.getNumIntensityClasses();
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//
//        for (int i = 1; i <= 6; i++) {
//            final String option = "NumberOfIntensityClasses_" + i;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            myriMatchParameters.setNumIntensityClasses(i);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("NumberOfIntensityClasses");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//        int localId = -1;
//        int localSelection = 0;
//        for (int option : resultsMap.keySet()) {
//            System.out.println(" option #NumberOfIntensityClasses " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
//            if (resultsMap.get(option) > localId) {
//                localId = resultsMap.get(option);
//                localSelection = option;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//        return selectedOption;
//    }
//
//    public int optimizeClassSizeMultiplier(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        MyriMatchParameters myriMatchParameters = (MyriMatchParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.myriMatch.getIndex());
//        Integer selectedOption = myriMatchParameters.getClassSizeMultiplier();
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//
//        for (int i = 1; i <= 5; i++) {
//            final String option = "classSizeMultiplier_" + i;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            myriMatchParameters.setClassSizeMultiplier(i);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("classSizeMultiplier");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//        int localId = -1;
//        int localSelection = 0;
//        for (int option : resultsMap.keySet()) {
//            System.out.println(" option #classSizeMultiplier " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
//            if (resultsMap.get(option) > localId) {
//                localId = resultsMap.get(option);
//                localSelection = option;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//        return selectedOption;
//    }
//
//    public int optimizeNumberOfBatches(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        MyriMatchParameters myriMatchParameters = (MyriMatchParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.myriMatch.getIndex());
//        Integer selectedOption = myriMatchParameters.getNumberOfBatches();
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//
//        for (int i = 30; i <= 80;) {
//            final String option = "NumberOfBatches_" + i;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            myriMatchParameters.setNumberOfBatches(i);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("NumberOfBatches");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//            i += 10;
//        }
//        int localId = -1;
//        int localSelection = 0;
//        for (int option : resultsMap.keySet()) {
//            System.out.println(" option #NumberOfBatches " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
//            if (resultsMap.get(option) > localId) {
//                localId = resultsMap.get(option);
//                localSelection = option;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//        return selectedOption;
//    }
//
//    public int optimizeMaxPeakCount(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        MyriMatchParameters myriMatchParameters = (MyriMatchParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.myriMatch.getIndex());
//        Integer selectedOption = myriMatchParameters.getMaxPeakCount();
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//
//        for (int i = 200; i <= 400;) {
//            final String option = "maxPeakCount_" + i;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            myriMatchParameters.setMaxPeakCount(i);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("maxPeakCount");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//            i += 50;
//        }
//        int localId = -1;
//        int localSelection = 0;
//        for (int option : resultsMap.keySet()) {
//            System.out.println(" option #maxPeakCount " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
//            if (resultsMap.get(option) > localId) {
//                localId = resultsMap.get(option);
//                localSelection = option;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//        return selectedOption;
//    }
//
//    public int optimizeMinimumPeaks(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        Integer selectedOption = xtandemParameters.getMinPeaksPerSpectrum();
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        int lastValue = 0;
//        for (int i = 5; i <= 100;) {
//            final String option = "minpeaksNum_" + i;
//            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//            xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setMinPeaksPerSpectrum(i);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("minpeaksNum");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//            if (lastValue > resultsMap.get(i)) {
//                break;
//            }
//            lastValue = resultsMap.get(i);
//            i += 10;
//
//        }
//        int localId = -1;
//        int localSelection = 0;
//        for (int option : resultsMap.keySet()) {
//            if (resultsMap.get(option) > localId) {
//                localId = resultsMap.get(option);
//                localSelection = option;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//        return selectedOption;
//    }
//
//    public double optimizeNoiseSuppression(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Double, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        double selectedOption = xtandemParameters.getMinPrecursorMass();
//        final String option = "noiseSupression_" + false;
//        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//        xtandemParameters.setUseNoiseSuppression(false);
//        final ParameterScoreModel paramScore = new ParameterScoreModel();
//        paramScore.setParamId("noiseSupression");
//
//        Future future = executor.submit(() -> {
//            resultsMap.put(0.0, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//            paramScore.setScore(resultsMap.get(0.0));
//            paramScore.setParamValue(option);
//            parameterScoreSet.add(paramScore);
//        });
//        while (!future.isDone()) {
//        }
//        xtandemParameters.setUseNoiseSuppression(true);
//        for (double j = 500; j < 1600;) {
//            final String suboption = "noiseSupression_" + true + "_" + j;
//            final String subupdatedName = Configurations.DEFAULT_RESULT_NAME + "_" + suboption + "_" + msFileName;
//            final double i = j;
//            xtandemParameters.setMinPrecursorMass(j);
//            final ParameterScoreModel paramScore2 = new ParameterScoreModel();
//            paramScore2.setParamId("noiseSupression");
//            future = executor.submit(() -> {
//                resultsMap.put(i, excuteSearch(optProtDataset, subupdatedName, suboption, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore2.setScore(resultsMap.get(i));
//                paramScore2.setParamValue(option);
//                parameterScoreSet.add(paramScore2);
//            });
//            while (!future.isDone()) {
//            }
//            j += 350;
//
//        }
//
//        int localId = -1;
//        double localSelection = 0;
//        for (double dRangeScore : resultsMap.keySet()) {
////            System.out.println(" option 1 " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//    }
////
//

    public boolean optimizeGenerateDecoyParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = sageParameters.getGenerateDecoys();

        for (int i = 0; i < 2; i++) {
            final String option = "generateDecoy_" + (i == 1);
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setGenerateDecoys(i == 1);
            final int j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("generateDecoy");
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
        }
        return selectedOption;

    }

    public boolean optimizeDeisotopParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = sageParameters.getDeisotope();

        for (int i = 0; i < 2; i++) {
            final String option = "Deisotope_" + (i == 1);
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setDeisotope(i == 1);
            final int j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("Deisotope");
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
        }
        return selectedOption;

    }

    public boolean optimizeChimericSpectraParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = sageParameters.getChimera();

        for (int i = 0; i < 2; i++) {
            final String option = "Chimera_" + (i == 1);
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setChimera(i == 1);
            final int j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("Chimera");
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
        }
        return selectedOption;

    }

    public boolean optimizeWideWindowParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = sageParameters.getWideWindow();

        for (int i = 0; i < 2; i++) {
            final String option = "WideWindow_" + (i == 1);
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setWideWindow(i == 1);
            final int j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("WideWindow");
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
        }
        return selectedOption;

    }

    public boolean optimizePredectRetentionTimeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = sageParameters.getPredictRt();

        for (int i = 0; i < 2; i++) {
            final String option = "PredictRt_" + (i == 1);
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setPredictRt(i == 1);
            final int j = i;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("PredictRt");
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setScore(resultsMap.get(j).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
        }
        return selectedOption;

    }

////
//
//    public boolean optimizeQuickAcetyl(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isProteinQuickAcetyl();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useQuickAcetyl = (i == 1);
//            final String option = "useQuickAcetyl_" + useQuickAcetyl;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setProteinQuickAcetyl(useQuickAcetyl);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("QuickAcetyl");
//            Future future = MainUtilities.getExecutorService().submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
////            System.out.println("at quick acetyle " + useQuickAcetyl + " " + resultsMap.get(j));
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//        return selectedOption;
//    }
//
//    public boolean optimizeQuickPyrolidone(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isQuickPyrolidone();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useQuickPyrolidone = (i == 1);
//            final String option = "useQuickPyrolidone_" + useQuickPyrolidone;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setQuickPyrolidone(useQuickPyrolidone);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("QuickPyrolidone");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//            System.out.println("at quick prolyien  " + useQuickPyrolidone + " " + resultsMap.get(j));
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//    }
////
//
//    public boolean optimizeStPBias(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isStpBias();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useStpBias = (i == 1);
//            final String option = "useStpBias_" + useStpBias;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setStpBias(useStpBias);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("StpBias");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//
////            System.out.println("at useStpBias  " + useStpBias + " " + resultsMap.get(j));
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//    }
////
//
//    public boolean optimizeUseSmartPlus3Model(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        MyriMatchParameters myriMatchParameters = (MyriMatchParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.myriMatch.getIndex());
//
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = myriMatchParameters.getUseSmartPlusThreeModel();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useSmartPlus3Model = (i == 1);
//            final String option = "useSmartPlus3Model_" + useSmartPlus3Model;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            myriMatchParameters.setUseSmartPlusThreeModel(useSmartPlus3Model);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("useSmartPlus3Model");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
////
////    public boolean optimizeComputeXCorr(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        MyriMatchParameters myriMatchParameters = (MyriMatchParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.myriMatch.getIndex());
//
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = myriMatchParameters.getComputeXCorr();
//
//        for (int i = 0; i < 2; i++) {
//            boolean computeXCorr = (i == 1);
//            final String option = "computeXCorr_" + computeXCorr;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            myriMatchParameters.setUseSmartPlusThreeModel(computeXCorr);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("computeXCorr");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
//
//
//    public boolean optimizeRefineUnanticipatedCleavage(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        xtandemParameters.setRefine(true);
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isRefineUnanticipatedCleavages();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useRefineUnanticipatedCleavages = (i == 1);
//            final String option = "useRefineUnanticipatedCleavages_" + useRefineUnanticipatedCleavages;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setRefineUnanticipatedCleavages(useRefineUnanticipatedCleavages);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("UnanticipatedCleavages");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
////
//
//    public boolean optimizeRefineSimiEnzymaticCleavage(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        xtandemParameters.setRefine(true);
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isRefineSemi();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useRefineSimiEnzymaticCleavage = (i == 1);
//            final String option = "useRefineSimiEnzymaticCleavage_" + useRefineSimiEnzymaticCleavage;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setRefineSemi(useRefineSimiEnzymaticCleavage);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("SimiEnzymaticCleavage");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
////
//
//    public boolean optimizePotintialModification(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        xtandemParameters.setRefine(true);
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isPotentialModificationsForFullRefinment();
//
//        for (int i = 0; i < 2; i++) {
//            boolean usePotintialModification = (i == 1);
//            final String option = "usePotintialModification_" + usePotintialModification;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setPotentialModificationsForFullRefinment(usePotintialModification);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("PotintialModification");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
////
//
//    public boolean optimizeRefinePointMutations(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        xtandemParameters.setRefine(true);
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isRefinePointMutations();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useRefinePointMutations = (i == 1);
//            final String option = "useRefinePointMutations_" + useRefinePointMutations;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setRefinePointMutations(useRefinePointMutations);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("PointMutations");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
////
////
//
//    public boolean optimizeRefineSnAPs(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        xtandemParameters.setRefine(true);
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isRefineSnaps();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useRefineSnAPs = (i == 1);
//            final String option = "useRefineSnAPs_" + useRefineSnAPs;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setRefineSnaps(useRefineSnAPs);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("SnAPs");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
//
//    public boolean optimizeRefineSpectrumSynthesis(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        int idRate = optProtDataset.getActiveIdentificationNum();
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
//        xtandemParameters.setRefine(true);
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        boolean selectedOption = xtandemParameters.isRefineSpectrumSynthesis();
//
//        for (int i = 0; i < 2; i++) {
//            boolean useRefineSpectrumSynthesis = (i == 1);
//            final String option = "useRefineSpectrumSynthesis_" + useRefineSpectrumSynthesis;
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//            xtandemParameters.setRefineSpectrumSynthesis(useRefineSpectrumSynthesis);
//            final int j = i;
//            final ParameterScoreModel paramScore = new ParameterScoreModel();
//            paramScore.setParamId("SpectrumSynthesis");
//
//            Future future = executor.submit(() -> {
//                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setScore(resultsMap.get(j));
//                paramScore.setParamValue(option);
//                parameterScoreSet.add(paramScore);
//            });
//            while (!future.isDone()) {
//            }
//        }
//
//        int localId = -1;
//        int localSelection = 0;
//        for (int dRangeScore : resultsMap.keySet()) {
//            if (resultsMap.get(dRangeScore) > localId) {
//                localId = resultsMap.get(dRangeScore);
//                localSelection = dRangeScore;
//
//            }
//        }
//        if (localId >= idRate) {
//            selectedOption = localSelection == 1;
//            optProtDataset.setActiveIdentificationNum(localId);
//        }
//
//        return selectedOption;
//
//    }
}
