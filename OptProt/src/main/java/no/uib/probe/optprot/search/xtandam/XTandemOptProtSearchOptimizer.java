package no.uib.probe.optprot.search.xtandam;

import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.io.IoUtil;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.parameters.identification.tool_specific.XtandemParameters;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;
import no.uib.probe.optprot.model.OptimisedSearchResults;
import no.uib.probe.optprot.model.ParameterScoreModel;
import no.uib.probe.optprot.model.RawScoreModel;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.search.DefaultOptProtSearchOptimizer;
import no.uib.probe.optprot.search.SearchExecuter;
import no.uib.probe.optprot.util.MainUtilities;
import no.uib.probe.optprot.util.SpectraFileUtilities;

/**
 *
 * @author yfa041
 */
public class XTandemOptProtSearchOptimizer extends DefaultOptProtSearchOptimizer {

    /**
     * The compomics PTM factory.
     */
    private final ModificationFactory ptmFactory = ModificationFactory.getInstance();

    private final SearchingSubDataset optProtDataset;
    private final SearchInputSetting searchInputSetting;
    private final File generatedIdentificationParametersFile;
    private final OptimisedSearchResults optimisedSearchResults;
    private final IdentificationParameters identificationParameters;
    private final Map<String, TreeSet<ParameterScoreModel>> parameterScoreMap;

    public XTandemOptProtSearchOptimizer(SearchingSubDataset optProtDataset, SearchInputSetting searchInputSetting, File generatedIdentificationParametersFile) throws IOException {

        this.optProtDataset = optProtDataset;
        this.searchInputSetting = searchInputSetting;
        this.generatedIdentificationParametersFile = generatedIdentificationParametersFile;
        this.identificationParameters = IdentificationParameters.getIdentificationParameters(generatedIdentificationParametersFile);
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

        parameterScoreMap.put("XtandemSpectrumDynamicRange", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemNumberOfPeaks", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemMinimumFragmentMz", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemMinimumPeaks", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemNoiseSuppression", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemParentIsotopExpansion", new TreeSet<>(Collections.reverseOrder()));

        parameterScoreMap.put("XtandemQuickAcetyl", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemQuickPyrolidone", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemStPBias", new TreeSet<>(Collections.reverseOrder()));

        parameterScoreMap.put("XtandemUseRefine", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemUnanticipatedCleavage", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemRefineSimiEnzymaticCleavage", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemPotintialModification", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemPointMutations", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemSnAPs", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemSpectrumSynthesis", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("XtandemRefVarPTM", new TreeSet<>(Collections.reverseOrder()));
    }

    private String digestionParameterOpt;
    private boolean simiEnzymaticCleavage = false;

    public void startProcess(List<String> paramOrder) throws IOException {
        digestionParameterOpt = identificationParameters.getSearchParameters().getDigestionParameters().getCleavageParameter().name();
        for (String param : paramOrder) {
            System.out.println("-------------------------------------------param " + param + "-------------------------------------------");
            if (param.equalsIgnoreCase("DigestionParameter_1") && searchInputSetting.isOptimizeDigestionParameter()) {

                optimisedSearchResults.setDigestionParameter("enzyme");
                String value = this.optimizeEnzymeParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("EnzymeParameter"));
                optimisedSearchResults.setEnzymeName(value);
                if (!value.equalsIgnoreCase(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName())) {
                    int nMissesCleavages = identificationParameters.getSearchParameters().getDigestionParameters().getnMissedCleavages(value);
                    identificationParameters.getSearchParameters().getDigestionParameters().clearEnzymes();
                    identificationParameters.getSearchParameters().getDigestionParameters().addEnzyme(EnzymeFactory.getInstance().getEnzyme(value));
                    identificationParameters.getSearchParameters().getDigestionParameters().setnMissedCleavages(value, nMissesCleavages);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                }

                value = this.optimizeSpecificityParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("SpecificityParameter"));
                if (!value.equalsIgnoreCase(identificationParameters.getSearchParameters().getDigestionParameters().getSpecificity(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName()).name())) {
                    identificationParameters.getSearchParameters().getDigestionParameters().setSpecificity(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName(), DigestionParameters.Specificity.valueOf(value));
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                }
                if (optProtDataset.getIdentificationRate() < 10) {
                    digestionParameterOpt = this.optimizeDigestionParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("DigestionParameter"));
                }

                continue;

            }

            if (param.equalsIgnoreCase("FragmentIonTypesParameter") && searchInputSetting.isOptimizeFragmentIonTypesParameter()) {
                String value = this.optimizeFragmentIonTypesParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("FragmentIonTypesParameter"));
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
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                }
                continue;

            }

//confusing param
            if (param.equalsIgnoreCase("DigestionParameter_2") && searchInputSetting.isOptimizeDigestionParameter()) {
                int value = this.optimizeMaxMissCleavagesParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("MaxMissCleavagesParameter"));
                if (value != identificationParameters.getSearchParameters().getDigestionParameters().getnMissedCleavages(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName())) {
                    identificationParameters.getSearchParameters().getDigestionParameters().setnMissedCleavages(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName(), value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                }
                continue;

            }

            if (param.equalsIgnoreCase("FragmentToleranceParameter") && searchInputSetting.isOptimizeFragmentToleranceParameter()) {
                double value = this.optimizeFragmentToleranceParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("FragmentToleranceParameter"));
                if (value != identificationParameters.getSearchParameters().getFragmentIonAccuracy()) {
                    identificationParameters.getSearchParameters().setFragmentIonAccuracy(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                }
                continue;
            }

            if (param.equalsIgnoreCase("PrecursorChargeParameter") && searchInputSetting.isOptimizePrecursorChargeParameter()) {

                int[] values = this.optimizePrecursorChargeParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("PrecursorChargeParameter"));
                if (values[1] != identificationParameters.getSearchParameters().getMaxChargeSearched()) {
                    identificationParameters.getSearchParameters().setMinChargeSearched(values[0]);
                    identificationParameters.getSearchParameters().setMaxChargeSearched(values[1]);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                }
                continue;
            }

            if (param.equalsIgnoreCase("XtandemAdvancedParameter") && searchInputSetting.isOptimizeXtandemAdvancedParameter()) {
                advancedParam = true;
                XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
                useRefinment = optimizeUseRefine(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemUseRefine"));
                if (useRefinment != xtandemParameters.isRefine()) {
                    xtandemParameters.setRefine(useRefinment);
                }
                if (!useRefinment) {
                    System.out.println("disable second stage");

                }
                continue;
            }

            if (useRefinment && param.equalsIgnoreCase("XtandemAdvancedParameter_A") && searchInputSetting.isOptimizeXtandemAdvancedParameter()) {
                advancedParam = true;
                XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());

                int ivalue = optimizeSpectrumPeaksNumber(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemNumberOfPeaks"));
                if (ivalue != xtandemParameters.getnPeaks()) {
                    xtandemParameters.setnPeaks(ivalue);
                }
                double dvalue = optimizeMinimumFragmentMz(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemMinimumFragmentMz"));
                if (dvalue != xtandemParameters.getMinFragmentMz()) {
                    xtandemParameters.setMinFragmentMz(dvalue);
                }
                ivalue = optimizeMinimumPeaks(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemMinimumPeaks"));
                if (ivalue != xtandemParameters.getMinPeaksPerSpectrum()) {
                    xtandemParameters.setMinPeaksPerSpectrum(ivalue);
                }
                boolean bvalue = optimizeParentIsotopExpansion(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemParentIsotopExpansion"));
                if (bvalue != xtandemParameters.getParentMonoisotopicMassIsotopeError()) {
                    xtandemParameters.setParentMonoisotopicMassIsotopeError(bvalue);//                  
                }
                bvalue = optimizeQuickAcetyl(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemQuickAcetyl"));
                if (bvalue) {
                    xtandemParameters.setProteinQuickAcetyl(bvalue);
                }
                bvalue = optimizeQuickPyrolidone(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemQuickPyrolidone"));
                if (bvalue) {
                    xtandemParameters.setQuickPyrolidone(bvalue);
                }

                bvalue = optimizeRefineUnanticipatedCleavage(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemUnanticipatedCleavage"));
                if (bvalue != xtandemParameters.isRefineUnanticipatedCleavages()) {
                    xtandemParameters.setRefineUnanticipatedCleavages(bvalue);
                }

                simiEnzymaticCleavage = optimizeRefineSimiEnzymaticCleavage(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemRefineSimiEnzymaticCleavage"));


                bvalue = optimizeRefineSpectrumSynthesis(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemSpectrumSynthesis"));
                if (bvalue != xtandemParameters.isRefineSpectrumSynthesis()) {
                    xtandemParameters.setRefineSpectrumSynthesis(bvalue);
                }

                IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);

                continue;
            }
            if (param.equalsIgnoreCase("PrecursorToleranceParameter") && searchInputSetting.isOptimizePrecursorToleranceParameter()) {

                double value = this.optimizePrecursorToleranceParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("PrecursorToleranceParameter"));
                if (value != identificationParameters.getSearchParameters().getPrecursorAccuracy()) {
                    identificationParameters.getSearchParameters().setPrecursorAccuracy(value);
                    if (value > 1) {
                        identificationParameters.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
                    } else {
                        identificationParameters.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.DA);
                    }
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                }
            }

            if (param.equalsIgnoreCase("ModificationParameter") && searchInputSetting.isOptimizeModificationParameter()) {
                Map<String, Set<String>> modificationsResults = this.optimizeModificationsParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("ModificationsParameter"));
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
                IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                MainUtilities.resetExecutorService();
                continue;
            }
            if (param.equalsIgnoreCase("XtandemAdvancedParameter_B") && searchInputSetting.isOptimizeXtandemAdvancedParameter()) {
                advancedParam = true;
                XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());

                double dvalue = optimizeSpectrumDynamicRange(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemSpectrumDynamicRange"));
                if (dvalue != xtandemParameters.getDynamicRange()) {
                    xtandemParameters.setDynamicRange(dvalue);
                }//                
//                
                boolean bvalue = optimizePotintialModification(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemPotintialModification"));
                if (bvalue) {
                    xtandemParameters.setPotentialModificationsForFullRefinment(bvalue);
                }
                bvalue = optimizeRefinePointMutations(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemPointMutations"));
                if (bvalue) {
                    xtandemParameters.setRefinePointMutations(bvalue);
                }

                Set<String> refVM = this.optimizeRefinVariableMod(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemRefVarPTM"));
                identificationParameters.getSearchParameters().getModificationParameters().clearRefinementModifications();
                for (String mod : refVM) {
                    identificationParameters.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(mod));
                }

                IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
//              
                continue;
            }
            if (param.equalsIgnoreCase("PrecursorToleranceParameter") && searchInputSetting.isOptimizePrecursorToleranceParameter()) {

                double value = this.optimizePrecursorToleranceParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("PrecursorToleranceParameter"));
                if (value != identificationParameters.getSearchParameters().getPrecursorAccuracy()) {
                    identificationParameters.getSearchParameters().setPrecursorAccuracy(value);
                    if (value > 1) {
                        identificationParameters.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
                    } else {
                        identificationParameters.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.DA);
                    }
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                }
            }
        }

        if (!digestionParameterOpt.equalsIgnoreCase(identificationParameters.getSearchParameters().getDigestionParameters().getCleavageParameter().name())) {
            optimisedSearchResults.setDigestionParameter(digestionParameterOpt);
            identificationParameters.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf(digestionParameterOpt));
            IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
        }
        if (simiEnzymaticCleavage) {
            XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            xtandemParameters.setRefineSemi(simiEnzymaticCleavage);
            IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);

        }

        for (String key
                : parameterScoreMap.keySet()) {
            System.out.println(key + "  " + parameterScoreMap.get(key));
        }

    }
    private boolean advancedParam = true;
    private boolean useRefinment = false;

    @Override
    public synchronized RawScoreModel excuteSearch(SearchingSubDataset optProtDataset, String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addSpectraList, SearchInputSetting optProtSearchSettings, File identificationParametersFile, boolean pairData) {
//        try {
        if (!optProtSearchSettings.getXTandemEnabledParameters().getParamsToOptimize().isEnabledParam(paramOption.split("_")[0])) {
            return new RawScoreModel();
        }
        if (!advancedParam && tempIdParam.getSearchParameters().getModificationParameters().getFixedModifications().size() == 1) {
            boolean terminalMod = ptmFactory.getModification(tempIdParam.getSearchParameters().getModificationParameters().getFixedModifications().get(0)).getModificationType().isCTerm() || ptmFactory.getModification(tempIdParam.getSearchParameters().getModificationParameters().getFixedModifications().get(0)).getModificationType().isNTerm();
            if (terminalMod) {
                return new RawScoreModel();
            }

        }
        if (defaultOutputFileName.contains("_resultsf_Carbamilation of protein N-term") || defaultOutputFileName.contains("resultsf_Acetylation of protein N-term")) {
            return new RawScoreModel();
        }

        SearchParameters searchParameters = tempIdParam.getSearchParameters();
        if (!advancedParam) {
            XtandemParameters xtandemParameters = (XtandemParameters) searchParameters.getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            xtandemParameters.setQuickPyrolidone(false);
            xtandemParameters.setProteinQuickAcetyl(false);
            xtandemParameters.setStpBias(false);
            xtandemParameters.setRefine(false);
        }

        Future<File> f = MainUtilities.getLongExecutorService().submit(() -> {
            File resultOutput = SearchExecuter.executeSearch(defaultOutputFileName, optProtSearchSettings, optProtDataset.getSubMsFile(), optProtDataset.getSubFastaFile(), tempIdParam, identificationParametersFile);
            return resultOutput;
        });
        File resultOutput = null;
        try {
            resultOutput = f.get();
        } catch (InterruptedException | ExecutionException ex) {
            Logger.getLogger(XTandemOptProtSearchOptimizer.class.getName()).log(Level.SEVERE, null, ex);
        }

        final List<SpectrumMatch> validatedMaches = SpectraFileUtilities.getValidatedIdentificationResults(resultOutput, optProtDataset.getSubMsFile(), Advocate.xtandem, tempIdParam);
        RawScoreModel rawScore = SpectraFileUtilities.getComparableRawScore(optProtDataset, validatedMaches, Advocate.xtandem, pairData);//(optProtDataset, resultOutput, optProtDataset.getSubMsFile(), Advocate.sage, tempIdParam, updateDataReference);

//        MainUtilities.deleteFolder(resultOutput);
        if (addSpectraList && rawScore.isSensitiveChange()) {
            rawScore.setSpectrumMatchResult(validatedMaches);
        }
        return (rawScore);

    }
//

    public double optimizeSpectrumDynamicRange(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("spectrumDR");
        Map<Double, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//         xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int spectraCounter = (int)Math.round(optProtDataset.getValidatedIdRefrenceData().length*1.01);
//        IdentificationParameters tempIdParam = oreginaltempIdParam;//IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        double selectedOption = xtandemParameters.getDynamicRange();
        for (double i = 60.0; i < 220;) {
            if (i == selectedOption) {
                i += 20;
            }
            final String option = "spectrumDR_" + i;

            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setDynamicRange(i);
            final double j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()&& scoreModel.getpValue()<0.05) {
                    if (scoreModel.getData().length <= spectraCounter) {
                        break;
                    }
                    spectraCounter = Math.max(spectraCounter, scoreModel.getData().length);
                    spectraCounter = (int)Math.round(spectraCounter*1.01);
                    resultsMap.put(j, scoreModel);

                } 
//                else if (!scoreModel.isSensitiveChange() && !scoreModel.isSameData() && scoreModel.getImprovmentScore() != -100) {
////                    break;
//                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

            i += 20;
        }
        xtandemParameters.setDynamicRange(selectedOption);
        TreeMap<RawScoreModel, Double> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());

        if (!resultsMap.isEmpty()) {
            for (double option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }

    public int optimizeSpectrumPeaksNumber(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("peaksNum");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        IdentificationParameters tempIdParam = oreginaltempIdParam;//IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        Integer selectedOption = xtandemParameters.getnPeaks();
        int spectraCounter = 0;
        for (int i = selectedOption; i <= 100;) {
            if (i == selectedOption) {
                i += 10;
            }

            final String option = "peaksNum_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setnPeaks(i);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    if (scoreModel.getData().length <= spectraCounter) {
                        break;
                    }
                    spectraCounter = Math.max(spectraCounter, scoreModel.getData().length);
                    resultsMap.put(j, scoreModel);

                } else if (!scoreModel.isSensitiveChange() && !scoreModel.isSameData() && scoreModel.getImprovmentScore() != -100) {
                    break;
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            i += 10;
        }
        xtandemParameters.setnPeaks(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);

        return selectedOption;
    }

    public double optimizeMinimumFragmentMz(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("minimumFragmentMz");

        Map<Double, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        IdentificationParameters tempIdParam = oreginaltempIdParam;
        XtandemParameters xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        int spectraCounter = 0;
        double selectedOption = xtandemParameters.getMinFragmentMz();
        for (double i = xtandemParameters.getMinFragmentMz(); i <= xtandemParameters.getMinFragmentMz() + 100.0;) {
            if (i == selectedOption) {
                i += 25;
            }
            final String option = "minimumFragmentMz_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setMinFragmentMz(i);
            final double j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    if (scoreModel.getData().length <= spectraCounter) {
                        break;
                    }
                    spectraCounter = Math.max(spectraCounter, scoreModel.getData().length);
                    resultsMap.put(j, scoreModel);

                } else if (!scoreModel.isSensitiveChange() && !scoreModel.isSameData() && scoreModel.getImprovmentScore() != -100) {
                    break;
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            i += 25;
        }
        TreeMap<RawScoreModel, Double> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        xtandemParameters.setMinFragmentMz(selectedOption);
        if (!resultsMap.isEmpty()) {
            for (double option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);

        return selectedOption;
    }

    public int optimizeMinimumPeaks(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("minpeaksNum");
        int spectraCounter = 0;
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//         xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());

        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        Integer selectedOption = xtandemParameters.getMinPeaksPerSpectrum();
        for (int i = xtandemParameters.getMinPeaksPerSpectrum(); i <= 20;) {
            if (i == selectedOption) {
                i += 1;
            }
            final String option = "minpeaksNum_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setMinPeaksPerSpectrum(i);
            final int j = i;
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange()) {
                    if (scoreModel.getData().length < spectraCounter) {
                        break;
                    }
                    spectraCounter = Math.max(spectraCounter, scoreModel.getData().length);
                    resultsMap.put(j, scoreModel);
                } else if (!scoreModel.isSensitiveChange() && !scoreModel.isSameData() && scoreModel.getImprovmentScore() != -100) {
                    break;
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            i += 1;

        }

        xtandemParameters.setMinPeaksPerSpectrum(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);

//         System.exit(0);
        return selectedOption;
    }
//

    public double optimizeNoiseSuppression(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("noiseSupression");
        Map<Double, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption1 = xtandemParameters.isUseNoiseSuppression();
        double selectedOption2 = xtandemParameters.getMinPrecursorMass();
        final String option = "noiseSupression_" + false;
        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
        xtandemParameters.setUseNoiseSuppression(false);

        if (selectedOption1 != false) {
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();

                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(0.0, scoreModel);

                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        xtandemParameters.setUseNoiseSuppression(true);
        for (double j = 500; j < 1600;) {
            final String suboption = "noiseSupression_" + true + "_" + j;
            final String subupdatedName = Configurations.DEFAULT_RESULT_NAME + "_" + suboption + "_" + msFileName;
            final double i = j;
            xtandemParameters.setMinPrecursorMass(j);
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, subupdatedName, suboption, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });

            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(i, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            j += 350;
        }
        xtandemParameters.setUseNoiseSuppression(selectedOption1);
        xtandemParameters.setMinPrecursorMass(selectedOption2);
        TreeMap<RawScoreModel, Double> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (double option2 : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option2), option2);
            }
            selectedOption2 = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption2 + "");
        parameterScoreSet.add(paramScore);

        return selectedOption2;
    }

    public boolean optimizeParentIsotopExpansion(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("parentMonoisotopicMassIsotopeError");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.getParentMonoisotopicMassIsotopeError();

        for (int i = 0; i < 2; i++) {
            if (((i == 1)) == selectedOption) {
                continue;
            }
            final String option = "parentMonoisotopicMassIsotopeError_" + (i == 1);
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setParentMonoisotopicMassIsotopeError(i == 1);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        xtandemParameters.setParentMonoisotopicMassIsotopeError(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }
    ////
    //

    public boolean optimizeQuickAcetyl(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("QuickAcetyl");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isProteinQuickAcetyl();

        for (int i = 0; i < 2; i++) {
            boolean useQuickAcetyl = (i == 1);
            if (useQuickAcetyl == selectedOption) {
                continue;
            }
            final String option = "useQuickAcetyl_" + useQuickAcetyl;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setProteinQuickAcetyl(useQuickAcetyl);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();

                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }
        xtandemParameters.setProteinQuickAcetyl(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }

        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }
//

    public boolean optimizeQuickPyrolidone(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("QuickPyrolidone");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isQuickPyrolidone();

        for (int i = 0; i < 2; i++) {
            boolean useQuickPyrolidone = (i == 1);
            if (useQuickPyrolidone == selectedOption) {
                continue;
            }
            final String option = "useQuickPyrolidone_" + useQuickPyrolidone;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setQuickPyrolidone(useQuickPyrolidone);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {

                RawScoreModel scoreModel = f.get();

                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        xtandemParameters.setQuickPyrolidone(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }
////
//

    public boolean optimizeStPBias(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("StpBias");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isStpBias();

        for (int i = 0; i < 2; i++) {
            boolean useStpBias = (i == 1);
            if (useStpBias == selectedOption) {
                continue;
            }
            final String option = "useStpBias_" + useStpBias;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setStpBias(useStpBias);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        xtandemParameters.setStpBias(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
//            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }
////

    public boolean optimizeUseRefine(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("useRefineuseRefine");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefine();

        for (int i = 0; i < 2; i++) {
            boolean useRefine = (i == 1);
            if (useRefine == selectedOption) {
                continue;
            }
            final String option = "useRefine_" + useRefine;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setRefine(useRefine);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange() || scoreModel.isSameData()) {
                    resultsMap.put(j, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        xtandemParameters.setRefine(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;

    }
////
//

    public Set<String> optimizeRefinVariableMod(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("refineVariableModifications");
        TreeMap<RawScoreModel, String> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());

        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        Map<String, RawScoreModel> twoDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        Map<String, RawScoreModel> threeDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());

        Map<String, RawScoreModel> fourDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());

        TreeSet<RawScoreModel> sorterSet = new TreeSet<>(Collections.reverseOrder());

        MainUtilities.cleanOutputFolder();
        paramScore.setParamId("refineVariableModifications");//           
        for (String vMod : oreginaltempIdParam.getSearchParameters().getModificationParameters().getVariableModifications()) {
            if (ptmFactory.getModification(vMod) == null) {
                continue;
            }
            oreginaltempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
            oreginaltempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod));

            final String option = vMod;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, searchInputSetting, generatedIdentificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                System.out.println("var ref " + vMod + "  " + scoreModel);
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(option, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }
        //2d refinment
        if (resultsMap.size() > 1) {
            for (String selectedRef : resultsMap.keySet()) {
                for (String vMod : resultsMap.keySet()) {
                    if (vMod.equals(selectedRef)) {
                        continue;
                    }
                    oreginaltempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                    oreginaltempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod));
                    oreginaltempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(selectedRef));
                    final String option = vMod + "_" + selectedRef;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;

                    Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                        RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, searchInputSetting, generatedIdentificationParametersFile, true);
                        return scoreModel;
                    });
                    try {
                        RawScoreModel scoreModel = f.get();
                        if (scoreModel.isSignificatChange()) {
                            sorterSet.add(scoreModel);
                            sorterSet.add(resultsMap.get(selectedRef));
                            sorterSet.add(resultsMap.get(vMod));
                            if (sorterSet.first() == scoreModel) {
                                twoDResultsMap.put(option, scoreModel);
                            }
                            sorterSet.clear();

                        }
                    } catch (ExecutionException | InterruptedException ex) {
                        ex.printStackTrace();
                    }
                }
            }
        }
        //3d refinment
        if (!twoDResultsMap.isEmpty() && resultsMap.size() > 2) {
            for (String selectedRef : resultsMap.keySet()) {
                for (String vMod : twoDResultsMap.keySet()) {
                    if (vMod.contains(selectedRef)) {
                        continue;
                    }
                    oreginaltempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                    oreginaltempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[0]));
                    oreginaltempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[1]));
                    oreginaltempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(selectedRef));
                    final String option = vMod + "_" + selectedRef;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;
                    Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                        RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, searchInputSetting, generatedIdentificationParametersFile, true);
                        return scoreModel;
                    });
                    try {
                        RawScoreModel scoreModel = f.get();
                        if (scoreModel.isSignificatChange()) {
                            sorterSet.add(scoreModel);
                            sorterSet.add(resultsMap.get(selectedRef));
                            sorterSet.add(twoDResultsMap.get(vMod));
                            if (sorterSet.first() == scoreModel) {
                                threeDResultsMap.put(option, scoreModel);
                            }
                            sorterSet.clear();

                        }
                    } catch (ExecutionException | InterruptedException ex) {
                        ex.printStackTrace();
                    }

                }
            }

        }
        //4d refinment
        if (!threeDResultsMap.isEmpty() && resultsMap.size() > 3) {
            for (String selectedRef : resultsMap.keySet()) {
                for (String vMod : threeDResultsMap.keySet()) {
                    if (vMod.equals(selectedRef)) {
                        continue;
                    }
                    oreginaltempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                    oreginaltempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[0]));
                    oreginaltempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[1]));
                    oreginaltempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[2]));
                    oreginaltempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(selectedRef));
                    final String option = vMod + "_" + selectedRef;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;
                    Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                        RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, searchInputSetting, generatedIdentificationParametersFile, true);
                        return scoreModel;
                    });
                    try {
                        RawScoreModel scoreModel = f.get();
                        if (scoreModel.isSignificatChange()) {
                            sorterSet.add(scoreModel);
                            sorterSet.add(resultsMap.get(selectedRef));
                            sorterSet.add(threeDResultsMap.get(vMod));
                            if (sorterSet.first() == scoreModel) {
                                fourDResultsMap.put(option, scoreModel);
                            }
                            sorterSet.clear();

                        }
                    } catch (ExecutionException | InterruptedException ex) {
                        ex.printStackTrace();
                    }
                }
            }

        }
        resultsMap.clear();
        resultsMap.putAll(resultsMap);
        resultsMap.putAll(twoDResultsMap);
        resultsMap.putAll(threeDResultsMap);
        resultsMap.putAll(fourDResultsMap);
        sortedResultsMap.clear();
        for (String key : resultsMap.keySet()) {
            sortedResultsMap.put(resultsMap.get(key), key);
        }

        Set<String> refinementVarModMap = new HashSet<>();
        if (!sortedResultsMap.isEmpty()) {
            String varMod = sortedResultsMap.firstEntry().getValue();
            if (varMod.contains("_")) {
                refinementVarModMap.addAll(Arrays.asList(varMod.split("_")));
            } else {
                refinementVarModMap.add(varMod);
            }
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(refinementVarModMap.toString());
        parameterScoreSet.add(paramScore);
        return refinementVarModMap;

    }

    public boolean optimizeRefineUnanticipatedCleavage(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("UnanticipatedCleavages");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        xtandemParameters.setRefine(true);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefineUnanticipatedCleavages();
        MainUtilities.resetExecutorService();
        for (int i = 0; i < 2; i++) {
            boolean useRefineUnanticipatedCleavages = (i == 1);
//            if (useRefineUnanticipatedCleavages == selectedOption) {
//                continue;
//            }
            final String option = "useRefineUnanticipatedCleavages_" + useRefineUnanticipatedCleavages;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setRefineUnanticipatedCleavages(useRefineUnanticipatedCleavages);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {

                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        xtandemParameters.setRefineUnanticipatedCleavages(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }

        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }
////
//

    public boolean optimizeRefineSimiEnzymaticCleavage(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("SimiEnzymaticCleavage");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefineSemi();
        MainUtilities.resetExecutorService();
        for (int i = 0; i < 2; i++) {
            boolean useRefineSimiEnzymaticCleavage = (i == 1);
            if (useRefineSimiEnzymaticCleavage == selectedOption) {
                continue;
            }
            final String option = "useRefineSimiEnzymaticCleavage_" + useRefineSimiEnzymaticCleavage;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setRefineSemi(useRefineSimiEnzymaticCleavage);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                System.out.println("simi inzematic thingy " + scoreModel + "  " + optProtDataset.getValidatedIdRefrenceData().length * 1.1);
                if (scoreModel.isSignificatChange() && scoreModel.getTotalNumber() > (optProtDataset.getValidatedIdRefrenceData().length * 1.1)) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setComments("Cause slow searching");
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        xtandemParameters.setRefineSemi(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);

        return selectedOption;

    }
////
//

    public boolean optimizePotintialModification(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("PotintialModification");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        xtandemParameters.setRefine(true);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isPotentialModificationsForFullRefinment();
        MainUtilities.resetExecutorService();
        for (int i = 0; i < 2; i++) {
            boolean usePotintialModification = (i == 1);
            if (usePotintialModification == selectedOption) {
                continue;
            }
            final String option = "usePotintialModification_" + usePotintialModification;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setPotentialModificationsForFullRefinment(usePotintialModification);
            final int j = i;
            MainUtilities.resetExecutorService();
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
              
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            MainUtilities.resetExecutorService();
        }
        xtandemParameters.setPotentialModificationsForFullRefinment(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;

    }
////
//

    public boolean optimizeRefinePointMutations(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("PointMutations");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        xtandemParameters.setRefine(true);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefinePointMutations();
        MainUtilities.resetExecutorService();
        for (int i = 0; i < 2; i++) {
            boolean useRefinePointMutations = (i == 1);
            if (useRefinePointMutations == selectedOption) {
                continue;
            }
            final String option = "useRefinePointMutations_" + useRefinePointMutations;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setRefinePointMutations(useRefinePointMutations);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        xtandemParameters.setRefinePointMutations(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }

    public boolean optimizeRefineSnAPs(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("SnAPs");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefineSnaps();

        for (int i = 0; i < 2; i++) {
            boolean useRefineSnAPs = (i == 1);
            if (useRefineSnAPs == selectedOption) {
                continue;
            }
            final String option = "useRefineSnAPs_" + useRefineSnAPs;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setRefineSnaps(useRefineSnAPs);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        xtandemParameters.setRefineSnaps(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;

    }
//

    public boolean optimizeRefineSpectrumSynthesis(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("SpectrumSynthesis");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefineSpectrumSynthesis();
        for (int i = 0; i < 2; i++) {
            boolean useRefineSpectrumSynthesis = (i == 1);
            if (useRefineSpectrumSynthesis == selectedOption) {
                continue;
            }
            final String option = "useRefineSpectrumSynthesis_" + useRefineSpectrumSynthesis;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setRefineSpectrumSynthesis(useRefineSpectrumSynthesis);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        xtandemParameters.setRefineSpectrumSynthesis(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;

    }

}
