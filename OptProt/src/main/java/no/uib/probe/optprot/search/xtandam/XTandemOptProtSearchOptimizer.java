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
        optProtDataset.setParameterScoreMap(parameterScoreMap);
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
        
        XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        xtandemParameters.setRefinePointMutations(false);
        xtandemParameters.setProteinPtmComplexity(6);
        
    }
    
    private String digestionParameterOpt;
    private boolean simiEnzymaticCleavage = false;
    private String enzymeSpecificityOpt = "specific";
    
    public void startProcess(List<String> paramOrder) throws IOException {
        digestionParameterOpt = identificationParameters.getSearchParameters().getDigestionParameters().getCleavageParameter().name();
        searchInputSetting.setDigestionParameterOpt(digestionParameterOpt);
        for (String param : paramOrder) {
            System.out.println("-------------------------------------------param " + param + "-------------------------------------------");
            if (param.equalsIgnoreCase("DigestionParameter_1") && searchInputSetting.isOptimizeDigestionParameter()) {
                String[] values = this.optimizeEnzymeParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("EnzymeParameter"));
                identificationParameters.getSearchParameters().getDigestionParameters().clearEnzymes();
                if (!values[0].equalsIgnoreCase("")) {
                    optimisedSearchResults.setEnzymeName(values[0]);
                    int nMissesCleavages = Integer.parseInt(values[2]);                 
                    identificationParameters.getSearchParameters().getDigestionParameters().addEnzyme(EnzymeFactory.getInstance().getEnzyme(values[0]));
                    enzymeSpecificityOpt = values[1];
                    identificationParameters.getSearchParameters().getDigestionParameters().setnMissedCleavages(values[0], nMissesCleavages);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                }
                MainUtilities.cleanOutputFolder();
                continue;
            }
            if (param.equalsIgnoreCase("DigestionTypeParameter") && searchInputSetting.isOptimizeDigestionParameter()) {
                digestionParameterOpt = this.optimizeDigestionCleavageParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("DigestionParameter"));
                searchInputSetting.setDigestionParameterOpt(digestionParameterOpt);
                MainUtilities.cleanOutputFolder();
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
                MainUtilities.cleanOutputFolder();
                continue;
                
            }

//confusing param
            if (param.equalsIgnoreCase("FragmentToleranceParameter") && searchInputSetting.isOptimizeFragmentToleranceParameter()) {
                double value = this.optimizeFragmentToleranceParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("FragmentToleranceParameter"));
                if (value != identificationParameters.getSearchParameters().getFragmentIonAccuracy()) {
                    identificationParameters.getSearchParameters().setFragmentIonAccuracy(value);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                }
                MainUtilities.cleanOutputFolder();
                continue;
            }
            
            if (param.equalsIgnoreCase("PrecursorChargeParameter") && searchInputSetting.isOptimizePrecursorChargeParameter()) {
                
                int[] values = this.optimizePrecursorChargeParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("PrecursorChargeParameter"));
                if (values[1] != identificationParameters.getSearchParameters().getMaxChargeSearched()) {
                    identificationParameters.getSearchParameters().setMinChargeSearched(values[0]);
                    identificationParameters.getSearchParameters().setMaxChargeSearched(values[1]);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                }
                MainUtilities.cleanOutputFolder();
                continue;
            }
            
            if (param.equalsIgnoreCase("XtandemAdvancedParameter") && searchInputSetting.isOptimizeXtandemAdvancedParameter()) {
                XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
                useRefinment = optimizeUseRefine(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemUseRefine"));
                if (useRefinment != xtandemParameters.isRefine()) {
                    xtandemParameters.setRefine(useRefinment);
                }
                if (!useRefinment) {
                    System.out.println("Error----->>> disable second stage");                    
                }
                boolean bvalue = optimizeQuickAcetyl(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemQuickAcetyl"));
                if (bvalue) {
                    xtandemParameters.setProteinQuickAcetyl(bvalue);
                }
                bvalue = optimizeQuickPyrolidone(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemQuickPyrolidone"));
                if (bvalue) {
                    xtandemParameters.setQuickPyrolidone(bvalue);
                }
                IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                
                MainUtilities.cleanOutputFolder();
                continue;
            }
            
            if (useRefinment && param.equalsIgnoreCase("XtandemAdvancedParameter_A") && searchInputSetting.isOptimizeXtandemAdvancedParameter()) {
//                advancedParam = true;
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
                
                MainUtilities.cleanOutputFolder();
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
                MainUtilities.cleanOutputFolder();
                continue;
            }
            if (param.equalsIgnoreCase("XtandemAdvancedParameter_B") && searchInputSetting.isOptimizeXtandemAdvancedParameter()) {
                
                XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
                
                double dvalue = optimizeSpectrumDynamicRange(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemSpectrumDynamicRange"));
                if (dvalue != xtandemParameters.getDynamicRange()) {
                    xtandemParameters.setDynamicRange(dvalue);
                }
                boolean bvalue = optimizePotintialModification(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemPotintialModification"));
                if (bvalue) {
                    xtandemParameters.setPotentialModificationsForFullRefinment(bvalue);
                    Set<String> refVM = this.optimizeRefinVariableMod(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemRefVarPTM"));
                    identificationParameters.getSearchParameters().getModificationParameters().clearRefinementModifications();
                    for (String mod : refVM) {
                        identificationParameters.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(mod));
                    }
                }
                bvalue = optimizeRefinePointMutations(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("XtandemPointMutations"));
                if (bvalue) {
                    xtandemParameters.setRefinePointMutations(bvalue);
                }
                
                IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
//              
                MainUtilities.cleanOutputFolder();
                continue;
            }
//            if (param.equalsIgnoreCase("PrecursorToleranceParameter") && searchInputSetting.isOptimizePrecursorToleranceParameter()) {
//
//                double value = this.optimizePrecursorToleranceParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("PrecursorToleranceParameter"));
//                if (value != identificationParameters.getSearchParameters().getPrecursorAccuracy()) {
//                    identificationParameters.getSearchParameters().setPrecursorAccuracy(value);
//                    if (value > 1) {
//                        identificationParameters.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
//                    } else {
//                        identificationParameters.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.DA);
//                    }
//                    IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
//                }
//            }
        }
        
        if (!digestionParameterOpt.equalsIgnoreCase(identificationParameters.getSearchParameters().getDigestionParameters().getCleavageParameter().name())) {
            optimisedSearchResults.setDigestionParameter(digestionParameterOpt);
            identificationParameters.getSearchParameters().getDigestionParameters().clearEnzymes();
            identificationParameters.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf(digestionParameterOpt));
            IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
        }
        if (simiEnzymaticCleavage) {
            XtandemParameters xtandemParameters = (XtandemParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            xtandemParameters.setRefineSemi(simiEnzymaticCleavage);
            IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
            
        }
        if (!enzymeSpecificityOpt.equalsIgnoreCase("specific")) {
            
            identificationParameters.getSearchParameters().getDigestionParameters().setSpecificity(identificationParameters.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName(), DigestionParameters.Specificity.valueOf(enzymeSpecificityOpt));
            IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
            
        }
        
        for (String key
                : parameterScoreMap.keySet()) {
            System.out.println(key + "  " + parameterScoreMap.get(key));
        }
        
    }
    private boolean useRefinment = false;
    
    @Override
    public synchronized RawScoreModel excuteSearch(SearchingSubDataset optProtDataset, String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addSpectraList, SearchInputSetting optProtSearchSettings, File identificationParametersFile, boolean pairData) {
        
        if (!optProtSearchSettings.getXTandemEnabledParameters().getParamsToOptimize().isEnabledParam(paramOption.split("_")[0])) {
            System.out.println(" param: " + paramOption + "  not supporten by xtandem");
            return new RawScoreModel(paramOption);
        }
        if (defaultOutputFileName.contains("_resultsf_Carbamilation of protein N-term") || defaultOutputFileName.contains("resultsf_Acetylation of protein N-term")) {
            System.out.println(" param: " + paramOption + "  not supporten by xtandem");
            return new RawScoreModel(paramOption);
        }
        
        if (paramOption.contains("Pyrolidone from")) {
            XtandemParameters xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            if (xtandemParameters.isQuickPyrolidone()) {
                System.out.println(" quick pyro applied : " + paramOption + " ");
                return new RawScoreModel(paramOption);
            }
        }
        if (paramOption.contains("Acetylation of protein N-term")) {
            XtandemParameters xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
            if (xtandemParameters.isProteinQuickAcetyl()) {
                System.out.println(" quick acetylationapplied : " + paramOption + " ");
                return new RawScoreModel(paramOption);
            }
        }
        if (paramOption.contains("_")) {
            paramOption = paramOption.split("_")[1];
        }
        
        Future<File> f = MainUtilities.getLongExecutorService().submit(() -> {
            File resultOutput = SearchExecuter.executeSearch(defaultOutputFileName, optProtSearchSettings, optProtDataset.getSubMsFile(), optProtDataset.getSubFastaFile(), tempIdParam, identificationParametersFile);
            return resultOutput;
        });
        File resultOutput = null;
        try {
            while (!f.isDone()) {
//                System.out.print("----------------");
            }
            resultOutput = f.get();
        } catch (InterruptedException | ExecutionException ex) {
            ex.printStackTrace();
        }
        
        final List<SpectrumMatch> validatedMaches = SpectraUtilities.getValidatedIdentificationResults(resultOutput, optProtDataset.getSubMsFile(), Advocate.xtandem, tempIdParam);
        RawScoreModel rawScore = SpectraUtilities.getComparableRawScore(optProtDataset, validatedMaches, Advocate.xtandem, addSpectraList, paramOption);//(optProtDataset, resultOutput, optProtDataset.getSubMsFile(), Advocate.sage, tempIdParam, updateDataReference);

        if (addSpectraList && rawScore.isSensitiveChange()) {
            rawScore.setSpectrumMatchResult(validatedMaches);
        }
        return (rawScore);
        
    }
    
    public double optimizeSpectrumDynamicRange(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("spectrumDR");
        Map<Double, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int spectraCounter = (int) Math.round(optProtDataset.getActiveIdentificationNum() * 1.01);
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
                if (scoreModel.isSignificatChange()) {
                    if (scoreModel.getSpectrumMatchResult().size() <= spectraCounter) {
                        break;
                    }
                    spectraCounter = Math.max(spectraCounter, scoreModel.getSpectrumMatchResult().size());
                    spectraCounter = (int) Math.round(spectraCounter * 1.01);
                    resultsMap.put(j, scoreModel);
                    
                }
                
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
            double impact = Math.round((double) (resultsMap.get(selectedOption).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
            
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }
    
    public int optimizeSpectrumPeaksNumber(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("peaksNum");
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        Integer selectedOption = xtandemParameters.getnPeaks();
        double topScore = 0;
        for (int i = 50; i <= 100;) {
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
                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j + "", scoreModel);
                    if (scoreModel.getFinalScore() <= topScore) {
                        break;
                    }
                    topScore = scoreModel.getFinalScore();
                    
                } else if (!scoreModel.isSensitiveChange() && !scoreModel.isSameData() && scoreModel.getImprovmentScore() != -100 && topScore > 0) {
                    break;
                }
                
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            i += 10;
        }
        xtandemParameters.setnPeaks(selectedOption);
        if (!resultsMap.isEmpty()) {
            selectedOption = Integer.valueOf(SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber()));
            double impact = Math.round((double) (resultsMap.get(selectedOption + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedOption + ""));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        
        return selectedOption;
    }
    
    public double optimizeMinimumFragmentMz(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("minimumFragmentMz");
        
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        IdentificationParameters tempIdParam = oreginaltempIdParam;
        XtandemParameters xtandemParameters = (XtandemParameters) tempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        int spectraCounter = optProtDataset.getActiveIdentificationNum();
        double selectedOption = xtandemParameters.getMinFragmentMz();
        for (double i = 100; i <= 300;) {
            if (i == selectedOption) {
                i += 50;
            }
            final String option = "minimumFragmentMz_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setMinFragmentMz(i);
            final double j = i;
            
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, optimisedSearchParameter, generatedIdentificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j + "", scoreModel);
                    
                }
                
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            i += 50;
        }
        
        if (!resultsMap.isEmpty()) {
            selectedOption = Double.parseDouble(SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber()));
            double impact = Math.round((double) (resultsMap.get(selectedOption + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedOption + ""));
        }
        
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        
        return selectedOption;
    }
    
    public int optimizeMinimumPeaks(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("minpeaksNum");
        double lastScore = 0.0;
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//         xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());

        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        Integer selectedOption = xtandemParameters.getMinPeaksPerSpectrum();
        for (int i = xtandemParameters.getMinPeaksPerSpectrum(); i <= xtandemParameters.getnPeaks(); i++) {
            if (i == selectedOption) {
                i += 1;
            }
            final String option = "minpeaksNum_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setMinPeaksPerSpectrum(i);
            final int j = i;
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange()) {
                    if (scoreModel.getFinalScore() <= lastScore) {
                        break;
                    }
                    lastScore = scoreModel.getFinalScore();
                    resultsMap.put(j, scoreModel);
                } else if (!scoreModel.isSensitiveChange() && scoreModel.getImprovmentScore() != -100) {
                    break;
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        if (!resultsMap.isEmpty()) {
            TreeSet<Integer> sorter = new TreeSet<>(resultsMap.keySet());
            selectedOption = sorter.last();
            double impact = Math.round((double) (resultsMap.get(selectedOption).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedOption));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        
        return selectedOption;
    }
    
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
            double impact = Math.round((double) (resultsMap.get(selectedOption2).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption2 + "");
        parameterScoreSet.add(paramScore);
        
        return selectedOption2;
    }
    
    public boolean optimizeParentIsotopExpansion(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        MainUtilities.cleanOutputFolder();
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
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
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
        xtandemParameters.setParentMonoisotopicMassIsotopeError(selectedOption);
        TreeMap<RawScoreModel, Integer> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue() == 1;
            double impact = Math.round((double) (resultsMap.get(1).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }
    
    public boolean optimizeQuickAcetyl(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("QuickAcetyl");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
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
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.getFinalScore() > 0 || scoreModel.isSameData() || scoreModel.isSensitiveChange()) {
                    resultsMap.put(j, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            
        }
        xtandemParameters.setProteinQuickAcetyl(selectedOption);
        if (!resultsMap.isEmpty()) {
            boolean use = selectedOption;
            for (int option : resultsMap.keySet()) {
                use = (option == 1);
                double impact = Math.round((double) (resultsMap.get(option).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                optProtDataset.setActiveScoreModel(resultsMap.get(option));
            }
            selectedOption = use;
            
        }
        
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }
    
    public boolean optimizeQuickPyrolidone(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("QuickPyrolidone");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
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
        if (!resultsMap.isEmpty()) {
            boolean use = selectedOption;
            for (int option : resultsMap.keySet()) {
                use = (option == 1);
                double impact = Math.round((double) (resultsMap.get(option).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                optProtDataset.setActiveScoreModel(resultsMap.get(option));
            }
            selectedOption = use;
            
        }
        
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }
    
    public boolean optimizeStPBias(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("StpBias");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
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
        if (!resultsMap.isEmpty()) {
            boolean use = selectedOption;
            for (int option : resultsMap.keySet()) {
                use = (option == 1);
                double impact = Math.round((double) (resultsMap.get(option).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                optProtDataset.setActiveScoreModel(resultsMap.get(option));
            }
            selectedOption = use;
            
        }
        
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }
////

    public boolean optimizeUseRefine(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("useRefineuseRefine");
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefine();
        boolean useRefine = true;
        final String option = "useRefine_" + useRefine;
        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
        xtandemParameters.setRefine(useRefine);
        Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, generatedIdentificationParametersFile, true);
            return scoreModel;
        });
        try {
            RawScoreModel scoreModel = f.get();
            double impact = Math.round((double) (scoreModel.getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(scoreModel);
            selectedOption = true;
            paramScore.setScore(scoreModel.getSpectrumMatchResult().size());
            paramScore.setParamValue(selectedOption + "");
            parameterScoreSet.add(paramScore);
        } catch (ExecutionException | InterruptedException ex) {
            ex.printStackTrace();
        }
        
        return selectedOption;
        
    }
    
    public Set<String> optimizeRefinVariableMod(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("refineVariableModifications");
        TreeMap<RawScoreModel, String> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        Map<String, RawScoreModel> twoDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        Map<String, RawScoreModel> threeDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        
        Map<String, RawScoreModel> fourDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        
        TreeSet<RawScoreModel> sorterSet = new TreeSet<>(Collections.reverseOrder());
        
        MainUtilities.cleanOutputFolder();
        paramScore.setParamId("refineVariableModifications");//    

        Set<String> potintialVariableMod = new LinkedHashSet<>(oreginaltempIdParam.getSearchParameters().getModificationParameters().getVariableModifications());
//        }
        for (String vMod : potintialVariableMod) {
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
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(refinementVarModMap.toString());
        parameterScoreSet.add(paramScore);
        return refinementVarModMap;
        
    }
    
    public boolean optimizeRefineUnanticipatedCleavage(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("UnanticipatedCleavages");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        XtandemParameters xtandemParameters = (XtandemParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.xtandem.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = xtandemParameters.isRefineUnanticipatedCleavages();
        MainUtilities.resetExecutorService();
        for (int i = 0; i < 2; i++) {
            boolean useRefineUnanticipatedCleavages = (i == 1);
            if (useRefineUnanticipatedCleavages == selectedOption) {
                continue;
            }
            final String option = "useRefineUnanticipatedCleavages_" + useRefineUnanticipatedCleavages;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            xtandemParameters.setRefineUnanticipatedCleavages(useRefineUnanticipatedCleavages);
            final int j = i;
            
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
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
        if (!resultsMap.isEmpty()) {
            boolean use = selectedOption;
            for (int option : resultsMap.keySet()) {
                use = (option == 1);
                double impact = Math.round((double) (resultsMap.get(option).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                optProtDataset.setActiveScoreModel(resultsMap.get(option));
            }
            selectedOption = use;
            
        }
        
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }
    
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
                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j, scoreModel);
                    paramScore.setComments("Cause slow searching");
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        xtandemParameters.setRefineSemi(selectedOption);
        
        if (!resultsMap.isEmpty()) {
            boolean use = selectedOption;
            for (int option : resultsMap.keySet()) {
                use = (option == 1);
                double impact = Math.round((double) (resultsMap.get(option).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                optProtDataset.setActiveScoreModel(resultsMap.get(option));
            }
            selectedOption = use;
            
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
        
    }
    
    public boolean optimizePotintialModification(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("PotintialModification");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
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
                
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
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
        if (!resultsMap.isEmpty()) {
            boolean use = selectedOption;
            for (int option : resultsMap.keySet()) {
                use = (option == 1);
                double impact = Math.round((double) (resultsMap.get(option).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                optProtDataset.setActiveScoreModel(resultsMap.get(option));
            }
            selectedOption = use;
            
        }
        
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        
        parameterScoreSet.add(paramScore);
        return selectedOption;
        
    }
    
    public boolean optimizeRefinePointMutations(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("PointMutations");
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
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
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                System.out.println(optProtDataset.getActiveIdentificationNum() + " - use point mutation " + scoreModel);
                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j, scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        xtandemParameters.setRefinePointMutations(selectedOption);
        if (!resultsMap.isEmpty()) {
            boolean use = selectedOption;
            for (int option : resultsMap.keySet()) {
                use = (option == 1);
                
                double impact = Math.round((double) (resultsMap.get(option).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                optProtDataset.setActiveScoreModel(resultsMap.get(option));
            }
            selectedOption = use;
            
        }
        
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
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
                
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
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
        if (!resultsMap.isEmpty()) {
            boolean use = selectedOption;
            for (int option : resultsMap.keySet()) {
                use = (option == 1);
                double impact = Math.round((double) (resultsMap.get(option).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                optProtDataset.setActiveScoreModel(resultsMap.get(option));
            }
            selectedOption = use;
            
        }
        
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
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
                
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
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
        xtandemParameters.setRefineSpectrumSynthesis(selectedOption);
        
        if (!resultsMap.isEmpty()) {
            boolean use = selectedOption;
            for (int option : resultsMap.keySet()) {
                use = (option == 1);
                
                double impact = Math.round((double) (resultsMap.get(option).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                optProtDataset.setActiveScoreModel(resultsMap.get(option));
            }
            selectedOption = use;
            
        }
        
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
        
    }
    
}
