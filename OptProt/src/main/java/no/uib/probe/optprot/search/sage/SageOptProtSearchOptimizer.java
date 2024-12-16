package no.uib.probe.optprot.search.sage;

import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.io.IoUtil;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
//import com.compomics.util.parameters.identification.tool_specific.MyriMatchParameters;
import com.compomics.util.parameters.identification.tool_specific.SageParameters;
import com.compomics.util.parameters.identification.tool_specific.XtandemParameters;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
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
     * The compomics PTM factory.
     */
    private final ModificationFactory ptmFactory = ModificationFactory.getInstance();

    private final SearchingSubDataset optProtDataset;
    private final SearchInputSetting searchInputSetting;
    private final File generatedIdentificationParametersFile;
    private final OptimisedSearchResults optimisedSearchResults;
    private final IdentificationParameters identificationParameters;
    private final Map<String, TreeSet<ParameterScoreModel>> parameterScoreMap;

    public Map<String, TreeSet<ParameterScoreModel>> getParameterScoreMap() {
        return parameterScoreMap;
    }

    public SageOptProtSearchOptimizer(SearchingSubDataset optProtDataset, SearchInputSetting searchInputSetting, File generatedIdentificationParametersFile) throws IOException {

        this.optProtDataset = optProtDataset;
        this.searchInputSetting = searchInputSetting;
        this.generatedIdentificationParametersFile = generatedIdentificationParametersFile;
        this.identificationParameters = IdentificationParameters.getIdentificationParameters(generatedIdentificationParametersFile);
        this.optimisedSearchResults = new OptimisedSearchResults();
        this.parameterScoreMap = new LinkedHashMap<>();
        optProtDataset.setParameterScoreMap(parameterScoreMap);
        MainUtilities.cleanOutputFolder();
       if(searchInputSetting.isOptimizeAllParameters()){
        SageParameters sageParameters = (SageParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        sageParameters.setMaxVariableMods(2);
        sageParameters.setNumPsmsPerSpectrum(1);
        sageParameters.setGenerateDecoys(false);
        sageParameters.setWideWindow(false);
        sageParameters.setPredictRt(true);}
//        System.out.println(" identificationParameters.getFastaParameters().getDecoyFlag() " + identificationParameters.getFastaParameters().getDecoyFlag() + "  oreginal size  " + optProtDataset.getOreginalDatasize() + "  total subsize " + optProtDataset.getTotalSpectraNumber());
//        sageParameters.setMinFragmentMz(150.0);
//        sageParameters.setMaxFragmentMz(1500.0);
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
        parameterScoreMap.put("SageChimericSpectraParameter", new TreeSet<>(Collections.reverseOrder()));//
        parameterScoreMap.put("SageWideWindowParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SagePredectRetentionTimeParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageNumberOfPeakParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageMinMatchedPeaksPeakParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageMaxFragmentChargeParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageNumberOfPeakParameter", new TreeSet<>(Collections.reverseOrder()));
        parameterScoreMap.put("SageBatchSizeParameter", new TreeSet<>(Collections.reverseOrder()));

    }

    private String digestionParameterOpt;
    private boolean simiEnzymaticCleavage = false;
    private String enzymeSpecificityOpt = "specific";

    public void startProcess(List<String> paramOrder) throws IOException {
        digestionParameterOpt = identificationParameters.getSearchParameters().getDigestionParameters().getCleavageParameter().name();
        searchInputSetting.setDigestionParameterOpt(digestionParameterOpt);
         MainUtilities.cleanOutputFolder();
//        if (!searchInputSetting.isOptimizeAllParameters()) {
            //run refrence search 
            runReferenceRun(optProtDataset, identificationParameters, searchInputSetting);
//        }

        for (String param : paramOrder) {
            System.out.println("-------------------------------------------param " + param + "-------------------------------------------");
            MainUtilities.cleanOutputFolder();
            if (param.equalsIgnoreCase("DigestionParameter_1") && searchInputSetting.isOptimizeDigestionParameter()) {
                String[] values = this.optimizeEnzymeParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("EnzymeParameter"));
                identificationParameters.getSearchParameters().getDigestionParameters().clearEnzymes();
                if (!values[0].equalsIgnoreCase("")) {
                    optimisedSearchResults.setEnzymeName(values[0]);
                    int nMissesCleavages = Integer.parseInt(values[2]);// identificationParameters.getSearchParameters().getDigestionParameters().getnMissedCleavages(value);                   
                    identificationParameters.getSearchParameters().getDigestionParameters().addEnzyme(EnzymeFactory.getInstance().getEnzyme(values[0]));
                    enzymeSpecificityOpt = values[1];
//                    identificationParameters.getSearchParameters().getDigestionParameters().setSpecificity(values[0], DigestionParameters.Specificity.valueOf(values[1]));
                    identificationParameters.getSearchParameters().getDigestionParameters().setnMissedCleavages(values[0], nMissesCleavages);
                    IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);
                }

                continue;
            }
            if (param.equalsIgnoreCase("DigestionTypeParameter") && searchInputSetting.isOptimizeDigestionParameter()) {
                digestionParameterOpt = this.optimizeDigestionCleavageParameter(optProtDataset, generatedIdentificationParametersFile, searchInputSetting, parameterScoreMap.get("DigestionParameter"));
                searchInputSetting.setDigestionParameterOpt(digestionParameterOpt);
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
            if (param.equalsIgnoreCase("SageAdvancedParameter_A") && searchInputSetting.isOptimizeSageAdvancedParameter()) {
                SageParameters sageParameters = (SageParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
                sageParameters.setMaxVariableMods(0);
                double[] dvalues = optimizeFragmentMzParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SageFragmentMzParameter"));
                if (dvalues[1] != sageParameters.getMaxFragmentMz() || dvalues[0] != sageParameters.getMinFragmentMz()) {
                    sageParameters.setMinFragmentMz(dvalues[0]);
                    sageParameters.setMaxFragmentMz(dvalues[1]);
                    System.out.println("optimizeFragmentMzParameter " + dvalues[0] + "--" + dvalues[1] + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
                }
                int value = optimizeIonMinIndexParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SageIonMinIndexParameter"));
                if (value != sageParameters.getMinIonIndex()) {
                    sageParameters.setMinIonIndex(value);
                    System.out.println("optimizeIonMinIndexParameter " + value + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
                }
                int[] values = optimizePeptideLengthParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SagePeptideLengthParameter"));
                if (values[1] != sageParameters.getMaxPeptideLength() || values[0] != sageParameters.getMinPeptideLength()) {
                    sageParameters.setMinPeptideLength(values[0]);
                    sageParameters.setMaxPeptideLength(values[1]);
                    System.out.println("peptide length " + values[0] + "--" + values[1] + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
                }
//                dvalues = optimizePeptideMassParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SagePeptideMassParameter"));
//                if (dvalues[1] != sageParameters.getMaxPeptideMass() || dvalues[0] != sageParameters.getMinPeptideMass()) {
//                    sageParameters.setMinPeptideMass(dvalues[0]);
//                    sageParameters.setMaxPeptideMass(dvalues[1]);
//                    System.out.println("optimizePeptideMassParameter " + dvalues[0] + "--" + dvalues[1] + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
//                }
//

                sageParameters.setMaxVariableMods(2);
                IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);

                continue;
            }

            if (param.equalsIgnoreCase("SageAdvancedParameter_B") && searchInputSetting.isOptimizeSageAdvancedParameter()) {
                SageParameters sageParameters = (SageParameters) identificationParameters.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
                
                
                /**maybe remove**/
//                  int value = optimizeIonMinIndexParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SageIonMinIndexParameter"));
//                if (value != sageParameters.getMinIonIndex()) {
//                    sageParameters.setMinIonIndex(value);
//                    System.out.println("optimizeIonMinIndexParameter " + value + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
//                }
            double[]  dvalues = optimizePeptideMassParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SagePeptideMassParameter"));
                if (dvalues[1] != sageParameters.getMaxPeptideMass() || dvalues[0] != sageParameters.getMinPeptideMass()) {
                    sageParameters.setMinPeptideMass(dvalues[0]);
                    sageParameters.setMaxPeptideMass(dvalues[1]);
                    System.out.println("optimizePeptideMassParameter " + dvalues[0] + "--" + dvalues[1] + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
                }
                /****until here**/
                
                
                int[] values = optimizeNumberOfPeakParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SageNumberOfPeakParameter"));
                if (values[1] != sageParameters.getMaxPeaks() || values[0] != sageParameters.getMinPeaks()) {
                    sageParameters.setMinPeaks(values[0]);
                    sageParameters.setMaxPeaks(values[1]);
                }

                boolean valueBoolean = optimizeDeisotopParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SageDeisotopParameter"));
                if (valueBoolean != sageParameters.getDeisotope()) {
                    sageParameters.setDeisotope(valueBoolean);
                }
                System.out.println("optimizeDeisotopParameter " + valueBoolean + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
                valueBoolean = optimizeChimericSpectraParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SageChimericSpectraParameter"));
                if (valueBoolean != sageParameters.getChimera()) {
                    sageParameters.setChimera(valueBoolean);
                }
                System.out.println("SageChimericSpectraParameter " + valueBoolean + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
                valueBoolean = optimizePredectRetentionTimeParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SagePredectRetentionTimeParameter"));
                if (valueBoolean != sageParameters.getPredictRt()) {
                    sageParameters.setPredictRt(valueBoolean);
                }
                System.out.println("SagePredectRetentionTimeParameter " + valueBoolean + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
//               

//                System.out.println("optimizeMaxFragmentChargeParameter " + value + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
               int  intvalue = optimizeMaxVariableModificationParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SageMaxVariableModificationParameter"));
                if (intvalue != sageParameters.getMaxVariableMods()) {
                    sageParameters.setMaxVariableMods(intvalue);
//                    System.out.println("optimizeMaxVariableModificationParameter " + value + "------------------------------------------------------------------------->>> 12 id rate " + optProtDataset.getActiveIdentificationNum());
//
                }
                intvalue = optimizeMinMatchedPeaksParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SageMinMatchedPeaksPeakParameter"));
                if (intvalue != sageParameters.getMinMatchedPeaks()) {
                    sageParameters.setMinMatchedPeaks(intvalue);
                }
//                      
                intvalue = optimizeMaxFragmentChargeParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SageMaxFragmentChargeParameter"));
                sageParameters.setMaxFragmentCharge(intvalue);

                //need adjustment
                sageParameters.setGenerateDecoys(true);
                valueBoolean = optimizeGenerateDecoyParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SageGenerateDecoyParameter"));
                if (valueBoolean != sageParameters.getGenerateDecoys()) {
                    sageParameters.setGenerateDecoys(valueBoolean);
                }
                valueBoolean = optimizeWideWindowParameter(optProtDataset, identificationParameters, searchInputSetting, parameterScoreMap.get("SageWideWindowParameter"));
                if (valueBoolean != sageParameters.getWideWindow()) {
                    sageParameters.setWideWindow(valueBoolean);
                }

//                System.out.println("widewindow " + sageParameters.getWideWindow());
                IdentificationParameters.saveIdentificationParameters(identificationParameters, generatedIdentificationParametersFile);

            }
//

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

    @Override
    public synchronized RawScoreModel excuteSearch(SearchingSubDataset optProtDataset, String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addSpectraList, SearchInputSetting optProtSearchSettings, File identificationParametersFile, boolean pairData) {
        if (!optProtSearchSettings.getSageEnabledParameters().getParamsToOptimize().isEnabledParam(paramOption.split("_")[0])) {
            System.out.println("param " + paramOption + " is not supported " + paramOption);
            return new RawScoreModel(paramOption);
        }
//        if (defaultOutputFileName.contains("_resultsf_Carbamilation of protein N-term") || defaultOutputFileName.contains("_resultsf_Acetylation of protein N-term") || defaultOutputFileName.contains("_resultsf_Pyrolidone from carbamidomethylated C")) {
//            System.out.println("param " + paramOption + " is not supported " + paramOption);
//            return new RawScoreModel();
//        }

        Future<File> f = MainUtilities.getLongExecutorService().submit(() -> {
            File resultOutput = SearchExecuter.executeSearch(defaultOutputFileName, optProtSearchSettings, optProtDataset.getSubMsFile(), optProtDataset.getSubFastaFile(), tempIdParam, identificationParametersFile);
            return resultOutput;
        });
        File resultOutput = null;
        try {
            while (!f.isDone()) {
            }
            resultOutput = f.get();
        } catch (InterruptedException | ExecutionException ex) {
            ex.printStackTrace();
        }

        List<SpectrumMatch> validatedMaches = SpectraUtilities.getValidatedIdentificationResults(resultOutput, optProtDataset.getSubMsFile(), Advocate.sage, tempIdParam);
//        System.out.println("at param name is " + paramOption + "  " + validatedMaches.size());
        if (paramOption.contains("_")) {
            paramOption = paramOption.split("_")[1];
        }
        RawScoreModel rawScore = SpectraUtilities.getComparableRawScore(optProtDataset, validatedMaches, Advocate.sage, addSpectraList, paramOption);//(optProtDataset, resultOutput, optProtDataset.getSubMsFile(), Advocate.sage, tempIdParam, updateDataReference);
        if (addSpectraList && rawScore.isSensitiveChange()) {
            rawScore.setSpectrumMatchResult(validatedMaches);
        }
        return (rawScore);
    }

    public int[] optimizePeptideLengthParameter(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("peptideLength");
        SageParameters sageParameter = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());

        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxPeptideLengthOption = sageParameter.getMaxPeptideLength();
        int selectedMinPeptideLengthOption = sageParameter.getMinPeptideLength();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());

        sageParameter.setMaxPeptideLength(35);
        int selectedV = selectedMinPeptideLengthOption;
        for (int i = 5; i <= 10; i++) {
            if (i == selectedMinPeptideLengthOption) {
                continue;
            }
            sageParameter.setMinPeptideLength(i);
            final String option = "minPeptideLength_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange() && scoreModel.getSpectrumMatchResult().size() >= optProtDataset.getCurrentScoreModel().getSpectrumMatchResult().size()) {
//                    System.out.println("add as selected score " + j);
                    resultsMap.put(j + "", scoreModel);
                    selectedV = j;
                } else {
                    break;
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        if (!resultsMap.isEmpty()) {
            String bestOption = selectedV + "";// SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber());//    
//            int[] scoreLengthValues = new int[resultsMap.size()];
//            String[] keys = new String[resultsMap.size()];
//            int i = 0;
//            for (String key : resultsMap.keySet()) {
//                scoreLengthValues[i] = resultsMap.get(key).getSpectrumMatchResult().size();
//                keys[i] = key;
//                System.out.println("at indix " + i + "  " + key + "  " + resultsMap.get(key).getFinalScore());
//                i++;
//            }
//            selectedV = Integer.parseInt(keys[SpectraUtilities.findMainDrop(scoreLengthValues)]);
            selectedMinPeptideLengthOption = selectedV;
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));

        }

        resultsMap.clear();
        sageParameter.setMinPeptideLength(selectedMinPeptideLengthOption);
        for (int i = 25; i < 35; i++) {
            if (i == selectedMaxPeptideLengthOption) {
                continue;
            }
            sageParameter.setMaxPeptideLength(i);
            final String option = "maxPeptideLength_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.getFinalScore() > 0) {
                    resultsMap.put(j + "", scoreModel);
                } else if (!scoreModel.isSameData()) {
                    break;
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        if (!resultsMap.isEmpty()) {
            selectedMaxPeptideLengthOption = Integer.parseInt(SpectraUtilities.compareScoresSet(resultsMap));
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedMaxPeptideLengthOption + ""));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedMinPeptideLengthOption + "-" + selectedMaxPeptideLengthOption + "");
        parameterScoreSet.add(paramScore);

        return new int[]{selectedMinPeptideLengthOption, selectedMaxPeptideLengthOption};
    }

    public double[] optimizePeptideMassParameter(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        SageParameters sageParameter = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("PeptideMass");
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedMaxPeptideMassOption = sageParameter.getMaxPeptideMass();
        double selectedMinPeptideMassOption = sageParameter.getMinPeptideMass();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        sageParameter.setMaxPeptideMass(selectedMaxPeptideMassOption);
        double lastScore = -100000.0;

        for (double i = 400; i <= 600;) {
            if (i == selectedMinPeptideMassOption) {
                i += 50.0;
            }
            sageParameter.setMinPeptideMass(i);
            final String option = "minPeptideMass_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            double j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange() && lastScore < scoreModel.getFinalScore()) {
                    lastScore = scoreModel.getFinalScore();
                    resultsMap.put(j + "", scoreModel);
                } else {
                    break;
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            i += 50.0;
        }

        if (!resultsMap.isEmpty()) {
            String bestOption = SpectraUtilities.compareScoresSet(resultsMap);//  
            selectedMinPeptideMassOption = Double.parseDouble(bestOption);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));
        }
        sageParameter.setMinPeptideMass(selectedMinPeptideMassOption);
        resultsMap.clear();
        for (double i = 4000; i <= 6000;) {
            if (i == selectedMaxPeptideMassOption) {
                i += 500;
            }
            sageParameter.setMaxPeptideMass(i);
            final String option = "maxPeptideMass_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            double j = i;
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
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

            i += 500;
        }

        if (!resultsMap.isEmpty()) {
            String bestOption = SpectraUtilities.compareScoresSet(resultsMap);//  
            selectedMaxPeptideMassOption = Double.parseDouble(bestOption);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedMinPeptideMassOption + " - " + selectedMaxPeptideMassOption + "");
        parameterScoreSet.add(paramScore);

        return new double[]{selectedMinPeptideMassOption, selectedMaxPeptideMassOption};
    }

    public int optimizeMaxVariableModificationParameter(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("maxVarPTMs");
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        int selectedOption = sageParameters.getMaxVariableMods();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());

        for (int i = 1; i <= 3; i++) {
            if (i == selectedOption) {
                continue;
            }
            final String option = "maxVarPTMs_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setMaxVariableMods(i);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange() || (scoreModel.isSensitiveChange() && j < selectedOption)) {
                    resultsMap.put(j + "", scoreModel);
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        if (!resultsMap.isEmpty()) {
            String bestOption = SpectraUtilities.compareScoresSet(resultsMap);//  
            selectedOption = Integer.parseInt(bestOption);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));
        }

        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);

        return selectedOption;
    }

    public double[] optimizeFragmentMzParameter(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("FragmentMz");
        SageParameters sageParameter = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedMaxFragmentMzOption = sageParameter.getMaxFragmentMz();
        double selectedMinFragmentMzOption = sageParameter.getMinFragmentMz();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        sageParameter.setMaxFragmentMz(selectedMaxFragmentMzOption);
        for (double i = 150; i <= 300;) {
            if (i == selectedMinFragmentMzOption) {
                i += 25.0;
                continue;
            }
            sageParameter.setMinFragmentMz(i);
            final String option = "minFragmentMz_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            double j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
//                System.out.println("min fz " + j + "  " + scoreModel);
                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j + "", scoreModel);
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            i += 25.0;
        }

        if (!resultsMap.isEmpty()) {

            String bestOption = SpectraUtilities.compareScoresSet(resultsMap);//  
            selectedMinFragmentMzOption = Double.parseDouble(bestOption);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));
        }
        sageParameter.setMinFragmentMz(selectedMinFragmentMzOption);
        resultsMap.clear();
        for (double i = 1500; i <= 3000;) {
            if (i == selectedMaxFragmentMzOption) {
                i += 250;
                continue;
            }
            sageParameter.setMaxFragmentMz(i);
            final String option = "maxFragmentMz_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            double j = i;
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();

                if (j < selectedMaxFragmentMzOption && (scoreModel.getFinalScore() > 0)) {
                    resultsMap.put(j + "", scoreModel);
                } else if (j > selectedMaxFragmentMzOption && !scoreModel.isSensitiveChange()) {
                    break;
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

            i += 250;
        }
        if (!resultsMap.isEmpty()) {
            String bestOption = SpectraUtilities.compareScoresSet(resultsMap);//  
            selectedMaxFragmentMzOption = Double.parseDouble(bestOption);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedMinFragmentMzOption + " - " + selectedMaxFragmentMzOption + "");
        parameterScoreSet.add(paramScore);
        return new double[]{selectedMinFragmentMzOption, selectedMaxFragmentMzOption};
    }

    public int optimizeIonMinIndexParameter(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("ionMinIndex");
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        int selectedOption = sageParameters.getMinIonIndex();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        for (int i = 0; i <= 6; i++) {
            if (i == selectedOption) {
                continue;
            }
            final String option = "ionMinIndex_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setMinIonIndex(i);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                System.out.println("ionindex : "+i+"  "+scoreModel);
                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j + "", scoreModel);
                } else {
                    break;
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        if (!resultsMap.isEmpty()) {
            String bestOption = SpectraUtilities.compareScoresSet(resultsMap);//  
            selectedOption = Integer.parseInt(bestOption);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }

    public boolean optimizeGenerateDecoyParameter(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("generateDecoy");
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = sageParameters.getGenerateDecoys();

        for (int i = 0; i < 2; i++) {
            boolean generateDecoy = (i == 1);
            if (generateDecoy == selectedOption) {
                continue;
            }

            final String option = "generateDecoy_" + generateDecoy;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setGenerateDecoys(generateDecoy);
            final int j = i;
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
//                System.out.println("generate decoy result " + generateDecoy + "  " + scoreModel);
                if (scoreModel.getFinalScore() > 0) {
                    resultsMap.put(j + "", scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        if (!resultsMap.isEmpty()) {
            int bestOption = Integer.parseInt(SpectraUtilities.compareScoresSet(resultsMap));
            selectedOption = bestOption == 1;
            double impact = Math.round((double) (resultsMap.get(bestOption + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption + ""));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;

    }

    public boolean optimizeDeisotopParameter(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("Deisotope");
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = sageParameters.getDeisotope();
        for (int i = 0; i < 2; i++) {
            if (selectedOption == (i == 1)) {
                continue;
            }
            final String option = "Deisotope_" + (i == 1);
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setDeisotope(i == 1);
            final int j = i;
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                System.out.println("deistop score "+scoreModel+"   "+selectedOption+"   "+(i == 1)+"   "+(selectedOption == (i == 1)));
                if (scoreModel.isSensitiveChange()) {
                    resultsMap.put(j + "", scoreModel);

                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        if (!resultsMap.isEmpty()) {
            String bestOption = SpectraUtilities.compareScoresSet(resultsMap);//  
            selectedOption = Integer.parseInt(bestOption) == 1;
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;

    }

    public boolean optimizeChimericSpectraParameter(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("Chimera");
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = sageParameters.getChimera();
        for (int i = 0; i < 2; i++) {
            if (selectedOption == (i == 1)) {
                continue;
            }
            final String option = "Chimera_" + (i == 1);
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setChimera(i == 1);
            final int j = i;
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j + "", scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        if (!resultsMap.isEmpty()) {
            String bestOption = SpectraUtilities.compareScoresSet(resultsMap);//  
            selectedOption = Integer.parseInt(bestOption) == 1;
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;

    }

    public boolean optimizeWideWindowParameter(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {

        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("WideWindow");
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = sageParameters.getWideWindow();

        for (int i = 0; i < 2; i++) {
            if (selectedOption == (i == 1)) {
                continue;
            }
            final String option = "WideWindow_" + (i == 1);
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setWideWindow(i == 1);
            final int j = i;
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
//                System.out.println("score model " + (option) + "  " + scoreModel);
                if (scoreModel.getFinalScore() > 1.5+optProtDataset.getComparisonsThreshold()) {
                    resultsMap.put(j + "", scoreModel);

                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        if (!resultsMap.isEmpty()) {
            String bestOption = SpectraUtilities.compareScoresSet(resultsMap);//  
            selectedOption = Integer.parseInt(bestOption) == 1;
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));
            paramScore.setComments("Slow processing");
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;

    }

    public boolean optimizePredectRetentionTimeParameter(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("PredictRt");
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        boolean selectedOption = sageParameters.getPredictRt();

        for (int i = 0; i < 2; i++) {
            if (selectedOption == (i == 1)) {
                continue;
            }
            final String option = "PredictRt_" + (i == 1);
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setPredictRt(i == 1);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(j + "", scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        if (!resultsMap.isEmpty()) {
            String bestOption = SpectraUtilities.compareScoresSet(resultsMap);//  
            selectedOption = Integer.parseInt(bestOption) == 1;
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;

    }

    public int[] optimizeNumberOfPeakParameter(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("#Peaks");
        SageParameters sageParameter = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxPeaksNumberOption = sageParameter.getMaxPeaks();
        int selectedMinPeaksNumberOption = sageParameter.getMinPeaks();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        sageParameter.setMaxPeaks(selectedMaxPeaksNumberOption);
        RawScoreModel prescoreModel = null;
        for (int i = 10; i <= 20; i++) {
            if (i == selectedMinPeaksNumberOption) {
                continue;
            }
            sageParameter.setMinPeaks(i);
            final String option = "minPeaks_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange() && scoreModel.getIdPSMNumber() >= optProtDataset.getCurrentScoreModel().getIdPSMNumber()) { //|| scoreModel.isSameData()
                    if (prescoreModel == null) {
                        prescoreModel = scoreModel;
                    } else if (prescoreModel.getFinalScore() < scoreModel.getFinalScore()) {
                        prescoreModel = scoreModel;
                    } else {
                        continue;
                    }
                    resultsMap.put(j + "", scoreModel);
                } else if (j > selectedMinPeaksNumberOption || scoreModel.getIdPSMNumber() < optProtDataset.getCurrentScoreModel().getIdPSMNumber()) {
                    break;
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        if (!resultsMap.isEmpty()) {
            String bestOption = SpectraUtilities.compareScoresSet(resultsMap);//            
            selectedMinPeaksNumberOption = Integer.parseInt(bestOption);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));

        }
        resultsMap.clear();
        sageParameter.setMinPeaks(selectedMinPeaksNumberOption);
        for (int i = 100; i <= 200;) {
            if (i == selectedMaxPeaksNumberOption) {
                i += 10;
                continue;
            }
            sageParameter.setMaxPeaks(i);
            final String option = "maxPeaks_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if ((scoreModel.getFinalScore() > 0)) {
                    resultsMap.put(j + "", scoreModel);
                } else if (j > selectedMaxPeaksNumberOption && !scoreModel.isSensitiveChange()) {
                    break;
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
            i += 10;
        }

        if (!resultsMap.isEmpty()) {

            selectedMaxPeaksNumberOption = Integer.parseInt(SpectraUtilities.compareScoresSet(resultsMap));
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedMaxPeaksNumberOption + ""));
            optProtDataset.setActiveScoreModel(resultsMap.get(selectedMaxPeaksNumberOption + ""));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedMinPeaksNumberOption + "-" + selectedMaxPeaksNumberOption + "");
        parameterScoreSet.add(paramScore);

        return new int[]{selectedMinPeaksNumberOption, selectedMaxPeaksNumberOption};
    }

    public int optimizeMinMatchedPeaksParameter(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("minMatchedPeaks");
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        int selectedOption = sageParameters.getMinMatchedPeaks();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        for (int i = 3; i <= 6; i++) {
            if (i == selectedOption) {
                continue;
            }
            final String option = "minMatchedPeaks_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setMinMatchedPeaks(i);
            final int j = i;
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange() && scoreModel.getIdPSMNumber() >= optProtDataset.getActiveIdentificationNum()) {
                    resultsMap.put(j + "", scoreModel);
                } else if (j > selectedOption) {
                    break;
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        if (!resultsMap.isEmpty()) {
            String bestOption = SpectraUtilities.compareScoresSet(resultsMap);//  
            selectedOption = Integer.parseInt(bestOption);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }

    public int optimizeMaxFragmentChargeParameter(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("maxFragmentCharge");
        SageParameters sageParameters = (SageParameters) oreginaltempIdParam.getSearchParameters().getAlgorithmSpecificParameters().get(Advocate.sage.getIndex());
        int selectedOption = 1;
        if (sageParameters.getMaxFragmentCharge() != null) {
            selectedOption = sageParameters.getMaxFragmentCharge();
        }
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        for (int i = 1; i <= 5; i++) {
            if (i == selectedOption) {
                continue;
            }
            final String option = "maxFragmentCharge_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            sageParameters.setMaxFragmentCharge(i);
            final int j = i;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, false, optimisedSearchParameter, generatedIdentificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.getFinalScore() > 0) {
                    resultsMap.put(j + "", scoreModel);
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        if (!resultsMap.isEmpty()) {
            String bestOption = SpectraUtilities.compareScoresSet(resultsMap);//  
            selectedOption = Integer.parseInt(bestOption);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestOption));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }
    
    public void runReferenceRun(SearchingSubDataset optProtDataset, IdentificationParameters oreginaltempIdParam, SearchInputSetting optimisedSearchParameter) throws IOException {

        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        final String option = "reference_run_default_";
        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

        Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, generatedIdentificationParametersFile, false);
            return scoreModel;
        });
        try {
            RawScoreModel scoreModel = f.get();
            System.out.println("reference run " + scoreModel);
            optProtDataset.setActiveScoreModel(scoreModel);
        } catch (ExecutionException | InterruptedException ex) {
            ex.printStackTrace();
        }
    }


}
