package no.uib.probe.optprot.search;

import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.biology.modifications.Modification;
import com.compomics.util.experiment.biology.modifications.ModificationCategory;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.io.IoUtil;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
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
import no.uib.probe.optprot.model.ParameterScoreModel;
import no.uib.probe.optprot.model.RawScoreModel;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.model.SortedPTMs;
import no.uib.probe.optprot.util.MainUtilities;
import no.uib.probe.optprot.util.SpectraUtilities;

/**
 *
 * @author yfa041
 */
public abstract class DefaultOptProtSearchOptimizer {

    /**
     * The compomics PTM factory.
     */
    private final ModificationFactory ptmFactory = ModificationFactory.getInstance();

    public String optimizeDigestionCleavageParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("CleavageParameter");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name();
        int idRate = optProtDataset.getActiveIdentificationNum();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        RawScoreModel oreginalScore = new RawScoreModel();
        oreginalScore.setTotalNumber(idRate);
        String[] cleavageParameters = new String[]{"wholeProtein", "unSpecific"};
        resultsMap.put(selectedOption, oreginalScore);
        int spectraCounter = optProtDataset.getActiveIdentificationNum();
        for (String cleavageParameter : cleavageParameters) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf(cleavageParameter));
            if (cleavageParameter.equalsIgnoreCase("enzyme")) {
                tempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(EnzymeFactory.getInstance().getEnzyme("Trypsin"));
            } else {
                tempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
            }
            final String option = cleavageParameter;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, optimisedSearchParameter, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.getFinalScore() > 1) {
                    spectraCounter = Math.max(spectraCounter, scoreModel.getSpectrumMatchResult().size());
                    resultsMap.put(option, f.get());

                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        int total = optProtDataset.getActiveIdentificationNum();
        if (!resultsMap.isEmpty()) {
            String bestScore = SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber());
            selectedOption = bestScore;
            if (selectedOption.equalsIgnoreCase("unSpecific")) {
                double impact = Math.round((double) (resultsMap.get(selectedOption).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                total = resultsMap.get(selectedOption).getSpectrumMatchResult().size();
                paramScore.setComments("Extremely slow processing");
            }
        }
        paramScore.setScore(total);
        paramScore.setParamValue(selectedOption);
        parameterScoreSet.add(paramScore);
        return selectedOption;

    }

    /**
     * Optimize digestion enzyme
     *
     * @param optProtDataset
     * @param identificationParametersFile
     * @param optimisedSearchParameter
     * @param parameterScoreSet
     * @return
     * @throws IOException
     */
    public String[] optimizeEnzymeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("Enzyme");
        String[] values = new String[3];
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        boolean optimizeOnlyMissedCleaveNum = false;
        boolean compareBetweenEnzymes = false;
        int missedClavageNumb = 0;
        if (optimisedSearchParameter.getDigestionParameterOpt().equalsIgnoreCase("enzyme") && (optProtDataset.getIdentificationRate() > 40)) {//            
            values[0] = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
            values[1] = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getSpecificity(values[0]).name();
            missedClavageNumb = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(values[0]);
            optimizeOnlyMissedCleaveNum = true;
        } else {
            oreginaltempIdParam.getSearchParameters().getDigestionParameters().setCleavageParameter(DigestionParameters.CleavageParameter.valueOf("enzyme"));
            compareBetweenEnzymes = true;
            values[1] = "specific";
        }
        Map<String, RawScoreModel> resultsMapI = Collections.synchronizedMap(new LinkedHashMap<>());
        Map<String, RawScoreModel> resultsMapII = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        //optimise enzyme  
        if (!optimizeOnlyMissedCleaveNum) {
            for (Enzyme enzyme : EnzymeFactory.getInstance().getEnzymes()) {
                if (enzyme.getName().replace(" ", "").equalsIgnoreCase("Trypsin(noPrule)")) {
                    continue;
                }
                oreginaltempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
                oreginaltempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(enzyme);
                oreginaltempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(enzyme.getName(), missedClavageNumb);
                final String option = enzyme.getName();
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, identificationParametersFile, false);
                    return scoreModel;
                });
                try {
                    RawScoreModel scoreModel = f.get();
                    if (scoreModel.isSignificatChange()) {
                        resultsMapI.put(option, scoreModel);
                    }
                    if (compareBetweenEnzymes) {
                        resultsMapII.put(option, scoreModel);
                    }
                } catch (InterruptedException | ExecutionException ex) {
                    ex.printStackTrace();
                }
            }
            if (compareBetweenEnzymes && !resultsMapII.isEmpty()) {
                String enzymeName = SpectraUtilities.compareScoresSet(resultsMapII, optProtDataset.getTotalSpectraNumber());
                values[0] = enzymeName;
            } else if (!resultsMapI.isEmpty()) {
                String enzymeName = SpectraUtilities.compareScoresSet(resultsMapI, optProtDataset.getTotalSpectraNumber());
                values[0] = enzymeName;
                optProtDataset.setActiveScoreModel(resultsMapI.get(enzymeName));
            } else {
                values[0] = "Trypsin";
            }

            //optimize specifty 
            oreginaltempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
            oreginaltempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(EnzymeFactory.getInstance().getEnzyme(values[0]));
            oreginaltempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(values[0], missedClavageNumb);
            resultsMapII.clear();
            resultsMapI.clear();
            for (int i = 0; i < DigestionParameters.Specificity.values().length; i++) {
                final String option = DigestionParameters.Specificity.getSpecificity(i).name();
                oreginaltempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(values[0], DigestionParameters.Specificity.getSpecificity(i));
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, identificationParametersFile, false);
                    return scoreModel;
                });
                try {
                    RawScoreModel scoreModel = f.get();
                    if (scoreModel.isSignificatChange()) {
                        resultsMapI.put(option, scoreModel);
                    }
                    if ((compareBetweenEnzymes && scoreModel.getFinalScore() > 0)) {
                        resultsMapII.put(option, scoreModel);
                    }
                } catch (ExecutionException | InterruptedException ex) {
                    ex.printStackTrace();
                }
            }
            if (compareBetweenEnzymes && !resultsMapII.isEmpty()) {
                String specifty = SpectraUtilities.compareScoresSet(resultsMapII, optProtDataset.getTotalSpectraNumber());
                values[1] = specifty;
            } else if (!resultsMapI.isEmpty()) {
                String specifty = SpectraUtilities.compareScoresSet(resultsMapI, optProtDataset.getTotalSpectraNumber());
                values[1] = specifty;
                double impact = Math.round((double) (resultsMapI.get(specifty).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                paramScore.setImpact(impact);
                optProtDataset.setActiveScoreModel(resultsMapI.get(specifty));

            } else {
                values[1] = "specific";
            }

        }
        oreginaltempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(values[0], DigestionParameters.Specificity.valueOf("specific"));
///number op missed cleavage
        resultsMapI.clear();
        for (int i = 0; i < 5; i++) {
            oreginaltempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(values[0], i);
            final String option = "missedCleavages_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                System.out.println("missed clavage numb " + i + "  " + scoreModel + "  " + optProtDataset.getActiveIdentificationNum() + "  " + scoreModel.getDataLengthFactor());

                if (scoreModel.isSignificatChange()) {
                    resultsMapI.put(i + "", scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        String numbOfMissedCleavage = "2";
        if (!resultsMapI.isEmpty()) {
            numbOfMissedCleavage = (SpectraUtilities.compareScoresSet(resultsMapI, optProtDataset.getTotalSpectraNumber()));
            double impact = Math.round((double) (resultsMapI.get(numbOfMissedCleavage).getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMapI.get(numbOfMissedCleavage));
        }
        values[2] = numbOfMissedCleavage;
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(Arrays.asList(values).toString());
        parameterScoreSet.add(paramScore);
        return values;

    }

    /**
     *
     * @param optProtDataset
     * @param identificationParametersFile
     * @param optimisedSearchParameter
     * @param parameterScoreSet
     * @return
     * @throws IOException
     */
//    public String optimizeSpecificityParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
//        final ParameterScoreModel paramScore = new ParameterScoreModel();
//        paramScore.setParamId("Specificity");
//        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//        if (!oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name().equalsIgnoreCase("enzyme")) {
//            return null;
//        }
//        String enzymeName = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
//        String selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getSpecificity(enzymeName).name();
//        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
//        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
//        for (int i = 0; i < DigestionParameters.Specificity.values().length; i++) {
//            final String option = DigestionParameters.Specificity.getSpecificity(i).name();
//            if (option.equalsIgnoreCase(selectedOption)) {
//                continue;
//            }
//            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
//            tempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(enzymeName, DigestionParameters.Specificity.getSpecificity(i));
//
//            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
//
//            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
//                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
//                return scoreModel;
//            });
//            try {
//                RawScoreModel scoreModel = f.get();
//                if (scoreModel.isSignificatChange()) {
//                    resultsMap.put(option, scoreModel);
//                }
//            } catch (ExecutionException | InterruptedException ex) {
//                ex.printStackTrace();
//            }
//        }
//
//        TreeMap<RawScoreModel, String> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
//        if (!resultsMap.isEmpty()) {
//            for (String option : resultsMap.keySet()) {
//                sortedResultsMap.put(resultsMap.get(option), option);
//            }
//            selectedOption = sortedResultsMap.firstEntry().getValue();
//            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
//        }
//        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
//        paramScore.setParamValue(selectedOption);
//        parameterScoreSet.add(paramScore);
//
//        return selectedOption;
//
//    }
    public String optimizeFragmentIonTypesParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("FragmentIons");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        String selectedOption = oreginaltempIdParam.getSearchParameters().getForwardIons() + "-" + oreginaltempIdParam.getSearchParameters().getRewindIons();
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SearchParameters tempSearchParameters = tempIdParam.getSearchParameters();
        ArrayList<Integer> selectedForwardIons = tempSearchParameters.getForwardIons();
        String[] forwardIons = new String[]{"b", "a", "c"};
        String[] rewindIons = new String[]{"y", "x", "z"};

        for (String forwardIon : forwardIons) {
            selectedForwardIons.clear();
            Integer forwardIonType = PeptideFragmentIon.getIonType(forwardIon);
            selectedForwardIons.add(forwardIonType);
            for (String rewindIon : rewindIons) {
                Integer rewindIonType = PeptideFragmentIon.getIonType(rewindIon);
                ArrayList<Integer> selectedRewindIons = new ArrayList<>();
                selectedRewindIons.add(rewindIonType);
                tempSearchParameters.setRewindIons(selectedRewindIons);
                String option = selectedForwardIons + "-" + selectedRewindIons;
                if (option.equalsIgnoreCase(selectedOption)) {
                    continue;
                }

                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
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
        }

        if (!resultsMap.isEmpty()) {
            String bestScore = SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber());
            selectedOption = bestScore;
            double impact = Math.round((double) (resultsMap.get(selectedOption + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestScore));
        }

        selectedOption = selectedOption.replace("[", "").replace("]", "");
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption);
        parameterScoreSet.add(paramScore);
        return selectedOption;
    }

    public Integer optimizeMaxMissCleavagesParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("missedCleavages");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        if (!oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name().equalsIgnoreCase("enzyme")) {
            return -1;
        }
        String enzymeName = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
        Integer selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(enzymeName);
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());

        for (int i = 0; i < 5; i++) {
            if (i == selectedOption) {
                continue;
            }
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(enzymeName, i);
            final String option = "missedCleavages_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(i + "", scoreModel);
                } else if (i > selectedOption && !scoreModel.isSensitiveChange()) {
                    break;
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }

        if (!resultsMap.isEmpty()) {
            String bestScore = SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber());
            selectedOption = Integer.valueOf(bestScore);
            double impact = Math.round((double) (resultsMap.get(selectedOption + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestScore));
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption + "");
        parameterScoreSet.add(paramScore);

        return selectedOption;

    }

    public double optimizeFragmentToleranceParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("fragmentAccuracy");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedOption = oreginaltempIdParam.getSearchParameters().getFragmentIonAccuracy();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double[] values = new double[]{0.01, 0.02, 0.05, 0.1, 0.2, 0.5};
        for (double i : values) {
            if (selectedOption == i) {
                continue;
            }
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setFragmentIonAccuracy(i);
            tempIdParam.getSearchParameters().setFragmentAccuracyType(SearchParameters.MassAccuracyType.DA);
            final String option = "fragmentAccuracy_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                System.out.println("fragment accurcy " + i + "  " + scoreModel + "  " + scoreModel.getTotalNumber());
                if (scoreModel.getFinalScore() > 1) {
                    resultsMap.put(i + "", scoreModel);
                } else if (i > selectedOption) {
                    break;
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

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

    public int[] optimizePrecursorChargeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("charge");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxChargeOption = oreginaltempIdParam.getSearchParameters().getMaxChargeSearched();
        int selectedMinChargeOption = oreginaltempIdParam.getSearchParameters().getMinChargeSearched();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        int spectraCounter = optProtDataset.getActiveIdentificationNum();
        for (int i = 1; i < 5; i++) {

            for (int j = 2; j <= 5; j++) {
                if (j <= i) {
                    continue;
                }
                tempIdParam.getSearchParameters().setMinChargeSearched(i);
                tempIdParam.getSearchParameters().setMaxChargeSearched(j);

                final String option = "charge-" + i + "," + j;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {

                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                    return scoreModel;
                });
                try {

                    RawScoreModel scoreModel = f.get();
                    if (scoreModel.isSignificatChange()) {
                        if (scoreModel.getSpectrumMatchResult().size() < spectraCounter) {
                            continue;
                        }
                        spectraCounter = Math.max(spectraCounter, scoreModel.getSpectrumMatchResult().size());
                        resultsMap.put(option, scoreModel);

                    }
                } catch (ExecutionException | InterruptedException ex) {
                    ex.printStackTrace();
                }
            }

        }

        if (!resultsMap.isEmpty()) {
            String bestScore = SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber());
            double impact = Math.round((double) (resultsMap.get(bestScore + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestScore));
            String[] topOption = bestScore.split("-")[1].split(",");
            selectedMinChargeOption = Integer.parseInt(topOption[0]);
            selectedMaxChargeOption = Integer.parseInt(topOption[1]);
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedMinChargeOption + "," + selectedMaxChargeOption);
        parameterScoreSet.add(paramScore);
        return new int[]{selectedMinChargeOption, selectedMaxChargeOption};

    }

    public int[] optimizeIsotopParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {

        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("isotop_");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxIsotopicCorrectionOption = oreginaltempIdParam.getSearchParameters().getMaxIsotopicCorrection();
        int selectedMinIsotopicCorrectioneOption = oreginaltempIdParam.getSearchParameters().getMinIsotopicCorrection();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);

        for (int i = -2; i < 2; i++) {
            for (int j = -1; j <= 2; j++) {
                if (j <= i) {
                    continue;
                }
                tempIdParam.getSearchParameters().setMinIsotopicCorrection(i);
                tempIdParam.getSearchParameters().setMaxIsotopicCorrection(j);
                final String option = "isotop_" + i + "," + j;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
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
        }

        if (!resultsMap.isEmpty()) {
            String bestScore = SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber());
            double impact = Math.round((double) (resultsMap.get(bestScore + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestScore));
            String[] topOption = bestScore.split("_")[1].split(",");
            selectedMinIsotopicCorrectioneOption = Integer.parseInt(topOption[0]);
            selectedMaxIsotopicCorrectionOption = Integer.parseInt(topOption[1]);
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedMinIsotopicCorrectioneOption + "," + selectedMaxIsotopicCorrectionOption);
        parameterScoreSet.add(paramScore);

        return new int[]{selectedMinIsotopicCorrectioneOption, selectedMaxIsotopicCorrectionOption};
    }

    public double optimizePrecursorToleranceParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("PrecursorAccuracy");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedOption = oreginaltempIdParam.getSearchParameters().getPrecursorAccuracy();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double[] iValues = new double[]{5, 10, 15, 20, 25};
        boolean toEnd = false;
        for (double i : iValues) {
            if (i == selectedOption) {
                continue;
            }
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setPrecursorAccuracy(i);
            tempIdParam.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
            final String option = "precursorAccuracy_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                return scoreModel;
            });
            try {

                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(i + "", scoreModel);
                } else if (i > selectedOption) {
                    toEnd = true;
                    break;

                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }

        if (!toEnd) {
            iValues = new double[]{0.1, 0.3, 0.5, 0.7, 0.9};
            if (!optProtDataset.isHighResolutionMassSpectrometers()) {
                for (double i : iValues) {
                    IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
                    tempIdParam.getSearchParameters().setPrecursorAccuracy(i);
                    tempIdParam.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.DA);
                    final String option = "precursorAccuracy_Da" + i;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
                    Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                        RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                        return scoreModel;
                    });
                    try {
                        RawScoreModel scoreModel = f.get();
                        if (scoreModel.isSignificatChange()) {
                            resultsMap.put(i + "", scoreModel);
                        }
                    } catch (ExecutionException | InterruptedException ex) {
                        ex.printStackTrace();
                    }

                }
            }
        }
        if (!resultsMap.isEmpty()) {
            String bestScore = SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber());
            double impact = Math.round((double) (resultsMap.get(bestScore + "").getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
            paramScore.setImpact(impact);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestScore));
            selectedOption = Double.parseDouble(bestScore);
        }

        paramScore.setScore(optProtDataset.getActiveIdentificationNum());

        parameterScoreSet.add(paramScore);
        if (selectedOption >= 5) {
            paramScore.setParamValue(selectedOption + "PPM");
            paramScore.setComments("High-Resolution Mass Spectrometers: Instruments like Orbitrap or Fourier Transform Ion Cyclotron Resonance (FT-ICR)");
        } else {
            paramScore.setParamValue(selectedOption + "Da");
            paramScore.setComments("Low-Resolution Mass Spectrometers: Quadrupole and ion trap mass spectrometers have lower mass accuracy");
        }
        return selectedOption;
    }

    //add common fixed
    //then run fixed mod
    //if not M included
    // add common variable
    //then potinatial ptms
    /**
     * This function responsible for selection best combination of
     * (fixed,variable,refinement ) modifications-PTMs
     *
     * @param optProtDataset input data-set
     * @param identificationParametersFile identification file
     * @param searchInputSetting
     * @param parameterScoreSet set to store the parameter information object
     * @return modification map
     * @throws IOException
     */
    public Map<String, Set<String>> optimizeModificationsParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting searchInputSetting, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);

        final ParameterScoreModel fixedModParamScore = new ParameterScoreModel();
        fixedModParamScore.setParamId("FixedModifications");

        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        ArrayList<String> selectedFixedModificationOption = new ArrayList<>(oreginaltempIdParam.getSearchParameters().getModificationParameters().getFixedModifications());
        ArrayList<String> selectedVariableModificationOption = new ArrayList<>();
        Map<String, Set<String>> modificationsResults = new HashMap<>();
        List<String> mods = new ArrayList<>();
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common));
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common_Biological));
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common_Artifact));
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);

        tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
        tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        TreeMap<RawScoreModel, String> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        Set<String> potintialVariableMod = new LinkedHashSet<>(ptmFactory.getModifications(ModificationCategory.Common));
        String commonFixedMod = "Carbamidomethylation of C";
        String commonVariableMod = "Oxidation of M";
        boolean addedCommonPTMs = false;
        Map<String, RawScoreModel> targtedFixedModificationScore = new TreeMap<>();
        //first stage common fixed modification
//        for ( : commonMods) {
        try {
            final String option = commonFixedMod;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "cvf_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
            tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(commonFixedMod));
            tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(commonFixedMod));
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            RawScoreModel scoreModel = f.get();
            System.out.println("common fixed " + scoreModel);
            if (scoreModel.isSensitiveChange()) {
                addedCommonPTMs = true;
                targtedFixedModificationScore.put("C", scoreModel);

                optProtDataset.setActiveScoreModel(scoreModel);
            }
        } catch (ExecutionException | InterruptedException ex) {
            ex.printStackTrace();
        }
        //second stage fixed modifications
        mods.remove(commonFixedMod);

        for (String modId : mods) {
            String modPattern = ptmFactory.getModification(modId).getPattern().toString();
            if (modPattern.equalsIgnoreCase("C") && addedCommonPTMs) {
                continue;
            }
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
            tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(modId));
            tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(modId));
            if (addedCommonPTMs) {
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(commonFixedMod));
                tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(commonFixedMod));
            }
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                resultsMap.put(modId, f.get());
                if (resultsMap.get(modId).isSignificatChange()) {
                    sortedResultsMap.put(resultsMap.get(modId), modId);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }
        for (RawScoreModel score : sortedResultsMap.keySet()) {
            String modificationId = sortedResultsMap.get(score);
            String modPattern = ptmFactory.getModification(modificationId).getPattern().toString();
            if (ptmFactory.getModification(modificationId).getModificationType().isNTerm()) {
                if (!targtedFixedModificationScore.containsKey("NTERM")) {
                    targtedFixedModificationScore.put("NTERM", resultsMap.get(modificationId));
                } else if (targtedFixedModificationScore.containsKey("") || modPattern.equalsIgnoreCase("")) {
                    continue;
                }
            }

            if (!targtedFixedModificationScore.containsKey(modPattern)) {
                targtedFixedModificationScore.put(modPattern, resultsMap.get(modificationId));
            }
        }

        resultsMap.clear();
        int modLimit = 0;
        for (String key : targtedFixedModificationScore.keySet()) {
            if (key.equals("C") && addedCommonPTMs) {
                continue;
            }
            RawScoreModel scoreModel = targtedFixedModificationScore.get(key);
            resultsMap.put(sortedResultsMap.get(scoreModel), scoreModel);
            modLimit++;
            if (modLimit > 3) {
                break;
            }
        }
        sortedResultsMap.clear();
        Set<String> clearSet = new LinkedHashSet<>(resultsMap.keySet());
//       
        for (String modId : clearSet) {
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "fAsm_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            if (addedCommonPTMs) {
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(commonFixedMod));
                tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(commonFixedMod));

            }
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel vScore = f.get();
                if (vScore.isSensitiveChange()) {
                    if (SpectraUtilities.isBetterScore(resultsMap.get(modId).getSpectrumMatchResult(), vScore.getSpectrumMatchResult(), optProtDataset.getTotalSpectraNumber()) > 0) {
                        resultsMap.remove(modId);
                    }
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }

        Map<String, RawScoreModel> oneDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>(resultsMap));
        Map<String, RawScoreModel> twoDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
        tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();

        int indexI = 0;
        for (String mod1Id : oneDResultsMap.keySet()) {
            int indexII = 0;
            for (String mod2Id : oneDResultsMap.keySet()) {
                if (indexII <= indexI) {
                    indexII++;
                    continue;
                }
                tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
                tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod1Id));
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id));
                tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(mod1Id));
                tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(mod2Id));
                if (addedCommonPTMs) {
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(commonFixedMod));
                    tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(commonFixedMod));

                }

                final String option = mod1Id + "_" + mod2Id;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;
                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, searchInputSetting, identificationParametersFile, false);
                    return scoreModel;
                });
                try {
                    RawScoreModel score2D = f.get();
                    Map<String, RawScoreModel> sortinglist = new LinkedHashMap<>();
                    sortinglist.put(option, score2D);
                    sortinglist.put(mod1Id, oneDResultsMap.get(mod1Id));
                    sortinglist.put(mod2Id, oneDResultsMap.get(mod2Id));
                    String bestScoreFixed = SpectraUtilities.compareScoresSet(sortinglist, optProtDataset.getTotalSpectraNumber());
                    System.out.println("best score 1d bestScoreFixed " + bestScoreFixed + "  ");
                    if (bestScoreFixed.equalsIgnoreCase(option)) {
                        twoDResultsMap.put(option, score2D);
                    }
                } catch (ExecutionException | InterruptedException ex) {
                    ex.printStackTrace();
                }
                indexII++;
            }
            indexI++;
        }
        Map<String, RawScoreModel> threeDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        if (!twoDResultsMap.isEmpty()) {
            for (String mod1Id : oneDResultsMap.keySet()) {
                for (String mod2Id : twoDResultsMap.keySet()) {
                    if (mod2Id.contains(mod1Id)) {
                        continue;
                    }
                    tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
                    if (addedCommonPTMs) {
                        tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(commonFixedMod));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(commonFixedMod));

                    }

                    tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod1Id));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[0]));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[1]));

                    tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(mod1Id));
                    tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(mod2Id.split("_")[0]));
                    tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(mod2Id.split("_")[1]));

                    final String option = mod1Id + "_" + mod2Id;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;
                    Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                        RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, searchInputSetting, identificationParametersFile, false);
                        return scoreModel;
                    });
                    try {
                        RawScoreModel score3D = f.get();
                        Map<String, RawScoreModel> sortinglist = new LinkedHashMap<>();
                        sortinglist.put(option, score3D);
                        sortinglist.put(mod1Id, oneDResultsMap.get(mod1Id));
                        sortinglist.put(mod2Id, twoDResultsMap.get(mod2Id));
                        String bestScoreFixed = SpectraUtilities.compareScoresSet(sortinglist, optProtDataset.getTotalSpectraNumber());
                        System.out.println("best score 1d bestScoreFixed " + bestScoreFixed + "  ");
                        if (bestScoreFixed.equalsIgnoreCase(option)) {
                            threeDResultsMap.put(option, score3D);
                        }
                    } catch (ExecutionException | InterruptedException ex) {
                        ex.printStackTrace();
                    }
                }
            }

        }
        Map<String, RawScoreModel> fourDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        if (!threeDResultsMap.isEmpty()) {
            for (String mod1Id : oneDResultsMap.keySet()) {
                for (String mod2Id : threeDResultsMap.keySet()) {
                    if (mod2Id.contains(mod1Id)) {
                        continue;
                    }
                    tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
                    tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
                    if (addedCommonPTMs) {
                        tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(commonFixedMod));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(commonFixedMod));

                    }
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod1Id));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[0]));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[1]));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[2]));

                    tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(mod1Id));
                    tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(mod2Id.split("_")[0]));
                    tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(mod2Id.split("_")[1]));
                    tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(mod2Id.split("_")[2]));

                    final String option = mod1Id + "_" + mod2Id;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;

                    Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                        RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, searchInputSetting, identificationParametersFile, false);
                        return scoreModel;
                    });
                    try {
                        RawScoreModel score4D = f.get();

                        Map<String, RawScoreModel> sortinglist = new LinkedHashMap<>();
                        sortinglist.put(option, score4D);
                        sortinglist.put(mod1Id, oneDResultsMap.get(mod1Id));
                        sortinglist.put(mod2Id, threeDResultsMap.get(mod2Id));
                        String bestScoreFixed = SpectraUtilities.compareScoresSet(sortinglist, optProtDataset.getTotalSpectraNumber());
                        System.out.println("best score 1d bestScoreFixed " + bestScoreFixed + "  ");

                        if (bestScoreFixed.equalsIgnoreCase(option)) {
                            fourDResultsMap.put(option, score4D);
                        }
                    } catch (ExecutionException | InterruptedException ex) {
                        ex.printStackTrace();
                    }

                }
            }

        }
        resultsMap.clear();
        resultsMap.putAll(oneDResultsMap);
        resultsMap.putAll(twoDResultsMap);
        resultsMap.putAll(threeDResultsMap);
        resultsMap.putAll(fourDResultsMap);
        sortedResultsMap.clear();
        for (String key : resultsMap.keySet()) {
            sortedResultsMap.put(resultsMap.get(key), key);
        }
        selectedFixedModificationOption.clear();
        if (addedCommonPTMs) {
            selectedFixedModificationOption.add(commonFixedMod);
        }
        if (!sortedResultsMap.isEmpty()) {
            String fixedMod = sortedResultsMap.firstEntry().getValue();
            if (fixedMod.contains("_")) {
                selectedFixedModificationOption.addAll(Arrays.asList(fixedMod.split("_")));
            } else {
                selectedFixedModificationOption.add(fixedMod);
            }
        }

        tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
        tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
        modificationsResults.put("fixedModifications", new HashSet<>(selectedFixedModificationOption));
        MainUtilities.cleanOutputFolder();
        //process refine fixed modifications
        for (String fixedMod : selectedFixedModificationOption) {
            tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(fixedMod));
            tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(fixedMod));
        }
        modificationsResults.put("refinmentFixedModifications", new HashSet<>(selectedFixedModificationOption));
        MainUtilities.cleanOutputFolder();
        MainUtilities.cleanOutputFolder();
        fixedModParamScore.setScore(optProtDataset.getActiveIdentificationNum());
        fixedModParamScore.setParamValue(selectedFixedModificationOption.toString());
        parameterScoreSet.add(fixedModParamScore);

        //test for full fixed before add common variable
        // stage 3 test for common variable modification 
        try {
            final String option = commonVariableMod + "_" + selectedFixedModificationOption;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "finalFixed_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(commonVariableMod));
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            RawScoreModel score4D = f.get();
            if (score4D.isSensitiveChange()) {
                addedCommonPTMs = true;
                double impact = Math.round((double) (score4D.getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                fixedModParamScore.setImpact(impact);
                optProtDataset.setActiveScoreModel(score4D);
                selectedVariableModificationOption.add(commonVariableMod);
            } else {
                addedCommonPTMs = false;
            }
        } catch (ExecutionException | InterruptedException ex) {
            ex.printStackTrace();
        }

        //process variable modifications
        final ParameterScoreModel variableModParamScore = new ParameterScoreModel();
        variableModParamScore.setParamId("VariableModifications");
        resultsMap.clear();
        oneDResultsMap.clear();
        twoDResultsMap.clear();
        threeDResultsMap.clear();
        fourDResultsMap.clear();
//        sorterSet.clear();
        TreeMap<RawScoreModel, String> sortingModificationMap = new TreeMap<>(Collections.reverseOrder());//
        Map<String, Set<String>> modifiedSpectrumMap = new LinkedHashMap<>();
//        Map<String, Set<String>> fullSpectrumMap = new LinkedHashMap<>();
        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
        //check oxidation as variable mod
        potintialVariableMod.removeAll(selectedFixedModificationOption);
        potintialVariableMod.removeAll(selectedVariableModificationOption);
        mods.removeAll(selectedFixedModificationOption);
        mods.removeAll(selectedVariableModificationOption);

        //init potintial variable mod
        List<SortedPTMs> tree = new ArrayList<>();

        for (String modificationId : mods) {
            Modification mod = ptmFactory.getModification(modificationId);
            if (mod.getModificationType().isNTerm() || mod.getModificationType().isCTerm()) {
                potintialVariableMod.add(modificationId);
                continue;
            }
            String modPattern = mod.getPattern().toString();
            tree.add(new SortedPTMs(modificationId, SpectraUtilities.isPotintialVariableModification(modPattern, optProtDataset.getCurrentScoreModel().getSpectrumMatchResult())));

        }
        Collections.sort(tree);
        Collections.reverse(tree);
        int limit = optProtDataset.getCurrentScoreModel().getSpectrumMatchResult().size() * 10 / 100;
        int i = 0;
        for (SortedPTMs k : tree) {
            if (k.getScore() < (limit)) {
                break;
            }
            potintialVariableMod.add(k.getName());
            i++;

        }
        for (String modId : potintialVariableMod) {
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));
            if (addedCommonPTMs) {
                tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(commonVariableMod));
            }

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange()) {
                    System.out.println("variable mod " + modId + "  " + scoreModel + "  total: " + optProtDataset.getActiveIdentificationNum() + "  " + tempIdParam.getSearchParameters().getModificationParameters().getFixedModifications() + "   " + tempIdParam.getSearchParameters().getModificationParameters().getVariableModifications());
                    sortingModificationMap.put(scoreModel, modId);
                    List<SpectrumMatch> spectraResults = scoreModel.getSpectrumMatchResult();
                    modifiedSpectrumMap.put(modId, SpectraUtilities.getModifiedSpectrumSet(spectraResults));
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }

//        //filter 1 d map to max 4 modifications
        int modcounter = 0;
        if (!sortingModificationMap.isEmpty()) {
            for (RawScoreModel rsm : sortingModificationMap.keySet()) {
                String variableMod = sortingModificationMap.get(rsm);
                oneDResultsMap.put(variableMod, rsm);
                modcounter++;
//                if ((modcounter >= 4 && addedCommonPTMs)||modcounter >= 5 ) {
//                    break;
//                }
            }
        }
        //start process 2d variable mod
//        sorterSet.clear();
//        TreeSet<RawScoreModel> sorterSet = new TreeSet<>(Collections.reverseOrder());

        indexI = 0;
        if (oneDResultsMap.size() > 1) {
            System.out.println("start 2d vm process " + oneDResultsMap.keySet());
            for (String variableModI : oneDResultsMap.keySet()) {
                int indexII = 0;
                for (String variableModII : oneDResultsMap.keySet()) {
                    if (indexII <= indexI) {
                        indexII++;
                        continue;
                    }

                    Set<String> intersectionSet = SpectraUtilities.getIntersectionSet(modifiedSpectrumMap.get(variableModI), modifiedSpectrumMap.get(variableModII));
                    Set<String> fullSet = new HashSet<>(oneDResultsMap.get(variableModI).getSpecTitles());
                    fullSet.addAll(oneDResultsMap.get(variableModII).getSpecTitles());
                    int effectiveSize = fullSet.size() - intersectionSet.size();
                    final String option = variableModI + "_" + variableModII;
                    System.out.println("option : " + option + "  " + effectiveSize + "  " + oneDResultsMap.get(variableModI).getSpecTitles().size() + "  " + oneDResultsMap.get(variableModII).getSpecTitles().size());
                    if ( effectiveSize > oneDResultsMap.get(variableModI).getSpecTitles().size() && effectiveSize > oneDResultsMap.get(variableModII).getSpecTitles().size()) {
                        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                        if (addedCommonPTMs) {
                            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(commonVariableMod));

                        }
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModI));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII));
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
                        Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                            return scoreModel;
                        });
                        try {
                            RawScoreModel vModScore = f.get();
                            if (vModScore.isSensitiveChange()) {
//                                sorterSet.add(vModScore);
//                                sorterSet.add(oneDResultsMap.get(variableModI));
//                                sorterSet.add(oneDResultsMap.get(variableModII));

                                Map<String, RawScoreModel> sortinglist = new LinkedHashMap<>();
                                sortinglist.put(option, vModScore);
                                sortinglist.put(variableModI, oneDResultsMap.get(variableModI));
                                sortinglist.put(variableModII, oneDResultsMap.get(variableModII));
                                String bestScoreVariable = SpectraUtilities.compareScoresSet(sortinglist, optProtDataset.getTotalSpectraNumber());
                                if (bestScoreVariable.equalsIgnoreCase(option)) {
//                                if (sorterSet.first() == vModScore && SpectraUtilities.isBetterScore(oneDResultsMap.get(variableModI).getSpectrumMatchResult(), vModScore.getSpectrumMatchResult(), optProtDataset.getTotalSpectraNumber()) > 0 && SpectraUtilities.isBetterScore(oneDResultsMap.get(variableModII).getSpectrumMatchResult(), vModScore.getSpectrumMatchResult(), optProtDataset.getTotalSpectraNumber()) > 0) {
                                    twoDResultsMap.put(option, vModScore);
                                    modifiedSpectrumMap.put(option, SpectraUtilities.getModifiedSpectrumSet(vModScore.getSpectrumMatchResult()));
                                }
                            }
//                            sorterSet.clear();
                        } catch (ExecutionException | InterruptedException ex) {
                            ex.printStackTrace();
                        }

                    }

                    indexII++;
                }
                indexI++;
            }
        }

        MainUtilities.cleanOutputFolder();
        //process 3 d 
        if (!twoDResultsMap.isEmpty() && oneDResultsMap.size() > 2) {
            for (String variableModI : oneDResultsMap.keySet()) {
                for (String variableModII : twoDResultsMap.keySet()) {
                    if (variableModII.contains(variableModI)) {
                        continue;
                    }

                    Set<String> intersectionSet = SpectraUtilities.getIntersectionSet(modifiedSpectrumMap.get(variableModI), modifiedSpectrumMap.get(variableModII));
                    Set<String> fullSet = new HashSet<>(oneDResultsMap.get(variableModI).getSpecTitles());
                    fullSet.addAll(twoDResultsMap.get(variableModII).getSpecTitles());
                    int effectiveSize = fullSet.size() - intersectionSet.size();
                    final String option = variableModI + "_" + variableModII;

                    if (effectiveSize > oneDResultsMap.get(variableModI).getSpecTitles().size() && effectiveSize > twoDResultsMap.get(variableModII).getSpecTitles().size()) {
                        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                        if (addedCommonPTMs) {
                            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(commonVariableMod));

                        }
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModI));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[0]));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[1]));
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
                        Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                            return scoreModel;
                        });
                        try {
                            RawScoreModel vModScore = f.get();
                            if (vModScore.isSignificatChange()) {
//                                sorterSet.add(vModScore);
//                                sorterSet.add(oneDResultsMap.get(variableModI));
//                                sorterSet.add(twoDResultsMap.get(variableModII));

                                Map<String, RawScoreModel> sortinglist = new LinkedHashMap<>();
                                sortinglist.put(option, vModScore);
                                sortinglist.put(variableModI, oneDResultsMap.get(variableModI));
                                sortinglist.put(variableModII, twoDResultsMap.get(variableModII));
                                String bestScoreVariable = SpectraUtilities.compareScoresSet(sortinglist, optProtDataset.getTotalSpectraNumber());
                                if (bestScoreVariable.equalsIgnoreCase(option)) {

//                                 System.out.println("3d scoring "+option+"  "+(sorterSet.first() == vModScore)+" first  "+SpectraUtilities.isBetterScore(oneDResultsMap.get(variableModI).getData(), vModScore.getData(), false)+"  final  "+SpectraUtilities.isBetterScore(twoDResultsMap.get(variableModII).getData(), vModScore.getData(), true));
//                                if (sorterSet.first() == vModScore && SpectraUtilities.isBetterScore(oneDResultsMap.get(variableModI).getSpectrumMatchResult(), vModScore.getSpectrumMatchResult(), optProtDataset.getTotalSpectraNumber()) > 0 && SpectraUtilities.isBetterScore(twoDResultsMap.get(variableModII).getSpectrumMatchResult(), vModScore.getSpectrumMatchResult(), optProtDataset.getTotalSpectraNumber()) > 0) {
                                    threeDResultsMap.put(option, vModScore);
                                    modifiedSpectrumMap.put(option, SpectraUtilities.getModifiedSpectrumSet(vModScore.getSpectrumMatchResult()));
                                }
                            }
//                            sorterSet.clear();
                        } catch (ExecutionException | InterruptedException ex) {
                            ex.printStackTrace();
                        }

                    }
                }
            }
        }
        if (!threeDResultsMap.isEmpty() && oneDResultsMap.size() > 3) {
            for (String variableModI : oneDResultsMap.keySet()) {
                for (String variableModII : threeDResultsMap.keySet()) {
                    if (variableModII.contains(variableModI)) {
                        continue;
                    }

                    Set<String> intersectionSet = SpectraUtilities.getIntersectionSet(modifiedSpectrumMap.get(variableModI), modifiedSpectrumMap.get(variableModII));
                    Set<String> fullSet = new HashSet<>(oneDResultsMap.get(variableModI).getSpecTitles());
                    fullSet.addAll(threeDResultsMap.get(variableModII).getSpecTitles());
                    int effectiveSize = fullSet.size() - intersectionSet.size();
                    final String option = variableModI + "_" + variableModII;

                    if (effectiveSize > oneDResultsMap.get(variableModI).getSpecTitles().size() && effectiveSize > threeDResultsMap.get(variableModII).getSpecTitles().size()) {
                        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                        if (addedCommonPTMs) {
                            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(commonVariableMod));

                        }
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModI));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[0]));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[1]));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[2]));
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
                        Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                            return scoreModel;
                        });
                        try {
                            RawScoreModel vModScore = f.get();
                            if (vModScore.isSensitiveChange()) {
//                                sorterSet.add(vModScore);
//                                sorterSet.add(oneDResultsMap.get(variableModI));
//                                sorterSet.add(threeDResultsMap.get(variableModII));

                                Map<String, RawScoreModel> sortinglist = new LinkedHashMap<>();
                                sortinglist.put(option, vModScore);
                                sortinglist.put(variableModI, oneDResultsMap.get(variableModI));
                                sortinglist.put(variableModII, threeDResultsMap.get(variableModII));
                                String bestScoreVariable = SpectraUtilities.compareScoresSet(sortinglist, optProtDataset.getTotalSpectraNumber());
                                if (bestScoreVariable.equalsIgnoreCase(option)) {

//                                if (sorterSet.first() == vModScore && SpectraUtilities.isBetterScore(oneDResultsMap.get(variableModI).getSpectrumMatchResult(), vModScore.getSpectrumMatchResult(), optProtDataset.getTotalSpectraNumber()) > 0 && SpectraUtilities.isBetterScore(threeDResultsMap.get(variableModII).getSpectrumMatchResult(), vModScore.getSpectrumMatchResult(), optProtDataset.getTotalSpectraNumber()) > 0) {
                                    fourDResultsMap.put(option, vModScore);
                                    modifiedSpectrumMap.put(option, SpectraUtilities.getModifiedSpectrumSet(vModScore.getSpectrumMatchResult()));
                                }
                            }
//                            sorterSet.clear();
                        } catch (ExecutionException | InterruptedException ex) {
                            ex.printStackTrace();
                        }

                    }
                }
            }
        }
        resultsMap.clear();
        resultsMap.putAll(oneDResultsMap);
        resultsMap.putAll(twoDResultsMap);
        resultsMap.putAll(threeDResultsMap);
        resultsMap.putAll(fourDResultsMap);
        System.out.println("oneDResultsMap " + oneDResultsMap + "  twoDResultsMap " + twoDResultsMap + "  threeDResultsMap" + threeDResultsMap + "   fourDResultsMap" + fourDResultsMap);
        if (!resultsMap.isEmpty()) {
            for (String key : resultsMap.keySet()) {
                System.out.println("reach finals -----------<<>> "+key);
                sortedResultsMap.put(resultsMap.get(key), key);
            }
            String bestCombMod = SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber());//sortedResultsMap.firstEntry().getValue(); //
            System.out.println("be svariable mod " + sortedResultsMap.firstEntry().getValue() + " vs " + bestCombMod);
//            System.out.println("best comp " + bestCombMod + "  " + resultsMap.keySet());
            if (bestCombMod.contains("_")) {
                selectedVariableModificationOption.addAll(Arrays.asList(bestCombMod.split("_")));
            } else {
                selectedVariableModificationOption.add(bestCombMod);
            }
//            System.out.println("best option is " + bestCombMod);
        }

//        sortedResultsMap.clear();
//        int index = 0;
//        for (String key : resultsMap.keySet()) {
//            int index2 = 0;
//            for (String key2 : resultsMap.keySet()) {
//                if (index2 <= index) {
//                    index2++;
//                    continue;
//                }
//                boolean better = SpectraUtilities.isBetterScore(resultsMap.get(key).getData(), resultsMap.get(key2).getData(), false) > 0;
//                System.out.println(key + " vs " + key2 + "  better : " + better);
//                index2++;
//            }
//            index++;
//        }
//        for (String key : resultsMap.keySet()) {
//            sortedResultsMap.put(resultsMap.get(key), key);
//
//        }
////        selectedVariableModificationOption.clear();
//        if (!sortedResultsMap.isEmpty()) {
//            String varMod = sortedResultsMap.firstEntry().getValue();
//            if (varMod.contains("_")) {
//                selectedVariableModificationOption.addAll(Arrays.asList(varMod.split("_")));
//            } else {
//                selectedVariableModificationOption.add(varMod);
//            }
//        }
//        String vtOption = "";
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        for (String mod : selectedVariableModificationOption) {
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod));
//            vtOption += mod + "_";
        }

        if (!selectedVariableModificationOption.isEmpty()) {
            final String finalvOption = selectedVariableModificationOption.get(0);
            final String fvupdatedName = Configurations.DEFAULT_RESULT_NAME + "fv_" + finalvOption + "_" + msFileName;
            Future<RawScoreModel> vFuture = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, fvupdatedName, finalvOption, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel variableScoreModel = vFuture.get();
                double impact = Math.round((double) (variableScoreModel.getSpectrumMatchResult().size() - optProtDataset.getActiveIdentificationNum()) * 100.0 / (double) optProtDataset.getActiveIdentificationNum()) / 100.0;
                variableModParamScore.setImpact(impact);
//                if (variableScoreModel.isSensitiveChange()) {
                optProtDataset.setActiveScoreModel(variableScoreModel);
//                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        variableModParamScore.setScore(optProtDataset.getActiveIdentificationNum());
        variableModParamScore.setParamValue(selectedVariableModificationOption.toString());
        parameterScoreSet.add(variableModParamScore);
        modificationsResults.put("variableModifications", new HashSet<>(selectedVariableModificationOption));
        optProtDataset.setPotintialVariableMod(potintialVariableMod);
        MainUtilities.cleanOutputFolder();
        return modificationsResults;

    }

    public abstract RawScoreModel excuteSearch(SearchingSubDataset optProtDataset, String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addPeptideMasses, SearchInputSetting optProtSearchSettings, File identificationParametersFile, boolean pairData);

}
