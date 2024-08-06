package no.uib.probe.optprot.search;

import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.biology.modifications.ModificationCategory;
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
//                System.out.println("is significant " + option + "  -->> " + scoreModel + "   " + optProtDataset.getValidatedIdRefrenceData().length);
                if (scoreModel.isSignificatChange()) {
                    if (scoreModel.getSizeEffect() >= 0.05) {
                        continue;
                    }
                    spectraCounter = Math.max(spectraCounter, scoreModel.getSpectrumMatchResult().size());

                    resultsMap.put(option, f.get());
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
//        System.out.println("# of validated before " + optProtDataset.getValidatedIdRefrenceData().length + "  ");
        TreeMap<RawScoreModel, String> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (String option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            if (!selectedOption.equalsIgnoreCase(sortedResultsMap.firstEntry().getValue())) {
                selectedOption = sortedResultsMap.firstEntry().getValue();
//                optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey().getData());

            }
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption);
        parameterScoreSet.add(paramScore);
//        System.out.println("# of validated after " + optProtDataset.getActiveIdentificationNum());

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
    public String[] optimizeEnzymeParameter(SearchingSubDataset optProtDataset,
            File identificationParametersFile, SearchInputSetting optimisedSearchParameter,
            TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("Enzyme");
        String[] values = new String[3];
//        if (optProtDataset.getIdentificationRate() > 30) {
//            paramScore.setScore(optProtDataset.getActiveIdentificationNum());
//            paramScore.setParamValue("Trypsin");
//            parameterScoreSet.add(paramScore);
//            return "Trypsin";
//        }

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
//        int nMissesCleavages = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(selectedOption);
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
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
                    if (scoreModel.isSignificatChange() || compareBetweenEnzymes) {
                        resultsMap.put(option, scoreModel);
                    }
                } catch (InterruptedException | ExecutionException ex) {
                    ex.printStackTrace();
                }
            }
            if (compareBetweenEnzymes) {
                String enzymeName = SpectraUtilities.compareScoresSet(resultsMap,optProtDataset.getTotalSpectraNumber());
                values[0] = enzymeName;
            } else if (!resultsMap.isEmpty()) {
                String enzymeName = SpectraUtilities.compareScoresSet(resultsMap,optProtDataset.getTotalSpectraNumber());
                values[0] = enzymeName;
                optProtDataset.setActiveScoreModel(resultsMap.get(enzymeName));
            } else {
                values[0] = "Trypsin";
            }

            //optimize specifty 
            oreginaltempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
            oreginaltempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(EnzymeFactory.getInstance().getEnzyme(values[0]));
            oreginaltempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(values[0], missedClavageNumb);
            resultsMap.clear();
            for (int i = 0; i < DigestionParameters.Specificity.values().length; i++) {
                final String option = DigestionParameters.Specificity.getSpecificity(i).name();
//                System.out.println("option " + option);
//            if (option.equalsIgnoreCase(selectedOption)) {
//                continue;
//            }
                oreginaltempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(values[0], DigestionParameters.Specificity.getSpecificity(i));

                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, oreginaltempIdParam, true, optimisedSearchParameter, identificationParametersFile, false);
                    return scoreModel;
                });
                try {
                    RawScoreModel scoreModel = f.get();
                    if (scoreModel.isSignificatChange() || (compareBetweenEnzymes && scoreModel.getFinalScore() > 0)) {
                        resultsMap.put(option, scoreModel);
                    }
                } catch (ExecutionException | InterruptedException ex) {
                    ex.printStackTrace();
                }
            }
            if (compareBetweenEnzymes && !resultsMap.isEmpty()) {
                String specifty = SpectraUtilities.compareScoresSet(resultsMap,optProtDataset.getTotalSpectraNumber());
                values[1] = specifty;
            } else if (!resultsMap.isEmpty()) {
                String specifty = SpectraUtilities.compareScoresSet(resultsMap,optProtDataset.getTotalSpectraNumber());
                values[1] = specifty;
                optProtDataset.setActiveScoreModel(resultsMap.get(specifty));

            } else {
                values[1] = "specific";
            }

        }
        int spectraCounter = optProtDataset.getActiveIdentificationNum();
//        System.out.println("at specifity " + values[1]);
        oreginaltempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(values[0], DigestionParameters.Specificity.valueOf(values[1]));
///number op missed cleavage
        resultsMap.clear();
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
//                    if (scoreModel.getSpectrumMatchResult().size() < spectraCounter && scoreModel.getSizeEffect() >= 0.05) {
//                        continue;
//                    }
                    spectraCounter = Math.max(spectraCounter, scoreModel.getSpectrumMatchResult().size());
                    resultsMap.put(i + "", scoreModel);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }
        }
        String numbOfMissedCleavage = "2";
        if (!resultsMap.isEmpty()) {
            numbOfMissedCleavage = (SpectraUtilities.compareScoresSet(resultsMap,optProtDataset.getTotalSpectraNumber()));
            optProtDataset.setActiveScoreModel(resultsMap.get(numbOfMissedCleavage));
        }
//        else {
////            values[0] = "";
//        }
        values[2] = numbOfMissedCleavage;
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(Arrays.asList(values).toString());
        parameterScoreSet.add(paramScore);
//        System.out.println("selected value " + values[2]);
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
    public String optimizeSpecificityParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("Specificity");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        if (!oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name().equalsIgnoreCase("enzyme")) {
            return null;
        }
        String enzymeName = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
        String selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getSpecificity(enzymeName).name();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        for (int i = 0; i < DigestionParameters.Specificity.values().length; i++) {
            final String option = DigestionParameters.Specificity.getSpecificity(i).name();
            if (option.equalsIgnoreCase(selectedOption)) {
                continue;
            }
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(enzymeName, DigestionParameters.Specificity.getSpecificity(i));

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

        TreeMap<RawScoreModel, String> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (String option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
        }
        paramScore.setScore(optProtDataset.getActiveIdentificationNum());
        paramScore.setParamValue(selectedOption);
        parameterScoreSet.add(paramScore);

        return selectedOption;

    }

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

        TreeMap<RawScoreModel, String> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (String option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
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
//        if (!resultsMap.isEmpty()) {
//            String best = SpectraUtilities.compareScoresSet(resultsMap);
////            System.out.println("best of missed clav " + best);
//
//        }

        TreeMap<RawScoreModel, String> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (String option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = Integer.valueOf(sortedResultsMap.firstEntry().getValue());
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
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
        double[] values = new double[]{0.02, 0.05, 0.1, 0.2, 0.5};
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
            selectedOption = Double.parseDouble(SpectraUtilities.compareScoresSet(resultsMap,optProtDataset.getTotalSpectraNumber()));
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
//        if (!resultsMap.isEmpty()) {
//            String best = SpectraUtilities.compareScoresSet(resultsMap);
////            System.out.println("best of charge_ " + best);
//
//        }

        TreeMap<RawScoreModel, String> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (String option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            String[] topOption = sortedResultsMap.firstEntry().getValue().split("-")[1].split(",");
            selectedMinChargeOption = Integer.parseInt(topOption[0]);
            selectedMaxChargeOption = Integer.parseInt(topOption[1]);
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
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
        TreeMap<RawScoreModel, String> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (String option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            String[] topOption = sortedResultsMap.firstEntry().getValue().split("_")[1].split(",");
            selectedMinIsotopicCorrectioneOption = Integer.parseInt(topOption[0]);
            selectedMaxIsotopicCorrectionOption = Integer.parseInt(topOption[1]);
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
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
        double[] iValues = new double[]{10, 15, 20, 25};
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
                } else if (i > selectedOption && !scoreModel.isSensitiveChange()) {
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
//        if (!resultsMap.isEmpty()) {
//            String best = SpectraUtilities.compareScoresSet(resultsMap);
//            System.out.println("best tolerance is " + best);
//        }

        TreeMap<RawScoreModel, Double> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (String option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), Double.valueOf(option));
            }
            selectedOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setActiveScoreModel(sortedResultsMap.firstEntry().getKey());
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
        TreeSet<String> commonMods = new TreeSet<>(potintialVariableMod);
        String commonFixedMod = "Carbamidomethylation of C";
        String commonVariableMod = "Oxidation of M";
        boolean addedCommonfixed = false;
        Map<String, RawScoreModel> targtedFixedModificationScore = new TreeMap<>();
        RawScoreModel tempRawScore = null;

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
            System.out.println("common ptm " + commonFixedMod + "  " + scoreModel);
            if (scoreModel.isSensitiveChange()) {      
                addedCommonfixed=true;
                targtedFixedModificationScore.put("C", scoreModel);
               optProtDataset.setActiveScoreModel(scoreModel);
            }
        } catch (ExecutionException | InterruptedException ex) {
            ex.printStackTrace();
        }       
        for (String modId : mods) {
            if (modId.equalsIgnoreCase(commonFixedMod)) {
                continue;
            }
            String modPattern = ptmFactory.getModification(modId).getPattern().toString();
            if (modPattern.equalsIgnoreCase("C")&&addedCommonfixed) {
                continue;
            }
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
            tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(modId));
            tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(modId));
            double spectraNum=optProtDataset.getActiveIdentificationNum();
            if (addedCommonfixed) {
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(commonFixedMod));
                tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(commonFixedMod));
            }
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                resultsMap.put(modId, f.get());
                System.out.println("Fixed mod " + modId + "  " + resultsMap.get(modId) + "  " +"  "+spectraNum );

                if (resultsMap.get(modId).getSpectrumMatchResult()!=null && resultsMap.get(modId).getSpectrumMatchResult().size()*1.1>spectraNum && !resultsMap.get(modId).isSameData()) {
                    System.out.println("potintial to add " + modId + "  " + resultsMap.get(modId));
                    potintialVariableMod.add(modId);
                    
                }
                if (resultsMap.get(modId).isSignificatChange()) {
                    sortedResultsMap.put(resultsMap.get(modId), modId);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }
//System.exit(0);
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
            if (key.equals("C") && commonFixedMod != null) {
                continue;
            }
            RawScoreModel scoreModel = targtedFixedModificationScore.get(key);
            resultsMap.put(sortedResultsMap.get(scoreModel), scoreModel);
            modLimit++;
            if (modLimit > 3) {
                break;
            }
        }
//        if (commonFixedMod != null && !resultsMap.containsKey(commonFixedMod)) {
//            resultsMap.put(commonFixedMod, targtedFixedModificationScore.get("C"));
//        }
        sortedResultsMap.clear();
        Set<String> clearSet = new LinkedHashSet<>(resultsMap.keySet());
        TreeSet<RawScoreModel> sorterSet = new TreeSet<>(Collections.reverseOrder());
        for (String modId : clearSet) {
            /*if exist better keep it as fixed mod*/
//            if (modId.equals("Carbamidomethylation of C")) {
//                continue;
//            }
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "fAsm_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            if (commonFixedMod != null) {
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(commonFixedMod));
                tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(commonFixedMod));

            }
            if (commonVariableMod != null) {
                tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(commonVariableMod));

            }
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel vScore = f.get();
                if (vScore.isSignificatChange()) {

//                    double[] compscores = SpectraUtilities.compareData(resultsMap.get(modId), vScore.getData(), true);
//                    if (compscores[3] > resultsMap.get(modId).getFinalScore() * 1.1) {
                    if (SpectraUtilities.isBetterScore(resultsMap.get(modId).getSpectrumMatchResult(), vScore.getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber()) > 0) {
                        resultsMap.remove(modId);
                        if (modId.equalsIgnoreCase("Oxidation of M")) {
                            selectedVariableModificationOption.add(modId);
                        }
                    }
                }

                sorterSet.clear();
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }

        Map<String, RawScoreModel> oneDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>(resultsMap));
        Map<String, RawScoreModel> twoDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
        tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
        if (commonVariableMod != null) {
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(commonVariableMod));

        }

        int indexI = 0;
        sorterSet.clear();
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
                if (commonFixedMod != null) {
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
                    sorterSet.add(score2D);
                    sorterSet.add(oneDResultsMap.get(mod1Id));
                    sorterSet.add(oneDResultsMap.get(mod2Id));
                    if (sorterSet.first() == score2D) {
                        twoDResultsMap.put(option, score2D);
                    }
                    sorterSet.clear();
                } catch (ExecutionException | InterruptedException ex) {
                    ex.printStackTrace();
                }
                indexII++;
            }
            indexI++;
        }
        sorterSet.clear();
        Map<String, RawScoreModel> threeDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        if (!twoDResultsMap.isEmpty()) {
            for (String mod1Id : oneDResultsMap.keySet()) {
                for (String mod2Id : twoDResultsMap.keySet()) {
                    if (mod2Id.contains(mod1Id)) {
                        continue;
                    }
                    tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
                    if (commonFixedMod != null) {
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
                        sorterSet.add(score3D);
                        sorterSet.add(oneDResultsMap.get(mod1Id));
                        sorterSet.add(twoDResultsMap.get(mod2Id));
                        if (sorterSet.first() == score3D) {
                            threeDResultsMap.put(option, score3D);
                        }
                        sorterSet.clear();
                    } catch (ExecutionException | InterruptedException ex) {
                        ex.printStackTrace();
                    }
                }
            }

        }
//        resultsMap.clear();
//
        Map<String, RawScoreModel> fourDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        if (!threeDResultsMap.isEmpty()) {
            for (String mod1Id : oneDResultsMap.keySet()) {
                for (String mod2Id : threeDResultsMap.keySet()) {
                    if (mod2Id.contains(mod1Id)) {
                        continue;
                    }
                    tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
                    tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
                    if (commonFixedMod != null) {
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
                        sorterSet.add(score4D);
                        sorterSet.add(oneDResultsMap.get(mod1Id));
                        sorterSet.add(threeDResultsMap.get(mod2Id));
                        if (sorterSet.first() == score4D) {
                            fourDResultsMap.put(option, score4D);
                        }
                        sorterSet.clear();
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
        if (commonFixedMod != null) {
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
//        System.out.println("final fixed mod " + selectedFixedModificationOption);
//        System.out.println("final potintial v mod " + potintialVariableMod);
        //update the optprot data 
        if (commonVariableMod != null) {
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(commonVariableMod));

        }
        if (true) {
            final String option = "Oxidation of M_" + commonFixedMod;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "finalFixed_" + option + "_" + msFileName;

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel score4D = f.get();
                optProtDataset.setActiveScoreModel(score4D);
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }

        //process variable modifications
        final ParameterScoreModel variableModParamScore = new ParameterScoreModel();
        variableModParamScore.setParamId("VariableModifications");
        resultsMap.clear();
        oneDResultsMap.clear();
        twoDResultsMap.clear();
        threeDResultsMap.clear();
        fourDResultsMap.clear();
        sorterSet.clear();
        TreeMap<RawScoreModel, String> sortingModificationMap = new TreeMap<>(Collections.reverseOrder());
        Map<String, Set<String>> modifiedSpectrumMap = new LinkedHashMap<>();
//        Map<String, Set<String>> fullSpectrumMap = new LinkedHashMap<>();
        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
        //check oxidation as variable mod

        for (String modId : potintialVariableMod) {  //mods
            if (selectedFixedModificationOption.contains(modId) || modId.equalsIgnoreCase(commonVariableMod)) {
                continue;
            }
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));

            if (commonVariableMod != null) {
                tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(commonVariableMod));

            }

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                System.out.println("variable mod " + modId + "  " + scoreModel + "  total: " + optProtDataset.getActiveIdentificationNum() + "  " + tempIdParam.getSearchParameters().getModificationParameters().getFixedModifications() + "   ");
//                System.out.println("xtandemParameters "+xtandemParameters.isRefine()+"   "+xtandemParameters.isQuickPyrolidone()+"  "+xtandemParameters.isProteinQuickAcetyl());
                if (scoreModel.isSensitiveChange()) {
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
                if (modcounter > 3) {
                    break;
                }
            }
        }
        //start process 2d variable mod
        sorterSet.clear();
        indexI = 0;
        if (oneDResultsMap.size() > 1) {
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
                    if (effectiveSize > oneDResultsMap.get(variableModI).getSpecTitles().size() && effectiveSize > oneDResultsMap.get(variableModII).getSpecTitles().size()) {
                        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                        if (commonVariableMod != null) {
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
                                sorterSet.add(vModScore);
                                sorterSet.add(oneDResultsMap.get(variableModI));
                                sorterSet.add(oneDResultsMap.get(variableModII));
                                if (sorterSet.first() == vModScore && SpectraUtilities.isBetterScore(oneDResultsMap.get(variableModI).getSpectrumMatchResult(), vModScore.getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber()) > 0 && SpectraUtilities.isBetterScore(oneDResultsMap.get(variableModII).getSpectrumMatchResult(), vModScore.getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber()) > 0) {
                                    twoDResultsMap.put(option, vModScore);
                                    modifiedSpectrumMap.put(option, SpectraUtilities.getModifiedSpectrumSet(vModScore.getSpectrumMatchResult()));
                                }
                            }
                            sorterSet.clear();
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
                        if (commonVariableMod != null) {
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
                                sorterSet.add(vModScore);
                                sorterSet.add(oneDResultsMap.get(variableModI));
                                sorterSet.add(twoDResultsMap.get(variableModII));

//                                 System.out.println("3d scoring "+option+"  "+(sorterSet.first() == vModScore)+" first  "+SpectraUtilities.isBetterScore(oneDResultsMap.get(variableModI).getData(), vModScore.getData(), false)+"  final  "+SpectraUtilities.isBetterScore(twoDResultsMap.get(variableModII).getData(), vModScore.getData(), true));
                                if (sorterSet.first() == vModScore && SpectraUtilities.isBetterScore(oneDResultsMap.get(variableModI).getSpectrumMatchResult(), vModScore.getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber()) > 0 && SpectraUtilities.isBetterScore(twoDResultsMap.get(variableModII).getSpectrumMatchResult(), vModScore.getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber())> 0) {
                                    threeDResultsMap.put(option, vModScore);
                                    modifiedSpectrumMap.put(option, SpectraUtilities.getModifiedSpectrumSet(vModScore.getSpectrumMatchResult()));
                                }
                            }
                            sorterSet.clear();
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
                        if (commonVariableMod != null) {
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
                                sorterSet.add(vModScore);
                                sorterSet.add(oneDResultsMap.get(variableModI));
                                sorterSet.add(threeDResultsMap.get(variableModII));
                                if (sorterSet.first() == vModScore && SpectraUtilities.isBetterScore(oneDResultsMap.get(variableModI).getSpectrumMatchResult(), vModScore.getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber()) > 0 && SpectraUtilities.isBetterScore(threeDResultsMap.get(variableModII).getSpectrumMatchResult(), vModScore.getSpectrumMatchResult(),optProtDataset.getTotalSpectraNumber())> 0) {
                                    fourDResultsMap.put(option, vModScore);
                                    modifiedSpectrumMap.put(option, SpectraUtilities.getModifiedSpectrumSet(vModScore.getSpectrumMatchResult()));
                                }
                            }
                            sorterSet.clear();
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
        sortedResultsMap.clear();
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
        for (String key : resultsMap.keySet()) {
            sortedResultsMap.put(resultsMap.get(key), key);
        }
        selectedVariableModificationOption.clear();
        if (!sortedResultsMap.isEmpty()) {
            String varMod = sortedResultsMap.firstEntry().getValue();
            if (varMod.contains("_")) {
                selectedVariableModificationOption.addAll(Arrays.asList(varMod.split("_")));
            } else {
                selectedVariableModificationOption.add(varMod);
            }
        }
        if (commonVariableMod != null) {
            selectedVariableModificationOption.add(commonVariableMod);
        }
        String vtOption = "";
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        for (String mod : selectedVariableModificationOption) {
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod));
            vtOption += mod + "_";
        }
        final String finalvOption = vtOption;
        final String fvupdatedName = Configurations.DEFAULT_RESULT_NAME + "fv_" + vtOption + "_" + msFileName;
        Future<RawScoreModel> vFuture = MainUtilities.getExecutorService().submit(() -> {
            RawScoreModel scoreModel = excuteSearch(optProtDataset, fvupdatedName, finalvOption, tempIdParam, false, searchInputSetting, identificationParametersFile, false);
            return scoreModel;
        });
        try {
            RawScoreModel variableScoreModel = vFuture.get();
            if (variableScoreModel.isSensitiveChange()) {
                optProtDataset.setActiveScoreModel(variableScoreModel);
            }
        } catch (ExecutionException | InterruptedException ex) {
            ex.printStackTrace();
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
