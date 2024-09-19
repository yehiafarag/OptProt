package no.uib.probe.optprot.search;

import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.experiment.biology.ions.impl.PeptideFragmentIon;
import com.compomics.util.experiment.biology.modifications.Modification;
import com.compomics.util.experiment.biology.modifications.ModificationCategory;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import com.compomics.util.experiment.identification.Advocate;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.io.IoUtil;
import com.compomics.util.parameters.identification.IdentificationParameters;
import com.compomics.util.parameters.identification.search.DigestionParameters;
import com.compomics.util.parameters.identification.search.SearchParameters;
import com.compomics.util.parameters.identification.tool_specific.DirecTagParameters;
import com.compomics.util.parameters.identification.tool_specific.NovorParameters;
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
import org.apache.batik.svggen.font.table.Table;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

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
        RawScoreModel oreginalScore = new RawScoreModel("CleavageParameter");
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
        double threshold = 1;
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
                    if (scoreModel.getFinalScore() > threshold) {
                        resultsMapI.put(option, scoreModel);
                    }
                    if (compareBetweenEnzymes && scoreModel.getFinalScore() > threshold) {
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
                    if (scoreModel.getFinalScore() > threshold) {
                        resultsMapI.put(option, scoreModel);
                    }
                    if ((compareBetweenEnzymes && scoreModel.getFinalScore() > threshold)) {
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
        double threshold = 1 + optProtDataset.getComparisonsThreshold();
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
                if ((i < selectedOption && scoreModel.isSensitiveChange()) || (scoreModel.getFinalScore() > threshold && scoreModel.getTotalNumber() >= optProtDataset.getCurrentScoreModel().getTotalNumber())) {
                    resultsMap.put(i + "", scoreModel);
                    threshold++;
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
        double threshold = 1;
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
                if (scoreModel.getFinalScore() > threshold) {
                    threshold += 0.5;
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
                        if (scoreModel.getFinalScore() > threshold) {
                            threshold += 0.5;
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

        final ParameterScoreModel fixedModParamScore = new ParameterScoreModel();
        fixedModParamScore.setParamId("FixedModifications");

        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        ArrayList<String> selectedFixedModificationOption = new ArrayList<>();
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
//        TreeMap<RawScoreModel, String> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        Set<String> potintialMods = new LinkedHashSet<>();//ptmFactory.getModifications(ModificationCategory.Common)
        String commonFixedMod = "Carbamidomethylation of C";
        potintialMods.add(commonFixedMod);
        String commonVariableMod = "Oxidation of M";
//        boolean addedCommonFixedPTMs = false;
        Map<String, RawScoreModel> targtedFixedModificationScore = new TreeMap<>();
        Map<String, RawScoreModel> fullFixedModificationScore = new LinkedHashMap<>();
        ArrayList<SortedPTMs> fullFixedModificationEffect = new ArrayList<>();
        //first stage common fixed modification
        String prefix = "f_";
        resultsMap.putAll(this.checkModificationsScores(selectedFixedModificationOption, selectedVariableModificationOption, potintialMods, true, msFileName, tempIdParam, optProtDataset, identificationParametersFile, searchInputSetting, prefix, false));
        if (!resultsMap.isEmpty()) {
            String bestMod = SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber());
            selectedFixedModificationOption.add(bestMod);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestMod));
            potintialMods.clear();
//            addedCommonFixedPTMs = true;
            targtedFixedModificationScore.put("C", resultsMap.get(bestMod));
            MainUtilities.cleanOutputFolder();
            resultsMap.clear();

        }
        //try variable coomon mod first
        // stage 2 test for common variable modification   
        final ParameterScoreModel variableModParamScore = new ParameterScoreModel();
        variableModParamScore.setParamId("VariableModifications");
//        if (!addedCommonVariablePTMs) {
        potintialMods.add(commonVariableMod);
        prefix = "v_";
        resultsMap.putAll(this.checkModificationsScores(selectedFixedModificationOption, selectedVariableModificationOption, potintialMods, false, msFileName, tempIdParam, optProtDataset, identificationParametersFile, searchInputSetting, prefix, false));
        if (!resultsMap.isEmpty()) {
            String bestMod = SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber());
            selectedVariableModificationOption.add(bestMod);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestMod));
            potintialMods.clear();
            MainUtilities.cleanOutputFolder();
            resultsMap.clear();
        }
        MainUtilities.cleanOutputFolder();
        //process variable modifications
        mods.removeAll(selectedFixedModificationOption);
        mods.removeAll(selectedVariableModificationOption);
        //second stage fixed modifications
        for (String modId : mods) {
            String modPattern = ptmFactory.getModification(modId).getPattern().toString();
            if (modPattern.equals("")) {
                modPattern = ptmFactory.getModification(modId).getModificationType().isNTerm() + "-" + ptmFactory.getModification(modId).getModificationType().isCTerm();
            }
            if (targtedFixedModificationScore.containsKey(modPattern)) {
                continue;
            }
            if (!selectedVariableModificationOption.isEmpty() && modPattern.equalsIgnoreCase("M")) {
                continue;
            }
            potintialMods.add(modId);
        }
        prefix = "f_";
        fullFixedModificationScore.putAll(this.checkModificationsScores(selectedFixedModificationOption, selectedVariableModificationOption, potintialMods, true, msFileName, tempIdParam, optProtDataset, identificationParametersFile, searchInputSetting, prefix, true));
        List<Double> scoresList = new ArrayList<>();
        for (String modId : fullFixedModificationScore.keySet()) {
            RawScoreModel scoreModel = fullFixedModificationScore.get(modId);
            if (scoreModel.isSensitiveChange()) {
                resultsMap.put(modId, scoreModel);
            }
            double avg1 = SpectraUtilities.getModificationFrequentScore(modId, scoreModel.getSpecTitles(), optProtDataset.getCurrentScoreModel().getSpecTitles());
            double b1 = avg1;
            if (avg1 > -1) {
                scoresList.add(b1);
            }
            fullFixedModificationEffect.add(new SortedPTMs(modId, b1, 0));

        }

        potintialMods.clear();
        potintialMods.addAll(resultsMap.keySet());
        prefix = "FAV_";
        Map<String, RawScoreModel> clearResults = this.checkModificationsScores(selectedFixedModificationOption, selectedVariableModificationOption, potintialMods, false, msFileName, tempIdParam, optProtDataset, identificationParametersFile, searchInputSetting, prefix, false);
        for (String modId : clearResults.keySet()) {
            RawScoreModel vScore = clearResults.get(modId);
            if (vScore.isSensitiveChange()) {
                if (SpectraUtilities.isBetterScore(resultsMap.get(modId).getSpectrumMatchResult(), vScore.getSpectrumMatchResult(), optProtDataset.getTotalSpectraNumber()) > 0) {
                    resultsMap.remove(modId);
                    System.out.println("remove from fixed " + modId + "   left in there " + resultsMap.keySet());
                }
            }
        }
//        System.exit(0);
        //2nd fixed modifications
        if (!resultsMap.isEmpty()) {
            String bestMod = SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber());
            selectedFixedModificationOption.add(bestMod);
            optProtDataset.setActiveScoreModel(resultsMap.get(bestMod));
            potintialMods.clear();
            System.out.println("check best " + bestMod);
            String modPattern = ptmFactory.getModification(bestMod).getPattern().toString();
            System.out.println("mod pattern " + modPattern);
            if (modPattern.equals("")) {
                modPattern = ptmFactory.getModification(bestMod).getModificationType().isNTerm() + "-" + ptmFactory.getModification(bestMod).getModificationType().isCTerm();
            }
            System.out.println("keep going");
            targtedFixedModificationScore.put(modPattern, resultsMap.get(bestMod));
            MainUtilities.cleanOutputFolder();
            resultsMap.remove(bestMod);
        }
        //3rd and 4th  fixed modification processing 
        boolean test = !resultsMap.isEmpty();
        int counter = 3;
        while (test) {
            for (String modId : resultsMap.keySet()) {
                String modPattern = ptmFactory.getModification(modId).getPattern().toString();
                if (modPattern.equals("")) {
                    modPattern = ptmFactory.getModification(modId).getModificationType().isNTerm() + "-" + ptmFactory.getModification(modId).getModificationType().isCTerm();
                }

                if (targtedFixedModificationScore.containsKey(modPattern)) {
                    continue;
                }
                System.out.println("mod pattern " + modPattern + "  " + targtedFixedModificationScore.keySet());
                potintialMods.add(modId);
            }
            resultsMap.clear();
            if (!potintialMods.isEmpty()) {
                System.out.println("---run fixed mod -->>" + potintialMods);
                prefix = "f_";
                resultsMap.putAll(this.checkModificationsScores(selectedFixedModificationOption, selectedVariableModificationOption, potintialMods, true, msFileName, tempIdParam, optProtDataset, identificationParametersFile, searchInputSetting, prefix, false));
                if (!resultsMap.isEmpty()) {
                    String bestMod = SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber());
                    selectedFixedModificationOption.add(bestMod);
                    optProtDataset.setActiveScoreModel(resultsMap.get(bestMod));
                    potintialMods.clear();
                    String modPattern = ptmFactory.getModification(bestMod).getPattern().toString();
                    targtedFixedModificationScore.put(modPattern, resultsMap.get(bestMod));
                    MainUtilities.cleanOutputFolder();
                    resultsMap.remove(bestMod);
                    test = !resultsMap.isEmpty() || counter <= 4;
                    counter++;
                } else {
                    test = false;
                }
            } else {
                break;
            }

        }
        modificationsResults.put("fixedModifications", new HashSet<>(selectedFixedModificationOption));
        modificationsResults.put("refinmentFixedModifications", new HashSet<>(selectedFixedModificationOption));

        MainUtilities.cleanOutputFolder();
        fixedModParamScore.setScore(optProtDataset.getActiveIdentificationNum());
        fixedModParamScore.setParamValue(selectedFixedModificationOption.toString());
        parameterScoreSet.add(fixedModParamScore);
        resultsMap.clear();
        potintialMods.clear();
        double[] scoreArr = scoresList.stream().mapToDouble(Double::doubleValue).toArray();
        DescriptiveStatistics scoreDS = new DescriptiveStatistics(scoreArr);
        double avgScore = scoreDS.getPercentile(40);
        System.out.println("data from scores q1    " + scoreDS.getPercentile(25) + "  median " + scoreDS.getPercentile(50) + "    q3 " + scoreDS.getPercentile(75) + "   avg " + avgScore);

        Collections.sort(fullFixedModificationEffect);
        int ptmcounter = 0;
        for (SortedPTMs modification : fullFixedModificationEffect) {
//            Modification mod = ptmFactory.getModification(modification.getName());
            //get modified intersection with avctive spectra 

//            String modPattern = mod.getPattern().toString();
//            if (modPattern.equals("")) {
//                modPattern = mod.getModificationType().isNTerm() + "-" + mod.getModificationType().isCTerm();
//            }
//            if (targtedFixedModificationScore.containsKey(modPattern)) {
//                continue;
//            }
//            if (mod.getModificationType().isCTerm() || mod.getModificationType().isNTerm()) {
//                potintialMods.add(modification.getName());
//                continue;
//            }
            if (modification.getScore() >= avgScore || modification.getScore() == 0) {
                potintialMods.add(modification.getName());
                System.out.println(ptmcounter + " mod name   " + modification.getName() + "  mod score " + modification.getScore() + "   " + avgScore);
                ptmcounter++;
            }

        }

        System.out.println("  " + " --> potintial " + potintialMods.size() + "  / " + (potintialMods.size() / 2) + "  vs " + "  " + "  " + potintialMods);
//        System.exit(0);
//        System.exit(0);
        resultsMap.clear();
        prefix = "v_";
        counter = 0;
        double thre = 0;
        Map<String, TreeMap<Double, String>> filterVMMap = new LinkedHashMap<>();
        while (counter < 4) {
            resultsMap.putAll(this.checkModificationsScores(selectedFixedModificationOption, selectedVariableModificationOption, potintialMods, false, msFileName, tempIdParam, optProtDataset, identificationParametersFile, searchInputSetting, counter + "" + prefix, false));
            filterVMMap.clear();
            if (!resultsMap.isEmpty()) {
                for (String modId : resultsMap.keySet()) {
                    Modification mod = ptmFactory.getModification(modId);
                    //get modified intersection with avctive spectra 
                    String modPattern = mod.getPattern().toString();
                    if (modPattern.equals("")) {
                        modPattern = mod.getModificationType().isNTerm() + "-" + mod.getModificationType().isCTerm();
                    }
                    if (!filterVMMap.containsKey(modPattern)) {
                        filterVMMap.put(modPattern, new TreeMap<>(Collections.reverseOrder()));
                    }
                    filterVMMap.get(modPattern).put(resultsMap.get(modId).getFinalScore(), modId);
                }
                boolean first;
                System.out.println("filtered vm are " + filterVMMap);
                Set<String> toKeep = new HashSet();
                for (String patteren : filterVMMap.keySet()) {
                    first = true;
                    TreeMap<Double, String> modPatMap = filterVMMap.get(patteren);
                    for (String mod : modPatMap.values()) {
                        if (first) {
                            first = false;
                            continue;
                        }
                        toKeep.add(mod);
                    }
                }
                for (String remove : toKeep) {
                    resultsMap.remove(remove);
                }
                String bestMod = SpectraUtilities.compareScoresSet(resultsMap, optProtDataset.getTotalSpectraNumber());
                System.out.println("------------------------------------------------------->>>> result map round " + counter + "  " + bestMod + "  " + resultsMap.get(bestMod).isSignificatChange());
                System.out.println("best vm score " + bestMod + "  " + resultsMap.get(bestMod) + "   " + optProtDataset.getCurrentScoreModel().getSpectrumMatchResult().size());

                if (resultsMap.get(bestMod).getFinalScore() > thre) {
                    selectedVariableModificationOption.add(bestMod);
                    optProtDataset.setActiveScoreModel(resultsMap.get(bestMod));
                    MainUtilities.cleanOutputFolder();
                    resultsMap.remove(bestMod);
                    thre += 0.5;
                } else {
                    resultsMap.clear();
                }
                if (resultsMap.isEmpty()) {
                    break;
                }

                TreeSet<SortedPTMs> sorePtms = new TreeSet<>(Collections.reverseOrder());
                for (String mod : resultsMap.keySet()) {
                    if (resultsMap.get(mod).isSignificatChange()) {
                        sorePtms.add(new SortedPTMs(mod, resultsMap.get(mod).getFinalScore(), 0));

                    }
                }
                resultsMap.clear();
                potintialMods.clear();
                int subCounter = selectedVariableModificationOption.size();
                for (SortedPTMs mod : sorePtms) {
                    if (subCounter > 4) {
                        break;
                    }
                    potintialMods.add(mod.getName());
                    subCounter++;
                    System.out.println(counter + "-->> mod " + mod.getName() + "  " + mod.getScore());

                }
                counter++;

            } else {
                break;
            }
        }
        variableModParamScore.setScore(optProtDataset.getActiveIdentificationNum());
        variableModParamScore.setParamValue(selectedVariableModificationOption.toString());
        System.out.println("final vm " + selectedVariableModificationOption);
        parameterScoreSet.add(variableModParamScore);
        modificationsResults.put("variableModifications", new HashSet<>(selectedVariableModificationOption));
//        optProtDataset.setPotintialVariableMod(selectedVariableModificationOption);
        MainUtilities.cleanOutputFolder();
        return modificationsResults;

    }

    private Map<String, RawScoreModel> checkModificationsScores(ArrayList<String> selectedFixedModificationOption, ArrayList<String> selectedVariableModificationOption, Set<String> modifications, boolean fixed, String msFileName, IdentificationParameters tempIdParam, SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting searchInputSetting, String prefix, boolean addAll) {
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        for (String modId : modifications) {
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + prefix + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
//            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));
            for (String fixedMod : selectedFixedModificationOption) {
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(fixedMod));
                tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(fixedMod));
            }
            for (String variableMod : selectedVariableModificationOption) {
                tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableMod));
            }
            if (fixed) {
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(modId));
                tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(modId));
            } else {
                tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));
            }
            System.out.println("modification " + modId + "  " + selectedFixedModificationOption + "   " + selectedVariableModificationOption);
            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSensitiveChange() || addAll) {
                    resultsMap.put(modId, scoreModel);//                 
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }
        tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
        tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        return resultsMap;

    }

    public abstract RawScoreModel excuteSearch(SearchingSubDataset optProtDataset, String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addPeptideMasses, SearchInputSetting optProtSearchSettings, File identificationParametersFile, boolean pairData);

}
