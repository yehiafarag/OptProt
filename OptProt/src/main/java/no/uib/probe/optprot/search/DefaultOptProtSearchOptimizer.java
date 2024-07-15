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
import no.uib.probe.optprot.util.SpectraFileUtilities;

/**
 *
 * @author yfa041
 */
public abstract class DefaultOptProtSearchOptimizer {

    /**
     * The compomics PTM factory.
     */
    private final ModificationFactory ptmFactory = ModificationFactory.getInstance();

    public DefaultOptProtSearchOptimizer() {
    }

    public String optimizeDigestionParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
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
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(option, f.get());;
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
            if (!selectedOption.equalsIgnoreCase(sortedResultsMap.firstEntry().getValue())) {
                selectedOption = sortedResultsMap.firstEntry().getValue();
                optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());

            }
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
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
    public String optimizeEnzymeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("Enzyme");

        if (optProtDataset.getIdentificationRate() > 30) {
            paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
            paramScore.setParamValue("Trypsin");
            parameterScoreSet.add(paramScore);
            return "Trypsin";
        }

        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        if (!oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name().equalsIgnoreCase("enzyme")) {
            return null;
        }
        String selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
        int nMissesCleavages = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(selectedOption);
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        //optimise enzyme  
        Enzyme[] enzymes = new Enzyme[EnzymeFactory.getInstance().getEnzymes().size() - 1];
//        enzymes[0] = EnzymeFactory.getInstance().getEnzyme("Trypsin");
        int i = 0;
        for (Enzyme enzyme : EnzymeFactory.getInstance().getEnzymes()) {
            if (enzyme.getName().equals("Trypsin")) {
                continue;
            }
            enzymes[i] = enzyme;
            i++;
        }
        for (Enzyme enzyme : enzymes) {
            if (enzyme.getName().replace(" ", "").equalsIgnoreCase("Trypsin(noPrule)") || enzyme.getName().equalsIgnoreCase(selectedOption)) {
                continue;
            }
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
            tempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(enzyme);
            tempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(enzyme.getName(), nMissesCleavages);
            final String option = enzyme.getName();
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
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
        TreeMap<RawScoreModel, String> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (String option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        paramScore.setParamValue(selectedOption);
        parameterScoreSet.add(paramScore);
        return selectedOption;

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
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
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
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }

        selectedOption = selectedOption.replace("[", "").replace("]", "");
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
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
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
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
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(i, scoreModel);
                } else if (i > selectedOption && !scoreModel.isSensitiveChange()) {
                    break;
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
            selectedOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
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
        Map<Double, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double[] values = new double[]{0.02, 0.05, 0.07, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5};
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
                    resultsMap.put(i, scoreModel);
                } else if (i > selectedOption && !scoreModel.isSensitiveChange()) {
                    break;
                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }
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

    public int[] optimizePrecursorChargeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        final ParameterScoreModel paramScore = new ParameterScoreModel();
        paramScore.setParamId("charge");
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxChargeOption = oreginaltempIdParam.getSearchParameters().getMaxChargeSearched();
        int selectedMinChargeOption = oreginaltempIdParam.getSearchParameters().getMinChargeSearched();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);

        for (int i = 1; i < 5; i++) {

            for (int j = 2; j <= 5; j++) {
                if (j <= i) {
                    continue;
                }
                tempIdParam.getSearchParameters().setMinChargeSearched(i);
                tempIdParam.getSearchParameters().setMaxChargeSearched(j);

                final String option = "charge_" + i + "," + j;
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
            selectedMinChargeOption = Integer.parseInt(topOption[0]);
            selectedMaxChargeOption = Integer.parseInt(topOption[1]);
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
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
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }
        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
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
        Map<Double, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
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
                    resultsMap.put(i, scoreModel);
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
                            resultsMap.put(i, scoreModel);
                        }
                    } catch (ExecutionException | InterruptedException ex) {
                        ex.printStackTrace();
                    }

                }
            }
        }
        TreeMap<RawScoreModel, Double> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());
        if (!resultsMap.isEmpty()) {
            for (double option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
        }

        paramScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);

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
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
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
        Set<String> potintialVariableMod = new LinkedHashSet<>();
        for (String modId : mods) {
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(modId));

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                resultsMap.put(modId, f.get());
                if (resultsMap.get(modId).getImprovmentScore() > -0.5) {
                    potintialVariableMod.add(modId);
                }
                if (resultsMap.get(modId).isSignificatChange()) {
                    sortedResultsMap.put(resultsMap.get(modId), modId);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }
        Map<String, RawScoreModel> targtedFixedModificationScore = new TreeMap<>();
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
            RawScoreModel scoreModel = targtedFixedModificationScore.get(key);
            resultsMap.put(sortedResultsMap.get(scoreModel), scoreModel);
            modLimit++;
            if (modLimit > 3) {
                break;
            }
        }
        sortedResultsMap.clear();
        Set<String> clearSet = new LinkedHashSet<>(resultsMap.keySet());
        TreeSet<RawScoreModel> sorterSet = new TreeSet<>(Collections.reverseOrder());
        for (String modId : clearSet) {
            /*if exist better keep it as fixed mod*/
            if (modId.equals("Carbamidomethylation of C")) {
                continue;
            }
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "fAsm_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, true);
                return scoreModel;
            });
            try {
                RawScoreModel vScore = f.get();
                if (vScore.isSignificatChange()) {
//                    double[] compscores = SpectraFileUtilities.compareData(resultsMap.get(modId).getData(), vScore.getData(), true);
//                    if (compscores[3] > resultsMap.get(modId).getFinalScore() * 1.1) {
                    if (SpectraFileUtilities.isBetterScore(resultsMap.get(modId).getData(), vScore.getData(), true)) {
                        resultsMap.remove(modId);
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
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod1Id));
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id));
                final String option = mod1Id + "_" + mod2Id;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;
                Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, searchInputSetting, identificationParametersFile, true);
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
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod1Id));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[0]));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[1]));
                    final String option = mod1Id + "_" + mod2Id;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;
                    Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                        RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, searchInputSetting, identificationParametersFile, true);
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
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod1Id));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[0]));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[1]));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[2]));
                    final String option = mod1Id + "_" + mod2Id;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;

                    Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                        RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, searchInputSetting, identificationParametersFile, true);
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
        if (!sortedResultsMap.isEmpty()) {
            String fixedMod = sortedResultsMap.firstEntry().getValue();
            if (fixedMod.contains("_")) {
                selectedFixedModificationOption.addAll(Arrays.asList(fixedMod.split("_")));
            } else {
                selectedFixedModificationOption.add(fixedMod);
            }
        }
        modificationsResults.put("fixedModifications", new HashSet<>(selectedFixedModificationOption));
        MainUtilities.cleanOutputFolder();
        tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
        tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();

        //process refine fixed modifications
        String ftoption = "";
        for (String fixedMod : selectedFixedModificationOption) {
            tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(fixedMod));
            ftoption += fixedMod + "_";
        }
        final String ftupdatedName = Configurations.DEFAULT_RESULT_NAME + "ff_" + ftoption + "_" + msFileName;
        for (String fixedModifications : resultsMap.keySet()) {
            for (String str : fixedModifications.split("_")) {
                tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(str));
            }
        }
        modificationsResults.put("refinmentFixedModifications", new HashSet<>(selectedFixedModificationOption));
        MainUtilities.cleanOutputFolder();
        RawScoreModel fixedModScore;

        final String finalFTOption = ftoption;
        Future<RawScoreModel> fFuture = MainUtilities.getExecutorService().submit(() -> {
            RawScoreModel scoreModel = excuteSearch(optProtDataset, ftupdatedName, finalFTOption, tempIdParam, false, searchInputSetting, identificationParametersFile, true);
            return scoreModel;
        });
        try {
            fixedModScore = fFuture.get();
            if (fixedModScore.isSignificatChange()) {
                optProtDataset.setValidatedIdRefrenceData(fixedModScore.getData());
            }
        } catch (ExecutionException | InterruptedException ex) {
            ex.printStackTrace();
        }
        MainUtilities.cleanOutputFolder();
        fixedModParamScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        fixedModParamScore.setParamValue(selectedFixedModificationOption.toString());
        parameterScoreSet.add(fixedModParamScore);

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
        for (String modId : potintialVariableMod) {  //mods
            if (selectedFixedModificationOption.contains(modId)) {
                continue;
            }
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));

            Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, false);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    sortingModificationMap.put(scoreModel, modId);
                    List<SpectrumMatch> spectraResults = scoreModel.getSpectrumMatchResult();
                    modifiedSpectrumMap.put(modId, SpectraFileUtilities.getModifiedSpectrumSet(spectraResults));
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

                    Set<String> intersectionSet = SpectraFileUtilities.getIntersectionSet(modifiedSpectrumMap.get(variableModI), modifiedSpectrumMap.get(variableModII));
                    Set<String> fullSet = new HashSet<>(oneDResultsMap.get(variableModI).getSpecTitles());
                    fullSet.addAll(oneDResultsMap.get(variableModII).getSpecTitles());
                    int effectiveSize = fullSet.size() - intersectionSet.size();
                    final String option = variableModI + "_" + variableModII;
                    if (effectiveSize > oneDResultsMap.get(variableModI).getSpecTitles().size() && effectiveSize > oneDResultsMap.get(variableModII).getSpecTitles().size()) {
                        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModI));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII));
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
                        Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, true);
                            return scoreModel;
                        });
                        try {
                            RawScoreModel vModScore = f.get();
                            if (vModScore.isSignificatChange()) {
                                sorterSet.add(vModScore);
                                sorterSet.add(oneDResultsMap.get(variableModI));
                                sorterSet.add(oneDResultsMap.get(variableModII));
                                if (sorterSet.first() == vModScore) {
                                    twoDResultsMap.put(option, vModScore);
                                    modifiedSpectrumMap.put(option, SpectraFileUtilities.getModifiedSpectrumSet(vModScore.getSpectrumMatchResult()));
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

                    Set<String> intersectionSet = SpectraFileUtilities.getIntersectionSet(modifiedSpectrumMap.get(variableModI), modifiedSpectrumMap.get(variableModII));
                    Set<String> fullSet = new HashSet<>(oneDResultsMap.get(variableModI).getSpecTitles());
                    fullSet.addAll(twoDResultsMap.get(variableModII).getSpecTitles());
                    int effectiveSize = fullSet.size() - intersectionSet.size();
                    final String option = variableModI + "_" + variableModII;

                    if (effectiveSize > oneDResultsMap.get(variableModI).getSpecTitles().size() && effectiveSize > twoDResultsMap.get(variableModII).getSpecTitles().size()) {
                        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModI));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[0]));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[1]));
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
                        Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, true);
                            return scoreModel;
                        });
                        try {
                            RawScoreModel vModScore = f.get();
                            if (vModScore.isSignificatChange()) {
                                sorterSet.add(vModScore);
                                sorterSet.add(oneDResultsMap.get(variableModI));
                                sorterSet.add(twoDResultsMap.get(variableModII));
                                if (sorterSet.first() == vModScore) {
                                    threeDResultsMap.put(option, vModScore);
                                    modifiedSpectrumMap.put(option, SpectraFileUtilities.getModifiedSpectrumSet(vModScore.getSpectrumMatchResult()));
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

                    Set<String> intersectionSet = SpectraFileUtilities.getIntersectionSet(modifiedSpectrumMap.get(variableModI), modifiedSpectrumMap.get(variableModII));
                    Set<String> fullSet = new HashSet<>(oneDResultsMap.get(variableModI).getSpecTitles());
                    fullSet.addAll(threeDResultsMap.get(variableModII).getSpecTitles());
                    int effectiveSize = fullSet.size() - intersectionSet.size();
                    final String option = variableModI + "_" + variableModII;

                    if (effectiveSize > oneDResultsMap.get(variableModI).getSpecTitles().size() && effectiveSize > threeDResultsMap.get(variableModII).getSpecTitles().size()) {
                        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModI));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[0]));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[1]));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[2]));
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
                        Future<RawScoreModel> f = MainUtilities.getExecutorService().submit(() -> {
                            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, searchInputSetting, identificationParametersFile, true);
                            return scoreModel;
                        });
                        try {
                            RawScoreModel vModScore = f.get();
                            if (vModScore.isSignificatChange()) {
                                sorterSet.add(vModScore);
                                sorterSet.add(oneDResultsMap.get(variableModI));
                                sorterSet.add(threeDResultsMap.get(variableModII));
                                if (sorterSet.first() == vModScore) {
                                    fourDResultsMap.put(option, vModScore);
                                    modifiedSpectrumMap.put(option, SpectraFileUtilities.getModifiedSpectrumSet(vModScore.getSpectrumMatchResult()));
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

        String vtOption = "";
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        for (String mod : selectedVariableModificationOption) {
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod));
            vtOption += mod + "_";
        }
        final String finalvOption = vtOption;
        final String fvupdatedName = Configurations.DEFAULT_RESULT_NAME + "fv_" + vtOption + "_" + msFileName;
        Future<RawScoreModel> vFuture = MainUtilities.getExecutorService().submit(() -> {
            RawScoreModel scoreModel = excuteSearch(optProtDataset, fvupdatedName, finalvOption, tempIdParam, false, searchInputSetting, identificationParametersFile, true);
            return scoreModel;
        });
        try {
            RawScoreModel variableScoreModel = vFuture.get();
            if (variableScoreModel.isSignificatChange()) {
                optProtDataset.setValidatedIdRefrenceData(variableScoreModel.getData());
            }
        } catch (ExecutionException | InterruptedException ex) {
            ex.printStackTrace();
        }

        variableModParamScore.setScore(optProtDataset.getValidatedIdRefrenceData().length);
        variableModParamScore.setParamValue(selectedVariableModificationOption.toString());
        parameterScoreSet.add(variableModParamScore);
        modificationsResults.put("variableModifications", new HashSet<>(selectedVariableModificationOption));

//        //test refine mod 
      
        MainUtilities.cleanOutputFolder();
        return modificationsResults;

    }

    public abstract RawScoreModel excuteSearch(SearchingSubDataset optProtDataset, String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addPeptideMasses, SearchInputSetting optProtSearchSettings, File identificationParametersFile, boolean pairData);

}
