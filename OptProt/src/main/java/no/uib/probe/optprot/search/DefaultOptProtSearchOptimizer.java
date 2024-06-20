/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
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
import no.uib.probe.optprot.model.StringSorter;
import no.uib.probe.optprot.util.MainUtilities;
import static no.uib.probe.optprot.util.MainUtilities.executor;
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

            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("CleavageParameter");

            Future<RawScoreModel> f = executor.submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(option, f.get());
                    paramScore.setScore(resultsMap.get(option).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
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
        return selectedOption;

    }

    public String optimizeEnzymeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {

        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        if (!oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name().equalsIgnoreCase("enzyme")) {
            return null;
        }
        String selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
        int nMissesCleavages = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(selectedOption);

//        int idRate = optProtDataset.getActiveIdentificationNum();
        Map<String, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        //optimise enzyme  
        Enzyme[] enzymes = new Enzyme[EnzymeFactory.getInstance().getEnzymes().size()];
        enzymes[0] = EnzymeFactory.getInstance().getEnzyme("Trypsin");
        int i = 1;
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
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("Enzyme");
            Future<RawScoreModel> f = executor.submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(option, scoreModel);
                    paramScore.setScore(scoreModel.getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
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

        return selectedOption;

    }

    public String optimizeSpecificityParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        if (!oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name().equalsIgnoreCase("enzyme")) {
            return null;
        }
        String enzymeName = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
        String selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getSpecificity(enzymeName).name();
//        int idRate = optProtDataset.getActiveIdentificationNum();
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
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("Specificity");

            Future<RawScoreModel> f = executor.submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(option, scoreModel);
                    paramScore.setScore(resultsMap.get(option).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
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
            System.out.println("updated id ref data " + selectedOption);
        }
        return selectedOption;

    }

    public String optimizeFragmentIonTypesParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
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
                final ParameterScoreModel paramScore = new ParameterScoreModel();
                paramScore.setParamId("FragmentIons");
                Future<RawScoreModel> f = executor.submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                    return scoreModel;
                });
                try {
                    RawScoreModel scoreModel = f.get();
                    if (scoreModel.isSignificatChange()) {
                        resultsMap.put(option, scoreModel);
                        paramScore.setScore(resultsMap.get(option).getTotalNumber());
                        paramScore.setParamValue(option);
                        parameterScoreSet.add(paramScore);
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
            System.out.println("updated id ref data " + selectedOption);
        }

        selectedOption = selectedOption.replace("[", "").replace("]", "");

        return selectedOption;
    }

    public Integer optimizeMaxMissCleavagesParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        if (!oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name().equalsIgnoreCase("enzyme")) {
            return -1;
        }
        String enzymeName = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
        Integer selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(enzymeName);
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());

        for (int i = 0; i < 8; i++) {
            if (i == selectedOption) {
                continue;
            }
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(enzymeName, i);
            final String option = "missedCleavages_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("missedCleavages");
            Future<RawScoreModel> f = executor.submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(i, scoreModel);
                    paramScore.setScore(resultsMap.get(i).getTotalNumber());
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
            selectedOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
            System.out.println("updated id ref data " + selectedOption);
        }

        return selectedOption;

    }

    public double optimizeFragmentToleranceParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
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
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("fragmentAccuracy");

            Future<RawScoreModel> f = executor.submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(i, scoreModel);
                    paramScore.setScore(resultsMap.get(i).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
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
            System.out.println("updated id ref data " + selectedOption);
        }
        return selectedOption;

    }

    public int[] optimizePrecursorChargeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxChargeOption = oreginaltempIdParam.getSearchParameters().getMaxChargeSearched();
        int selectedMinChargeOption = oreginaltempIdParam.getSearchParameters().getMinChargeSearched();
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());

        for (int i = 2; i < 6; i++) {
            if (i == selectedMaxChargeOption) {
                continue;
            }

            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setMaxChargeSearched(i);
            if (i <= selectedMinChargeOption) {
//                selectedMinChargeOption = i;
                tempIdParam.getSearchParameters().setMinChargeSearched(i - 1);
            } else {
                tempIdParam.getSearchParameters().setMinChargeSearched(selectedMinChargeOption);
            }
            final String option = "maxCharge_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("maxCharge");

            Future<RawScoreModel> f = executor.submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(i, scoreModel);
                    paramScore.setScore(resultsMap.get(i).getTotalNumber());
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
            selectedMaxChargeOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
            System.out.println("updated id ref data selectedMaxCharge" + selectedMaxChargeOption);
        }
        resultsMap.clear();

        for (int i = selectedMaxChargeOption; i > 0; i--) {
            if (i == selectedMinChargeOption) {
                continue;
            }
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setMaxChargeSearched(selectedMaxChargeOption);
            tempIdParam.getSearchParameters().setMinChargeSearched(i);
            final String option = "minCharge_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("minCharge");

            Future<RawScoreModel> f = executor.submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(i, scoreModel);
                    paramScore.setScore(resultsMap.get(i).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }
        sortedResultsMap.clear();
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedMinChargeOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
            System.out.println("updated id ref data MinCharge " + selectedMinChargeOption);
        }
        MainUtilities.cleanOutputFolder();
        return new int[]{selectedMinChargeOption, selectedMaxChargeOption};

    }

    public int[] optimizeIsotopParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxIsotopicCorrectionOption = oreginaltempIdParam.getSearchParameters().getMaxIsotopicCorrection();
        int selectedMinIsotopicCorrectioneOption = oreginaltempIdParam.getSearchParameters().getMinIsotopicCorrection();
        Map<Integer, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        for (int i = -8; i <= 8; i++) {
            if (i == selectedMinIsotopicCorrectioneOption) {
                continue;
            }
            if (i > selectedMaxIsotopicCorrectionOption) {
                tempIdParam.getSearchParameters().setMaxIsotopicCorrection(i + 1);
            } else {
                tempIdParam.getSearchParameters().setMaxIsotopicCorrection(selectedMaxIsotopicCorrectionOption);
            }
            tempIdParam.getSearchParameters().setMinIsotopicCorrection(i);
            final String option = "minIsotop_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("minIsotop");
            Future<RawScoreModel> f = executor.submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(i, scoreModel);
                    paramScore.setScore(resultsMap.get(i).getTotalNumber());
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
            selectedMinIsotopicCorrectioneOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
            System.out.println("updated id ref data minIstop " + selectedMinIsotopicCorrectioneOption);
        }

        resultsMap.clear();

        tempIdParam.getSearchParameters().setMinIsotopicCorrection(selectedMinIsotopicCorrectioneOption);
        for (int i = selectedMinIsotopicCorrectioneOption; i < 9; i++) {
            if (i == selectedMaxIsotopicCorrectionOption) {
                continue;
            }
            tempIdParam.getSearchParameters().setMaxIsotopicCorrection(i);
            final String option = "maxIsotop_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("maxIsotop");

            Future<RawScoreModel> f = executor.submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(i, scoreModel);
                    paramScore.setScore(resultsMap.get(i).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }
        sortedResultsMap.clear();
        if (!resultsMap.isEmpty()) {
            for (int option : resultsMap.keySet()) {
                sortedResultsMap.put(resultsMap.get(option), option);
            }
            selectedMaxIsotopicCorrectionOption = sortedResultsMap.firstEntry().getValue();
            optProtDataset.setValidatedIdRefrenceData(sortedResultsMap.firstEntry().getKey().getData());
            System.out.println("updated id ref data MaxIsotop" + selectedMaxIsotopicCorrectionOption);
        }

        return new int[]{selectedMinIsotopicCorrectioneOption, selectedMaxIsotopicCorrectionOption};
    }

    public double optimizePrecursorToleranceParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedOption = oreginaltempIdParam.getSearchParameters().getPrecursorAccuracy();
        Map<Double, RawScoreModel> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double[] iValues = new double[]{5, 10, 15, 20};
        for (double i : iValues) {
            if (i == selectedOption) {
                continue;
            }
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setPrecursorAccuracy(i);
            tempIdParam.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
            final String option = "precursorAccuracy_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            final ParameterScoreModel paramScore = new ParameterScoreModel();
            paramScore.setParamId("PrecursorAccuracy_PPM");
            System.out.println("at-------------------------------------------------- value to search " + i);
            Future<RawScoreModel> f = executor.submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option + "_noupdate", tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                return scoreModel;
            });
            try {

                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    resultsMap.put(i, scoreModel);
                    paramScore.setScore(resultsMap.get(i).getTotalNumber());
                    paramScore.setParamValue(option);
                    parameterScoreSet.add(paramScore);
                    paramScore.setComments("High-Resolution Mass Spectrometers: Instruments like Orbitrap or Fourier Transform Ion Cyclotron Resonance (FT-ICR)");
                    System.out.println("value " + i + " ppm " + scoreModel);

                }
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }

        iValues = new double[]{0.1, 0.3, 0.5, 0.7, 0.9};
        if (!optProtDataset.isHighResolutionMassSpectrometers()) {
            for (double i : iValues) {
                IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
                tempIdParam.getSearchParameters().setPrecursorAccuracy(i);
                tempIdParam.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.DA);
                final String option = "precursorAccuracy_Da" + i;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
                final ParameterScoreModel paramScore = new ParameterScoreModel();
                paramScore.setParamId("PrecursorAccuracy_Da");

                Future<RawScoreModel> f = executor.submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                    return scoreModel;
                });
                try {
                    RawScoreModel scoreModel = f.get();
                    if (scoreModel.isSignificatChange()) {
                        resultsMap.put(i, scoreModel);
                        paramScore.setScore(resultsMap.get(i).getTotalNumber());
                        paramScore.setParamValue(option);
                        parameterScoreSet.add(paramScore);
                        paramScore.setComments("Low-Resolution Mass Spectrometers: Quadrupole and ion trap mass spectrometers have lower mass accuracy");
                        System.out.println("value " + i + " kd " + scoreModel);
                    }
                } catch (ExecutionException | InterruptedException ex) {
                    ex.printStackTrace();
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

        return selectedOption;
    }

    public Map<String, Set<String>> optimizeModificationsParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter, TreeSet<ParameterScoreModel> parameterScoreSet) throws IOException {
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
        for (String modId : mods) {
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(modId));

            Future<RawScoreModel> f = executor.submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                return scoreModel;
            });
            try {
                resultsMap.put(modId, f.get());
            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }

        TreeMap<RawScoreModel, String> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());

        for (String modificationId : resultsMap.keySet()) {
            if (resultsMap.get(modificationId).isSignificatChange()) {
                sortedResultsMap.put(resultsMap.get(modificationId), modificationId);
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
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "fAsm_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));

            Future<RawScoreModel> f = executor.submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                return scoreModel;
            });
            try {
                RawScoreModel vScore = f.get();
                if (vScore.isSignificatChange()) {
                    boolean compscores = SpectraFileUtilities.compare2RawScoresData(vScore, resultsMap.get(modId));
                    if (compscores) {
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
                Future<RawScoreModel> f = executor.submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
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
                    Future<RawScoreModel> f = executor.submit(() -> {
                        RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
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

                    Future<RawScoreModel> f = executor.submit(() -> {
                        RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
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
        final ParameterScoreModel fixedModParamScore = new ParameterScoreModel();
        fixedModParamScore.setParamId("FixedModifications");
        final String finalFTOption = ftoption;
        Future<RawScoreModel> fFuture = executor.submit(() -> {
            RawScoreModel scoreModel = excuteSearch(optProtDataset, ftupdatedName, finalFTOption, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
            return scoreModel;
        });
        try {
            fixedModScore = fFuture.get();
            if (fixedModScore.isSignificatChange()) {
                fixedModParamScore.setScore(fixedModScore.getTotalNumber());
                fixedModParamScore.setParamValue(ftoption);
                parameterScoreSet.add(fixedModParamScore);
                optProtDataset.setValidatedIdRefrenceData(fixedModScore.getData());
            }
        } catch (ExecutionException | InterruptedException ex) {
            ex.printStackTrace();
        }
        MainUtilities.cleanOutputFolder();
        //process variable modifications
        resultsMap.clear();
        oneDResultsMap.clear();
        twoDResultsMap.clear();
        threeDResultsMap.clear();
        fourDResultsMap.clear();
        sorterSet.clear();
        TreeMap<RawScoreModel, String> sortingModificationMap = new TreeMap<>(Collections.reverseOrder());
        Map<String, Set<String>> modifiedSpectrumMap = new LinkedHashMap<>();
        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
        Map<String, Integer> uniqModCompI = new LinkedHashMap<>();
        for (String modId : mods) {
            if (selectedFixedModificationOption.contains(modId)) {
                continue;
            }
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));
            Future<RawScoreModel> f = executor.submit(() -> {
                RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, optimisedSearchParameter, identificationParametersFile);
                return scoreModel;
            });
            try {
                RawScoreModel scoreModel = f.get();
                if (scoreModel.isSignificatChange()) {
                    sortingModificationMap.put(scoreModel, modId);
                    List<SpectrumMatch> spectraResults = scoreModel.getSpectrumMatchResult();
                    modifiedSpectrumMap.put(modId, SpectraFileUtilities.getModifiedSpectrumSet(ptmFactory.getModification(modId), spectraResults));

                }

            } catch (ExecutionException | InterruptedException ex) {
                ex.printStackTrace();
            }

        }
        int counter = 0;
        for (RawScoreModel model : sortingModificationMap.keySet()) {
            String modId = sortingModificationMap.get(model);
            uniqModCompI.put(modId, modifiedSpectrumMap.get(modId).size());
            if (counter >= 3) {
                break;
            }
            counter++;
        }

        sorterSet.clear();
        Map<String, Integer> uniqModCompII = new LinkedHashMap<>();
//        if (uniqModCompI.size() > 1) {
            for (String variableModI : uniqModCompI.keySet()) {
                for (String variableModII : uniqModCompI.keySet()) {
                    StringSorter nameSortingTree = new StringSorter();
                    nameSortingTree.add(variableModI);
                    nameSortingTree.add(variableModII);
                    if (variableModII.equalsIgnoreCase(variableModI)|| uniqModCompII.containsKey(nameSortingTree.toString())) {                       
                        continue;
                    }
                    Set<String> filterset = new HashSet<>(modifiedSpectrumMap.get(variableModI));
                    filterset.addAll(modifiedSpectrumMap.get(variableModII));
                    if (filterset.size() > modifiedSpectrumMap.get(variableModI).size() && filterset.size() > modifiedSpectrumMap.get(variableModII).size()) {
                        uniqModCompII.put(nameSortingTree.toString(), filterset.size());
                    }
                }
            }
//        }

        Map<String, Integer> uniqModCompIII = new LinkedHashMap<>();
//        if (modifiedSpectrumMap.size() > 2) {
            for (String variableModI : uniqModCompI.keySet()) {
                for (String variableModII : uniqModCompII.keySet()) {
                    StringSorter nameSortingTree = new StringSorter();
                    nameSortingTree.add(variableModI);
                    nameSortingTree.add(variableModII.split("_")[0]);
                    nameSortingTree.add(variableModII.split("_")[1]);
                    if (variableModII.contains(variableModI) || uniqModCompIII.containsKey(nameSortingTree.toString())) {
                        continue;
                    }
                    Set<String> filterset = new HashSet<>(modifiedSpectrumMap.get(variableModI));
                    filterset.addAll(modifiedSpectrumMap.get(variableModII.split("_")[0]));
                    filterset.addAll(modifiedSpectrumMap.get(variableModII.split("_")[1]));
                    if (filterset.size() > modifiedSpectrumMap.get(variableModI).size() && filterset.size() > uniqModCompII.get(variableModII)) {
                        uniqModCompIII.put(nameSortingTree.toString(), filterset.size());
                    }
                }
            }

//        }

        Map<String, Integer> uniqModCompIV = new LinkedHashMap<>();
//        if (modifiedSpectrumMap.size() > 3) {
            for (String variableModI : uniqModCompI.keySet()) {
                for (String variableModII : uniqModCompIII.keySet()) {
                    StringSorter nameSortingTree = new StringSorter();
                    nameSortingTree.add(variableModI);
                    nameSortingTree.add(variableModII.split("_")[0]);
                    nameSortingTree.add(variableModII.split("_")[1]);
                    nameSortingTree.add(variableModII.split("_")[2]);
                    if (variableModII.contains(variableModI) || uniqModCompIV.containsKey(nameSortingTree.toString())) {
                        continue;
                    }
                    Set<String> filterset = new HashSet<>(modifiedSpectrumMap.get(variableModI));
                    filterset.addAll(modifiedSpectrumMap.get(variableModII.split("_")[0]));
                    filterset.addAll(modifiedSpectrumMap.get(variableModII.split("_")[1]));
                    filterset.addAll(modifiedSpectrumMap.get(variableModII.split("_")[2]));
                    if (filterset.size() > modifiedSpectrumMap.get(variableModI).size() && filterset.size() > uniqModCompIII.get(variableModII)) {
                        uniqModCompIV.put(nameSortingTree.toString(), filterset.size());
                    }
                }
            }

//        }

        uniqModCompI.putAll(uniqModCompII);
        uniqModCompI.putAll(uniqModCompIII);
        uniqModCompI.putAll(uniqModCompIV);
        TreeMap<Integer, Set<String>> sortingVarModMap = new TreeMap<>(Collections.reverseOrder());
        for (String key : uniqModCompI.keySet()) {
            int count = uniqModCompI.get(key);
            if (!sortingVarModMap.containsKey(count)) {
                sortingVarModMap.put(count, new LinkedHashSet<>());
            }
            sortingVarModMap.get(count).add(key);
        }

        if (!sortingVarModMap.isEmpty()) {

            for (String modsm : sortingVarModMap.firstEntry().getValue()) {
                for (RawScoreModel model : sortingModificationMap.keySet()) {
                    for (String mod : modsm.split("_")) {
                        if (sortingModificationMap.get(model).equalsIgnoreCase(mod)) {
                            oneDResultsMap.put(mod, model);
                        }
                    }

                }

                break;

            }
            double startingrange = sortingVarModMap.firstEntry().getValue().iterator().next().split("_").length;
            for (int score : sortingVarModMap.keySet()) {

                for (String values : sortingVarModMap.get(score)) {
                    if (values.split("_").length == startingrange) {
                        System.out.println("# mod: " + startingrange + " score: " + score + "  values " + values + "   ");
                        startingrange--;
                    }
                    if (startingrange == 0) {
                        break;
                    }
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

                    Set<String> filterset = new HashSet<>(modifiedSpectrumMap.get(variableModI));
                    filterset.addAll(modifiedSpectrumMap.get(variableModII));
                    final String option = variableModI + "_" + variableModII;
                    if (filterset.size() > modifiedSpectrumMap.get(variableModI).size() && filterset.size() > modifiedSpectrumMap.get(variableModII).size()) {
                        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModI));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII));
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
                        Future<RawScoreModel> f = executor.submit(() -> {
                            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, optimisedSearchParameter, identificationParametersFile);
                            return scoreModel;
                        });
                        try {
                            RawScoreModel vModScore = f.get();
                            if (vModScore.isSignificatChange()) {
//                                change sorter to other way
                                sorterSet.add(vModScore);
                                sorterSet.add(oneDResultsMap.get(variableModI));
                                sorterSet.add(oneDResultsMap.get(variableModII));
                                if (sorterSet.first() == vModScore) {
                                    twoDResultsMap.put(option, vModScore);
                                    modifiedSpectrumMap.put(option, SpectraFileUtilities.getModifiedSpectrumSet(ptmFactory.getModification(variableModI), vModScore.getSpectrumMatchResult()));
                                    modifiedSpectrumMap.get(option).addAll(SpectraFileUtilities.getModifiedSpectrumSet(ptmFactory.getModification(variableModII), vModScore.getSpectrumMatchResult()));
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

                    Set<String> filterset = new HashSet<>(modifiedSpectrumMap.get(variableModI));
                    filterset.addAll(modifiedSpectrumMap.get(variableModII));
                    final String option = variableModI + "_" + variableModII;

                    if (filterset.size() > modifiedSpectrumMap.get(variableModI).size() && filterset.size() > modifiedSpectrumMap.get(variableModII).size()) {
                        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModI));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[0]));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[1]));
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
                        Future<RawScoreModel> f = executor.submit(() -> {
                            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, optimisedSearchParameter, identificationParametersFile);
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
                                    modifiedSpectrumMap.put(option, SpectraFileUtilities.getModifiedSpectrumSet(ptmFactory.getModification(variableModI), vModScore.getSpectrumMatchResult()));
                                    modifiedSpectrumMap.get(option).addAll(SpectraFileUtilities.getModifiedSpectrumSet(ptmFactory.getModification(variableModII.split("_")[0]), vModScore.getSpectrumMatchResult()));
                                    modifiedSpectrumMap.get(option).addAll(SpectraFileUtilities.getModifiedSpectrumSet(ptmFactory.getModification(variableModII.split("_")[1]), vModScore.getSpectrumMatchResult()));
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

                    Set<String> filterset = new HashSet<>(modifiedSpectrumMap.get(variableModI));
                    filterset.addAll(modifiedSpectrumMap.get(variableModII));
                    final String option = variableModI + "_" + variableModII;

                    if (filterset.size() > modifiedSpectrumMap.get(variableModI).size() && filterset.size() > modifiedSpectrumMap.get(variableModII).size()) {
                        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModI));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[0]));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[1]));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII.split("_")[2]));
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
                        Future<RawScoreModel> f = executor.submit(() -> {
                            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, true, optimisedSearchParameter, identificationParametersFile);
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
                                    modifiedSpectrumMap.put(option, SpectraFileUtilities.getModifiedSpectrumSet(ptmFactory.getModification(variableModI), vModScore.getSpectrumMatchResult()));
                                    modifiedSpectrumMap.get(option).addAll(SpectraFileUtilities.getModifiedSpectrumSet(ptmFactory.getModification(variableModII.split("_")[0]), vModScore.getSpectrumMatchResult()));
                                    modifiedSpectrumMap.get(option).addAll(SpectraFileUtilities.getModifiedSpectrumSet(ptmFactory.getModification(variableModII.split("_")[1]), vModScore.getSpectrumMatchResult()));
                                    modifiedSpectrumMap.get(option).addAll(SpectraFileUtilities.getModifiedSpectrumSet(ptmFactory.getModification(variableModII.split("_")[2]), vModScore.getSpectrumMatchResult()));
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
        final ParameterScoreModel variable = new ParameterScoreModel();
        variable.setParamId("VariableModifications");
        String vtOption = "";
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        for (String mod : selectedVariableModificationOption) {
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod));
            vtOption += mod + "_";
        }
        final String finalvOption = vtOption;
        final String fvupdatedName = Configurations.DEFAULT_RESULT_NAME + "fv_" + vtOption + "_" + msFileName;
        Future<RawScoreModel> vFuture = executor.submit(() -> {
            RawScoreModel scoreModel = excuteSearch(optProtDataset, fvupdatedName, finalvOption, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
            return scoreModel;
        });
        try {
            RawScoreModel variableScoreModel = vFuture.get();
            if (variableScoreModel.isSignificatChange()) {
                variable.setScore(variableScoreModel.getTotalNumber());
                variable.setParamValue(ftoption);
                parameterScoreSet.add(variable);
                optProtDataset.setValidatedIdRefrenceData(variableScoreModel.getData());
            }
        } catch (ExecutionException | InterruptedException ex) {
            ex.printStackTrace();
        }
        modificationsResults.put("variableModifications", new HashSet<>(selectedVariableModificationOption));

//        //test refine mod 
        MainUtilities.cleanOutputFolder();

        if (optimisedSearchParameter.getSelectedSearchEngine().getIndex() == Advocate.xtandem.getIndex()) {
            resultsMap.clear();
            oneDResultsMap.clear();
            twoDResultsMap.clear();
            threeDResultsMap.clear();
            fourDResultsMap.clear();
            sorterSet.clear();

//           
            for (String vMod : resultsMap.keySet()) {
                if (ptmFactory.getModification(vMod) == null) {
                    continue;
                }
                tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod));

                final String option = vMod;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;
                final ParameterScoreModel paramScore = new ParameterScoreModel();
                paramScore.setParamId("refineVariableModifications");

                Future<RawScoreModel> f = executor.submit(() -> {
                    RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                    return scoreModel;
                });
                try {
                    RawScoreModel scoreModel = f.get();
                    if (scoreModel.isSignificatChange()) {
                        oneDResultsMap.put(option, scoreModel);
                        paramScore.setScore(oneDResultsMap.get(option).getTotalNumber());
                        paramScore.setParamValue(option);
                        parameterScoreSet.add(paramScore);
                    }
                } catch (ExecutionException | InterruptedException ex) {
                    ex.printStackTrace();
                }

            }
            //2d refinment
            resultsMap.clear();
            if (oneDResultsMap.size() > 1) {
                for (String selectedRef : oneDResultsMap.keySet()) {
                    for (String vMod : oneDResultsMap.keySet()) {
                        if (vMod.equals(selectedRef)) {
                            continue;
                        }
                        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(selectedRef));
                        final String option = vMod + "_" + selectedRef;
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;
                        final ParameterScoreModel paramScore = new ParameterScoreModel();
                        paramScore.setParamId("refinmentVariableModifications");
                        Future<RawScoreModel> f = executor.submit(() -> {
                            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                            return scoreModel;
                        });
                        try {
                            RawScoreModel scoreModel = f.get();
                            if (scoreModel.isSignificatChange()) {
                                paramScore.setScore(scoreModel.getTotalNumber());
                                paramScore.setParamValue(option);
                                parameterScoreSet.add(paramScore);
                                sorterSet.add(scoreModel);
                                sorterSet.add(oneDResultsMap.get(selectedRef));
                                sorterSet.add(oneDResultsMap.get(vMod));
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
            if (!twoDResultsMap.isEmpty() && oneDResultsMap.size() > 2) {
                for (String selectedRef : oneDResultsMap.keySet()) {
                    for (String vMod : twoDResultsMap.keySet()) {
                        if (vMod.contains(selectedRef)) {
                            continue;
                        }
                        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[0]));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[1]));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(selectedRef));
                        final String option = vMod + "_" + selectedRef;
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;
                        final ParameterScoreModel paramScore = new ParameterScoreModel();
                        paramScore.setParamId("refinmentVariableModifications");
                        Future<RawScoreModel> f = executor.submit(() -> {
                            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                            return scoreModel;
                        });
                        try {
                            RawScoreModel scoreModel = f.get();
                            if (scoreModel.isSignificatChange()) {
                                paramScore.setScore(scoreModel.getTotalNumber());
                                paramScore.setParamValue(option);
                                parameterScoreSet.add(paramScore);
                                sorterSet.add(scoreModel);
                                sorterSet.add(oneDResultsMap.get(selectedRef));
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
            if (!threeDResultsMap.isEmpty() && oneDResultsMap.size() > 3) {
                for (String selectedRef : oneDResultsMap.keySet()) {
                    for (String vMod : threeDResultsMap.keySet()) {
                        if (vMod.equals(selectedRef)) {
                            continue;
                        }
                        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[0]));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[1]));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[2]));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(selectedRef));
                        final String option = vMod + "_" + selectedRef;
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;
                        final ParameterScoreModel paramScore = new ParameterScoreModel();
                        paramScore.setParamId("refinmentVariableModifications");
                        Future<RawScoreModel> f = executor.submit(() -> {
                            RawScoreModel scoreModel = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile);
                            return scoreModel;
                        });
                        try {
                            RawScoreModel scoreModel = f.get();
                            if (scoreModel.isSignificatChange()) {
                                paramScore.setScore(scoreModel.getTotalNumber());
                                paramScore.setParamValue(option);
                                parameterScoreSet.add(paramScore);
                                sorterSet.add(scoreModel);
                                sorterSet.add(oneDResultsMap.get(selectedRef));
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
            resultsMap.putAll(oneDResultsMap);
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
            modificationsResults.put("refinmentVariableModifications", refinementVarModMap);
        }
        MainUtilities.cleanOutputFolder();
        System.out.println("---->>> final  mod list ----- " + modificationsResults);

        return modificationsResults;

    }

    public abstract RawScoreModel excuteSearch(SearchingSubDataset optProtDataset, String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addPeptideMasses, SearchInputSetting optProtSearchSettings, File identificationParametersFile);

}
