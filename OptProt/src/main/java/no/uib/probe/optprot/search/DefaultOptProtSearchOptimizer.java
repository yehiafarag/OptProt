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
import com.compomics.util.experiment.identification.matches.ModificationMatch;
import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import com.compomics.util.experiment.identification.spectrum_assumptions.PeptideAssumption;
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
import java.util.concurrent.Future;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;
import no.uib.probe.optprot.model.ParameterScoreModel;
import no.uib.probe.optprot.model.SearchInputSetting;
import no.uib.probe.optprot.util.MainUtilities;
import static no.uib.probe.optprot.util.MainUtilities.executor;
import org.apache.commons.math3.util.Precision;

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

    public String optimizeDigestionParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name();
        int idRate = optProtDataset.getActiveIdentificationNum();
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        String[] cleavageParameters = new String[]{"wholeProtein", "unSpecific"};
        resultsMap.put(selectedOption, idRate);
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

            Future f = executor.submit(() -> {
                resultsMap.put(option, excuteSearch(optProtDataset, updatedName, "DigestionParameter", option + "_noupdate", tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!f.isDone()) {
            }

        }
//       
        System.out.println("at selected option " + selectedOption + "  " + idRate);
        int localId = -1;
        String localSelection = "";
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;

            }
        }

        double factor = (localId - idRate) * 100.0 / idRate;//* optProtDataset.getDataEpsilon();
        factor = Math.round(factor);
        if (factor >= optProtDataset.getAcceptedIDRatioThreshold() && resultsMap.get(localSelection) > resultsMap.get(selectedOption)) {
            selectedOption = localSelection;
            optProtDataset.setTempIdentificationNum(localId);
        }
        return selectedOption;

    }

    public String optimizeEnzymeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {

        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        if (!oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name().equalsIgnoreCase("enzyme")) {
            return null;
        }
        String selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
        int nMissesCleavages = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(selectedOption);

        int idRate = optProtDataset.getActiveIdentificationNum();
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
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
            if (enzyme.getName().equalsIgnoreCase(selectedOption) && !selectedOption.equalsIgnoreCase("Trypsin")) {
                continue;
            }
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
            tempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(enzyme);
            tempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(enzyme.getName(), nMissesCleavages);
            final String option = enzyme.getName();
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

            Future f = executor.submit(() -> {
                resultsMap.put(option, excuteSearch(optProtDataset, updatedName, "EnzymeParameter", option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                while (resultsMap.get(option) == null) {

                }

            });
            while (!f.isDone()) {
            }
            if (enzyme.getName().equals("Trypsin") && resultsMap.get(option) >= (optProtDataset.getActiveIdentificationNum())) {
                break;
            }

        }
        int localId = -1;
        String localSelection = "";
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;

            }
        }
        double factor = (localId - idRate) * 100.0 / idRate;//* optProtDataset.getDataEpsilon();
        factor = Math.round(factor);
        if (factor >= optProtDataset.getAcceptedIDRatioThreshold() && resultsMap.get(localSelection) > resultsMap.get(selectedOption)) {
            selectedOption = localSelection;
            optProtDataset.setActiveIdentificationNum(localId);
        }
        if (localSelection.equalsIgnoreCase(selectedOption)) {
            optProtDataset.setActiveIdentificationNum(Math.max(localId, idRate));
        }
        return selectedOption;

    }

    public String optimizeSpecificityParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        if (!oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name().equalsIgnoreCase("enzyme")) {
            return null;
        }
        String enzymeName = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
        String selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getSpecificity(enzymeName).name();
        int idRate = optProtDataset.getActiveIdentificationNum();
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());

        for (int i = 0; i < DigestionParameters.Specificity.values().length; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(enzymeName, DigestionParameters.Specificity.getSpecificity(i));
            final String option = DigestionParameters.Specificity.getSpecificity(i).name();
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

            Future f = executor.submit(() -> {
                resultsMap.put(option, excuteSearch(optProtDataset, updatedName, "SpecificityParameter", option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!f.isDone()) {
            }
            break;
        }

        int localId = -1;
        String localSelection = "";
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;

            }
        }
        double factor = (localId - idRate) * 100.0 / idRate;//* optProtDataset.getDataEpsilon();
        factor = Math.round(factor);
        if (factor >= optProtDataset.getAcceptedIDRatioThreshold() && resultsMap.get(localSelection) > resultsMap.get(selectedOption)) {
            selectedOption = localSelection;
            optProtDataset.setActiveIdentificationNum(localId);
        }
        if (localSelection.equalsIgnoreCase(selectedOption)) {
            optProtDataset.setActiveIdentificationNum(Math.max(localId, idRate));
        }

        return selectedOption;

    }

    public Integer optimizeMaxMissCleavagesParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        if (!oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name().equalsIgnoreCase("enzyme")) {
            return -1;
        }
        String enzymeName = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
        Integer selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(enzymeName);
        int idRate = optProtDataset.getActiveIdentificationNum();
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());

        for (int i = 0; i < 8; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(enzymeName, i);
            final String option = "missedCleavages_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, "MaxMissCleavagesParameter", option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
        }
        int localId = -1;
        int localSelection = 0;
        for (int option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;

            }
        }
        double factor = (localId - idRate) * 100.0 / idRate;
        if (factor >= optProtDataset.getAcceptedIDRatioThreshold() && resultsMap.get(localSelection) > resultsMap.get(selectedOption)) {
            selectedOption = localSelection;
            optProtDataset.setActiveIdentificationNum(localId);
        }
        if (localSelection == (selectedOption)) {
            optProtDataset.setActiveIdentificationNum(Math.max(localId, idRate));
        }
        return selectedOption;

    }

    public String optimizeFragmentIonTypesParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        String selectedOption = oreginaltempIdParam.getSearchParameters().getForwardIons() + "-" + oreginaltempIdParam.getSearchParameters().getRewindIons();
        int idRate = optProtDataset.getActiveIdentificationNum();
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
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;

                Future f = executor.submit(() -> {
                    resultsMap.put(option, excuteSearch(optProtDataset, updatedName, "FragmentIonTypesParameter", option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                    MainUtilities.cleanOutputFolder();
                });
                while (!f.isDone()) {
                }

            }
        }
        int localId = -1;
        String localSelection = "";
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;

            }
        }
        double factor = (localId - idRate) * 100.0 / idRate;
        factor = Math.round(factor);
        if (factor >= optProtDataset.getAcceptedIDRatioThreshold() && resultsMap.get(localSelection) > resultsMap.get(selectedOption)) {
            selectedOption = localSelection;
            optProtDataset.setActiveIdentificationNum(localId);
        }
        if (localSelection.equalsIgnoreCase(selectedOption)) {
            optProtDataset.setActiveIdentificationNum(Math.max(localId, idRate));
        }
        selectedOption = selectedOption.replace("[", "").replace("]", "");
        return selectedOption;
    }

    public double optimizePrecursorToleranceParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedOption = oreginaltempIdParam.getSearchParameters().getPrecursorAccuracy();
        int idRate = optProtDataset.getActiveIdentificationNum();
        Map<Double, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double[] iValues = new double[]{5, 10, 15, 20};
        for (double i : iValues) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setPrecursorAccuracy(i);
            tempIdParam.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
            final String option = "precursorAccuracy_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            final double j = i;

            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, "PrecursorToleranceParameter", option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
//                paramScore.setComments("High-Resolution Mass Spectrometers: Instruments like Orbitrap or Fourier Transform Ion Cyclotron Resonance (FT-ICR)");
            });
            while (!future.isDone()) {
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
                final double j = i;

                Future future = executor.submit(() -> {
                    resultsMap.put(j, excuteSearch(optProtDataset, updatedName, "PrecursorToleranceParameter", option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                });
                while (!future.isDone()) {
                }
//                paramScore.setComments("Low-Resolution Mass Spectrometers: Quadrupole and ion trap mass spectrometers have lower mass accuracy");
            }
        }

        int localId = -1;
        double localSelection = 0;
        for (double option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;

            }
        }
        double factor = (localId - idRate) * 100.0 / idRate;//* optProtDataset.getDataEpsilon();;
        factor = Math.round(factor);
        if (factor >= optProtDataset.getAcceptedIDRatioThreshold()) {
            selectedOption = localSelection;
            optProtDataset.setActiveIdentificationNum(localId);
        }
        return selectedOption;
    }

    public double optimizeFragmentToleranceParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedOption = oreginaltempIdParam.getSearchParameters().getFragmentIonAccuracy();
        int idRate = optProtDataset.getActiveIdentificationNum();
        Map<Double, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        double[] values = new double[]{0.02, 0.05, 0.07, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5};
        for (double i : values) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setFragmentIonAccuracy(i);
            tempIdParam.getSearchParameters().setFragmentAccuracyType(SearchParameters.MassAccuracyType.DA);
            final String option = "fragmentAccuracy_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            double j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, "FragmentToleranceParameter", option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
        }
        int localId = -1;
        double localSelection = 0;
        for (double option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;
            }
        }
        double factor = (localId - idRate) * 100.0 / idRate;
        factor = Math.round(factor);
        if (factor >= optProtDataset.getAcceptedIDRatioThreshold() && resultsMap.get(localSelection) > resultsMap.get(selectedOption)) {
            selectedOption = localSelection;
            optProtDataset.setActiveIdentificationNum(localId);
        }
        if (localSelection == (selectedOption)) {
            optProtDataset.setActiveIdentificationNum(Math.max(localId, idRate));
        }
        return selectedOption;

    }

    public int[] optimizePrecursorChargeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxChargeOption = oreginaltempIdParam.getSearchParameters().getMaxChargeSearched();
        int selectedMinChargeOption = oreginaltempIdParam.getSearchParameters().getMinChargeSearched();
        int idRate = optProtDataset.getActiveIdentificationNum();
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());

        for (int i = 2; i < 6; i++) {

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
            int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, "PrecursorChargeParameter", option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());

            });
            while (!future.isDone()) {
            }

        }
        int localId = -1;
        int localSelection = 0;
        for (int option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;

            }
        }
        double factor = (localId - idRate) * 100.0 / idRate;//* optProtDataset.getDataEpsilon();;
        factor = Math.round(factor);
        if (factor >= optProtDataset.getAcceptedIDRatioThreshold() && resultsMap.get(localSelection) > resultsMap.get(selectedMaxChargeOption)) {
            selectedMaxChargeOption = localSelection;
            optProtDataset.setActiveIdentificationNum(localId);

        }
        if (localSelection == (selectedMaxChargeOption)) {
            optProtDataset.setActiveIdentificationNum(Math.max(localId, idRate));
        }

        resultsMap.clear();

        for (int i = selectedMaxChargeOption; i > 0; i--) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setMaxChargeSearched(selectedMaxChargeOption);
            tempIdParam.getSearchParameters().setMinChargeSearched(i);
            final String option = "minCharge_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;

            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, "PrecursorChargeParameter", option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
        }
        localId = -1;
        localSelection = 0;
        for (int option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;
            }
        }
        factor = (localId - optProtDataset.getActiveIdentificationNum()) * 100.0 / idRate;//* optProtDataset.getDataEpsilon();;
        factor = Math.round(factor);
        if (factor >= optProtDataset.getAcceptedIDRatioThreshold() && resultsMap.get(localSelection) > resultsMap.get(selectedMinChargeOption)) {
            selectedMinChargeOption = localSelection;
            optProtDataset.setActiveIdentificationNum(localId);
        }
        if (localSelection == (selectedMinChargeOption)) {
            optProtDataset.setActiveIdentificationNum(Math.max(localId, idRate));
        }
        MainUtilities.cleanOutputFolder();
        return new int[]{selectedMinChargeOption, selectedMaxChargeOption};

    }

    public int[] optimizeIsotopParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxIsotopicCorrectionOption = oreginaltempIdParam.getSearchParameters().getMaxIsotopicCorrection();
        int selectedMinIsotopicCorrectioneOption = oreginaltempIdParam.getSearchParameters().getMinIsotopicCorrection();
        int idRate = optProtDataset.getActiveIdentificationNum();
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        for (int i = -8; i <= 8; i++) {
            if (i > selectedMaxIsotopicCorrectionOption) {
                tempIdParam.getSearchParameters().setMaxIsotopicCorrection(i + 1);
            } else {
                tempIdParam.getSearchParameters().setMaxIsotopicCorrection(selectedMaxIsotopicCorrectionOption);
            }
            tempIdParam.getSearchParameters().setMinIsotopicCorrection(i);
            final String option = "minIsotop_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, "IsotopParameter", option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }

        }

        int localId = -1;
        int localSelection = 0;
        for (int option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;

            }
        }
        double factor = (localId - idRate) * 100.0 / idRate;//* optProtDataset.getDataEpsilon();;
        factor = Math.round(factor);
        if (factor >= optProtDataset.getAcceptedIDRatioThreshold() && resultsMap.get(localSelection) > resultsMap.get(selectedMinIsotopicCorrectioneOption)) {
            selectedMinIsotopicCorrectioneOption = localSelection;
            optProtDataset.setActiveIdentificationNum(localId);
            idRate = localId;
        }
        resultsMap.clear();
        tempIdParam.getSearchParameters().setMinIsotopicCorrection(selectedMinIsotopicCorrectioneOption);
        for (int i = selectedMinIsotopicCorrectioneOption; i < 9; i++) {
            tempIdParam.getSearchParameters().setMaxIsotopicCorrection(i);
            final String option = "maxIsotop_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;

            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, "IsotopParameter", option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
        }
        localId = -1;
        localSelection = 0;
        for (int option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;

            }
        }
        factor = (localId - idRate) * 100.0 / idRate;//* optProtDataset.getDataEpsilon();;
        factor = Math.round(factor);
        if (factor >= optProtDataset.getAcceptedIDRatioThreshold() && resultsMap.get(localSelection) > resultsMap.get(selectedMaxIsotopicCorrectionOption)) {
            selectedMaxIsotopicCorrectionOption = localSelection;
            optProtDataset.setActiveIdentificationNum(localId);
        }
        return new int[]{selectedMinIsotopicCorrectioneOption, selectedMaxIsotopicCorrectionOption};
    }

    public Map<String, Set<String>> optimizeModificationsParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        ArrayList<String> selectedFixedModificationOption = new ArrayList<>(oreginaltempIdParam.getSearchParameters().getModificationParameters().getFixedModifications());
        ArrayList<String> selectedVariableModificationOption = new ArrayList<>();
        int idRate = optProtDataset.getActiveIdentificationNum();
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
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        for (String modId : mods) {
            final String option =  modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(modId));
            Future future = executor.submit(() -> {
                resultsMap.put(modId, excuteSearch(optProtDataset, updatedName, "ModificationsParameter", option + "_noupdate"+"_fixed", tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }

        }

        TreeMap<Integer, Set<String>> sortedResultsMap = new TreeMap<>(Collections.reverseOrder());

        for (String modificationId : resultsMap.keySet()) {

            if (!sortedResultsMap.containsKey(resultsMap.get(modificationId))) {
                sortedResultsMap.put(resultsMap.get(modificationId), new HashSet<>());
            }
            sortedResultsMap.get(resultsMap.get(modificationId)).add(modificationId);
        }

        Map<String, Map<String, Integer>> targtedFixedModificationScore = new TreeMap<>();
        for (int score : sortedResultsMap.keySet()) {
            for (String modificationId : sortedResultsMap.get(score)) {
                String modPattern = ptmFactory.getModification(modificationId).getPattern().toString();
                if (!targtedFixedModificationScore.containsKey(modPattern)) {
                    targtedFixedModificationScore.put(modPattern, new TreeMap<>());
                    targtedFixedModificationScore.get(modPattern).put(modificationId, resultsMap.get(modificationId));
                    System.out.println("at mod is " + modificationId + "  target " + modPattern + "  score " + score);
                }

            }
        }

        sortedResultsMap.clear();
        int modLimit = 0;
        for (String key : targtedFixedModificationScore.keySet()) {
            for (String modId : targtedFixedModificationScore.get(key).keySet()) {
                int score = targtedFixedModificationScore.get(key).get(modId);
                if (score < idRate) {
                    continue;
                }
                if (!sortedResultsMap.containsKey(score)) {
                    sortedResultsMap.put(score, new LinkedHashSet<>());
                }
                sortedResultsMap.get(score).add(modId);
                modLimit++;
                if (modLimit >= 5) {
                    break;
                }
            }
            if (modLimit >= 5) {
                break;
            }
        }
        resultsMap.clear();//
        for (int score : sortedResultsMap.keySet()) {
            for (String mod : sortedResultsMap.get(score)) {
                resultsMap.put(mod, score);
            }
        }

        Map<String, Integer> oneDResultsMap = new TreeMap<>(resultsMap);

        Set<String> clearSet = new LinkedHashSet<>(oneDResultsMap.keySet());
        for (String modId : clearSet) {
            final String option =  modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "m_" + option + "_" + msFileName;

            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));

            Future future = executor.submit(() -> {
                int vmodScore = excuteSearch(optProtDataset, updatedName, "ModificationsParameter", option + "_noupdate"+"_variable", tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                double factor = (vmodScore - oneDResultsMap.get(modId)) * 100.0 / oneDResultsMap.get(modId);//* optProtDataset.getDataEpsilon();;
                factor = Math.round(factor);
                if (factor >= optProtDataset.getAcceptedIDRatioThreshold()) {
                    oneDResultsMap.remove(modId);
                }

            });
            while (!future.isDone()) {
            }
        }

        Map<String, Integer> twoDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        int indexI = 0;
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
                Future future = executor.submit(() -> {
                    int score = excuteSearch(optProtDataset, updatedName, "ModificationsParameter", option + "_noupdate"+"_fixed_Mods", tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                    double factor = (score - idRate) * 100.0 / idRate;//* optProtDataset.getDataEpsilon();;
                    factor = Math.round(factor);
                    if (score > oneDResultsMap.get(mod1Id) && score > oneDResultsMap.get(mod2Id) && factor >= optProtDataset.getAcceptedIDRatioThreshold()) {
                        twoDResultsMap.put(option, score);
                    }
                });
                while (!future.isDone()) {
                }

                indexII++;
            }
            indexI++;
        }
        //filter 2 d results

        Map<String, Integer> threeDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        if (oneDResultsMap.size() > 2) {
            for (String mod1Id : oneDResultsMap.keySet()) {
                for (String mod2Id : twoDResultsMap.keySet()) {
                    if (mod2Id.contains(mod1Id)) {
                        continue;
                    }
                    tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod1Id));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[0]));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[1]));
                    final String option =  mod1Id + "_" + mod2Id;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;

                    Future future = executor.submit(() -> {
                        int score = excuteSearch(optProtDataset, updatedName, "ModificationsParameter", option + "_noupdate"+"_fixed_Mods", tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                        double factor = (score - idRate) * 100.0 / idRate;//* optProtDataset.getDataEpsilon();;
                        factor = Math.round(factor);
                        if (score > oneDResultsMap.get(mod1Id) && score > twoDResultsMap.get(mod2Id) && factor >= optProtDataset.getAcceptedIDRatioThreshold()) {
                            threeDResultsMap.put(option, score);
                        }

                    });
                    while (!future.isDone()) {
                    }

                }
            }

        }
        resultsMap.clear();

        if (oneDResultsMap.size() > 3) {
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
                    final String option =  mod1Id + "_" + mod2Id;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;

                    Future future = executor.submit(() -> {
                        int score = excuteSearch(optProtDataset, updatedName, "ModificationsParameter", option + "_noupdate"+"_fixed_Mods", tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                        double factor = (score - idRate) * 100.0 / idRate;//* optProtDataset.getDataEpsilon();;
                        factor = Math.round(factor);
                        if (score > oneDResultsMap.get(mod1Id) && score > threeDResultsMap.get(mod2Id) && factor >= optProtDataset.getAcceptedIDRatioThreshold()) {
                            resultsMap.put(option, score);
                        }
                    });
                    while (!future.isDone()) {
                    }

                }
            }

        }
        resultsMap.putAll(oneDResultsMap);
        resultsMap.putAll(twoDResultsMap);
        resultsMap.putAll(threeDResultsMap);
        TreeMap<Integer, Set<String>> orderedFixedModificationScore = new TreeMap<>(Collections.reverseOrder());
        for (String modificationId : resultsMap.keySet()) {
            if (resultsMap.get(modificationId) < idRate) {
                continue;
            }
            if (!orderedFixedModificationScore.containsKey(resultsMap.get(modificationId))) {
                orderedFixedModificationScore.put(resultsMap.get(modificationId), new LinkedHashSet<>());
            }
            orderedFixedModificationScore.get(resultsMap.get(modificationId)).add(modificationId);
        }
        if (!orderedFixedModificationScore.isEmpty() && idRate < orderedFixedModificationScore.firstKey()) {
            selectedFixedModificationOption.clear();
            Set<String> fmods = orderedFixedModificationScore.get(orderedFixedModificationScore.firstKey());

            for (var fm : fmods) {
                if (fm.contains("_")) {
                    selectedFixedModificationOption.addAll(Arrays.asList(fm.split("_")));
                } else {
                    selectedFixedModificationOption.add(fm);
                }
                break;
            }
            optProtDataset.setActiveIdentificationNum(orderedFixedModificationScore.firstKey());
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
        MainUtilities.cleanOutputFolder();
        final int fixedModScore = excuteSearch(optProtDataset, ftupdatedName, "ModificationsParameter", ftoption+"_refine_fixed", tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
        MainUtilities.cleanOutputFolder();

        //process variable modifications
        resultsMap.clear();
        Map<String, Set<String>> modifiedSpectrumMap = new LinkedHashMap<>();
        Map<String, Set<String>> fullSpectrumMap = new LinkedHashMap<>();
        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
        for (String modId : mods) {
            if (selectedFixedModificationOption.contains(modId)) {
                continue;
            }
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));

            Future future = executor.submit(() -> {
                List<SpectrumMatch> results = (excuteSearch(optProtDataset, updatedName, "ModificationsParameter", option + "_noupdate"+"_variable" , tempIdParam, false, optimisedSearchParameter, identificationParametersFile));
                int score = results.size();

                double factor = (score - fixedModScore) * 100.0 / fixedModScore;//* optProtDataset.getDataEpsilon();
                factor = Math.round(factor);
                if (factor >= optProtDataset.getAcceptedIDRatioThreshold()) {
                    modifiedSpectrumMap.put(modId, new HashSet<>());
                    fullSpectrumMap.put(modId, new HashSet<>());
                    //count the modification in the peptide
                    for (SpectrumMatch sm : results) {
                        if (sm == null) {
                            continue;
                        }
                        PeptideAssumption peptideAssm = sm.getAllPeptideAssumptions().toList().get(0);
                        if (peptideAssm.getPeptide().getVariableModifications().length > 0) {
                            for (ModificationMatch mm : peptideAssm.getPeptide().getVariableModifications()) {
                                if (mm.getModification().endsWith(ptmFactory.getModification(modId).getPattern().toString()) && Precision.equals(ptmFactory.getModification(modId).getMass(), Double.parseDouble(mm.getModification().split("@")[0]), 0.01)) {
                                    modifiedSpectrumMap.get(modId).add(sm.getSpectrumTitle());
                                }
                            }

                        }
                        fullSpectrumMap.get(modId).add(sm.getSpectrumTitle());
                    }
                    resultsMap.put(modId, score);
                }
            }
            );
            while (!future.isDone()) {
            }

        }
        oneDResultsMap.clear();
        TreeMap<Integer, List<String>> sorting1DMap = new TreeMap<>(Collections.reverseOrder());
        for (String mod : resultsMap.keySet()) {
            if (!modifiedSpectrumMap.get(mod).isEmpty()) {
                oneDResultsMap.put(mod, resultsMap.get(mod));
                if (!sorting1DMap.containsKey(resultsMap.get(mod))) {
                    sorting1DMap.put(resultsMap.get(mod), new ArrayList<>());
                }
                sorting1DMap.get(resultsMap.get(mod)).add(mod);
            }
        }
        //filter 0ne d map to max 4 modifications
        oneDResultsMap.clear();
        int modcounter = 0;
        if (!sorting1DMap.isEmpty()) {
            for (List<String> variableMods : sorting1DMap.values()) {
                for (String str : variableMods) {
                    oneDResultsMap.put(str, resultsMap.get(str));
                    modcounter++;
                }
                if (modcounter >= 4) {
                    break;
                }
            }
        }
        //start process 2d variable mod
        twoDResultsMap.clear();
        if (!sorting1DMap.isEmpty()) {
            for (String variableModI : sorting1DMap.firstEntry().getValue()) {
                for (String variableModII : oneDResultsMap.keySet()) {
                    if (variableModI.equalsIgnoreCase(variableModII)) {
                        continue;
                    }
                    Set<String> filterset = new HashSet<>(modifiedSpectrumMap.get(variableModI));
                    filterset.addAll(modifiedSpectrumMap.get(variableModII));
                    final String option =  variableModI + "_" + variableModII;
                    Set<String> fullset = new HashSet<>(fullSpectrumMap.get(variableModI));
                    fullset.addAll(fullSpectrumMap.get(variableModII));
                    modifiedSpectrumMap.put(option, new HashSet<>());
                    int intersectionCounter = 0;
                    for (String spec : filterset) {
                        if (modifiedSpectrumMap.get(variableModI).contains(spec) && modifiedSpectrumMap.get(variableModII).contains(spec)) {
                            intersectionCounter++;
                        } else {
                            modifiedSpectrumMap.get(option).add(spec);
                        }
                    }
                    if ((fullset.size() - intersectionCounter) > oneDResultsMap.get(variableModI)) {
                        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModI));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(variableModII));
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;

                        Future future = executor.submit(() -> {
                            int score = excuteSearch(optProtDataset, updatedName, "ModificationsParameter", option + "_noupdate"+"_variable_Mods", tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                            double factor = (score - Math.max(oneDResultsMap.get(variableModI), oneDResultsMap.get(variableModII))) * 100.0 / Math.max(oneDResultsMap.get(variableModI), oneDResultsMap.get(variableModII));// * optProtDataset.getDataEpsilon();
//                            System.out.println("at 2d score " + variableModI + ": " + oneDResultsMap.get(variableModI) + "  " + variableModII + ": " + oneDResultsMap.get(variableModII) + "  score " + score);
                            factor = Math.round(factor);
                            if (factor >= optProtDataset.getAcceptedIDRatioThreshold()) {
//                            if (!Precision.equals(oneDResultsMap.get(variableModI) * 0.01, score * 0.01, optProtDataset.getDataEpsilon()) && score > oneDResultsMap.get(variableModI) && !Precision.equals(oneDResultsMap.get(variableModII) * 0.01, score * 0.01, optProtDataset.getDataEpsilon()) && score > oneDResultsMap.get(variableModII)) {
                                twoDResultsMap.put(option, score);
                            }

                        });
                        while (!future.isDone()) {
                        }

                    }
                }
            }
        }
        MainUtilities.cleanOutputFolder();
        resultsMap.clear();
        threeDResultsMap.clear();
        if (!twoDResultsMap.isEmpty()) {
            for (String mod1Id : oneDResultsMap.keySet()) {
                for (String mod2Id : twoDResultsMap.keySet()) {
                    if (mod2Id.contains(mod1Id)) {
                        continue;
                    }
                    Set<String> filterset = new HashSet<>(modifiedSpectrumMap.get(mod1Id));
                    filterset.addAll(modifiedSpectrumMap.get(mod2Id.split("_")[0]));
                    filterset.addAll(modifiedSpectrumMap.get(mod2Id.split("_")[1]));
                    final String option =  mod1Id + "_" + mod2Id;
                    Set<String> fullset = new HashSet<>(fullSpectrumMap.get(mod1Id));
                    fullset.addAll(fullSpectrumMap.get(mod2Id.split("_")[0]));
                    fullset.addAll(fullSpectrumMap.get(mod2Id.split("_")[1]));
                    modifiedSpectrumMap.put(option, new HashSet<>());
                    int intersectionCounter = 0;
                    for (String spec : filterset) {
                        if ((modifiedSpectrumMap.get(mod1Id).contains(spec) && modifiedSpectrumMap.get(mod2Id.split("_")[0]).contains(spec)) || (modifiedSpectrumMap.get(mod1Id).contains(spec) && modifiedSpectrumMap.get(mod2Id.split("_")[1]).contains(spec)) || (modifiedSpectrumMap.get(mod2Id.split("_")[1]).contains(spec) && modifiedSpectrumMap.get(mod2Id.split("_")[0]).contains(spec))) {
                            intersectionCounter++;
                        } else {
                            modifiedSpectrumMap.get(option).add(spec);
                        }
                    }
                    if ((fullset.size() - intersectionCounter) > oneDResultsMap.get(mod1Id) && (fullset.size() - intersectionCounter) > twoDResultsMap.get(mod2Id)) {
                        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod1Id));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod2Id.split("_")[0]));
                        tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod2Id.split("_")[1]));
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;

                        Future future = executor.submit(() -> {
                            int score = excuteSearch(optProtDataset, updatedName, "ModificationsParameter", option + "_noupdate"+"_variable_Mods", tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                            double factor = (score - Math.max(oneDResultsMap.get(mod1Id), twoDResultsMap.get(mod2Id))) * 100.0 / Math.max(oneDResultsMap.get(mod1Id), twoDResultsMap.get(mod2Id));// * optProtDataset.getDataEpsilon();
                            factor = Math.round(factor);
                            if (factor >= optProtDataset.getAcceptedIDRatioThreshold()) {
                                threeDResultsMap.put(option, score);
                            }
                        });
                        while (!future.isDone()) {
                        }

                    }
                }
            }
            if (!threeDResultsMap.isEmpty()) {
                for (String mod1Id : oneDResultsMap.keySet()) {
                    for (String mod2Id : threeDResultsMap.keySet()) {
                        if (mod2Id.contains(mod1Id)) {
                            continue;
                        }
                        Set<String> filterset = new HashSet<>(modifiedSpectrumMap.get(mod1Id));
                        filterset.addAll(modifiedSpectrumMap.get(mod2Id.split("_")[0]));
                        filterset.addAll(modifiedSpectrumMap.get(mod2Id.split("_")[1]));
                        filterset.addAll(modifiedSpectrumMap.get(mod2Id.split("_")[2]));
                        final String option = mod1Id + "_" + mod2Id;
                        Set<String> fullset = new HashSet<>(fullSpectrumMap.get(mod1Id));
                        fullset.addAll(fullSpectrumMap.get(mod2Id.split("_")[0]));
                        fullset.addAll(fullSpectrumMap.get(mod2Id.split("_")[1]));
                        fullset.addAll(fullSpectrumMap.get(mod2Id.split("_")[2]));
                        modifiedSpectrumMap.put(option, new HashSet<>());
                        int intersectionCounter = 0;
                        for (String spec : filterset) {
                            if ((modifiedSpectrumMap.get(mod1Id).contains(spec) && modifiedSpectrumMap.get(mod2Id.split("_")[0]).contains(spec)) || (modifiedSpectrumMap.get(mod1Id).contains(spec) && modifiedSpectrumMap.get(mod2Id.split("_")[1]).contains(spec)) || (modifiedSpectrumMap.get(mod2Id.split("_")[1]).contains(spec) && modifiedSpectrumMap.get(mod2Id.split("_")[0]).contains(spec)) || (modifiedSpectrumMap.get(mod1Id).contains(spec) && modifiedSpectrumMap.get(mod2Id.split("_")[2]).contains(spec)) || (modifiedSpectrumMap.get(mod2Id.split("_")[1]).contains(spec) && modifiedSpectrumMap.get(mod2Id.split("_")[2]).contains(spec)) || (modifiedSpectrumMap.get(mod2Id.split("_")[2]).contains(spec) && modifiedSpectrumMap.get(mod2Id.split("_")[0]).contains(spec))) {
                                intersectionCounter++;
                            } else {
                                modifiedSpectrumMap.get(option).add(spec);
                            }
                        }
                        if ((fullset.size() - intersectionCounter) > oneDResultsMap.get(mod1Id) && (fullset.size() - intersectionCounter) > threeDResultsMap.get(mod2Id)) {
                            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod1Id));
                            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod2Id.split("_")[0]));
                            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod2Id.split("_")[1]));
                            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod2Id.split("_")[2]));
                            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
                            Future future = executor.submit(() -> {
                                int score = excuteSearch(optProtDataset, updatedName, "ModificationsParameter", option + "_noupdate"+"_variable_Mods" , tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                                if (score > oneDResultsMap.get(mod1Id) && score > threeDResultsMap.get(mod2Id)) {
                                    resultsMap.put(option, score);
                                }
                            });
                            while (!future.isDone()) {
                            }

                        }
                    }
                }

            }
        }

        System.out.println("final I   " + oneDResultsMap);//
        System.out.println("final II  " + twoDResultsMap);//
        System.out.println("final III " + threeDResultsMap);//
        System.out.println("final IV  " + resultsMap);//

        resultsMap.putAll(oneDResultsMap);
        resultsMap.putAll(twoDResultsMap);
        resultsMap.putAll(threeDResultsMap);

        TreeMap<Integer, Set<String>> orderedVariableModificationScore = new TreeMap<>(Collections.reverseOrder());
        for (String modificationId : resultsMap.keySet()) {
            if (resultsMap.get(modificationId) < fixedModScore) {
                continue;
            }

            if (!orderedVariableModificationScore.containsKey(resultsMap.get(modificationId))) {
                orderedVariableModificationScore.put(resultsMap.get(modificationId), new LinkedHashSet<>());
            }
            orderedVariableModificationScore.get(resultsMap.get(modificationId)).add(modificationId);
        }
//
        int varModScore;
        if (!orderedVariableModificationScore.isEmpty() && fixedModScore < orderedVariableModificationScore.firstKey()) {
            varModScore = orderedVariableModificationScore.firstKey();
            selectedVariableModificationOption.clear();
            Set<String> fmods = orderedVariableModificationScore.get(orderedVariableModificationScore.firstKey());
            for (var fm : fmods) {
                if (fm.contains("_")) {
                    selectedVariableModificationOption.addAll(Arrays.asList(fm.split("_")));
                } else {
                    selectedVariableModificationOption.add(fm);
                }
                break;
            }
        } else {
            varModScore = fixedModScore;
        }
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        for (String selectedVM : selectedVariableModificationOption) {
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(selectedVM));
        }

        //  maybe run one last run 
        final String vftupdatedName = Configurations.DEFAULT_RESULT_NAME + "fv_" + ftoption + "_" + msFileName;
        System.out.println("selectedVariableModificationOption " + selectedVariableModificationOption);
        int k = excuteSearch(optProtDataset, vftupdatedName,"ModificationsParameter", ftoption, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
        optProtDataset.setActiveIdentificationNum(k);

        //test refine mod 
        MainUtilities.cleanOutputFolder();
        if (optimisedSearchParameter.getSelectedSearchEngine().getIndex() == Advocate.xtandem.getIndex()) {
            Map<String, Integer> refinementVarModMap = Collections.synchronizedMap(new LinkedHashMap<>());
            for (String vMod : resultsMap.keySet()) {
                if (ptmFactory.getModification(vMod) == null) {
                    continue;
                }
                tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod));

                final String option =vMod;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;
                final ParameterScoreModel paramScore = new ParameterScoreModel();
                paramScore.setParamId("refineVariableModifications");
                Future future = executor.submit(() -> {
                    int score = excuteSearch(optProtDataset, updatedName,"ModificationsParameter", option+"_refine_variable", tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                    if (score > varModScore) {
                        refinementVarModMap.put(option, score);
                    }
                });
                while (!future.isDone()) {
                }

            }
            //2d refinment
            resultsMap.clear();
            if (refinementVarModMap.size() > 1) {
                for (String selectedRef : refinementVarModMap.keySet()) {
                    for (String vMod : refinementVarModMap.keySet()) {
                        if (vMod.equals(selectedRef)) {
                            continue;
                        }
                        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(selectedRef));
                        final String option =vMod + "_" + selectedRef;
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName; 
                        Future future = executor.submit(() -> {
                            int score = excuteSearch(optProtDataset, updatedName, "ModificationsParameter",option+ "_refine_variable_", tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                            if (score > refinementVarModMap.get(vMod)) {
                                resultsMap.put(option, score);
                            }
                        });
                        while (!future.isDone()) {
                        }

                    }
                }
            }
            //3d refinment
            Map<String, Integer> threeDrefinementVarModMap = Collections.synchronizedMap(new LinkedHashMap<>());
            if (!resultsMap.isEmpty()) {
                for (String selectedRef : refinementVarModMap.keySet()) {
                    for (String vMod : resultsMap.keySet()) {
                        if (vMod.equals(selectedRef)) {
                            continue;
                        }
                        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[0]));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[1]));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(selectedRef));
                        final String option = vMod + "_" + selectedRef;
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;
                        Future future = executor.submit(() -> {
                            int score = excuteSearch(optProtDataset, updatedName,"ModificationsParameter", option+"_refine_variable", tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                            if (score > resultsMap.get(vMod)) {
                                threeDrefinementVarModMap.put(option, score);
                            }
                        });
                        while (!future.isDone()) {
                        }

                    }
                }

            }
            //4d refinment
            Map<String, Integer> fourDrefinementVarModMap = Collections.synchronizedMap(new LinkedHashMap<>());
            if (!threeDResultsMap.isEmpty()) {
                for (String selectedRef : refinementVarModMap.keySet()) {
                    for (String vMod : threeDResultsMap.keySet()) {
                        if (vMod.equals(selectedRef)) {
                            continue;
                        }
                        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[0]));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[1]));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod.split("_")[2]));
                        tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(selectedRef));
                        final String option =vMod + "_" + selectedRef;
                        final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;
                       
                        Future future = executor.submit(() -> {
                            int score = excuteSearch(optProtDataset, updatedName,"ModificationsParameter", option+"_refine_variable", tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                            if (score > threeDResultsMap.get(vMod)) {
                                fourDrefinementVarModMap.put(option, score);
                            }
                        });
                        while (!future.isDone()) {
                        }

                    }
                }

            }
            refinementVarModMap.putAll(resultsMap);
            refinementVarModMap.putAll(threeDResultsMap);
            refinementVarModMap.putAll(fourDrefinementVarModMap);
            modificationsResults.put("refinmentVariableModifications", new HashSet<>(refinementVarModMap.keySet()));
        }
        modificationsResults.put("variableModifications", new HashSet<>(selectedVariableModificationOption));
        MainUtilities.cleanOutputFolder();
        return modificationsResults;

    }

    public abstract List<SpectrumMatch> excuteSearch(SearchingSubDataset optProtDataset, String defaultOutputFileName, String paramId, String paramOption, IdentificationParameters tempIdParam, boolean addPeptideMasses, SearchInputSetting optProtSearchSettings, File identificationParametersFile);

}
