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
import java.util.concurrent.Future;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.dataset.model.SearchingSubDataset;
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
        int idRate = optProtDataset.getIdentificationNum();
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        String[] cleavageParameters = new String[]{"enzyme", "wholeProtein", "unSpecific"};
        // Creates a pool of 3 threads
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
                resultsMap.put(option, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!f.isDone()) {
            }
            
            if (resultsMap.get(option) >= optProtDataset.getIdentificationNum()) {
                break;
            }
        }
//        
        for (String option : resultsMap.keySet()) {
            if (resultsMap.get(option) > idRate) {
                idRate = resultsMap.get(option);
                selectedOption = option;
                if (optProtDataset.getIdentificationNum() < idRate) {
                    optProtDataset.setIdentificationNum(idRate);
                }
            }
        }
        return selectedOption;
        
    }
    
    public String optimizeEnzymeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        if (!oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name().equalsIgnoreCase("enzyme")) {
            return null;
        }
        String selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
        int idRate = optProtDataset.getIdentificationNum();
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
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().clearEnzymes();
            tempIdParam.getSearchParameters().getDigestionParameters().addEnzyme(enzyme);
            final String option = enzyme.getName();
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            Future f = executor.submit(() -> {
                resultsMap.put(option, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!f.isDone()) {
            }
            if (enzyme.getName().equals("Trypsin") && resultsMap.get(option) >= (optProtDataset.getIdentificationNum())) {
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
        if (localId >= idRate) {
            selectedOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
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
        int idRate = optProtDataset.getIdentificationNum();
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        
        for (int i = 0; i < DigestionParameters.Specificity.values().length; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setSpecificity(enzymeName, DigestionParameters.Specificity.getSpecificity(i));
            final String option = DigestionParameters.Specificity.getSpecificity(i).name();
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            Future f = executor.submit(() -> {
                resultsMap.put(option, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
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
        if (localId >= idRate) {
            selectedOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return selectedOption;
        
    }
    
    public Integer optimizeMaxMissCleavagesParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        if (!oreginaltempIdParam.getSearchParameters().getDigestionParameters().getCleavageParameter().name().equalsIgnoreCase("enzyme")) {
            return null;
        }
        String enzymeName = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getEnzymes().get(0).getName();
        Integer selectedOption = oreginaltempIdParam.getSearchParameters().getDigestionParameters().getnMissedCleavages(enzymeName);
        int idRate = optProtDataset.getIdentificationNum();
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        
        for (int i = 0; i < 8; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getDigestionParameters().setnMissedCleavages(enzymeName, i);
            final String option = "missedCleavages_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            final int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
            System.out.println("# of missed clev " + i + " #  " + resultsMap.get(i) + "  " + j + "  " + idRate);
        }
        int localId = -1;
        int localSelection = 0;
        for (int option : resultsMap.keySet()) {
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
        }
        return selectedOption;
        
    }
    
    public String optimizeFragmentIonTypesParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        String selectedOption = oreginaltempIdParam.getSearchParameters().getForwardIons() + "-" + oreginaltempIdParam.getSearchParameters().getRewindIons();
        int idRate = optProtDataset.getIdentificationNum();
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        SearchParameters tempSearchParameters = tempIdParam.getSearchParameters();
        ArrayList<Integer> selectedForwardIons = tempSearchParameters.getForwardIons();
        String[] forwardIons = new String[]{"a", "b", "c"};
        String[] rewindIons = new String[]{"x", "y", "z"};
        
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
                    resultsMap.put(option, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
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
        if (localId >= idRate) {
            selectedOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
        }
        
        selectedOption = selectedOption.replace("[", "").replace("]", "");
        return selectedOption;
    }
    
    public double optimizePrecursorToleranceParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedOption = oreginaltempIdParam.getSearchParameters().getPrecursorAccuracy();
        int idRate = optProtDataset.getIdentificationNum();
        Map<Double, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        
        for (double i = 5; i < 30;) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setPrecursorAccuracy(i);
            tempIdParam.getSearchParameters().setPrecursorAccuracyType(SearchParameters.MassAccuracyType.PPM);
            final String option = "precursorAccuracy_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            final double j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
            i += 5;
        }
        int localId = -1;
        double localSelection = 0;
        for (double option : resultsMap.keySet()) {
//            System.out.println(" option 1 " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
        }
        return selectedOption;
    }
    
    public double optimizeFragmentToleranceParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        double selectedOption = oreginaltempIdParam.getSearchParameters().getPrecursorAccuracy();
        int idRate = optProtDataset.getIdentificationNum();
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
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
        }
        int localId = -1;
        double localSelection = 0;
        for (double option : resultsMap.keySet()) {
//            System.out.println(" option 2 " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;
                
            }
        }
        if (localId >= idRate) {
            selectedOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
        }
        return selectedOption;
        
    }
    
    public int[] optimizePrecursorChargeParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxChargeOption = oreginaltempIdParam.getSearchParameters().getMaxChargeSearched();
        int selectedMinChargeOption = oreginaltempIdParam.getSearchParameters().getMinChargeSearched();
        int idRate = optProtDataset.getIdentificationNum();
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        
        for (int i = 1; i < 6; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setMaxChargeSearched(i);
            if (i < selectedMinChargeOption) {
                selectedMinChargeOption = i;
                tempIdParam.getSearchParameters().setMinChargeSearched(i);
            }
            final String option = "maxCharge_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
        }
        int localId = -1;
        int localSelection = 0;
        for (int option : resultsMap.keySet()) {
            System.out.println(" max charge " + option + "  " + resultsMap.get(option) + " > " + localId + "  ");
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;
                
            }
        }
        if (localId >= idRate) {
            selectedMaxChargeOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
            idRate = localId;
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
                
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
        }
        
        localId = -1;
        localSelection = 0;
        for (int option : resultsMap.keySet()) {
//            System.out.println(" min charge " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;
                
            }
        }
        if (localId >= idRate) {
            selectedMinChargeOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
        }
        MainUtilities.cleanOutputFolder();
        return new int[]{selectedMinChargeOption, selectedMaxChargeOption};
        
    }
    
    public int[] optimizeIsotopParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        int selectedMaxIsotopicCorrectionOption = oreginaltempIdParam.getSearchParameters().getMaxIsotopicCorrection();
        int selectedMinIsotopicCorrectioneOption = oreginaltempIdParam.getSearchParameters().getMinIsotopicCorrection();
        int idRate = optProtDataset.getIdentificationNum();
        Map<Integer, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        for (int i = 7; i > -8; i--) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            if (i > selectedMaxIsotopicCorrectionOption) {
                tempIdParam.getSearchParameters().setMaxIsotopicCorrection(i);
            }
            tempIdParam.getSearchParameters().setMinChargeSearched(i);
            final String option = "minIsotop_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
        }
        
        int localId = -1;
        int localSelection = 0;
        for (int option : resultsMap.keySet()) {
//            System.out.println(" min isotopic  " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;
                
            }
        }
        if (localId >= idRate) {
            selectedMinIsotopicCorrectioneOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
        }
        resultsMap.clear();
        
        idRate = -1;
        for (int i = selectedMinIsotopicCorrectioneOption; i < 9; i++) {
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().setMaxIsotopicCorrection(i);
            tempIdParam.getSearchParameters().setMinIsotopicCorrection(selectedMinIsotopicCorrectioneOption);
            final String option = "maxIsotop_" + i;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "_" + option + "_" + msFileName;
            int j = i;
            Future future = executor.submit(() -> {
                resultsMap.put(j, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
        }
        localId = -1;
        localSelection = 0;
        for (int option : resultsMap.keySet()) {
//            System.out.println(" max isotopic  " + option + "  " + resultsMap.get(option) + " > " + localId + "  " + idRate);
            if (resultsMap.get(option) > localId) {
                localId = resultsMap.get(option);
                localSelection = option;
                
            }
        }
        if (localId >= idRate) {
            selectedMaxIsotopicCorrectionOption = localSelection;
            optProtDataset.setIdentificationNum(localId);
        }
        
        return new int[]{selectedMinIsotopicCorrectioneOption, selectedMaxIsotopicCorrectionOption};
    }
    
    public Map<String, Set<String>> optimizeModificationsParameter(SearchingSubDataset optProtDataset, File identificationParametersFile, SearchInputSetting optimisedSearchParameter) throws IOException {
        IdentificationParameters oreginaltempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
        String msFileName = IoUtil.removeExtension(optProtDataset.getSubMsFile().getName());
        ArrayList<String> selectedFixedModificationOption = new ArrayList<>(oreginaltempIdParam.getSearchParameters().getModificationParameters().getFixedModifications());
        ArrayList<String> selectedVariableModificationOption = new ArrayList<>(oreginaltempIdParam.getSearchParameters().getModificationParameters().getVariableModifications());
        int idRate = optProtDataset.getIdentificationNum();
        Map<String, Set<String>> modificationsResults = new HashMap<>();
        
        List<String> mods = new ArrayList<>();
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common));
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common_Biological));
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common_Artifact));
        
        Map<String, Integer> resultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        for (String modId : mods) {
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;
            
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
            tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
            tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(modId));
            
            Future future = executor.submit(() -> {
                resultsMap.put(modId, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
            });
            while (!future.isDone()) {
            }
        }
        System.out.println("top id number I " + optProtDataset.getIdentificationNum());
        Map<String, Map<String, Integer>> targtedFixedModificationScore = new TreeMap<>();
        for (String modificationId : resultsMap.keySet()) {
            if (resultsMap.get(modificationId) < idRate) {
                continue;
            }
            
            String modPattern = ptmFactory.getModification(modificationId).getPattern().toString();
            if (!targtedFixedModificationScore.containsKey(modPattern)) {
                targtedFixedModificationScore.put(modPattern, new TreeMap<>());
            }
            targtedFixedModificationScore.get(modPattern).put(modificationId, resultsMap.get(modificationId));
        }
        System.out.println("top id number II " + optProtDataset.getIdentificationNum());
        Map<String, Integer> oneDResultsMap = new TreeMap<>();
        for (String aa : targtedFixedModificationScore.keySet()) {
            Map<String, Integer> modMap = targtedFixedModificationScore.get(aa);
            int selectedoptionNum = -1;
            String selectedMod = "";
            for (String mod : modMap.keySet()) {
                if (selectedoptionNum < modMap.get(mod)) {
                    selectedoptionNum = modMap.get(mod);
                    selectedMod = mod;
                }
            }
            oneDResultsMap.put(selectedMod, selectedoptionNum);
            
        }
        
        Set<String> clearSet = new LinkedHashSet<>(oneDResultsMap.keySet());
        for (String modId : clearSet) {
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "m_" + option + "_" + msFileName;
            
            IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
            tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
            tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));
            
            Future future = executor.submit(() -> {
                int vmodScore = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
//                double epsilon = 0.00005 * Configurations.REFINED_MS_SIZE;
                if (!Precision.equals(vmodScore * 0.01, oneDResultsMap.get(modId) * 0.01, optProtDataset.getDataEpsilon()) && vmodScore > oneDResultsMap.get(modId)) {
                    oneDResultsMap.remove(modId);
                }
            });
            while (!future.isDone()) {
            }
        }
        resultsMap.clear();
        Map<String, Integer> twoDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        int indexI = 0;
        for (String mod1Id : oneDResultsMap.keySet()) {
            int indexII = 0;
            for (String mod2Id : oneDResultsMap.keySet()) {
                if (indexII <= indexI) {
                    indexII++;
                    continue;
                }
                IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
                tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
                tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod1Id));
                tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id));
                final String option = mod1Id + "_" + mod2Id;
                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;
                Future future = executor.submit(() -> {
                    twoDResultsMap.put(option, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                });
                while (!future.isDone()) {
                }
                indexII++;
            }
            indexI++;
        }
        
        if (oneDResultsMap.size() > 2) {
            for (String mod1Id : oneDResultsMap.keySet()) {
                for (String mod2Id : twoDResultsMap.keySet()) {
                    if (mod2Id.contains(mod1Id)) {
                        continue;
                    }
                    IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
                    tempIdParam.getSearchParameters().getModificationParameters().clearFixedModifications();
                    tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
                    tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                    tempIdParam.getSearchParameters().getModificationParameters().getRefinementFixedModifications().clear();
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod1Id));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[0]));
                    tempIdParam.getSearchParameters().getModificationParameters().addFixedModification(ptmFactory.getModification(mod2Id.split("_")[1]));
                    final String option = mod1Id + "_" + mod2Id;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "f_" + option + "_" + msFileName;
                    Future future = executor.submit(() -> {
                        resultsMap.put(option, excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size());
                    });
                    while (!future.isDone()) {
                    }
                }
            }
            
        }
        
        resultsMap.putAll(twoDResultsMap);
        resultsMap.putAll(oneDResultsMap);
        
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
        }
        
        System.out.println("final fixed mod id " + selectedFixedModificationOption);//
        modificationsResults.put("fixedModifications", new HashSet<>(selectedFixedModificationOption));
        MainUtilities.cleanOutputFolder();
        
        IdentificationParameters tempIdParam = IdentificationParameters.getIdentificationParameters(identificationParametersFile);
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
//        System.out.println("excute fm search before : " + excuteSearch(optProtDataset, ftupdatedName, ftoption, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size() + "/" + idRate);
        MainUtilities.cleanOutputFolder();
        for (String fixedModifications : resultsMap.keySet()) {
            if (resultsMap.get(fixedModifications) > idRate) {
                for (String str : fixedModifications.split("_")) {
                    tempIdParam.getSearchParameters().getModificationParameters().addRefinementFixedModification(ptmFactory.getModification(str));
                }
            }
        }
        int fixedModScore;
        if (!selectedFixedModificationOption.isEmpty()) {
            fixedModScore = excuteSearch(optProtDataset, ftupdatedName, ftoption, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
        } else {
            fixedModScore = idRate;
        }
//        System.out.println("excute fm search after : " + fixedModScore + "/" + idRate);

        optProtDataset.setIdentificationNum(fixedModScore);
        MainUtilities.cleanOutputFolder();
        //process variable modifications
        resultsMap.clear();
        Map<String, Set<String>> modifiedSpectrumMap = new LinkedHashMap<>();
        Map<String, Set<String>> fullSpectrumMap = new LinkedHashMap<>();
        for (String modId : mods) {
            if (selectedFixedModificationOption.contains(modId)) {
                continue;
            }
            final String option = modId;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
            tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
            tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(modId));
            
            Future future = executor.submit(() -> {
                List<SpectrumMatch> results = (excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile));
                int score = results.size();
//                double epsilon = 0.00005 * Configurations.REFINED_MS_SIZE;
                if (!Precision.equals(fixedModScore * 0.01, score * 0.01, optProtDataset.getDataEpsilon()) && score > fixedModScore) {
                    modifiedSpectrumMap.put(modId, new HashSet<>());
                    fullSpectrumMap.put(modId, new HashSet<>());
                    //count the modification in the peptide
                    for (SpectrumMatch sm : results) {
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

//                    
            }
            );
            while (!future.isDone()) {
            }
            
        }
        oneDResultsMap.clear();
        TreeMap<Integer, List<String>> sorting1DMap = new TreeMap<>();
        for (String mod : resultsMap.keySet()) {
            if (!modifiedSpectrumMap.get(mod).isEmpty()) {
                oneDResultsMap.put(mod, resultsMap.get(mod));
                if (!sorting1DMap.containsKey(resultsMap.get(mod))) {
                    sorting1DMap.put(resultsMap.get(mod), new ArrayList<>());
                }
                sorting1DMap.get(resultsMap.get(mod)).add(mod);
            }
        }

        //start process 2d variable mod
        twoDResultsMap.clear();
        for (String variableModI : sorting1DMap.lastEntry().getValue()) {
            
            for (String variableModII : oneDResultsMap.keySet()) {
                if (variableModI.equalsIgnoreCase(variableModII)) {
                    continue;
                }
                Set<String> filterset = new HashSet<>(modifiedSpectrumMap.get(variableModI));
                filterset.addAll(modifiedSpectrumMap.get(variableModII));
                final String option = variableModI + "_" + variableModII;
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
                        int score = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                        if (!Precision.equals(oneDResultsMap.get(variableModI) * 0.01, score * 0.01, optProtDataset.getDataEpsilon()) && score > oneDResultsMap.get(variableModI) && !Precision.equals(oneDResultsMap.get(variableModII) * 0.01, score * 0.01, optProtDataset.getDataEpsilon()) && score > oneDResultsMap.get(variableModII)) {
                            System.out.println("at -------------------------------------------------->>  2d score is  " + option + " score:  " + score + "  " + oneDResultsMap.get(variableModI) + "  " + oneDResultsMap.get(variableModII) + "   " + optProtDataset.getIdentificationNum());
                            twoDResultsMap.put(option, score);
                        }
                        
                    });
                    while (!future.isDone()) {
                    }
                }
            }
        }
        MainUtilities.cleanOutputFolder();

//
//        indexI = 0;
//        for (String variableModI : modifiedSpectrumMap.keySet()) {
//            System.out.println("-------------------before------------------------------->> to include v mod " + variableModI + "  #spectra " + modifiedSpectrumMap.get(variableModI).size() + "   " + optProtDataset.getIdentificationNum());
//            tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(variableModI));
//            int indexII = 0;
//            Set<String> filterset = new HashSet<>(modifiedSpectrumMap.get(variableModI));
//            Set<String> toClean = new HashSet<>();
//            for (String variableModII : modifiedSpectrumMap.keySet()) {
//                if (indexII <= indexI) {
//                    indexII++;
//                    continue;
//                }
//                int intersectionCounter = 0;
//                for (String spec : filterset) {
//                    if (modifiedSpectrumMap.get(variableModI).contains(spec) && modifiedSpectrumMap.get(variableModII).contains(spec)) {
////                        modifiedSpectrumMap.get(variableModII).remove(spec);
//                        toClean.add(spec);
//                        intersectionCounter++;
//                    }
//
//                }
//                System.out.println("intersection counter between (" + variableModI + ")and (" + variableModII + ") is " + intersectionCounter + "   " + modifiedSpectrumMap.get(variableModII).size());
//                indexII++;
//            }
//            indexI++;
////            modifiedSpectrumMap.get(variableModI).removeAll(toClean);
//            System.out.println("-----------------after--------------------------------->> to include v mod " + variableModI + "  #spectra " + modifiedSpectrumMap.get(variableModI).size());
//        }
//        MainUtilities.cleanOutputFolder();
//
//        for (String mod : oneDResultsMap.keySet()) {
//            if (oneDResultsMap.get(mod) > optProtDataset.getIdentificationNum()) {
//                optProtDataset.setIdentificationNum(idRate);
//            }
//        }
//
//        indexI = 0;
//        twoDResultsMap.clear();
        resultsMap.clear();
//        double epsilon = 0.00005 * Configurations.REFINED_MS_SIZE;
//        for (String mod1Id : oneDResultsMap.keySet()) {
//            int indexII = 0;
//            for (String mod2Id : oneDResultsMap.keySet()) {
//                if (indexII <= indexI) {
//                    indexII++;
//                    continue;
//                }
//                tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
////                tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
//                tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod1Id));
//                tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod2Id));
//                final String option = mod1Id + "_" + mod2Id;
//                final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
//                Future future = executor.submit(() -> {
//                    int score = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
////                    if (score > optProtDataset.getIdentificationNum()) {
//                    System.out.println("at -------------------------------------------------->>  2d value is  " + option + " score:  " + score + "  " + oneDResultsMap.get(mod1Id) + "  " + oneDResultsMap.get(mod2Id) + "   " + optProtDataset.getIdentificationNum());
//
//                    if (!Precision.equals(oneDResultsMap.get(mod1Id) * 0.01, score * 0.01, optProtDataset.getDataEpsilon()) && score > oneDResultsMap.get(mod1Id) && !Precision.equals(oneDResultsMap.get(mod2Id) * 0.01, score * 0.01, optProtDataset.getDataEpsilon()) && score > oneDResultsMap.get(mod2Id)) {
//
////                    if (score > oneDResultsMap.get(mod1Id) && score > oneDResultsMap.get(mod2Id)) {
//                        System.out.println("at -------------------------------------------------->>  2d score is  " + option + " score:  " + score + "  " + oneDResultsMap.get(mod1Id) + "  " + oneDResultsMap.get(mod2Id) + "   " + optProtDataset.getIdentificationNum());
//                        twoDResultsMap.put(option, score);
////                        optProtDataset.setIdentificationNum(score);
//                    }
//                });
//                while (!future.isDone()) {
//                }
//                indexII++;
//            }
//            indexI++;
//        }

        Map<String, Integer> threeDResultsMap = Collections.synchronizedMap(new LinkedHashMap<>());
        
        if (!twoDResultsMap.isEmpty()) {
            for (String mod1Id : oneDResultsMap.keySet()) {
                for (String mod2Id : twoDResultsMap.keySet()) {
                    if (mod2Id.contains(mod1Id)) {
                        continue;
                    }
                    
                    Set<String> filterset = new HashSet<>(modifiedSpectrumMap.get(mod1Id));
                    filterset.addAll(modifiedSpectrumMap.get(mod2Id.split("_")[0]));
                    filterset.addAll(modifiedSpectrumMap.get(mod2Id.split("_")[1]));
                    final String option = mod1Id + "_" + mod2Id;
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
                            int score = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                            if (!Precision.equals(oneDResultsMap.get(mod1Id) * 0.01, score * 0.01, optProtDataset.getDataEpsilon()) && score > oneDResultsMap.get(mod1Id) && !Precision.equals(twoDResultsMap.get(mod2Id) * 0.01, score * 0.01, optProtDataset.getDataEpsilon()) && score > twoDResultsMap.get(mod2Id)) {
                                System.out.println("at-------------------------------------------------->>  3d score is  " + option + " score:  " + score + "  " + oneDResultsMap.get(mod1Id) + "  " + twoDResultsMap.get(mod2Id));
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
//                        tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
                            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod1Id));
                            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod2Id.split("_")[0]));
                            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod2Id.split("_")[1]));
                            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(mod2Id.split("_")[2]));
                            
                            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "v_" + option + "_" + msFileName;
                            Future future = executor.submit(() -> {
                                int score = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                                
                                if (score > oneDResultsMap.get(mod1Id) && score > threeDResultsMap.get(mod2Id)) {
                                    System.out.println("at-------------------------------------------------->>  4d score is  " + option + " score:  " + score + " mod1Id " + oneDResultsMap.get(mod1Id) + "  " + threeDResultsMap.get(mod2Id));
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
//            System.out.println("Variable mods " + modificationId + " score " + resultsMap.get(modificationId));

        }
        int varModScore = orderedVariableModificationScore.firstKey();
//

        if (!orderedVariableModificationScore.isEmpty() && fixedModScore < orderedVariableModificationScore.firstKey()) {
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
        }
        tempIdParam.getSearchParameters().getModificationParameters().clearVariableModifications();
        for (String selectedVM : selectedVariableModificationOption) {
            tempIdParam.getSearchParameters().getModificationParameters().addVariableModification(ptmFactory.getModification(selectedVM));
        }
        //test refine mod 
        MainUtilities.cleanOutputFolder();
        Map<String, Integer> refinementVarModMap = Collections.synchronizedMap(new LinkedHashMap<>());
        for (String vMod : resultsMap.keySet()) {
            
            tempIdParam.getSearchParameters().getModificationParameters().clearRefinementModifications();
            tempIdParam.getSearchParameters().getModificationParameters().addRefinementVariableModification(ptmFactory.getModification(vMod));
            
            final String option = vMod;
            final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;
            Future future = executor.submit(() -> {
                int score = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                System.out.println("at-------------------------------------------------->>  refinment score is  " + option + " score:  " + score + " vmscore  " + varModScore);
                
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
                    final String option = vMod + "_" + selectedRef;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;
                    Future future = executor.submit(() -> {
                        int score = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                        if (score > refinementVarModMap.get(vMod)) {
                            System.out.println("at-------------------------------------------------->>  refinment score is  " + option + " score:  " + score + " vmscore  " + refinementVarModMap.get(vMod));
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
                        int score = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                        if (score > resultsMap.get(vMod)) {
                            System.out.println("at-------------------------------------------------->>  refinment score is  " + option + " score:  " + score + " vmscore  " + resultsMap.get(vMod));
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
                    final String option = vMod + "_" + selectedRef;
                    final String updatedName = Configurations.DEFAULT_RESULT_NAME + "rv_" + option + "_" + msFileName;
                    Future future = executor.submit(() -> {
                        int score = excuteSearch(optProtDataset, updatedName, option, tempIdParam, false, optimisedSearchParameter, identificationParametersFile).size();
                        if (score > threeDResultsMap.get(vMod)) {
                            System.out.println("at-------------------------------------------------->>  refinment score is  " + option + " score:  " + score + " vmscore  " + threeDResultsMap.get(vMod));
                            fourDrefinementVarModMap.put(option, score);
                        }
                    });
                    while (!future.isDone()) {
                    }
                    
                }
            }
            
        }
        refinementVarModMap.putAll(resultsMap);
        refinementVarModMap.putAll(resultsMap);
        refinementVarModMap.putAll(threeDResultsMap);
        refinementVarModMap.putAll(fourDrefinementVarModMap);
        System.out.println("refinementVarModMap " + refinementVarModMap);
        
        modificationsResults.put("variableModifications", new HashSet<>(selectedVariableModificationOption));
        modificationsResults.put("refinmentVariableModifications", new HashSet<>(refinementVarModMap.keySet()));
        MainUtilities.cleanOutputFolder();
        
        return modificationsResults;
        
    }
    
    public abstract ArrayList<SpectrumMatch> excuteSearch(SearchingSubDataset optProtDataset, String defaultOutputFileName, String paramOption, IdentificationParameters tempIdParam, boolean addPeptideMasses, SearchInputSetting optProtSearchSettings, File identificationParametersFile);
    
}
