/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.configurations;

import com.compomics.util.experiment.biology.enzymes.Enzyme;
import com.compomics.util.experiment.biology.enzymes.EnzymeFactory;
import com.compomics.util.experiment.biology.modifications.ModificationCategory;
import com.compomics.util.experiment.biology.modifications.ModificationFactory;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author yfa041
 */
public class SearchEngineParameterConfigurations {

    private final Map<String, Boolean> paramMap;

    public SearchEngineParameterConfigurations() {
        paramMap = new HashMap<>();
        paramMap.put("enzyme", true);
        paramMap.put("wholeProtein", true);
        paramMap.put("unSpecific", true);
        for (Enzyme enzyme : EnzymeFactory.getInstance().getEnzymes()) {
            paramMap.put(enzyme.getName(), true);
        }
        paramMap.put("specific", true);
        paramMap.put("semiSpecific", true);
        paramMap.put("specificNTermOnly", true);
        paramMap.put("specificCTermOnl", true);
        paramMap.put("missedCleavages", true);
        paramMap.put("[0]-[3]", true);
        paramMap.put("[0]-[4]", true);
        paramMap.put("[0]-[5]", true);

        paramMap.put("[1]-[3]", true);
        paramMap.put("[1]-[4]", true);
        paramMap.put("[1]-[5]", true);

        paramMap.put("[2]-[3]", true);
        paramMap.put("[2]-[4]", true);
        paramMap.put("[2]-[5]", true);
        paramMap.put("fragmentAccuracy", true);
        paramMap.put("precursorAccuracy", true);
        paramMap.put("maxCharge", true);
        paramMap.put("minCharge", true);
        paramMap.put("maxIsotop", true);
        paramMap.put("minIsotop", true);
        paramMap.put("reference", true);
        ModificationFactory ptmFactory = ModificationFactory.getInstance();
        List<String> mods = new ArrayList<>(ptmFactory.getModifications(ModificationCategory.Common_Biological));
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common));
        mods.addAll(ptmFactory.getModifications(ModificationCategory.Common_Artifact));

        paramMap.put("spectrumDR", true);
        paramMap.put("peaksNum", true);
        paramMap.put("minimumFragmentMz", true);
        paramMap.put("minpeaksNum", true);
        paramMap.put("noiseSupression", true);
        paramMap.put("parentMonoisotopicMassIsotopeError", true);
        paramMap.put("useQuickAcetyl", true);
        paramMap.put("useStpBias", true);
        paramMap.put("useQuickPyrolidone", true);
        paramMap.put("useRefine", true);
        paramMap.put("useRefineUnanticipatedCleavages", true);
        paramMap.put("useRefineSimiEnzymaticCleavage", true);
        paramMap.put("usePotintialModification", true);
        paramMap.put("useRefinePointMutations", true);
        paramMap.put("useRefineSnAPs", true);
        paramMap.put("useRefineSpectrumSynthesis", true);
        paramMap.put("", false);
        paramMap.put(null, false);
        paramMap.put("maxPrecursorMass", true);
        paramMap.put("minPrecursorMass", true);

        paramMap.put("minPeptideLength", true);
        paramMap.put("maxPeptideLength", true);
        paramMap.put("maxVarPTMs", true);
        paramMap.put("fragmentationMethod", true);
        paramMap.put("enzymatricTerminals", true);
        paramMap.put("useSmartPlus3Model", true);
        paramMap.put("computeXCorr", true);
        paramMap.put("TICCutoff", true);
        paramMap.put("NumberOfIntensityClasses", true);
        paramMap.put("classSizeMultiplier", true);
        paramMap.put("NumberOfBatches", true);
        paramMap.put("maxPeakCount", true);

        for (String mod : mods) {
            paramMap.put(mod, true);
        }
    }

    public void disableMinIsotop() {
        paramMap.replace("minIsotop", false);
    }

    public void disableMaxIsotop() {
        paramMap.replace("maxIsotop", false);
    }

    public void disableMinCharge() {
        paramMap.replace("minCharge", false);
    }

    public void disableMissedCleavages() {
        paramMap.replace("missedCleavages", false);
    }

    public boolean isEnabledParam(String param) {
        Boolean b = paramMap.get(param);
        if (b == null) {
            b = false;
        }
        return b;
    }

    public void disableSpecificEnzyme() {
        paramMap.replace("specificEnzyme", false);
    }

    public void disableSemiSpecificEnzyme() {
        paramMap.replace("semiSpecificEnzyme", false);
    }

    public void disableSpecificNTermOnlyEnzyme() {
        paramMap.replace("specificNTermOnly", false);
    }

    public boolean isSpecificCTermOnlyEnzyme() {
        return paramMap.get("specificCTermOnly");
    }

    public void disableSpecificCTermOnlyEnzyme() {
        paramMap.replace("specificCTermOnly ", false);
    }

    public void disableEnzyme() {
        this.paramMap.replace("enzyme", false);
    }

    public void setWholeProtein(boolean wholeProtein) {
        this.paramMap.replace("wholeProtein", false);
    }

    public void disableUnSpecific(boolean unSpecific) {
        paramMap.replace("unSpecific", false);
    }

    public void disableParam(String paramName) {
        if (paramMap.containsKey(paramName)) {
            paramMap.replace(paramName, false);
        } else {
            System.out.println("param " + paramName + " not supported");
        }
    }

    public void disableEnzyme(Enzyme enzyme) {
        paramMap.replace(enzyme.getName(), false);
    }

    public void disableFragmentTypes() {
        paramMap.put("[0]-[3]", false);
        paramMap.put("[0]-[4]", false);
        paramMap.put("[0]-[5]", false);

        paramMap.put("[1]-[3]", false);
        paramMap.put("[1]-[4]", false);
        paramMap.put("[1]-[5]", false);

        paramMap.put("[2]-[3]", false);
        paramMap.put("[2]-[4]", false);
        paramMap.put("[2]-[5]", false);
    }

}
