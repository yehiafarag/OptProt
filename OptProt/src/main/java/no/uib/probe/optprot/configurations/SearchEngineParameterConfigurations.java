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
        paramMap.put("specificEnzyme", true);
        paramMap.put("semiSpecificEnzyme", true);
        paramMap.put("specificNTermOnlyEnzyme", true);
        paramMap.put("specificCTermOnlyEnzyme", true);
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
        return paramMap.get(param);
    }

    public void disableSpecificEnzyme() {
        paramMap.replace("specificEnzyme", false);
    }

    public void disableSemiSpecificEnzyme() {
        paramMap.replace("semiSpecificEnzyme", false);
    }

    public void disableSpecificNTermOnlyEnzyme() {
        paramMap.replace("specificNTermOnlyEnzyme", false);
    }

    public boolean isSpecificCTermOnlyEnzyme() {
        return paramMap.get("specificCTermOnlyEnzyme");
    }

    public void disableSpecificCTermOnlyEnzyme() {
        paramMap.replace("specificCTermOnlyEnzyme ", false);
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

    public void disableEnzyme(Enzyme enzyme) {
        paramMap.replace(enzyme.getName(), false);
    }

}
