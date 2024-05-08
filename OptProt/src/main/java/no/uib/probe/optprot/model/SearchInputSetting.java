/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.model;

import com.compomics.util.experiment.identification.Advocate;
import java.io.File;
import no.uib.probe.optprot.configurations.Configurations;
import no.uib.probe.optprot.search.myrimatch.MyriMatchEnabledParameters;
import no.uib.probe.optprot.search.xtandam.OptProtXtandemAdvancedSearchParameter;
import no.uib.probe.optprot.search.xtandam.XTandemEnabledParameters;

/**
 *
 * @author yfa041
 */
public class SearchInputSetting {

    public SearchInputSetting() {
    }

    private boolean runOmssa = false;
    private boolean runXTandem = false;
    private boolean runMsgf = false;
    private boolean runMsAmanda = false;
    private boolean runMyriMatch = false;
    private boolean runComet = false;
    private boolean runTide = false;
    private boolean runAndromeda = false;
    private boolean runMetaMorpheus = false;
    private boolean runSage = false;
    private boolean runNovor = false;
    private boolean runDirecTag = false;
    private File omssaFolder = null;
    private final File xTandemFolder = new File(Configurations.XTANDEM_FOLDER);
    private File msgfFolder = null;
    private File msAmandaFolder = null;
    private final File myriMatchFolder = new File(Configurations.MYRIMATCH_FOLDER);
    private final File cometFolder = new File(Configurations.COMET_FOLDER);
    private File tideFolder = null;
    private File tideIndexLocation = null;
    private File andromedaFolder = null;
    private File metaMorpheusFolder = null;
    private File sageFolder;
    private final File novorFolder = new File(Configurations.NOVOR_FOLDER);
    private final File direcTagFolder = new File(Configurations.DIRECTAG_FOLDER);
    private File makeblastdbFolder = null;

    //parameter to optimise
    private boolean optimizeDigestionParameter;
    private boolean optimizeEnzymeParameter;
    private boolean optimizeSpecificityParameter;
    private boolean optimizeMaxMissCleavagesParameter;
    private boolean optimizeFragmentIonTypesParameter;
    private boolean optimizePrecursorToleranceParameter;
    private boolean optimizeFragmentToleranceParameter;
    private boolean optimizeIsotopsParameter;
    private boolean optimizeModificationParameter;
    private boolean recalibrateSpectraParameter;
    private boolean optimizeXtandemAdvancedParameter;
    private Advocate selectedSearchEngine = Advocate.xtandem;

    public Advocate getSelectedSearchEngine() {
        return selectedSearchEngine;
    }

    public void setSelectedSearchEngine(Advocate searchEngine) {
        this.selectedSearchEngine = searchEngine;
        if (searchEngine == Advocate.xtandem) {
            this.setRunXTandem(true);
        } else if (searchEngine == Advocate.myriMatch) {
            this.setRunMyriMatch(true);
        } else if (searchEngine == Advocate.direcTag) {
            this.setRunDirecTag(true);
        }

    }

    private final OptProtXtandemAdvancedSearchParameter xtandemOptProtAdvancedSearchParameters = new OptProtXtandemAdvancedSearchParameter();

    public boolean isOptimizeXtandemAdvancedParameter() {
        return optimizeXtandemAdvancedParameter;
    }

    public void setOptimizeXtandemAdvancedParameter(boolean optimizeXtandemAdvancedParameter) {
        this.optimizeXtandemAdvancedParameter = optimizeXtandemAdvancedParameter;
    }

    public boolean isRecalibrateSpectraParameter() {
        return recalibrateSpectraParameter;
    }

    public void setRecalibrateSpectraParameter(boolean recalibrateSpectraParameter) {
        this.recalibrateSpectraParameter = recalibrateSpectraParameter;
    }

    private XTandemEnabledParameters XTandemEnabledParameters = new XTandemEnabledParameters();
    private MyriMatchEnabledParameters MyriMatchEnabledParameters = new MyriMatchEnabledParameters();

    public XTandemEnabledParameters getXTandemEnabledParameters() {
        return XTandemEnabledParameters;
    }

    public boolean isOptimizeModificationParameter() {
        return optimizeModificationParameter;
    }

    public void setOptimizeModificationParameter(boolean optimizeModificationParameter) {
        this.optimizeModificationParameter = optimizeModificationParameter;
    }
    private boolean optimizePrecursorChargeParameter;

    public boolean isOptimizeFragmentToleranceParameter() {
        return optimizeFragmentToleranceParameter;
    }

    public void setOptimizeFragmentToleranceParameter(boolean optimizeFragmentToleranceParameter) {
        this.optimizeFragmentToleranceParameter = optimizeFragmentToleranceParameter;
    }

    public boolean isOptimizePrecursorToleranceParameter() {
        return optimizePrecursorToleranceParameter;
    }

    public void setOptimizePrecursorToleranceParameter(boolean optimizePrecursorToleranceParameter) {
        this.optimizePrecursorToleranceParameter = optimizePrecursorToleranceParameter;
    }

    public boolean isOptimizeDigestionParameter() {
        return optimizeDigestionParameter;
    }

    public void setOptimizeDigestionParameter(boolean optimizeDigestionParameter) {
        this.optimizeDigestionParameter = optimizeDigestionParameter;
    }

    public boolean isOptimizeEnzymeParameter() {
        return optimizeEnzymeParameter;
    }

    public void setOptimizeEnzymeParameter(boolean optimizeEnzymeParameter) {
        this.optimizeEnzymeParameter = optimizeEnzymeParameter;
    }

    public boolean isOptimizeSpecificityParameter() {
        return optimizeSpecificityParameter;
    }

    public void setOptimizeSpecificityParameter(boolean optimizeSpecificityParameter) {
        this.optimizeSpecificityParameter = optimizeSpecificityParameter;
    }

    public boolean isOptimizeMaxMissCleavagesParameter() {
        return optimizeMaxMissCleavagesParameter;
    }

    public void setOptimizeMaxMissCleavagesParameter(boolean optimizeMaxMissCleavagesParameter) {
        this.optimizeMaxMissCleavagesParameter = optimizeMaxMissCleavagesParameter;
    }

    public boolean isOptimizeFragmentIonTypesParameter() {
        return optimizeFragmentIonTypesParameter;
    }

    public void setOptimizeFragmentIonTypesParameter(boolean optimizeFragmentIonTypesParameter) {
        this.optimizeFragmentIonTypesParameter = optimizeFragmentIonTypesParameter;
    }

    public boolean isRunOmssa() {
        return runOmssa;
    }

    public void setRunOmssa(boolean runOmssa) {
        resetActiveSearchEngines();
        this.runOmssa = runOmssa;
    }

    public boolean isRunXTandem() {
        return runXTandem;
    }

    public void setRunXTandem(boolean runXTandem) {
        resetActiveSearchEngines();
        this.runXTandem = runXTandem;
    }

    public boolean isRunMsgf() {
        return runMsgf;
    }

    public void setRunMsgf(boolean runMsgf) {
        resetActiveSearchEngines();
        this.runMsgf = runMsgf;
    }

    public boolean isRunMsAmanda() {
        return runMsAmanda;
    }

    public void setRunMsAmanda(boolean runMsAmanda) {
        resetActiveSearchEngines();
        this.runMsAmanda = runMsAmanda;
    }

    public boolean isRunMyriMatch() {
        return runMyriMatch;
    }

    public void setRunMyriMatch(boolean runMyriMatch) {
        resetActiveSearchEngines();
        this.runMyriMatch = runMyriMatch;
    }

    public boolean isRunComet() {
        return runComet;
    }

    public void setRunComet(boolean runComet) {
        resetActiveSearchEngines();
        this.runComet = runComet;
    }

    public boolean isRunTide() {
        return runTide;
    }

    public void setRunTide(boolean runTide) {
        resetActiveSearchEngines();
        this.runTide = runTide;
    }

    public boolean isRunAndromeda() {
        return runAndromeda;
    }

    public void setRunAndromeda(boolean runAndromeda) {
        this.runAndromeda = runAndromeda;
    }

    public boolean isRunMetaMorpheus() {
        return runMetaMorpheus;
    }

    public void setRunMetaMorpheus(boolean runMetaMorpheus) {
        resetActiveSearchEngines();
        this.runMetaMorpheus = runMetaMorpheus;
    }

    public boolean isRunSage() {
        return runSage;
    }

    public void setRunSage(boolean runSage) {
        resetActiveSearchEngines();
        this.runSage = runSage;
    }

    public boolean isRunNovor() {
        return runNovor;
    }

    public void setRunNovor(boolean runNovor) {
        resetActiveSearchEngines();
        this.runNovor = runNovor;
    }

    public boolean isRunDirecTag() {
        return runDirecTag;
    }

    public void setRunDirecTag(boolean runDirecTag) {
        resetActiveSearchEngines();
        this.runDirecTag = runDirecTag;
    }

    public File getOmssaFolder() {
        return omssaFolder;
    }

    public void setOmssaFolder(File omssaFolder) {

        this.omssaFolder = omssaFolder;
    }

    public File getxTandemFolder() {
        return xTandemFolder;
    }

    public File getMsgfFolder() {
        return msgfFolder;
    }

    public void setMsgfFolder(File msgfFolder) {
        this.msgfFolder = msgfFolder;
    }

    public File getMsAmandaFolder() {
        return msAmandaFolder;
    }

    public void setMsAmandaFolder(File msAmandaFolder) {
        this.msAmandaFolder = msAmandaFolder;
    }

    public File getMyriMatchFolder() {
        return myriMatchFolder;
    }

    public File getCometFolder() {
        return cometFolder;
    }

    public File getTideFolder() {
        return tideFolder;
    }

    public void setTideFolder(File tideFolder) {
        this.tideFolder = tideFolder;
    }

    public File getTideIndexLocation() {
        return tideIndexLocation;
    }

    public void setTideIndexLocation(File tideIndexLocation) {
        this.tideIndexLocation = tideIndexLocation;
    }

    public File getAndromedaFolder() {
        return andromedaFolder;
    }

    public void setAndromedaFolder(File andromedaFolder) {
        this.andromedaFolder = andromedaFolder;
    }

    public File getMetaMorpheusFolder() {
        return metaMorpheusFolder;
    }

    public void setMetaMorpheusFolder(File metaMorpheusFolder) {
        this.metaMorpheusFolder = metaMorpheusFolder;
    }

    public File getSageFolder() {
        return sageFolder;
    }

    public void setSageFolder(File sageFolder) {
        this.sageFolder = sageFolder;
    }

    public File getNovorFolder() {
        return novorFolder;
    }

    public File getDirecTagFolder() {
        return direcTagFolder;
    }

    public File getMakeblastdbFolder() {
        return makeblastdbFolder;
    }

    public void setMakeblastdbFolder(File makeblastdbFolder) {
        this.makeblastdbFolder = makeblastdbFolder;
    }

    private void resetActiveSearchEngines() {
        runOmssa = false;
        runXTandem = false;
        runMsgf = false;
        runMsAmanda = false;
        runMyriMatch = false;
        runComet = false;
        runTide = false;
        runAndromeda = false;
        runMetaMorpheus = false;
        runSage = false;
        runNovor = false;
        runDirecTag = false;

    }

    public boolean isOptimizePrecursorChargeParameter() {
        return optimizePrecursorChargeParameter;
    }

    public void setOptimizePrecursorChargeParameter(boolean optimizePrecursorChargeParameter) {
        this.optimizePrecursorChargeParameter = optimizePrecursorChargeParameter;
    }

    public boolean isOptimizeIsotopsParameter() {
        return optimizeIsotopsParameter;
    }

    public void setOptimizeIsotopsParameter(boolean optimizeIsotopsParameter) {
        this.optimizeIsotopsParameter = optimizeIsotopsParameter;
    }

    public OptProtXtandemAdvancedSearchParameter getXtandemOptProtAdvancedSearchParameters() {
        return xtandemOptProtAdvancedSearchParameters;
    }

    public MyriMatchEnabledParameters getMyriMatchEnabledParameters() {
        return MyriMatchEnabledParameters;
    }

}
