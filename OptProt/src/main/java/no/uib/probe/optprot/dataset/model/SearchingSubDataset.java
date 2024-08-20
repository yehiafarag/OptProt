/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.dataset.model;

import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import no.uib.probe.optprot.model.ParameterScoreModel;
import no.uib.probe.optprot.model.RawScoreModel;
import org.jgrapht.util.DoublyLinkedList;

/**
 *
 * @author yfa041
 */
public class SearchingSubDataset {

    private Map<String, TreeSet<ParameterScoreModel>> parameterScoreMap;

    public Map<String, TreeSet<ParameterScoreModel>> getParameterScoreMap() {
        return parameterScoreMap;
    }

    public void setParameterScoreMap(Map<String, TreeSet<ParameterScoreModel>> parameterScoreMap) {
        this.parameterScoreMap = parameterScoreMap;
    }

    public int getIdentificationNum() {
        return identificationNum;
    }

    public void setIdentificationNum(int identificationNum) {
        this.identificationNum = identificationNum;
    }

    public File getSubDataFolder() {
        return subDataFolder;
    }

    public void setSubDataFolder(File subDataFolder) {
        this.subDataFolder = subDataFolder;
    }
    private File subDataFolder;
    private File subMsFile;
    private File subFastaFile;
    private File oreginalFastaFile;
    private boolean fullDataSpectaInput;
    private Set<String> potintialVariableMod;
    private final Map<String, Double> fullSpectraScore = new HashMap<>();

    public void setSpectraTitiles(String[] titiles) {
        for (String titile : titiles) {
            fullSpectraScore.put(titile, 0.0);
        }
    }

    public boolean isFullDataSpectaInput() {
        return fullDataSpectaInput;
    }

    public void setFullDataSpectaInput(boolean fullDataSpectaInput) {
        this.fullDataSpectaInput = fullDataSpectaInput;
    }

    public File getOreginalFastaFile() {
        return oreginalFastaFile;
    }

    public void setOreginalFastaFile(File oreginalFastaFile) {
        this.oreginalFastaFile = oreginalFastaFile;
    }
    private File searchSettingsFile;
    private double pValueThresholds = -1;

    public File getSearchSettingsFile() {
        return searchSettingsFile;
    }

    public void setSearchSettingsFile(File searchSettingsFile) {
        this.searchSettingsFile = searchSettingsFile;
    }
    private int defaultSettingIdentificationNum;
    private int totalSpectraNumber;
    private int tempIdentificationNum;
    private double processDelay;
    private boolean highResolutionMassSpectrometers = true;
    private double acceptedIDRatioThreshold = -1.0;
    private RawScoreModel currentScoreModel;

    public void setActiveScoreModel(RawScoreModel scoreModel) {
        this.currentScoreModel = scoreModel;
        this.updateValidatedIdRefrenceData(currentScoreModel.getSpectrumMatchResult());
    }

    public void updateValidatedIdRefrenceData(List<SpectrumMatch> validatedIdRefrenceData) {
        this.resetSpectraScoreMap();
        this.identificationNum = validatedIdRefrenceData.size();
        for (SpectrumMatch sm : validatedIdRefrenceData) {
            fullSpectraScore.replace(sm.getSpectrumTitle(), sm.getBestPeptideAssumption().getRawScore());
        }

    }

    private String datasetId;

    public void setAcceptedIDRatioThreshold(double acceptedIDRatioThreshold) {
        this.acceptedIDRatioThreshold = acceptedIDRatioThreshold;
    }

    public double getAcceptedIDRatioThreshold() {
        if (acceptedIDRatioThreshold == -1.0) {
            double d = this.getActiveIdentificationNum() * 100.0 / this.getTotalSpectraNumber();
            if (d >= 20) {
                acceptedIDRatioThreshold = 5;
            } else if (d <= 4.0) {
                acceptedIDRatioThreshold = 1;
            } else {
                acceptedIDRatioThreshold = d * 5 / 20;
            }
            acceptedIDRatioThreshold = Math.round(acceptedIDRatioThreshold);

        }
        return acceptedIDRatioThreshold;
    }

    public double getProcessDelay() {
        return processDelay;
    }

    public void setProcessDelay(double processDelay) {
        this.processDelay = processDelay;
    }

    public int getTempIdentificationNum() {
        return tempIdentificationNum;
    }

    public void setTempIdentificationNum(int tempIdentificationNum) {
        this.tempIdentificationNum = tempIdentificationNum;
    }

    private int oreginalDatasize;
    private int identificationNum;

    private double comparisonsThreshold = 0.0;

    public int getOreginalDatasize() {
        return oreginalDatasize;
    }

    public void setOreginalDatasize(int oreginalDatasize) {
        this.oreginalDatasize = oreginalDatasize;
    }

    public int getTotalSpectraNumber() {
        return totalSpectraNumber;
    }

    public void setTotalSpectraNumber(int totalSpectraNumber) {
        this.totalSpectraNumber = totalSpectraNumber;
    }

    public synchronized double getIdentificationRate() {
        return identificationNum * 100.0 / totalSpectraNumber;
    }

    public int getDefaultSettingIdentificationNum() {
        return defaultSettingIdentificationNum;
    }

    public synchronized int getActiveIdentificationNum() {
        return identificationNum;
    }

    public void setDefaultSettingIdentificationNum(int defaultSettingIdentificationNum) {
        this.defaultSettingIdentificationNum = defaultSettingIdentificationNum;
    }

    public File getSubMsFile() {
        return subMsFile;
    }

    public void setSubMsFile(File subMsFile) {
        this.subMsFile = subMsFile;
    }

    public File getSubFastaFile() {
        return subFastaFile;
    }

    public void setSubFastaFile(File subFastaFile) {
        this.subFastaFile = subFastaFile;
    }

    public boolean isHighResolutionMassSpectrometers() {
        return highResolutionMassSpectrometers;
    }

    public void setHighResolutionMassSpectrometers(boolean highResolutionMassSpectrometers) {
        this.highResolutionMassSpectrometers = highResolutionMassSpectrometers;
    }

    public String getDatasetId() {
        return datasetId;
    }

    public void setDatasetId(String datasetId) {
        this.datasetId = datasetId;
    }

    public double getpValueThresholds() {
        if (pValueThresholds == -1) {
            double res = this.identificationNum * 100 / getTotalSpectraNumber();
            if (res <= 5) {
                pValueThresholds = 0.1;
            } else {
                pValueThresholds = 0.05;
            }
        }

        return pValueThresholds;
    }

    public Set<String> getPotintialVariableMod() {
        return potintialVariableMod;
    }

    public void setPotintialVariableMod(Set<String> potintialVariableMod) {
        this.potintialVariableMod = potintialVariableMod;
    }

    public double getComparisonsThreshold() {
        return comparisonsThreshold;
    }

    public void setComparisonsThreshold(double comparisonsThreshold) {
        this.comparisonsThreshold = comparisonsThreshold;
    }

    public RawScoreModel getCurrentScoreModel() {
        return currentScoreModel;
    }

    public void setCurrentScoreModel(RawScoreModel currentScoreModel) {
        this.currentScoreModel = currentScoreModel;
    }

    private void resetSpectraScoreMap() {
        for (String titile : fullSpectraScore.keySet()) {
            fullSpectraScore.replace(titile, 0.0);
        }

    }

    public Map<String, Double> getFullSpectraScore() {
        return fullSpectraScore;
    }

}
