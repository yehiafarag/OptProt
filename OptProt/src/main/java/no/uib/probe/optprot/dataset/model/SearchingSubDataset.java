/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.dataset.model;

import java.io.File;

/**
 *
 * @author yfa041
 */
public class SearchingSubDataset {

    public File getSubDataFolder() {
        return subDataFolder;
    }

    public void setSubDataFolder(File subDataFolder) {
        this.subDataFolder = subDataFolder;
    }
    private File subDataFolder;
    private File subMsFile;
    private File subFastaFile;
    private File searchSettingsFile;

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
    
    private double[] validatedIdRefrenceData;

    public double[] getValidatedIdRefrenceData() {
        return validatedIdRefrenceData;
    }

    public void setValidatedIdRefrenceData(double[] validatedIdRefrenceData) {
        this.validatedIdRefrenceData = validatedIdRefrenceData;
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
            acceptedIDRatioThreshold=Math.round(acceptedIDRatioThreshold);

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

    public synchronized void setActiveIdentificationNum(int activeIdentificationNum) {
        this.activeIdentificationNum = activeIdentificationNum;
    }
    private int oreginalDatasize;
    private int userReferenceIdentificationNum;
    private int activeIdentificationNum;

    public int getOreginalDatasize() {
        return oreginalDatasize;
    }

    public void setOreginalDatasize(int oreginalDatasize) {
        this.oreginalDatasize = oreginalDatasize;
    }

//    public double getDataEpsilon() {
//        double d =(double)(defaultSettingIdentificationNum *100.0) / (double)(totalSpectraNumber);
//        double r= totalSpectraNumber*100.0/oreginalDatasize;
//        System.out.println("r: "+r+"  "+totalSpectraNumber+"  "+oreginalDatasize+"  --   "+defaultSettingIdentificationNum);
//        d=d/r;
////        d= 1.0/d;
////        d=0.000001;
////        d=100.0/totalSpectraNumber;
//       return d;//(double)totalSpectraNumber  / (double) oreginalDatasize;
//    }
    public int getTotalSpectraNumber() {
        return totalSpectraNumber;
    }

    public void setTotalSpectraNumber(int totalSpectraNumber) {
        this.totalSpectraNumber = totalSpectraNumber;
    }

    public synchronized double getIdentificationRate() {
        return activeIdentificationNum * 100.0 / totalSpectraNumber;
    }

    public int getDefaultSettingIdentificationNum() {
        return defaultSettingIdentificationNum;
    }

    public synchronized int getActiveIdentificationNum() {
        return activeIdentificationNum;
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

    public int getUserReferenceIdentificationNum() {
        return userReferenceIdentificationNum;
    }

    public void setUserReferenceIdentificationNum(int userReferenceIdentificationNum) {
        this.userReferenceIdentificationNum = userReferenceIdentificationNum;
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

}
