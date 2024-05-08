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

    private File subMsFile;
    private File subFastaFile;
    private int identificationNum;
    private int totalSpectraNumber;
    private int oreginalDatasize;

    public int getOreginalDatasize() {
        return oreginalDatasize;
    }

    public void setOreginalDatasize(int oreginalDatasize) {
        this.oreginalDatasize = oreginalDatasize;
    }

    public double getDataEpsilon() {
        System.out.println("no.uib.probe.optprot.dataset.model.SearchingSubDataset.getDataEpsilon() "+oreginalDatasize+"  /  "+totalSpectraNumber);
        return (double)totalSpectraNumber  / (double) oreginalDatasize;
    }

    public int getTotalSpectraNumber() {
        return totalSpectraNumber;
    }

    public void setTotalSpectraNumber(int totalSpectraNumber) {
        this.totalSpectraNumber = totalSpectraNumber;
    }

    public double getIdentificationRate() {
        return identificationNum * 100.0 / totalSpectraNumber;
    }

    public int getIdentificationNum() {
        return identificationNum;
    }

    public void setIdentificationNum(int identificationNum) {
        this.identificationNum = identificationNum;
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

}
