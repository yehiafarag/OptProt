/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.model;

import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import java.util.List;

/**
 *
 * @author yfa041
 */
public class RawScoreModel implements Comparable<RawScoreModel> {

    private int totalNumber = 0;
    private double tTestStat;
    private double pValue;
    private boolean significatChange;
    private List<SpectrumMatch> spectrumMatchResult;
    private double[] data;
    private boolean restrictedComparison;
    private double improvmentScore;

    public double getImprovmentScore() {
        return improvmentScore;
    }

    public void setImprovmentScore(double improvmentScore) {
        this.improvmentScore = improvmentScore;
    }

    public void setRestrictedComparison(boolean restrictedComparison) {
        this.restrictedComparison = restrictedComparison;
    }

    public double[] getData() {
        return data;
    }

    public void setData(double[] data) {
        this.data = data;
    }

    public List<SpectrumMatch> getSpectrumMatchResult() {
        return spectrumMatchResult;
    }

    public void setSpectrumMatchResult(List<SpectrumMatch> spectrumMatchResult) {
        this.spectrumMatchResult = spectrumMatchResult;
    }

    public int getTotalNumber() {
        return totalNumber;
    }

    public void setTotalNumber(int totalNumber) {
        this.totalNumber = totalNumber;
    }

    public double gettTestStat() {
        return tTestStat;
    }

    public void settTestStat(double tTestStat) {
        this.tTestStat = tTestStat;
    }

    public double getpValue() {
        return pValue;
    }

    public void setpValue(double pValue) {
        this.pValue = pValue;
    }

    public boolean isSignificatChange() {
        return significatChange;
    }

    public void setSignificatChange(boolean significatChange) {
        this.significatChange = significatChange;
    }

    @Override
    public int compareTo(RawScoreModel rs) {
        return Double.compare(improvmentScore, rs.getImprovmentScore());
//        int sigComparison = Boolean.compare(this.isSignificatChange(), rs.isSignificatChange());
//        if (sigComparison != 0) {
//            return sigComparison;
//        }
//        int tstatComparison = Double.compare(this.tTestStat, rs.tTestStat);
//        if (tstatComparison != 0) {
//            return tstatComparison;
//        }
//        System.out.println("up to total number "+tTestStat+"  "+rs.tTestStat);
//        return Integer.compare(this.totalNumber, rs.getTotalNumber());
        }

    @Override
    public String toString() {
        return "Param accepted: " + significatChange + " improvmentScore: " + improvmentScore + "  #Spectra: " + totalNumber;
    }

}
