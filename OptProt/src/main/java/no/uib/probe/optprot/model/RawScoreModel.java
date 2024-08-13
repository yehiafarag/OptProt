/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.model;

import com.compomics.util.experiment.identification.matches.SpectrumMatch;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import no.uib.probe.optprot.util.ScoreComparison;
import no.uib.probe.optprot.util.SpectraUtilities;

/**
 *
 * @author yfa041
 */
public class RawScoreModel implements Comparable<RawScoreModel> {

    private int totalNumber = 0;
    private double tTestStat;
    private double pValue;
    private boolean sameData;
    private int fullSpectraSize;

    public boolean isSameData() {
        return sameData;
    }

    public void setSameData(boolean sameData) {
        this.sameData = sameData;
    }

    public void setSensitiveChange(boolean sensitiveChange) {
        this.sensitiveChange = sensitiveChange;
    }
    private boolean significatChange;
    private boolean sensitiveChange;

    public boolean isSensitiveChange() {
        return sensitiveChange;
    }
    private List<SpectrumMatch> spectrumMatchResult;
    private double improvmentScore;
    private double finalScore;
    private double sizeEffect;
    private double dataLengthFactor;

    public double getSizeEffect() {
        return sizeEffect;
    }

    public void setSizeEffect(double sizeEffect) {
        this.sizeEffect = Math.round(sizeEffect * 100.0) / 100.0;
    }

    public double getFinalScore() {
        return finalScore;
    }

    public void setFinalScore(double finalScore) {
        this.finalScore = Math.round(finalScore * 1000.0) / 1000.0;
    }
    private Set<String> specTitles;

    public Set<String> getSpecTitles() {
        if (specTitles == null && spectrumMatchResult != null) {
            specTitles = new HashSet<>();
            for (SpectrumMatch sm : spectrumMatchResult) {
                specTitles.add(sm.getSpectrumTitle());
            }
        }
        return specTitles;
    }

    public double getImprovmentScore() {
        return improvmentScore;
    }

    public void setImprovmentScore(double improvmentScore) {
        this.improvmentScore = improvmentScore;
    }

    public List<SpectrumMatch> getSpectrumMatchResult() {
        return spectrumMatchResult;
    }

    public void setSpectrumMatchResult(List<SpectrumMatch> spectrumMatchResult) {
        this.spectrumMatchResult = spectrumMatchResult;
        this.totalNumber = spectrumMatchResult.size();
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
        if (significatChange) {
            sensitiveChange = true;
        }

    }

    public RawScoreModel() {
        this.sc = new ScoreComparison();
    }
    private final ScoreComparison sc;

    @Override
    public int compareTo(RawScoreModel rs) {
//        double compscore = SpectraUtilities.isBetterScore(spectrumMatchResult, rs.spectrumMatchResult, fullSpectraSize);
//        if (compscore == 0) {
//            return 0;
//        } else if (compscore > 0) {
//            return 1;
//        } else {
//            return -1;
//        }
        return Double.compare(finalScore, rs.getFinalScore());

    }

    @Override
    public String toString() {
        return "Param accepted: " + significatChange + "  final score: " + finalScore + " improvmentScore: " + improvmentScore + " tstat" + tTestStat + "  pvalue " + pValue + "  #Spectra: " + totalNumber + "  senstive improvment " + sensitiveChange + "  same data " + sameData + "  size effect " + sizeEffect;
    }

    public double getDataLengthFactor() {
        return dataLengthFactor;
    }

    public void setDataLengthFactor(double dataLengthFactor) {
        this.dataLengthFactor = dataLengthFactor;
    }

    public void setFullSpectraSize(int fullSpectraSize) {
        this.fullSpectraSize = fullSpectraSize;
    }

}
