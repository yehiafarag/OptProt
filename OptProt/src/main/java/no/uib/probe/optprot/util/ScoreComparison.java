/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.util;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.TTest;

public class ScoreComparison {

    public double calculateMean(double[] scores) {
        double sum = 0;
        for (double score : scores) {
            sum += score;
        }
        return sum / scores.length;
    }

    public double percentageImprovement(double mean1, double mean2) {
        return ((mean2 - mean1) / mean1) * 100.0;
    }

    public double calculateStandardDeviation(double[] scores) {
        double mean = calculateMean(scores);
        double sum = 0;
        for (double score : scores) {
            sum += Math.pow(score - mean, 2);
        }
        return Math.sqrt(sum / scores.length);
    }

    public double performTStatTest(double[] s1, double[] s2, boolean switchData) {
        double[] scores1, scores2;
        if (switchData) {
            scores1 = s2;
            scores2 = s1;
        } else {
            scores1 = s1;
            scores2 = s2;
        }

        if (scores1.length < 2 || scores2.length < 2) {
            if (scores1.length > scores2.length) {
                return -1;
            } else {
                return 1;
            }
        }
        TTest tTest = new TTest();
        if (scores1.length == scores2.length) {
            return tTest.pairedT(scores2, scores1);
        } else {
            return tTest.t(scores2, scores1);
        }
    }

    public double dataSizeEffect(double dataSizeFrom, double dataSizeTo) {
        dataSizeFrom++;
        dataSizeTo++;
        double effect = (dataSizeTo - dataSizeFrom) * 100.0 / Math.min(dataSizeFrom, dataSizeTo);

// 10.0 - (Math.min(dataSizeFrom, dataSizeTo) * 10.0 / Math.max(dataSizeFrom, dataSizeTo));
//        System.out.println("effect is " + dataSizeFrom + "  " + dataSizeTo + "  " + effect);
//        if (dataSizeFrom > dataSizeTo) {
//            effect *= -1.0;
//        }
        return logScaleNormalize(effect, Math.E);

    }

    public double performTTest(double[] s1, double[] s2, boolean switchData) {

        double[] scores1, scores2;
        if (switchData) {
            scores1 = s2;
            scores2 = s1;
        } else {
            scores1 = s1;
            scores2 = s2;
        }
        if (scores1.length < 2 || scores2.length < 2) {
            return 1;
        }
        TTest tTest = new TTest();
        return tTest.tTest(scores2, scores1);
    }

    public double normalize(double value, double min, double max) {
        return (value - min) / (max - min);
    }

    /**
     * Applies log scale normalization to an array of data.
     *
     * @param data The array of data to be normalized.
     * @param base The base of the logarithm to use.
     * @return The log scale normalized array of data.
     */
    public static double logScaleNormalize(double data, double base) {
        double sign = 1;
        if (data < 0) {
            sign = -1;
        }
        data = Math.abs(data);
        double normalizedData = (Math.log(data + 1) / Math.log(base)) * sign;
        if (Double.isNaN(normalizedData)) {
            normalizedData = 0;
        }
        return normalizedData;
    }

    public double calculateScore(double[] scores1, double[] scores2) {

        DescriptiveStatistics ds1;
        DescriptiveStatistics ds2 ;
        
         boolean switchData = false;
        if (scores2.length >= scores1.length) {
            ds1 = new DescriptiveStatistics(scores1);
            ds2 = new DescriptiveStatistics(scores2);
        } else {
            switchData = true;
            ds1 = new DescriptiveStatistics(scores2);
            ds2 = new DescriptiveStatistics(scores1);
        }

        double median1 = ds1.getPercentile(50);
        double median2 = ds2.getPercentile(50);

        double meanDifference = median2 - median1;
        double percentageImprovement = percentageImprovement(median1, median2);
        double tStat = performTStatTest(scores1, scores2,switchData);

        double stdDev1 = ds1.getStandardDeviation();//calculateStandardDeviation(scores1);
        double stdDev2 = ds2.getStandardDeviation();//calculateStandardDeviation(scores2);

        double stdDevDifference = stdDev1 - stdDev2;
        double pValue = performTTest(scores1, scores2,switchData);

        // Assuming the possible range for these metrics is between -200 and 200 for normalization purposes.
        // Adjust these ranges based on your specific use case.
        double normMeanDiff = logScaleNormalize(meanDifference, Math.E);
        double normPercentageImprovement = logScaleNormalize(percentageImprovement, Math.E);
        double normStdDevDiff = logScaleNormalize(stdDevDifference, Math.E);
        double normTStat = logScaleNormalize(tStat, Math.E);

        double normPValue = logScaleNormalize((1.0 - pValue), Math.E); // Inverting p-value for interpretation
        if (tStat < 0 && normPValue > 0) {
            normPValue *= -1.0;
        }
        double finalScore = (normMeanDiff + normPercentageImprovement + normStdDevDiff + normPValue + normTStat) / 5.0;

        if (percentageImprovement < 0 && finalScore > 0) {
//            finalScore *= -1.0;
            System.out.println("-----------------------------------2------------------------------>>> flip score case");
        }
        return finalScore;
    }

    public double[] calculateFinalScore(double[] scores1, double[] scores2) {

        DescriptiveStatistics ds1;
        DescriptiveStatistics ds2;
        boolean switchData = false;
        if (scores2.length >= scores1.length) {
            ds1 = new DescriptiveStatistics(scores1);
            ds2 = new DescriptiveStatistics(scores2);
        } else {
            switchData = true;
            ds1 = new DescriptiveStatistics(scores2);
            ds2 = new DescriptiveStatistics(scores1);
        }

        double mean1 = ds1.getMean();
        double mean2 = ds2.getMean();

        double meanDifference = mean2 - mean1;
        double percentageImprovement = percentageImprovement(mean1, mean2);

        if (scores1.length == scores2.length && mean1 == mean2 && percentageImprovement == 0) {
            return new double[]{0, 0, 0, 0, 0, 0};
        }

        double tStat = performTStatTest(scores1, scores2, switchData);

        double stdDev1 = ds1.getStandardDeviation();//calculateStandardDeviation(scores1);
        double stdDev2 = ds2.getStandardDeviation();//calculateStandardDeviation(scores2);

        double stdDevDifference = stdDev1 - stdDev2;
        double pValue = performTTest(scores1, scores2,switchData);

        // Assuming the possible range for these metrics is between -200 and 200 for normalization purposes.
        // Adjust these ranges based on your specific use case.
        double normMeanDiff = logScaleNormalize(meanDifference, Math.E);
        double normPercentageImprovement = logScaleNormalize(percentageImprovement, Math.E);
        double normStdDevDiff = logScaleNormalize(stdDevDifference, Math.E);
        double normTStat = logScaleNormalize(tStat, Math.E);
        double normPValue = logScaleNormalize((1.0 - pValue), Math.E); // Inverting p-value for interpretation
        if (tStat < 0) {
            normPValue *= -1.0;
        }
        double dataSizeEffect = 0;
        double finalScore;
        if (Double.isInfinite(Math.abs(normTStat))) {
            normTStat = normPercentageImprovement;
        }
        if (scores1.length == scores2.length) {
            finalScore = (normMeanDiff + normPercentageImprovement + normStdDevDiff + normPValue + normTStat) / 5.0;
        } else {
            dataSizeEffect = calculateCohensD(scores1, scores2,switchData);//
            if (Double.isNaN(dataSizeEffect) || Double.isInfinite(dataSizeEffect)) {
                dataSizeEffect = dataSizeEffect(scores1.length, scores2.length);
            }
            if (ds1.getN()> ds2.getN() && dataSizeEffect > 0) {
                dataSizeEffect *= -1.0;
            }
            finalScore = (normMeanDiff + normPercentageImprovement + normStdDevDiff + normPValue + normTStat + dataSizeEffect) / 6.0;
        }
        System.out.println("final score parts (meanDifference: " + scores1.length + " to " + scores2.length + " MD " + normMeanDiff + " , percentageImprovement: " + normPercentageImprovement + " , tStat: " + normTStat + " , stdDevDifference: " + normStdDevDiff + " , pValue: " + normPValue + " , dataSizeEffect: " + dataSizeEffect + ") " + finalScore + "  " + Double.isInfinite(Math.abs(normTStat)));

        if (percentageImprovement < 0 && finalScore > 0) {
//            finalScore *= -1.0;
            System.out.println("----------------------------------------------------------------->>> flip score case");
        }
        if (switchData) {
            System.out.println("data switch happened");
            finalScore *= -1.0;
        }

        return new double[]{percentageImprovement, tStat, pValue, dataSizeEffect, finalScore, 0};
    }

    public double[] calculateFinalScore(Double[] scores1, Double[] scores2) {
        double[] s1 = new double[scores1.length];
        for (int i = 0; i < s1.length; i++) {
            s1[i] = scores1[i];
        }
        double[] s2 = new double[scores2.length];
        for (int i = 0; i < s2.length; i++) {
            s2[i] = scores2[i];
        }
        return this.calculateFinalScore(s1, s2);
    }

    public static void main(String[] args) {
        double[] scores2 = {91.5, 70, 25, 12};
        double[] scores1 = {91.0, 89.0, 95.0, 87.0, 93.5, 90.0};

        ScoreComparison sc = new ScoreComparison();
        double finalScore = sc.calculateFinalScore(scores1, scores2)[4];

        System.out.println("Final Score: " + finalScore + "  " + sc.dataSizeEffect(scores1.length, scores2.length));

//         Interpretation
        if (finalScore > 0.75) {
            System.out.println("Significant positive enhancement in scores.");
        } else if (finalScore > 0.5) {
            System.out.println("Moderate positive enhancement in scores.");
        } else if (finalScore > 0.25) {
            System.out.println("Slight positive enhancement in scores.");
        } else {
            System.out.println("No significant enhancement or decline in scores.");
        }

    }

    public double calculateCohensD(double[] s1, double[] s2,boolean switchData) {
        double[] from,  to;
        if (switchData) {
            from = s2;
            to = s1;
        } else {
            from = s1;
            to = s2;
        }
        double mean1 = Arrays.stream(from).average().orElse(0.0);
        double mean2 = Arrays.stream(to).average().orElse(0.0);
        double variance1 = Arrays.stream(from).map(x -> Math.pow(x - mean1, 2)).sum() / (from.length - 1);
        double variance2 = Arrays.stream(to).map(x -> Math.pow(x - mean2, 2)).sum() / (to.length - 1);
        double pooledStdDev = Math.sqrt(((from.length - 1) * variance1 + (to.length - 1) * variance2) / (from.length + to.length - 2));
        double normCohenD = ((mean2 - mean1) / pooledStdDev);
        return logScaleNormalize(normCohenD, Math.E);
    }

    public List<Double> getElementsLargerThanParallel(double[] array, double specificValue) {
        // Use parallel streams to filter and collect the elements
        return Arrays.stream(array)
                .parallel()
                .filter(value -> value > specificValue)
                .boxed()
                .collect(Collectors.toList());
    }

    public double calculateRankBiserialCorrelation(double[] from, double[] to) {
        int n1 = from.length;
        int n2 = to.length;
        double[] combined = new double[n1 + n2];
        System.arraycopy(from, 0, combined, 0, n1);
        System.arraycopy(to, 0, combined, n1, n2);
        Arrays.sort(combined);

        double rankSum1 = 0;
        for (double value : from) {
            rankSum1 += Arrays.binarySearch(combined, value) + 1;  // Adding 1 for rank (1-based index)
        }

        double u1 = rankSum1 - n1 * (n1 + 1) / 2.0;
        double u2 = n1 * n2 - u1;
        double u = Math.min(u1, u2);
        double value = 1 - (2 * u) / (n1 * n2);
        if (from.length > to.length) {
            value *= -1.0;
        }
        return value;
    }

    public List<Double> removeOutliers(DescriptiveStatistics stats, double zScoreThreshold) {
        // Compute mean and standard deviation

        double mean = stats.getMean();
        double standardDeviation = stats.getStandardDeviation();

        // Filter out outliers
        return Arrays.stream(stats.getSortedValues())
                .filter(value -> Math.abs((value - mean) / standardDeviation) <= zScoreThreshold)
                .boxed()
                .collect(Collectors.toList());
    }

}
