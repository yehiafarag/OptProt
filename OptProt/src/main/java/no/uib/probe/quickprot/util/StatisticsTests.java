/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.quickprot.util;

import java.util.Arrays;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.inference.TTest;

public class StatisticsTests {

    public static void main(String[] args) {
        modeMedian();
//        // Example data
//        double[] sample1 = {45.2, 46.1, 47.3, 50.0, 48.5};
//        double[] sample2 = {52.0, 49.5, 53.2, 51.3, 50.4};
//
//        double z = zTest(sample1, sample2);
//        double pValue = pValueForZTest(z);
//
//        System.out.printf("Z-test statistic: %.4f%n", z);
//        System.out.printf("P-value: %.4f%n", pValue);
//        StatisticsTests.ratioComparisons(args);
    }

    public static double independentZTest(DescriptiveStatistics referenceSamples, DescriptiveStatistics improvedSample) {
//        int n1 = sample1.length;
//        int n2 = sample2.length;
//
//        double mean1 = mean(sample1);
//        double mean2 = mean(sample2);
//
//        double std1 = stdDev(sample1);
//        double std2 = stdDev(sample2);
//
//        double z = (mean1 - mean2) / Math.sqrt((std1 * std1 / n1) + (std2 * std2 / n2));
//        return z;
        double mean1 = referenceSamples.getMean();
        double mean2 = improvedSample.getMean();
        double variance1 = referenceSamples.getVariance();
        double variance2 = improvedSample.getVariance();

        double standardError = Math.sqrt((variance2 / improvedSample.getN()) + (variance1 / referenceSamples.getN()));
        double zScore = (mean2 - mean1) / standardError;
        return zScore;

    }

    public static double pairedZTest(DescriptiveStatistics referenceSamples, DescriptiveStatistics improvedSample) {
        double meanDifference = calculatePairedMeanDifference(referenceSamples.getValues(), improvedSample.getValues());
        double standardDeviation = calculatePairedStandardDeviation(referenceSamples.getValues(), improvedSample.getValues(), meanDifference);
        double standardError = standardDeviation / Math.sqrt(improvedSample.getN());
        double zScore = meanDifference / standardError;
        return zScore;
    }

    private static double calculatePairedMeanDifference(double[] refrenceSample, double[] improvedSample) {
        double sum = 0;
        for (int i = 0; i < improvedSample.length; i++) {
            sum += improvedSample[i] - refrenceSample[i];
        }
        return sum / improvedSample.length;
    }

    private static double calculatePairedStandardDeviation(double[] refrenceSample, double[] improvedSample, double meanDifference) {
        double sum = 0;
        for (int i = 0; i < improvedSample.length; i++) {
            sum += Math.pow((improvedSample[i] - refrenceSample[i]) - meanDifference, 2);
        }
        return Math.sqrt(sum / (improvedSample.length - 1));
    }

    public static double mean(double[] sample) {
        return Arrays.stream(sample).average().orElse(Double.NaN);
    }

    public static double stdDev(double[] sample) {
        double mean = mean(sample);
        double sum = 0;
        for (double num : sample) {
            sum += Math.pow(num - mean, 2);
        }
        return Math.sqrt(sum / (sample.length - 1));
    }

    public static double pValueForZTest(double z) {
        return 2 * (1 - cumulativeProbability(Math.abs(z)));
    }

    public static double cumulativeProbability(double z) {
        return 0.5 * (1 + erf(z / Math.sqrt(2)));
    }

    // Approximation of the error function
    public static double erf(double x) {
        double t = 1.0 / (1.0 + 0.5 * Math.abs(x));
        double tau = t * Math.exp(-x * x - 1.26551223 + t * (1.00002368
                + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806
                + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587
                + t * (-0.82215223 + t * 0.17087277)))))))));
        return x >= 0 ? 1 - tau : tau - 1;
    }

    public static double pValueForTTest(double[] fromData, double[] toData) {
        if (toData.length < 2 || fromData.length < 2) {
            return 100.0;
        }
        TTest tTest = new TTest();
        double pValue;
        if (toData.length == fromData.length) {
            pValue = tTest.pairedTTest(fromData, toData);
        } else {
            pValue = tTest.tTest(fromData, toData);//          
        }

        return pValue;
    }

    public static double tTest(double[] fromData, double[] toData) {

        if (toData.length < 2 || fromData.length < 2) {
            return -100.0;
        }

        // Paired t-test
        TTest tTest = new TTest();
        double tStatistic;
        if (toData.length > fromData.length) {
            double[] updatedData = new double[toData.length];
            System.arraycopy(fromData, 0, updatedData, 0, fromData.length);
            fromData = updatedData;
        } else if (toData.length < fromData.length) {
            double[] updatedData = new double[fromData.length];
            System.arraycopy(toData, 0, updatedData, 0, toData.length);
            toData = updatedData;
        }
        tStatistic = tTest.pairedT(fromData, toData);
        return tStatistic;
    }

    public static boolean isSignificatChange(double tStatistic, double pValue, int data1Size, int data2Size) {
        if (tStatistic < 0 && pValue <= 0.05) {//
            System.out.println("-------------------------------------------->>> the t-test  " + tStatistic + "   pvalue " + pValue + " case 0  " + "FALSE");
            return false;
        }
        if ((Double.isNaN(tStatistic) || Double.isNaN(pValue)) && data1Size == data2Size) {
            System.out.println("-------------------------------------------->>> the t-test  " + tStatistic + "   pvalue " + pValue + " case 1  " + "FALSE");
            return false;
        }
        if (tStatistic > 0 && data1Size * 1.02 < data2Size) {
            System.out.println("-------------------------------------------->>> the t-test  " + tStatistic + "   pvalue " + pValue + "  case 2  FALSE");
            return false;
        }
        if (tStatistic < 0 && pValue > 0.05 && data1Size < data2Size * 1.1) {
            System.out.println("-------------------------------------------->>> the t-test  " + tStatistic + "   pvalue " + pValue + "  case 3  FALSE");
            return false;

        }

//        
//        
//       else  else if ((tStatistic < 0 && pValueForZTest <= 0.05) || (((Double.isNaN(tStatistic) || Double.isNaN(pValueForZTest)) && resultScores.length != referenceScore.length))) {
//            System.out.println("-------------------------------------------->>> the t-test  " + tStatistic + "   pvalue " + pValueForZTest + "  case 4  FALSE");
//            return false;
//        } else if (resultScores.length * 1.02 > referenceScore.length && tStatistic > 0) {
//            System.out.println("-------------------------------------------->>> the t-test  " + tStatistic + "   pvalue " + pValueForZTest + "  case 5  TRUE");
////            Configurations.VALIDATED_ID_REF_DATA = resultScores;
//            return true;
//        }
        System.out.println("-------------------------------------------->>> the t-test  " + tStatistic + "   pvalue " + pValue + "  case 6  TRUE");
        System.out.println("else new condition tStatistic " + tStatistic + " Double.isNaN(tStatistic) " + Double.isNaN(tStatistic) + "  resultScores.length " + data1Size + "   " + data2Size);
        return true;

    }

    public static void ratioComparisons(String[] args) {
        // Example dataset: Revenue, Cost, and Employees for different departments
        double[] scores1 = {1000, 2600, 1900, 150};     // Revenue for Marketing, Sales, HR, IT
        double[] scores2 = {750, 2250, 1750, 75};        // Corresponding costs
        int[] spectraNumInQuartiles = {2750, 2750, 2750, 2750};
        String[] quartilesIds = {"Q1", "Q2", "Q3", "Q4"};
        System.out.println("Comparing Ratios in Same Category:");

        for (int i = 0; i < scores1.length; i++) {
            // Calculate Revenue-to-Cost Ratio
            double idToTagRatio = scores1[i] / scores2[i];

            // Calculate Revenue-to-Employee Ratio
            double idToTotal = scores1[i] / spectraNumInQuartiles[i];

            // Print the ratios for each department
            System.out.printf("\nQuartile: %s\n", quartilesIds[i]);
            System.out.printf("Identified spectra-to-Confedent Tag Ratio: %.2f\n", idToTagRatio);
            System.out.printf("identified spectra to total Ratio: %.2f\n", idToTotal);

            // Comparison (e.g., Ratio of Ratios)
            double ratioOfRatios = idToTagRatio / idToTotal;
            System.out.printf("Ratio of Ratios (identified vs total): %.2f\n", ratioOfRatios);

            // Alternatively, you can find the absolute difference between the ratios
            double difference = idToTagRatio - idToTotal;
            System.out.printf("Difference Between Ratios: %.2f\n", difference);
        }
    }

    public static double calculatePearsonCorrelationTest(double[] values1, double[] values2) {
        PearsonsCorrelation test = new PearsonsCorrelation();
        double correlation = test.correlation(values1, values2);
        correlation = Math.floor(correlation * 100.0) / 100.0;
        System.out.println("correlation is " + correlation);
        return correlation;

    }

    public static void modeMedian() {
        // Sample data
        double[] sample1 = {1.83, 0.50, 1.62, 2.48, 2.01, 1.32};
        double[] sample2 = {0.88, 0.65, 0.80, 1.15, 1.01, 0.75};

        // Combine samples to find the overall median
        double[] combined = new double[sample1.length + sample2.length];
        System.arraycopy(sample1, 0, combined, 0, sample1.length);
        System.arraycopy(sample2, 0, combined, sample1.length, sample2.length);
        double median = findMedian(combined);

        // Create contingency table
        long[][] contingencyTable = new long[2][2];
        for (double v : sample1) {
            if (v > median) {
                contingencyTable[0][0]++;
            } else {
                contingencyTable[1][0]++;
            }
        }
        for (double v : sample2) {
            if (v > median) {
                contingencyTable[0][1]++;
            } else {
                contingencyTable[1][1]++;
            }
        }

        // Perform chi-square test
        ChiSquareTest chiSquareTest = new ChiSquareTest();
        double pValue = chiSquareTest.chiSquareTest(contingencyTable);

        // Output the result
        System.out.println("P-value: " + pValue + "  value " + chiSquareTest.chiSquare(contingencyTable));
    }

    private static double findMedian(double[] data) {
        java.util.Arrays.sort(data);
        int middle = data.length / 2;
        if (data.length % 2 == 0) {
            return (data[middle - 1] + data[middle]) / 2.0;
        } else {
            return data[middle];
        }
    }

    public static double WilcoxonSignedRankTest(double[] before, double[] after) {

        double[] differences = new double[before.length];
        for (int i = 0; i < before.length; i++) {
            differences[i] = after[i] - before[i];
        }
        double[] absDifferences = Arrays.stream(differences).map(Math::abs).toArray();
        int[] ranks = rank(absDifferences);
        double positiveRankSum = 0;
        double negativeRankSum = 0;
        for (int i = 0; i < differences.length; i++) {
            if (differences[i] > 0) {
                positiveRankSum += ranks[i];
            } else if (differences[i] < 0) {
                negativeRankSum += ranks[i];
            }
        }
        double W = Math.min(positiveRankSum, negativeRankSum);
        return W;
    }

    private static int[] rank(double[] values) {
        int n = values.length;
        int[] ranks = new int[n];
        Double[] sortedValues = Arrays.stream(values).boxed().toArray(Double[]::new);
        Arrays.sort(sortedValues);

        for (int i = 0; i < n; i++) {
            ranks[i] = Arrays.asList(sortedValues).indexOf(values[i]) + 1;
        }
        return ranks;
    }

}
