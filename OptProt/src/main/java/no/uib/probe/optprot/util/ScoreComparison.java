/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.util;

import org.apache.commons.math3.stat.inference.TTest;

public class ScoreComparison {

    public double calculateMean(double[] scores) {
        double sum = 0;
        for (double score : scores) {
            sum += score;
        }
        return sum / scores.length;
    }

    public double percentageImprovement(double[] fromData, double[] toData) {
        double mean1 = calculateMean(fromData);
        double mean2 = calculateMean(toData);
        return ((mean2 - mean1) / mean1) * 100;
    }

    public double calculateStandardDeviation(double[] scores) {
        double mean = calculateMean(scores);
        double sum = 0;
        for (double score : scores) {
            sum += Math.pow(score - mean, 2);
        }
        return Math.sqrt(sum / scores.length);
    }

    public double performTStatTest(double[] scores1, double[] scores2) {
        if (scores1.length < 2 || scores2.length < 2) {
            return -100.0;
        }
        TTest tTest = new TTest();
        if (scores1.length == scores2.length) {
            return tTest.pairedT(scores2, scores1);
        } else {
            return tTest.t(scores2, scores1);
        }
    }

    public double performTTest(double[] scores1, double[] scores2) {
        if (scores1.length < 2 || scores2.length < 2) {
            return 1;
        }
        TTest tTest = new TTest();
        return tTest.tTest(scores2, scores1);
    }

    public double normalize(double value, double min, double max) {
        return (value - min) / (max - min);
    }

    public double calculateFinalScore(double[] scores1, double[] scores2) {
        double mean1 = calculateMean(scores1);
        double mean2 = calculateMean(scores2);

        double meanDifference = mean2 - mean1;

        double percentageImprovement = percentageImprovement(scores1, scores2);

        double stdDev1 = calculateStandardDeviation(scores1);
        double stdDev2 = calculateStandardDeviation(scores2);
        double stdDevDifference = Math.abs(stdDev2 - stdDev1);
        double pValue = performTTest(scores1, scores2);

        // Assuming the possible range for these metrics is between -200 and 200 for normalization purposes.
        // Adjust these ranges based on your specific use case.
        double normMeanDiff = normalize(meanDifference, 0, 100);

        double normPercentageImprovement = normalize(percentageImprovement, 0, 100);
        double normStdDevDiff = normalize(stdDevDifference, 0, 100);
        double normPValue = (1 - pValue) * 100; // Inverting p-value for interpretation

        if (meanDifference == 0) {
            System.out.println("-----------------------unsupport param was here---------------------- ");
        }
        System.out.println("mean  " + meanDifference + "  " + percentageImprovement + "  sdv " + stdDevDifference + "   meanDifference  " + meanDifference);
        // Combining normalized results into a final score
        double finalScore = (normMeanDiff + normPercentageImprovement + normStdDevDiff + normPValue) / 4;

        return finalScore;
    }

    public static void main(String[] args) {
        double[] scores1 = {91.0, 89.0, 95.0, 87.0, 93.5, 90.0};
        double[] scores2 = {91.5, 91.0, 95.5, 88.0, 94.5, 91.0};

        ScoreComparison sc = new ScoreComparison();
        double finalScore = sc.calculateFinalScore(scores1, scores2);

        System.out.println("Final Score: " + finalScore);

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
}
