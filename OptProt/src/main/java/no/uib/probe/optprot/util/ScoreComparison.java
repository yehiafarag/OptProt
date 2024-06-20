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

    public double percentageImprovement(double[] scores1, double[] scores2) {
        double mean1 = calculateMean(scores1);
        double mean2 = calculateMean(scores2);
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

    public double performTTest(double[] scores1, double[] scores2) {
        TTest tTest = new TTest();
        return tTest.tTest(scores1, scores2);
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
        double stdDevDifference = stdDev2 - stdDev1;
        double pValue = performTTest(scores1, scores2);

        // Assuming the possible range for these metrics is between -200 and 200 for normalization purposes.
        // Adjust these ranges based on your specific use case.
        double normMeanDiff = normalize(meanDifference, 0, 100);

        double normPercentageImprovement = normalize(percentageImprovement, 0, 100);
        double normStdDevDiff = normalize(stdDevDifference, 0, 100);
        double normPValue = 1 - pValue; // Inverting p-value for interpretation
        System.out.println("mean 1   " + meanDifference + "  " + percentageImprovement + "  sdv " + stdDevDifference + "   " + normMeanDiff + "  norm " + normPercentageImprovement);
        // Combining normalized results into a final score
        double finalScore = (normMeanDiff + normPercentageImprovement + normStdDevDiff + normPValue) / 4;

        return finalScore;
    }

    public static void main(String[] args) {
       double[] scores1 = {91.0, 89.0, 95.0, 87.0, 93.5, 90.0};
       double[] scores2 = {0.5, 0.0, 0.0, 0.0, 88.5};
        

        ScoreComparison sc = new ScoreComparison();
        double finalScore = sc.calculateFinalScore(scores1, scores2);

        System.out.println("Final Score: " + finalScore);
        System.out.println("Final Score: " + finalScore);

        // Interpretation
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
