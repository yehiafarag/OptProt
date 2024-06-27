/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.util;

import java.util.Arrays;
import org.apache.commons.math3.stat.inference.TTest;

public class StatisticsTests {

    public static void main(String[] args) {
        // Example data
        double[] sample1 = {45.2, 46.1, 47.3, 50.0, 48.5};
        double[] sample2 = {52.0, 49.5, 53.2, 51.3, 50.4};

        double z = zTest(sample1, sample2);
        double pValue = pValueForZTest(z);

        System.out.printf("Z-test statistic: %.4f%n", z);
        System.out.printf("P-value: %.4f%n", pValue);
    }

    public static double zTest(double[] sample1, double[] sample2) {
        int n1 = sample1.length;
        int n2 = sample2.length;

        double mean1 = mean(sample1);
        double mean2 = mean(sample2);

        double std1 = stdDev(sample1);
        double std2 = stdDev(sample2);

        double z = (mean1 - mean2) / Math.sqrt((std1 * std1 / n1) + (std2 * std2 / n2));
        return z;
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
}
