/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package no.uib.probe.optprot.util;
import java.util.Arrays;

public class ZTest {
    
    public static void main(String[] args) {
        // Example data
        double[] sample1 = {45.2, 46.1, 47.3, 50.0, 48.5};
        double[] sample2 = {52.0, 49.5, 53.2, 51.3, 50.4};
        
        double z = zTest(sample1, sample2);
        double pValue = pValue(z);
        
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
    
    public static double pValue(double z) {
        return 2 * (1 - cumulativeProbability(Math.abs(z)));
    }
    
    public static double cumulativeProbability(double z) {
        return 0.5 * (1 + erf(z / Math.sqrt(2)));
    }
    
    // Approximation of the error function
    public static double erf(double x) {
        double t = 1.0 / (1.0 + 0.5 * Math.abs(x));
        double tau = t * Math.exp(-x*x - 1.26551223 + t * (1.00002368 +
                        t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 +
                        t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
                        t * (-0.82215223 + t * 0.17087277)))))))));
        return x >= 0 ? 1 - tau : tau - 1;
    }
}
