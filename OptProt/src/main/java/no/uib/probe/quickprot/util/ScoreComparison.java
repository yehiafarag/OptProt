package no.uib.probe.quickprot.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;

public class ScoreComparison {

    private double calculateMean(double[] scores) {
        double sum = 0;
        for (double score : scores) {
            sum += score;
        }
        return sum / scores.length;
    }

    private double percentageImprovementIndependenet(double referenceMedian, double improvedMedian) { //(improvedMedian - referenceMedian) * 100.0) / referenceMedian)
        if (referenceMedian == 0) {
            return (improvedMedian * 100.0);
        }
        return (((improvedMedian - referenceMedian) * 100.0) / referenceMedian);//
    }

    private double percentageImprovementPaired(double[] referenceData, double[] improvedData) {
        if (improvedData.length == 0) {
            return 0;
        }
        double[] dArr = new double[referenceData.length];
        for (int i = 0; i < referenceData.length; i++) {
            double improvment = (improvedData[i] - referenceData[i]) / referenceData[i];
            improvment = improvment * 100.0;
            dArr[i] = improvment;
        }
        DescriptiveStatistics ds = new DescriptiveStatistics(dArr);
        if (ds.getPercentile(50.0) == 0.0) {//|| 
            return ds.getMean();
        }

        return ds.getPercentile(50.0);
    }
//rank-based Cohen's d

    public double medianBasedEffectSize(DescriptiveStatistics referenceSamples, DescriptiveStatistics improvedSample) {
// Step 1: Calculate the median of each group
        double medianGroup1 = referenceSamples.getPercentile(50);
        double medianGroup2 = improvedSample.getPercentile(50);

        // Step 2: Calculate the MAD for each group
        double madGroup1 = calculateMAD(referenceSamples);
        double madGroup2 = calculateMAD(improvedSample);

        // Step 3: Calculate the pooled MAD
        double pooledMAD = (madGroup1 + madGroup2) / 2;
        // Step 4: Calculate the median-based effect size
        double dMedian = (medianGroup2 - medianGroup1) / pooledMAD;
        if (Double.isInfinite(dMedian)) {
            dMedian = 0;
        }

        return dMedian;
    }

    private double calculateMAD(DescriptiveStatistics samples) {
        // Step 1: Calculate the median of the dataset
        double medianValue = samples.getPercentile(50);

        // Step 2: Calculate the absolute deviations from the median
        double[] absoluteDeviations = Arrays.stream(samples.getValues())
                .map(value -> Math.abs(value - medianValue))
                .toArray();

        // Step 3: Calculate the median of the absolute deviations
        double mad = new DescriptiveStatistics(absoluteDeviations).getPercentile(50);

        return mad;
    }

    public double calculateCohensD(DescriptiveStatistics referenceSamples, DescriptiveStatistics improvedSample) {

        double mean1 = referenceSamples.getMean();
        double mean2 = improvedSample.getMean();
        double variance1 = referenceSamples.getVariance();
        double variance2 = improvedSample.getVariance();
        double pooledStdDev = Math.sqrt(((referenceSamples.getN() - 1) * variance1 + (improvedSample.getN() - 1) * variance2) / (referenceSamples.getN() + improvedSample.getN() - 2));
        double normCohenD = ((mean2 - mean1) / pooledStdDev);
        return normCohenD;
    }

    public double calculateStandardDeviation(double[] scores) {
        double mean = calculateMean(scores);
        double sum = 0;
        for (double score : scores) {
            sum += Math.pow(score - mean, 2);
        }
        return Math.sqrt(sum / scores.length);
    }
    // Create a normal distribution with mean 0 and standard deviation 1

    public double dataSizeEffect(double dataSizeFrom, double dataSizeTo) {
        dataSizeFrom++;
        dataSizeTo++;
        double effect = (dataSizeTo - dataSizeFrom) * 100.0 / Math.min(dataSizeFrom, dataSizeTo);
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

    /**
     * Applies log scale normalization to an array of data.
     *
     * @param value The array of data to be normalized.
     * @param base The base of the logarithm to use.
     * @return The log scale normalized array of data.
     */
    public static double logScaleNormalize(double value, double base) {

        double sign = 1;
        if (value < 0) {
            sign = -1;
        }
        value = Math.abs(value);
        double normalizedData = (Math.log(value + 1) / Math.log(base)) * sign;
        if (Double.isNaN(normalizedData)) {
            normalizedData = 0;
        }
        return normalizedData;
    }

    /**
     * Scale improvement score based on z-score range.
     *
     * @param score The improvement score.
     * @return The normalized scare.
     */
    public static double scoreScaling(double score, double scoreMin, double scoreMax) {
        double minValue = -1.0;
        double maxValue = 1.0;

        score = Math.min(score, scoreMax);
        score = Math.max(score, scoreMin);

        double normalizedScore = ((score - scoreMin) / (scoreMax - scoreMin)) * (maxValue - minValue);
        normalizedScore += minValue;

        normalizedScore = Math.max(normalizedScore, minValue);
        normalizedScore = Math.min(normalizedScore, maxValue);

        return normalizedScore;
    }

    public double calculateScore(double[] referenceInputData, double[] testInputData, boolean paired) {

        double[] referenceScores = referenceInputData;
        double[] testScores = testInputData;
        if (referenceScores.length < 2 && testScores.length < 2) {
            return 0;//new double[]{0, 0, 0, 0, 0, 0};
        }
        System.out.println("start comparison "+paired);
        DescriptiveStatistics referenceData = new DescriptiveStatistics(referenceScores);
        DescriptiveStatistics improvedData = new DescriptiveStatistics(testScores);
        if (Double.isNaN(referenceData.getPercentile(50)) && Double.isNaN(improvedData.getPercentile(50))) {
            System.out.println("both nan " + referenceData.getGeometricMean() + "   " + referenceData.getN() + "  " + improvedData.getN());
            return 0;
        }
        double improvmentScorePersentage;
        if (paired) {
            improvmentScorePersentage = percentageImprovementPaired(referenceScores, testScores);

            return improvmentScorePersentage;
        } else {
            if (Double.isNaN(improvedData.getGeometricMean()) || Double.isNaN(referenceData.getGeometricMean())) {//|| (improvedData.getN() <= 4 && referenceData.getN() > 20)) {
                if (referenceData.getN() > improvedData.getN()) {
                     System.out.println("return -100 ------------2");
                    return -100;
                } else if (referenceData.getN() < improvedData.getN()) {
                    return 100;
                } else {
                    return 0;
                }

            }

            int df = StatisticsTests.calculateDegreesOfFreedom((int) referenceData.getN(), (int) improvedData.getN());
            if (df == 0) {
                return 0;
            }
            double zRefTest = StatisticsTests.independentZTest(referenceData, improvedData);
            // Significance level (alpha)
            double alpha = 0.05;
            // Degrees of freedom
            // Calculate the critical value for a two-tailed test
            double criticalValue = StatisticsTests.getCriticalValue(alpha, df);

            if (!(Math.abs(zRefTest) > Math.abs(criticalValue)) || Double.isNaN(zRefTest)) {//|| (improvedData.getN() <= 4 && referenceData.getN() > 20)) {
                System.out.println("unpatired data un balancd is here we need a way to compare or ignor ------improvedData.getGeometricMean() " + improvedData.getGeometricMean() + " referenceData.getGeometricMean() " + referenceData.getGeometricMean());

//                double prec = ((double) Math.abs(referenceData.getN() - improvedData.getN())) / (double) Math.max(referenceData.getN(), improvedData.getN());
//                prec *= 100.0;
                double referenceMedian = referenceData.getGeometricMean();// referenceData.getPercentile(50);
                double improvedMedian = improvedData.getGeometricMean();//improvedData.getPercentile(50);
                improvmentScorePersentage = percentageImprovementIndependenet(referenceMedian, improvedMedian);
//                if (referenceData.getN() > improvedData.getN()) {
//                    System.out.println("improvment pers " + improvmentScorePersentage + "  " + (-1.0 * prec));
//                    return -1.0 * prec;
//                } else {
//                    System.out.println("improvment pers " + improvmentScorePersentage + "  " + (prec));
//                    return prec;
//                }
                return improvmentScorePersentage;
            } else {
                double referenceMedian = referenceData.getGeometricMean();// referenceData.getPercentile(50);
                double improvedMedian = improvedData.getGeometricMean();//improvedData.getPercentile(50);
                improvmentScorePersentage = percentageImprovementIndependenet(referenceMedian, improvedMedian);
                return improvmentScorePersentage;

            }

        }

    }

    public double calculateScore1(double[] referenceInputData, double[] testInputData, boolean paired) {

        double[] referenceScores = referenceInputData;
        double[] testScores = testInputData;
        if (referenceScores.length < 2 && testScores.length < 2) {
            return 0;//new double[]{0, 0, 0, 0, 0, 0};
        }

        List<Double> valuesToScore = new ArrayList<>();
        DescriptiveStatistics referenceData = new DescriptiveStatistics(referenceScores);
        DescriptiveStatistics improvedData = new DescriptiveStatistics(testScores);
        if (Double.isNaN(referenceData.getPercentile(50)) && Double.isNaN(improvedData.getPercentile(50))) {
            System.out.println("both nan " + referenceData.getGeometricMean() + "   " + referenceData.getN() + "  " + improvedData.getN());
            return 0;//new double[]{0, 0, 0, 0, 0, 0};
        }

        double improvmentScorePersentage;
//        double comparisonScore;
//        double sizeEffect;
//        double cohensD = 0;
        if (paired) {
            improvmentScorePersentage = percentageImprovementPaired(referenceScores, testScores);
            return improvmentScorePersentage;
//            if (improvmentScorePersentage == 0) {
//                return 0;
//            }
//            comparisonScore = wilcoxonSignedRankTestPair(referenceData.getValues(), improvedData.getValues());
//            if (improvmentScorePersentage < 0) {
//                comparisonScore *= -1.0;
//            }
//            sizeEffect = medianBasedEffectSize(referenceData, improvedData);
//            cohensD = 0;//calculateCohensD(referenceData, improvedData);

        } else {
            int df = StatisticsTests.calculateDegreesOfFreedom((int) referenceData.getN(), (int) improvedData.getN());
            if (df == 0) {
                return 0;
            }
            double zRefTest = StatisticsTests.independentZTest(referenceData, improvedData);
            // Significance level (alpha)
            double alpha = 0.05;
            // Degrees of freedom
            // Calculate the critical value for a two-tailed test
            double criticalValue = StatisticsTests.getCriticalValue(alpha, df);

//            if ((Math.abs(zRefTest) > Math.abs(criticalValue)) != (Math.abs(mwTest) > Math.abs(criticalValue))) {
//            }
            // Calculate the critical value for a two-tailed test
//            if (!(Math.abs(zRefTest) > Math.abs(criticalValue))) {
//                double criticalValue2 = StatisticsTests.getCriticalValue(0.5, df);
//                System.out.println("update Critical value: " + criticalValue + "  diffrent is sig  " + (Math.abs(zRefTest) > Math.abs(criticalValue2)) + "  " + Math.abs(zRefTest) + " > " + Math.abs(criticalValue2));
//            }
//            if (!(Math.abs(zRefTest) > Math.abs(criticalValue)) || Double.isNaN(zRefTest)) {
//                System.out.println("not sig diffrenet " + referenceData.getN() + "  " + improvedData.getN() + "   ");
//                double prec = (Math.abs(referenceData.getN() - improvedData.getN())) / Math.max(referenceData.getN(), improvedData.getN());
//                if (referenceData.getN() > improvedData.getN()) {
//                    return -1.0 * prec;
//                } else {
//                    return prec;
//                }
//            }
//            if (Double.isNaN(referenceData.getGeometricMean()) || (referenceData.getN() <= 4 && improvedData.getN() > 20)) {
//                System.out.println("ref data was nan " + referenceData.getGeometricMean() + "   " + referenceData.getN() + "  " + improvedData.getN());
//                improvmentScorePersentage = 100.0;
//                return improvmentScorePersentage;//logScaleNormalize(improvmentScorePersentage, 10);
//            } else
            if (!(Math.abs(zRefTest) > Math.abs(criticalValue)) || Double.isNaN(zRefTest) || Double.isNaN(improvedData.getGeometricMean()) || Double.isNaN(referenceData.getGeometricMean())) {//|| (improvedData.getN() <= 4 && referenceData.getN() > 20)) {

                System.out.println("unpatired data un balancd is here ------");
                double prec = ((double) Math.abs(referenceData.getN() - improvedData.getN())) / (double) Math.max(referenceData.getN(), improvedData.getN());
                prec *= 100.0;
                double referenceMedian = referenceData.getGeometricMean();// referenceData.getPercentile(50);
                double improvedMedian = improvedData.getGeometricMean();//improvedData.getPercentile(50);
                improvmentScorePersentage = percentageImprovementIndependenet(referenceMedian, improvedMedian);
                if (referenceData.getN() > improvedData.getN()) {
                    System.out.println("improvment pers " + improvmentScorePersentage + "  " + (-1.0 * prec));
                    return -1.0 * prec;
                } else {
                    System.out.println("improvment pers " + improvmentScorePersentage + "  " + (prec));
                    return prec;
                }
            } else {
                double referenceMedian = referenceData.getGeometricMean();// referenceData.getPercentile(50);
                double improvedMedian = improvedData.getGeometricMean();//improvedData.getPercentile(50);
//                if (referenceMedian == 0 || improvedMedian == 0) {
////                    System.out.println("swich to mean before " + referenceMedian + "  " + improvedMedian);
//                    referenceMedian = referenceData.getMean();
//                    improvedMedian = improvedData.getMean();
////                    System.out.println("swich to mean after " + referenceMedian + "  " + improvedMedian);
//                }
//                comparisonScore = mannWhitneyTestIndependent(referenceData.getValues(), improvedData.getValues());
                improvmentScorePersentage = percentageImprovementIndependenet(referenceMedian, improvedMedian);
                return improvmentScorePersentage;
//                if (improvmentScorePersentage < 0) {
//                    comparisonScore *= -1.0;
//                }
//                cohensD = calculateCohensD(referenceData, improvedData);
//                sizeEffect = medianBasedEffectSize(referenceData, improvedData);
            }

        }

//        if (!paired && (improvmentScorePersentage > 0 != cohensD > 0)) {
//            System.out.println("final scores " + improvmentScorePersentage + "   SE  " + "  cohn d " + cohensD + "  ref: " + referenceData.getN() + "  target  " + improvedData.getN());
//            double referenceMedian = referenceData.getMean();
//            double improvedMedian = improvedData.getMean();
//            double improvmentScorePersentage2 = percentageImprovementIndependenet(referenceMedian, improvedMedian);
//            System.out.println("updated use avg  scores " + improvmentScorePersentage2);
//
////            for (double d : referenceInputData) {
////                System.out.print(d + ",");
////            }
////            System.out.println("");
////            System.out.println("vs");
////            for (double d : improvedInputData) {
////                System.out.print(d + ",");
////            }
////            System.out.println("");
////            System.exit(0);
//        }
//        if (cohensD < 0 && ((comparisonScore > 0) || (improvmentScorePersentage > 0) || (sizeEffect > 0))) {
//            System.out.println("------------------------------------------------------>>>>> cohen d negative " + cohensD + "  " + comparisonScore + "  " + improvmentScorePersentage + "  " + sizeEffect);
//            if (comparisonScore > 0) {
//                comparisonScore *= -1.0;
//            }
//            if (improvmentScorePersentage > 0) {
//                improvmentScorePersentage *= -1.0;
//            }
//            if (sizeEffect > 0) {
//                sizeEffect *= -1.0;
//            }
//        }
//        double normalisedComparisonScore = logScaleNormalize(comparisonScore, 10);
//        double normalisedImprovmentScore = improvmentScorePersentage;//logScaleNormalize(improvmentScorePersentage, 10);//scoreScaling(percentageImprovement,-30,30);
//        double normalisedSizeEffectScore = logScaleNormalize(sizeEffect, 10);
//        if (referenceData.getN() > 1 && improvedData.getN() > 1) {
////            valuesToScore.add(normalisedComparisonScore);
//            valuesToScore.add(normalisedImprovmentScore);
////            valuesToScore.add(normalisedSizeEffectScore);
////            valuesToScore.add(normalisedCohensDScore);
//        } else {
//            referenceData.addValue(0);
//            valuesToScore.add(0.0);
//        }
//
//        double finalScore = 0;
//        for (double ss : valuesToScore) {
//            finalScore += ss;
//        }
//        finalScore = finalScore / (double) valuesToScore.size();
//        if(comparisonScore<0.05)
//            finalScore=0;
//        if (flip) {
//            finalScore *= -1.0;
//        }
//        return improvmentScorePersentage;//finalScore;
    }

    public double calculateCohensD(double[] s1, double[] s2, boolean switchData) {
        double[] from, to;
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
        return (normCohenD);
    }

    public double oneToManyScore(double fromScore, DescriptiveStatistics stats) {
        double oldAverage = stats.getMean();
        // Step 2: Calculate the improvement score (percentage improvement)
        double improvementPercentage;
        if (fromScore == 0) {
            improvementPercentage = 1;
        } else {
            improvementPercentage = ((oldAverage - fromScore) / fromScore) * 100;
        }
        // Step 3: Compare the size of the arrays (if there are more or fewer old scores)
        int newSize = 1; // Single new score vs. old scores array
        int oldSize = (int) stats.getN();
        // Optional: Adjust the improvement calculation based on size
        double sizeFactor = (double) oldSize / newSize;
        double adjustedImprovement = improvementPercentage * sizeFactor;
        return adjustedImprovement;
    }

    public double mannWhitneyTestIndependent(double[] sample1, double[] sample2) {

        MannWhitneyUTest mannWhitneyUTest = new MannWhitneyUTest();
        if (sample1.length == 0 || sample2.length == 0) {
            return Double.NaN;
        }

//        double pValue = mannWhitneyUTest.mannWhitneyUTest(sample1, sample2);
//        if (pValue <= 0.05) {
//            int n1 = sample1.length;
//            int n2 = sample2.length;
//            double N = n1 + n2;
        double U = mannWhitneyUTest.mannWhitneyU(sample1, sample2);
////            Calculate rank  -biserial correlation
//            double rankBiserial = 1 - (2 * U) / (n1 * n2);
//
//            // Calculate Z-score (approximation)
//            double meanU = n1 * n2 / 2.0;
//            double stdU = Math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0);
//            double Z = (U - meanU) / stdU;
//
//            // Calculate effect size r
//            double r = Z / Math.sqrt(N);
//            System.out.println("at r is " + r);
        return U;
//            
//            return r;
//        }
//        return 0;
    }

    public double wilcoxonSignedRankTestPair(double[] referenceData, double[] improvedData) {

        WilcoxonSignedRankTest wilcoxonSignedRankTest = new WilcoxonSignedRankTest();
        double pValue = wilcoxonSignedRankTest.wilcoxonSignedRankTest(referenceData, improvedData, false);

        return pValue;
    }

    public static void main(String[] args) {
        double[] scores1 = {24.181788868689083, 26.381088683588953, 26.932822348455282, 28.183098374208907, 26.502771452416063, 26.313842050319227, 23.30745235741934, 24.101928241008064, 25.981828068689527, 23.034502191750313, 23.16842976349253};
        double[] scores2 = {23.85216786366755, 24.052747807215674, 22.680980745877108, 23.174512846898736, 31.306301299262387, 27.555997926642448, 23.74433648517002, 29.399078749253825, 26.16724952495946, 24.702858615186784};
//        double[] scores3 = {91.0, 89.0, 95.0, 87.0, 93.5, 90.0, 95.0, 89.0, 95.0, 87.0, 93.5, 90.0, 89.0, 95.0, 87.0, 93.5, 90.0, 20, 05};
//        double[] beforeWeights = {85, 78, 92, 70, 65, 90, 72, 88, 76, 82};
//        double[] afterWeights = {80, 75, 88, 68, 63, 85, 70, 84, 73, 78};
//
        DescriptiveStatistics ds1 = new DescriptiveStatistics(scores1);
        DescriptiveStatistics ds2 = new DescriptiveStatistics(scores2);
        ScoreComparison sc = new ScoreComparison();
        double cohenD = (sc.calculateCohensD(ds1, ds2));
        double medianImpro = sc.percentageImprovementIndependenet(ds1.getPercentile(50), ds2.getPercentile(50));
        double meanImpro = sc.percentageImprovementIndependenet(ds1.getGeometricMean(), ds2.getGeometricMean());
        double mean2Impro = sc.percentageImprovementIndependenet(ds1.getMean(), ds2.getMean());
        double normZ = StatisticsTests.independentZTest(ds1, ds2);
        System.out.println("improvment " + medianImpro + "   mean " + meanImpro + "  " + cohenD + "  " + normZ + "  " + mean2Impro);
//       
//        double normImpro = sc.logScaleNormalize(;
//        
//        sc.calculateCohensD(ds1, ds2);
//        double dMedian = (sc.medianBasedEffectSize(ds1, ds2));
//        double finalScore = (cohenD + normZ + normImpro) / 3.0;
//        double finalScore2 = (dMedian + normZ + normImpro) / 3.0;
//
////        double normZ2 = sc.mannWhitneyTestIndependent(ds1.getValues(), ds2.getValues());
//        double normZ3 = sc.wilcoxonSignedRankTestPair(scores3, scores2);
//        double normZ4 = StatisticsTests.WilcoxonSignedRankTest(scores3, scores2);
////        System.out.println("Final Score: " + normImpro + "  " + normZ + "  " + cohenD + "   " + dMedian + "------------>> " + finalScore + "  vs " + finalScore2);
////        System.out.println("z Score: " + "  " + normZ3 + "  vs " + normZ4);
//
////         Interpretation
//        if (finalScore > 0.75) {

//        } else if (finalScore > 0.5) {
//            System.out.println("Moderate positive enhancement in scores.");
//        } else if (finalScore > 0.25) {
//            System.out.println("Slight positive enhancement in scores.");
//        } else {
//            System.out.println("No significant enhancement or decline in scores.");
//        }
//
    }

}
