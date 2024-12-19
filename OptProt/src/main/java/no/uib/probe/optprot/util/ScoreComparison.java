package no.uib.probe.optprot.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;

public class ScoreComparison {

    private double calculateMean(double[] scores) {
        double sum = 0;
        for (double score : scores) {
            sum += score;
        }
        return sum / scores.length;
    }

    private double percentageImprovementIndependenet(double referenceMedian, double improvedMedian) {
        if (referenceMedian == 0) {
            return (improvedMedian * 100.0);
        }
        return (((improvedMedian - referenceMedian) * 100.0) / referenceMedian);//
    }

    private double percentageImprovementPaired(double[] referenceData, double[] improvedData) {
        if (improvedData.length == 0) {
            return 0;
        }
        double avgMedian = 0;
        List<Double> values = new ArrayList<>();
        for (int i = 0; i < referenceData.length; i++) {
            double improvment = (improvedData[i] - referenceData[i]) / referenceData[i];
            improvment = improvment * 100.0;
            values.add(improvment);

        }
        Collections.sort(values);
        if (values.size() % 2 == 0) {
            int mediIndex;
            if (values.size() < 3) {
                mediIndex = 0;
            } else {
                mediIndex = values.size() / 2;
            }

            System.out.println("at data size " + values.size() + "  " + mediIndex);
            avgMedian = (values.get(mediIndex) + (values.get(mediIndex + 1))) / 2.0;

        } else {
            int mediIndex = values.size() / 2;
            avgMedian = values.get(mediIndex);
        }

        return avgMedian;
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

    public double calculateScore(double[] referenceScores, double[] improvedScores, boolean paired) {
        if (referenceScores.length < 2 && improvedScores.length < 2) {
            return 0;//new double[]{0, 0, 0, 0, 0, 0};
        }

        List<Double> valuesToScore = new ArrayList<>();
        DescriptiveStatistics referenceData = new DescriptiveStatistics(referenceScores);
        DescriptiveStatistics improvedData = new DescriptiveStatistics(improvedScores);

        double improvmentScorePersentage;
//        double comparisonScore;
        double sizeEffect = 0;
        if (paired) {
//            comparisonScore = wilcoxonSignedRankTestPair(referenceData.getValues(), improvedData.getValues());
            improvmentScorePersentage = percentageImprovementPaired(referenceScores, improvedScores);
//            System.out.println("paired improvment score is " + improvmentScorePersentage+"  "+comparisonScore+"  ");
            if (improvmentScorePersentage == 0) {
                return 0;
            }
        } else {
            double referenceMedian = referenceData.getPercentile(50);
            double improvedMedian = improvedData.getPercentile(50);
            if (Double.isNaN(referenceData.getMean())) {
                referenceData = new DescriptiveStatistics(new double[improvedScores.length]);
//                comparisonScore = wilcoxonSignedRankTestPair(referenceData.getValues(), improvedData.getValues());
                improvmentScorePersentage = 100;
            } else if (Double.isNaN(improvedData.getMean())) {
                improvedData = new DescriptiveStatistics(new double[referenceScores.length]);
//                comparisonScore = wilcoxonSignedRankTestPair(referenceData.getValues(), improvedData.getValues());
                improvmentScorePersentage = -100;
//                System.out.println("persent in case of zere improved data " + comparisonScore + "  " + improvmentScorePersentage);
            } else {
//                comparisonScore = mannWhitneyTestIndependent(referenceData.getValues(), improvedData.getValues());
                improvmentScorePersentage = percentageImprovementIndependenet(referenceMedian, improvedMedian);
            }
//             System.out.println("independent improvment score is " + improvmentScorePersentage+"  "+comparisonScore+"  ");

        }
        sizeEffect = medianBasedEffectSize(referenceData, improvedData);

        System.out.println("final size effect is " + sizeEffect + "   " + referenceData.getN() + "   " + improvedData.getN());
        if (referenceData.getN() > improvedData.getN() && sizeEffect > 0.0) {
            sizeEffect *= -1.0;
        }
//        double normalisedComparisonScore = logScaleNormalize(comparisonScore, 10);
        double normalisedImprovmentScore = logScaleNormalize(improvmentScorePersentage, 10);//scoreScaling(percentageImprovement,-30,30);
//        if (Double.isInfinite(Math.abs(normalisedComparisonScore))) {
//            normalisedComparisonScore = normalisedImprovmentScore;
//        }
        if (referenceData.getN() > 1 && improvedData.getN() > 1) {
//            valuesToScore.add(normalisedComparisonScore);
            valuesToScore.add(normalisedImprovmentScore);
        } else {
            referenceData.addValue(0);
            valuesToScore.add(0.0);
        }
        if (!paired) {
            sizeEffect = logScaleNormalize(sizeEffect, 10);
            valuesToScore.add(sizeEffect);
        }
        double finalScore = 0;
        for (double ss : valuesToScore) {
            finalScore += ss;
        }
        finalScore = finalScore / (double) valuesToScore.size();
        return finalScore;
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
        double pValue = mannWhitneyUTest.mannWhitneyUTest(sample1, sample2);
        if (pValue <= 0.05) {
            return mannWhitneyUTest.mannWhitneyU(sample1, sample2);
        }
        return 0.0;
    }

    public double wilcoxonSignedRankTestPair(double[] sample1, double[] sample2) {

        WilcoxonSignedRankTest wilcoxonSignedRankTest = new WilcoxonSignedRankTest();
        double pValue = wilcoxonSignedRankTest.wilcoxonSignedRankTest(sample1, sample2, false);
        if (pValue <= 0.05) {
            return wilcoxonSignedRankTest.wilcoxonSignedRank(sample1, sample2);
        }
        return 0.0;
    }

    public static void main(String[] args) {
        double[] scores1 = {300.0, 78.0, 84.0, 76.0, 82.5, 40.0, 20.0};
        double[] scores2 = {91.0, 89.0, 95.0, 87.0, 93.5, 90.0, 95.0};//, 89.0, 95.0, 87.0, 93.5, 90.0, 89.0, 95.0, 87.0, 93.5, 90.0};
        double[] scores3 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        DescriptiveStatistics ds1 = new DescriptiveStatistics(scores3);
        DescriptiveStatistics ds2 = new DescriptiveStatistics(scores2);
        ScoreComparison sc = new ScoreComparison();
        double cohenD = sc.logScaleNormalize(sc.calculateCohensD(ds1, ds2), 2);
        double normImpro = sc.logScaleNormalize(sc.percentageImprovementIndependenet(ds1.getPercentile(50), ds2.getPercentile(50)), 10);
        double normZ = sc.logScaleNormalize(StatisticsTests.independentZTest(ds1, ds2), 10);
        sc.calculateCohensD(ds1, ds2);
        double dMedian = sc.logScaleNormalize(sc.medianBasedEffectSize(ds1, ds2), 10);
        double finalScore = (cohenD + normZ + normImpro) / 3.0;
        double finalScore2 = (dMedian + normZ + normImpro) / 3.0;

        double normZ2 = sc.logScaleNormalize(sc.mannWhitneyTestIndependent(ds1.getValues(), ds2.getValues()), 10);

        double normZ3 = sc.wilcoxonSignedRankTestPair(ds1.getValues(), ds2.getValues());

//        System.out.println("Final Score: " + normImpro + "  " + normZ + "  " + cohenD + "   " + dMedian + "------------>> " + finalScore + "  vs " + finalScore2);
        System.out.println("z Score: " + "  " + normZ3 + "   " + "------------>> " + "  vs ");

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
