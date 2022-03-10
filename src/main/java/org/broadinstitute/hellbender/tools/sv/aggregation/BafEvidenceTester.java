package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.mllib.clustering.GaussianMixture;
import org.apache.spark.mllib.clustering.GaussianMixtureModel;
import org.apache.spark.mllib.linalg.DenseMatrix;
import org.apache.spark.mllib.linalg.DenseVector;
import org.apache.spark.mllib.linalg.Matrix;
import org.apache.spark.mllib.linalg.Vector;
import org.apache.spark.mllib.stat.distribution.MultivariateGaussian;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.*;
import java.util.stream.Collectors;

public class BafEvidenceTester {

    private final SAMSequenceDictionary dictionary;

    public BafEvidenceTester(final SAMSequenceDictionary dictionary) {
        this.dictionary = dictionary;
    }

    public TTestResult calculateLogLikelihood(final SVCallRecord record, final List<BafEvidence> evidence,
                                         final Set<String> carrierSamples, final Set<String> backgroundSamples) {
        if (evidence == null || evidence.isEmpty()) {
            return null;
        }
        final List<BafEvidence> beforeBaf = new ArrayList<>(evidence.size());
        final List<BafEvidence> innerBaf = new ArrayList<>(evidence.size());
        final List<BafEvidence> afterBaf = new ArrayList<>(evidence.size());
        for (final BafEvidence baf : evidence) {
            if (baf.getStart() < record.getPositionA()) {
                beforeBaf.add(baf);
            } else if (baf.getStart() >= record.getPositionB()) {
                afterBaf.add(baf);
            } else {
                innerBaf.add(baf);
            }
        }
        final Map<String, Integer> beforeHetCounts = countSampleHets(beforeBaf);
        final Map<String, Integer> innerHetCounts = countSampleHets(innerBaf);
        final Map<String, Integer> afterHetCounts = countSampleHets(afterBaf);

        final Set<String> allSamples = Sets.union(carrierSamples, backgroundSamples);
        final Map<String, SampleStats> sampleStats = allSamples.stream().collect(Collectors.toMap(s -> s, s -> new SampleStats()));
        for (final BafEvidence baf : evidence) {
            if (baf.getStart() < record.getPositionA()) {
                sampleStats.get(baf.getSample()).beforeHetCount++;
            } else if (baf.getStart() >= record.getPositionB()) {
                sampleStats.get(baf.getSample()).afterHetCount++;
            } else {
                sampleStats.get(baf.getSample()).innerHetCount++;
            }
        }

        final double length = record.getLength();
        final double threshold = Math.min(50. / length, 0.001);
        int totalInnerCount = 0;
        int nullInnerCount = 0;
        final List<Double> nullRatios = new ArrayList<>();
        for (final Map.Entry<String, SampleStats> entry : sampleStats.entrySet()) {
            final String sample = entry.getKey();
            final SampleStats stats = entry.getValue();
            if (isRegionOfHomozygosity(stats.beforeHetCount, stats.innerHetCount, stats.afterHetCount, length, threshold)) {
                stats.isROH = true;
            } else {
                stats.deletionRatio = deletionTestStatistic(stats.beforeHetCount, stats.innerHetCount, stats.afterHetCount, length, threshold);
                if (!carrierSamples.contains(sample)) {
                    nullRatios.add(stats.deletionRatio);
                    nullInnerCount += stats.innerHetCount;
                }
            }
            totalInnerCount += stats.innerHetCount;
        }
        final double nullMean = nullInnerCount / (1.0 + nullRatios.size());
        final TTestResult result;
        if (nullRatios.size() > 10) {
            // Bayesian GMM?
            final double[] weights = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
            final Vector mean1 = new DenseVector(new double[]{0.});
            final Vector mean2 = new DenseVector(new double[]{0.5});
            final Vector mean3 = new DenseVector(new double[]{1.});
            final Matrix cov1 = new DenseMatrix(1, 1, new double[]{0.5});
            final Matrix cov2 = new DenseMatrix(1, 1, new double[]{0.5});
            final Matrix cov3 = new DenseMatrix(1, 1, new double[]{0.5});
            final MultivariateGaussian[] gaussians = {
                    new MultivariateGaussian(mean1, cov1),
                    new MultivariateGaussian(mean2, cov2),
                    new MultivariateGaussian(mean3, cov3),
            };
            final JavaSparkContext spark = new JavaSparkContext();
            final List<Vector> nullRatioVectors = nullRatios.stream().map(r -> new DenseVector(new double[]{r})).collect(Collectors.toList());
            final JavaRDD<Vector> data = spark.parallelize(nullRatioVectors);
            final GaussianMixtureModel gmm = new GaussianMixture().setK(3).run(data);
            result = tTest(sampleStats, nullRatios, totalInnerCount, gmm);
        } else {
            result = null;
        }
        return result;
    }

    private TTestResult tTest(final Map<String, SampleStats> sampleStats, final List<Double> nullRatios,
                              final int totalInnerCount, final GaussianMixtureModel gmm) {
        if (nullRatios.size() <= 10 || totalInnerCount < 10) {
            return null;
        }
        final double minNullRatio = nullRatios.stream().mapToDouble(d -> d).min().getAsDouble();
        final double maxNullRatio = nullRatios.stream().mapToDouble(d -> d).max().getAsDouble();
        if (maxNullRatio - minNullRatio < 0.0001) {
            return null;
        }
        final List<SampleStats> statsList = sampleStats.values().stream().filter(s -> !s.isROH).collect(Collectors.toList());
        if (statsList.isEmpty() || statsList.size() > nullRatios.size()) {
            return null;
        }
        final List<double[]> posterior = statsList.stream()
                .map(s -> new DenseVector(new double[]{s.deletionRatio})).map(gmm::predictSoft).collect(Collectors.toList());
        final double mean = nullRatios.stream().mapToDouble(Double::new).map(d -> Math.pow(10., -d)).sum() / nullRatios.size();
        final double logLik = posterior.stream()
                .mapToDouble(arr -> Math.log(MathUtils.arrayMax(arr)))
                .sum();
        return new TTestResult(mean, logLik);
    }

    private static final class SampleStats {
        public int beforeHetCount = 0;
        public int innerHetCount = 0;
        public int afterHetCount = 0;
        public Double deletionRatio = null;
        public boolean isROH = false;
    }

    public static final class TTestResult {
        public final double mean;
        public final double logLik;
        public TTestResult(final double mean, final double logLik) {
            this.mean = mean;
            this.logLik = logLik;
        }
    }

    protected static boolean isRegionOfHomozygosity(final int beforeHetCount, final int innerHetCount, final int afterHetCount, final double length, final double threshold) {
        return innerHetCount < threshold * length &&
            (beforeHetCount  < threshold * length || afterHetCount < threshold * length);
    }

    protected static double deletionTestStatistic(final int beforeHetCount, final int innerHetCount, final int afterHetCount, final double length, final double threshold) {
        final int flankHetCount = Math.min(beforeHetCount, afterHetCount);
        return Math.log10(innerHetCount + length * threshold) - Math.log10(flankHetCount + length * threshold);
    }

    protected static Map<String, Integer> countSampleHets(final Collection<BafEvidence> evidence) {
        return evidence.stream()
                .collect(Collectors.groupingBy(BafEvidence::getSample,
                        Collectors.collectingAndThen(Collectors.toList(), list -> list.stream().mapToInt(o -> 1).sum())));
    }

    public static final class BafEvidenceTestResult {
        private final double logLik;
        private final double snpRatio;
        private final double ksStat;

        public BafEvidenceTestResult(final double logLik, final double snpRatio, final double ksStat) {
            this.logLik = logLik;
            this.snpRatio = snpRatio;
            this.ksStat = ksStat;
        }

        public double getLogLikelihood() {
            return logLik;
        }

        public double getSnpRatio() {
            return snpRatio;
        }

        public double getKsStat() {
            return ksStat;
        }
    }
}
