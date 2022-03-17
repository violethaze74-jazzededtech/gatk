package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;

import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

public class BafEvidenceTester {

    private final SAMSequenceDictionary dictionary;
    private final Median median;

    public BafEvidenceTester(final SAMSequenceDictionary dictionary) {
        this.dictionary = dictionary;
        this.median = new Median();
    }

    private static void print(final SVCallRecord record, final String sample, final SampleStats stats, final Set<String> carrierSamples, final PrintStream printStream) {
        if (record.getContigA().equals("chrX") || record.getContigA().equals("chrY")) {
            return;
        }
        final List<String> builder = new ArrayList<>();
        builder.add(record.getId());
        builder.add(record.getContigA());
        builder.add(String.valueOf(record.getPositionA()));
        builder.add(String.valueOf(record.getPositionB()));
        builder.add(String.valueOf(record.getLength()));
        builder.add(record.getType().name());
        builder.add(String.valueOf(sample));
        builder.add(String.valueOf(carrierSamples.contains(sample)));
        builder.add(String.valueOf(stats.deletionRatio));
        builder.add(String.valueOf(stats.beforeHetCount));
        builder.add(String.valueOf(stats.innerHetCount));
        builder.add(String.valueOf(stats.afterHetCount));
        printStream.println(String.join("\t", builder));
    }

    private static void print2(final SVCallRecord record, final double v1, final double v2, final double v3, final PrintStream printStream) {
        if (record.getContigA().equals("chrX") || record.getContigA().equals("chrY")) {
            return;
        }
        final List<String> builder = new ArrayList<>();
        builder.add(record.getId());
        builder.add(record.getContigA());
        builder.add(String.valueOf(record.getPositionA()));
        builder.add(String.valueOf(record.getPositionB()));
        builder.add(String.valueOf(record.getLength()));
        builder.add(record.getType().name());
        builder.add(String.valueOf(v1));
        builder.add(String.valueOf(v2));
        builder.add(String.valueOf(v3));
        printStream.println(String.join("\t", builder));
    }

    public TestResult calculateLogLikelihood(final SVCallRecord record, final List<BafEvidence> evidence, final Set<String> excludedSamples) {
        if (evidence == null || evidence.isEmpty()) {
            return new TestResult(null, null);
        }
        final Set<String> carrierSamples = Sets.difference(record.getCarrierSampleSet(), excludedSamples);
        final Set<String> allSamples = Sets.difference(record.getAllSamples(), excludedSamples);
        final List<BafEvidence> filteredEvidence = evidence.stream().filter(baf -> allSamples.contains(baf.getSample())).collect(Collectors.toList());

        final List<BafEvidence> innerBaf = new ArrayList<>(filteredEvidence.size());
        final Map<String, SampleStats> sampleStats = allSamples.stream().collect(Collectors.toMap(s -> s, s -> new SampleStats()));
        for (final BafEvidence baf : filteredEvidence) {
            if (baf.getStart() < record.getPositionA()) {
                sampleStats.get(baf.getSample()).beforeHetCount++;
            } else if (baf.getStart() >= record.getPositionB()) {
                sampleStats.get(baf.getSample()).afterHetCount++;
            } else {
                sampleStats.get(baf.getSample()).innerHetCount++;
                innerBaf.add(baf);
            }
        }

        final Double deletionStatistic = calculateDeletionTestStatistic(record, sampleStats, carrierSamples);
        final Double duplicationStatistic = calculateDuplicationTestStatistic(innerBaf, carrierSamples);
        return new TestResult(deletionStatistic, duplicationStatistic);

        /*
        final double[] nullBaf = innerBaf.stream().filter(baf -> !carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).map(d -> Math.min(d, 1.0 - d)).toArray();
        final double[] carrierBaf = innerBaf.stream().filter(baf -> carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).map(d -> Math.min(d, 1.0 - d)).toArray();
        // Kolmogorov-Smirnov test
        if (nullBaf.length < 30 || carrierBaf.length < 30) {
            return null;
        }
        final double meanNull = new Median().evaluate(nullBaf);
        final double meanCarrier = new Median().evaluate(carrierBaf);

        final double ksP;
        if (nullBaf.length == 1) {
            final EmpiricalDistribution d = new EmpiricalDistribution();
            //d.load(carrierBaf);
            //ksP = 1.0 - d.cumulativeProbability(nullBaf[0]);
            return null;
        } else if (carrierBaf.length == 1) {
            final EmpiricalDistribution d = new EmpiricalDistribution();
            //d.load(nullBaf);
            //ksP = d.cumulativeProbability(carrierBaf[0]);
            return null;
        } else {
            final TTest ttest = new TTest();
            final double p = ttest.tTest(nullBaf, carrierBaf);
            ksP = meanCarrier < meanNull ? p : 1. - p;
            //final KolmogorovSmirnovTest ksTest = new KolmogorovSmirnovTest(); // TODO set random seed
            //ksP = 0; //ksTest.kolmogorovSmirnovTest(nullBaf, carrierBaf, false);
            //ksP = ksTest.approximateP(ksTest.kolmogorovSmirnovStatistic(nullBaf, carrierBaf), nullBaf.length, carrierBaf.length);
        }

        // Het count ratio stats (deletions)
        final Median median = new Median();
        final double[] nullRatiosArr = nullRatios.stream().mapToDouble(d -> d).toArray();
        final double[] carrierRatiosArr = carrierSamples.stream().map(sampleStats::get).filter(s -> !s.isROH).mapToDouble(s -> s.deletionRatio).toArray();
        final double medianCarrierRatio = median.evaluate(carrierRatiosArr);
        final double medianBackgroundRatio = median.evaluate(nullRatiosArr);
        final double medianBackgroundAbsoluteDeviation = median.evaluate(nullRatios.stream().mapToDouble(d -> Math.abs(medianBackgroundRatio - d)).toArray());
        final double hetCountRatioZ = (medianCarrierRatio - medianBackgroundRatio) / Math.max(0.001, medianBackgroundAbsoluteDeviation);

        print2(record, meanNull, meanCarrier, ksP, printStream);

        return new TestResult(hetCountRatioZ, meanCarrier);

         */

    }

    private Double calculateDuplicationTestStatistic(final List<BafEvidence> evidence, final Set<String> carrierSamples) {
        final double[] nullBaf = evidence.stream().filter(baf -> !carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).map(d -> Math.min(d, 1.0 - d)).toArray();
        final double[] carrierBaf = evidence.stream().filter(baf -> carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).map(d -> Math.min(d, 1.0 - d)).toArray();
        if (nullBaf.length < 30 || carrierBaf.length < 30) {
            return null;
        }
        final double meanNull = median.evaluate(nullBaf);
        final double meanCarrier = median.evaluate(carrierBaf);
        return meanNull - meanCarrier;
    }

    private Double calculateDeletionTestStatistic(final SVCallRecord record,
                                                  final Map<String, SampleStats> sampleStats,
                                                  final Set<String> carrierSamples) {
        final double length = record.getLength();
        final double threshold = Math.min(50. / length, 0.0005);
        int totalInnerCount = 0;
        final List<Double> nullRatios = new ArrayList<>();
        for (final Map.Entry<String, SampleStats> entry : sampleStats.entrySet()) {
            final String sample = entry.getKey();
            final SampleStats stats = entry.getValue();
            if (!isRegionOfHomozygosity(stats.beforeHetCount, stats.innerHetCount, stats.afterHetCount, length, threshold)) {
                stats.deletionRatio = calculateDeletionRatio(stats.beforeHetCount, stats.innerHetCount, stats.afterHetCount, length, threshold);
                if (!carrierSamples.contains(sample)) {
                    nullRatios.add(stats.deletionRatio);
                }
            }
            totalInnerCount += stats.innerHetCount;
        }

        if (nullRatios.size() <= 10 || totalInnerCount < 10) {
            return null;
        }
        final double minNullRatio = nullRatios.stream().mapToDouble(d -> d).min().getAsDouble();
        final double maxNullRatio = nullRatios.stream().mapToDouble(d -> d).max().getAsDouble();
        if (maxNullRatio - minNullRatio < 0.0001) {
            return null;
        }
        final List<SampleStats> carrierStats = carrierSamples.stream().map(sampleStats::get).filter(s -> s.deletionRatio != null).collect(Collectors.toList());
        if (carrierStats.isEmpty() || carrierStats.size() > nullRatios.size()) {
            return null;
        }
        final double[] carrierRatiosArr = carrierSamples.stream().map(sampleStats::get).filter(s -> s.deletionRatio != null).mapToDouble(s -> s.deletionRatio).toArray();
        return median.evaluate(carrierRatiosArr);
    }

    private static final class SampleStats {
        public int beforeHetCount = 0;
        public int innerHetCount = 0;
        public int afterHetCount = 0;
        public Double deletionRatio = null;
    }

    public static final class TestResult {
        private final Double delStat;
        private final Double dupStat;
        public TestResult(final Double delStat, final Double dupStat) {
            this.delStat = delStat;
            this.dupStat = dupStat;
        }

        public Double getDelStat() {
            return delStat;
        }

        public Double getDupStat() {
            return dupStat;
        }
    }

    protected static boolean isRegionOfHomozygosity(final int beforeHetCount, final int innerHetCount, final int afterHetCount, final double length, final double threshold) {
        return innerHetCount < threshold * length &&
            (beforeHetCount  < threshold * length || afterHetCount < threshold * length);
    }

    protected static double calculateDeletionRatio(final int beforeHetCount, final int innerHetCount, final int afterHetCount, final double length, final double threshold) {
        final int flankHetCount = Math.min(beforeHetCount, afterHetCount);
        return Math.log10(innerHetCount + length * threshold) - Math.log10(flankHetCount + length * threshold);
    }
}
