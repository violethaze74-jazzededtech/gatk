package org.broadinstitute.hellbender.tools.sv.aggregation;

import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.MannWhitneyU;
import picard.util.MathUtil;

import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

public class BafEvidenceTester {

    private final SAMSequenceDictionary dictionary;

    public BafEvidenceTester(final SAMSequenceDictionary dictionary) {
        this.dictionary = dictionary;
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
        builder.add(String.valueOf(stats.isROH));
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

    public TestResult calculateLogLikelihood(final SVCallRecord record, final List<BafEvidence> evidence, final Set<String> excludedSamples, final PrintStream printStream) {
        final Set<String> carrierSamples = Sets.difference(record.getCarrierSampleSet(), excludedSamples);
        final Set<String> allSamples = Sets.difference(record.getAllSamples(), excludedSamples);
        if (evidence == null || evidence.isEmpty()) {
            return null;
        }
        final List<BafEvidence> filteredEvidence = evidence.stream().filter(baf -> allSamples.contains(baf.getSample())).collect(Collectors.toList());
        final List<BafEvidence> beforeBaf = new ArrayList<>(filteredEvidence.size());
        final List<BafEvidence> innerBaf = new ArrayList<>(filteredEvidence.size());
        final List<BafEvidence> afterBaf = new ArrayList<>(filteredEvidence.size());
        for (final BafEvidence baf : filteredEvidence) {
            if (baf.getStart() < record.getPositionA()) {
                beforeBaf.add(baf);
            } else if (baf.getStart() >= record.getPositionB()) {
                afterBaf.add(baf);
            } else {
                innerBaf.add(baf);
            }
        }

        final Map<String, SampleStats> sampleStats = allSamples.stream().collect(Collectors.toMap(s -> s, s -> new SampleStats()));
        for (final BafEvidence baf : filteredEvidence) {
            if (baf.getStart() < record.getPositionA()) {
                sampleStats.get(baf.getSample()).beforeHetCount++;
            } else if (baf.getStart() >= record.getPositionB()) {
                sampleStats.get(baf.getSample()).afterHetCount++;
            } else {
                sampleStats.get(baf.getSample()).innerHetCount++;
            }
        }

        final double length = record.getLength();
        final double threshold = Math.min(50. / length, 0.0005);
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

        //Filter
        /*
        if (nullRatios.size() < 10) {
            return null;
        }
        if (nullRatios.size() <= 10 || totalInnerCount < 10) {
            return null;
        }
        final double minNullRatio = nullRatios.stream().mapToDouble(d -> d).min().getAsDouble();
        final double maxNullRatio = nullRatios.stream().mapToDouble(d -> d).max().getAsDouble();
        if (maxNullRatio - minNullRatio < 0.0001) {
            return null;
        }
        final List<SampleStats> carrierStats = carrierSamples.stream().map(sampleStats::get).filter(s -> !s.isROH).collect(Collectors.toList());
        if (carrierStats.isEmpty() || carrierStats.size() > nullRatios.size()) {
            return null;
        }
         */

        final double[] nullBaf = innerBaf.stream().filter(baf -> !carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).map(d -> Math.min(d, 1.0 - d)).toArray();
        final double[] carrierBaf = innerBaf.stream().filter(baf -> carrierSamples.contains(baf.getSample())).mapToDouble(BafEvidence::getValue).map(d -> Math.min(d, 1.0 - d)).toArray();
        // Mann-Whitney test
        if (nullBaf.length == 0 || carrierBaf.length == 0) {
            return null;
        }
        if (record.getId().equals("batch1_manta_chr22_chr22_485")) {
            int x = 0;
        }
        final MannWhitneyU mwu = new MannWhitneyU();
        final double meanNull = MathUtil.mean(nullBaf);
        final double meanCarrier = MathUtil.mean(carrierBaf);
        final double mannWhitneyP = mwu.test(nullBaf, carrierBaf, MannWhitneyU.TestType.FIRST_DOMINATES).getP();
        if (Double.isNaN(mannWhitneyP)) {
            return null; //TODO
        }

        // Het count ratio stats (deletions)
        final Median median = new Median();
        final double[] nullRatiosArr = nullRatios.stream().mapToDouble(d -> d).toArray();
        final double[] carrierRatiosArr = carrierSamples.stream().map(sampleStats::get).filter(s -> !s.isROH).mapToDouble(s -> s.deletionRatio).toArray();
        final double medianCarrierRatio = median.evaluate(carrierRatiosArr);
        final double medianBackgroundRatio = median.evaluate(nullRatiosArr);
        final double medianBackgroundAbsoluteDeviation = median.evaluate(nullRatios.stream().mapToDouble(d -> Math.abs(medianBackgroundRatio - d)).toArray());
        final double hetCountRatioZ = (medianCarrierRatio - medianBackgroundRatio) / Math.max(0.001, medianBackgroundAbsoluteDeviation);

        print2(record, meanNull, meanCarrier, mannWhitneyP, printStream);

        return new TestResult(hetCountRatioZ, meanCarrier);

    }

    /*
    private TTestResult tTest(final Map<String, SampleStats> sampleStats, final List<Double> nullRatios,
                              final int totalInnerCount, final GaussianMixtureModel gmm) {
        final List<double[]> posterior = statsList.stream()
                .map(s -> new DenseVector(new double[]{s.deletionRatio})).map(gmm::predictSoft).collect(Collectors.toList());
        final double mean = nullRatios.stream().mapToDouble(Double::new).map(d -> Math.pow(10., -d)).sum() / nullRatios.size();
        final double logLik = posterior.stream()
                .mapToDouble(arr -> Math.log(MathUtils.arrayMax(arr)))
                .sum();
        return new TTestResult(mean, logLik);
    }*/

    private static final class SampleStats {
        public int beforeHetCount = 0;
        public int innerHetCount = 0;
        public int afterHetCount = 0;
        public Double deletionRatio = null;
        public boolean isROH = false;
    }

    public static final class TestResult {
        public final double hetCountRatioZ;
        public final double mannWhitneyP;
        public TestResult(final double hetCountRatioZ, final double mannWhitneyP) {
            this.hetCountRatioZ = hetCountRatioZ;
            this.mannWhitneyP = mannWhitneyP;
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
