package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

import java.util.Arrays;

/**
 * This class has just one fundamental responsibility: calculating genotype likelihoods through the formula:
 *
 * Prob(reads | genotype) = product_{all reads} [[sum_{alleles in genotype} Prob(read | allele)]/ploidy]
 *
 * Note that this applies to non-somatic variant calling, where ploidy is a known integer and genotypes are given by the
 * number of copies of each allele.
 *
 * There is some unfortunate mixing-in of getting cached GenotypeAlleleCounts objects that will hopefully be improved
 * in future refactoring.
 *
 * COMPUTATIONAL NOTE
 * In the multiallelic calculation we accumulate the likelihood contribution of each read one allele at a time.  That is,
 * for genotype containing alleles A, B, C, we first fill an array with the likelihood contributions from allele A, then
 * we make a second pass and add the contributions from allele B, then allele C.  Traversing all the reads in each
 * allele row of the likelihoods array in this manner is cache-friendly and makes an enormous difference in runtime.
 *
 * The difference in performance comes from the fact that we index likelihoods first by allele, then by read.  Because of this,
 * likelihoods of consecutive reads with the same allele are adjacent in memory while likelihoods of consecutive alleles with the same read
 * are not.  In the former case looking up new likelihoods almost always results in a cache hit since many reads of the same allele
 * are loaded on the same cache page.
 *
 * If the cache-friendliness of this class is broken, it will show up as a severe regression in the runtime of its unit tests
 * for larger ploidies and allele counts.
 */
public class GenotypeLikelihoodCalculator {
    final int genotypeCount;

    final int alleleCount;

    final int ploidy;

    public GenotypeLikelihoodCalculator(final int ploidy, final int alleleCount) {
        genotypeCount = GenotypeIndexCalculator.genotypeCount(ploidy, alleleCount);
        this.alleleCount = alleleCount;
        this.ploidy = ploidy;
    }

    /**
     * Returns the number of possible genotypes given ploidy and the maximum allele index.
     * @return never {@code null}.
     */
    public int genotypeCount()  {
        return genotypeCount;
    }

    /**
     * Calculate the log10AlleleLikelihoods given the list of alleles and the likelihood map.
     *
     * @param log10AlleleLikelihoods the likelihood matrix all alleles vs all reads.
     *
     * @throws IllegalArgumentException if {@code alleleList} is {@code null} or {@code log10AlleleLikelihoods} is {@code null}
     *     or the alleleList size does not match the allele-count of this calculator, or there are missing allele vs
     *     read combinations in {@code log10AlleleLikelihoods}.
     *
     * @return never {@code null}.
     */
    public <EVIDENCE, A extends Allele> GenotypeLikelihoods log10GenotypeLikelihoods(final LikelihoodMatrix<EVIDENCE, A> log10AlleleLikelihoods) {
        final double[] log10GenotypeLikelihoods = computeLog10GenotypeLikelihoods(log10AlleleLikelihoods);
        return GenotypeLikelihoods.fromLog10Likelihoods(log10GenotypeLikelihoods);
    }

    /**
     * Compute the genotype log10 likelihoods as an array in the canonical genotype order.  That is, result[i] = Pr(reads | ith genotype)
     *
     * @param log10AlleleLikelihoods   log 10 likelihood matrix indexed by allele, then read
     * @return the log 10 likelihood of each genotype as an array
     */
    <EVIDENCE, A extends Allele> double[] computeLog10GenotypeLikelihoods(final LikelihoodMatrix<EVIDENCE, A> log10AlleleLikelihoods) {
        Utils.nonNull(log10AlleleLikelihoods);
        Utils.validateArg(log10AlleleLikelihoods.numberOfAlleles() == alleleCount, "mismatch between allele list and alleleCount");
        final int readCount = log10AlleleLikelihoods.evidenceCount();

        final double[][] log10LikelihoodsByAlleleAndRead = log10AlleleLikelihoods.asRealMatrix().getData();

        final boolean triallelicGenotypesPossible = alleleCount > 2 && ploidy > 2;
        final double[] perReadBuffer = triallelicGenotypesPossible ? new double[readCount] : null;

        // non-log space log10AlleleLikelihoods for multiallelic computation requires rescaling for stability when we
        // exponentiate away the log, and we store the scaling factor to bring back later
        final Pair<double[][], Double> rescaledNonLogLikelihoodsAndCorrection = !triallelicGenotypesPossible ? null :
                rescaledNonLogLikelihoods(log10AlleleLikelihoods);

        final double[] result = new double[genotypeCount];

        for (final GenotypeAlleleCounts gac : GenotypeAlleleCounts.iterable(ploidy, alleleCount)) {
            final int componentCount = gac.distinctAlleleCount();
            final int genotypeIndex = gac.index();
            if (componentCount == 1) {
                // homozygous case: log P(reads|AAAAA. . .) = sum_{reads} log P(read|A)
                final int allele = gac.alleleIndexAt(0);
                result[genotypeIndex] = MathUtils.sum(log10LikelihoodsByAlleleAndRead[allele]);
            } else if (componentCount == 2) {
                // biallelic het case: log P(reads | nA copies of A, nB copies of B) = sum_{reads} (log[(nA * P(read | A) + nB * P(read | B))] -log(ploidy))
                final double[] log10ReadLks1 = log10LikelihoodsByAlleleAndRead[gac.alleleIndexAt(0)];
                final int count1 = gac.alleleCountAt(0);
                final double log10Count1 = MathUtils.log10(count1);
                final double[] log10ReadLks2  = log10LikelihoodsByAlleleAndRead[gac.alleleIndexAt(1)];
                final double log10Count2 = MathUtils.log10(ploidy - count1);

                // note: if you are reading the multiallelic case below and have gotten paranoid about cache efficiency,
                // here the log10 likelihood matrix rows for *both* alleles are in the cache at once
                result[genotypeIndex] = new IndexRange(0, readCount).sum(r -> MathUtils.approximateLog10SumLog10(log10ReadLks1[r] + log10Count1, log10ReadLks2[r] + log10Count2))
                        - readCount * MathUtils.log10(ploidy);
            } else {
                // the multiallelic case is conceptually the same as the biallelic case but done in non-log space
                // We implement in a cache-friendly way by summing nA * P(read|A) over all alleles for each read, but iterating over reads as the inner loop
                Arrays.fill(perReadBuffer,0, readCount, 0);
                final double[][] rescaledNonLogLikelihoods = rescaledNonLogLikelihoodsAndCorrection.getLeft();
                final double log10Rescaling = rescaledNonLogLikelihoodsAndCorrection.getRight();
                gac.forEachAlleleIndexAndCount((a, f) -> new IndexRange(0, readCount).forEach(r -> perReadBuffer[r] += f * rescaledNonLogLikelihoods[a][r]));
                result[genotypeIndex] = new IndexRange(0, readCount).sum(r -> FastMath.log10(perReadBuffer[r])) - readCount * MathUtils.log10(ploidy) + log10Rescaling;
            }
        }
        return result;
    }


    /**
     * Given an input log10 log10Likelihoods matrix, subtract off the maximum of each read column so that each column's maximum is zero for numerical
     * stability.  (This is akin to dividing each read column by its maximum in non-log space).  Then exponentiate to enter non-log space, mutating
     * the log10Likelihoods matrix in-place.  Finally, record the sum of all log-10 subtractions, which is the total amount in log10 space
     * that must later be added to the overall likelihood, which is a sum over all reads (product in npon-log space).
     * @param log10Likelihoods    and input log-10 likelihoods matrix
     */
    private <EVIDENCE, A extends Allele> Pair<double[][], Double> rescaledNonLogLikelihoods(final LikelihoodMatrix<EVIDENCE,A> log10Likelihoods) {
        final double[][] log10LikelihoodsByAlleleAndRead = log10Likelihoods.asRealMatrix().getData();

        final int readCount = log10Likelihoods.evidenceCount();
        final double[] perReadMaxima = new double[readCount];
        Arrays.fill(perReadMaxima, 0, readCount, Double.NEGATIVE_INFINITY);

        // find the maximum log-likelihood over all alleles for each read
        // note how we traverse by read for cache-friendliness
        for (int a = 0; a < alleleCount; a++) {
            for (int r = 0; r < readCount; r++) {
                perReadMaxima[r] = FastMath.max(perReadMaxima[r], log10LikelihoodsByAlleleAndRead[a][r]);
            }
        }

        // subtract these maxima
        for (int a = 0; a < alleleCount; a++) {
            for (int r = 0; r < readCount; r++) {
                log10LikelihoodsByAlleleAndRead[a][r] -= perReadMaxima[r];
            }
        }

        final double scaleFactor = MathUtils.sum(perReadMaxima, 0, readCount);

        // switch to non-log now that we have rescaled for numerical stability
        new IndexRange(0, alleleCount).forEach(a -> MathUtils.applyToArrayInPlace(log10LikelihoodsByAlleleAndRead[a], x -> Math.pow(10.0, x)));

        // note that the variable name is now wrong
        return ImmutablePair.of(log10LikelihoodsByAlleleAndRead, scaleFactor);
    }

    /**
     * Returns the ploidy for this genotype likelihood calculator.
     * @return 0 or greater.
     */
    public int ploidy() {
        return ploidy;
    }

    /**
     * Returns the total number of alleles for this genotype calculator.
     * @return the number of alleles considered by this calculator.
     */
    public int alleleCount() {
        return alleleCount;
    }

    /**
     * Composes a genotype index map given a allele index recoding such that result[i] is the index of the old
     * genotype corresponding to the ith new genotype.
     *
     * @param newToOldAlleleMap allele recoding such that newToOldAlleleMap[i] is the index of the old allele
     *                               corresponding to the ith new allele
     *
     * @throws IllegalArgumentException if this calculator cannot handle the recoding provided. This is
     * the case when either {@code newToOldAlleleMap}'s length or any of its element (+ 1 as they are 0-based) is larger
     * this calculator's {@link #alleleCount()}. Also if any {@code oldToNewAllelesIndexMap} element is negative.
     *
     * @return never {@code null}.
     */
    public int[] newToOldGenotypeMap(final int[] newToOldAlleleMap, final GenotypesCache glCalcs) {
        Utils.nonNull(newToOldAlleleMap);
        final int newAlleleCount = newToOldAlleleMap.length;
        Utils.validateArg(newAlleleCount <= alleleCount,
                () -> String.format("New allele count %d exceeds old allele count %d.", newAlleleCount, alleleCount));
                ;
        final int[] result = new int[GenotypeIndexCalculator.genotypeCount(ploidy, newAlleleCount)];
        for (final GenotypeAlleleCounts newGAC : GenotypeAlleleCounts.iterable(ploidy, newAlleleCount)) {
            result[newGAC.index()] = GenotypeIndexCalculator.alleleCountsToIndex(newGAC, newToOldAlleleMap);
        }

        return result;
    }
}
