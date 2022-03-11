package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.lang3.tuple.Pair;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;

import java.util.Arrays;
import java.util.Iterator;

public class GenotypeLikelihoodCalculator implements Iterable<GenotypeAlleleCounts> {
    /**
     * Genotype table for this calculator.
     *
     * <p>It is ensure that it contains all the genotypes for this calculator ploidy and allele count, maybe more. For
     * that reason you must use {@link #genotypeCount} when iterating through this array and not relay on its length.</p>
     */
    private final GenotypeAlleleCounts[] genotypeAlleleCounts;

    final int genotypeCount;

    final int alleleCount;

    final int ploidy;

    /**
     * Cache of the last genotype-allele-count requested using {@link #genotypeAlleleCountsAt(int)}, when it
     * goes beyond the maximum genotype-allele-count static capacity. Check on that method documentation for details.
     */
    private GenotypeAlleleCounts lastOverheadCounts;

    private static final int INITIAL_READ_CAPACITY = 10;

    /**
     * How many reads the calculator supports.
     *
     * This figure is increased dynamically by {@code ensureReadCapacity}.
     */
    private int readCapacity = INITIAL_READ_CAPACITY;

    /**
     * Buffer field used as a temporary container when one value per read is stored.
     *
     * In the multiallelic calculation we accumulate the likelihood contribution of each read one allele at a time.  That is,
     * for genotype containing alleles A, B, C, we first fill the buffer with the likelihood contributions from allele A, then
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
    private double[] perReadBuffer = new double[INITIAL_READ_CAPACITY];


    public GenotypeLikelihoodCalculator(final int ploidy, final int alleleCount,
                                        final GenotypeAlleleCounts[][] genotypeTableByPloidy) {
        genotypeAlleleCounts = genotypeTableByPloidy[ploidy];
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
     * Returns the genotype associated to a particular likelihood index.
     *
     * <p>If {@code index} is larger than {@link GenotypeLikelihoodCalculators#MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR},
     *  this method will reconstruct that genotype-allele-count iteratively from the largest strongly referenced count available.
     *  or the last requested index genotype.
     *  </p>
     *
     * <p> Therefore if you are iterating through all genotype-allele-counts you should do sequentially and incrementally, to
     * avoid a large efficiency drop </p>.
     *
     * @param index query likelihood-index.
     * @return never {@code null}.
     */
    public GenotypeAlleleCounts genotypeAlleleCountsAt(final int index) {
        Utils.validateArg(index >= 0 && index < genotypeCount, () -> "invalid likelihood index: " + index + " >= " + genotypeCount
                    + " (genotype count for nalleles = " + alleleCount + " and ploidy " + ploidy);
        if (index < GenotypeLikelihoodCalculators.MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR) {
            return genotypeAlleleCounts[index];
        } else if (lastOverheadCounts == null || lastOverheadCounts.index() > index) {
            final GenotypeAlleleCounts result = genotypeAlleleCounts[GenotypeLikelihoodCalculators.MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR - 1].copy();
            result.increase(index - GenotypeLikelihoodCalculators.MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR + 1);
            lastOverheadCounts = result;
            return result.copy();
        } else {
            lastOverheadCounts.increase(index - lastOverheadCounts.index());
            return lastOverheadCounts.copy();
        }
    }

    /**
     * Calculate the likelihoods given the list of alleles and the likelihood map.
     *
     * @param likelihoods the likelihood matrix all alleles vs all reads.
     *
     * @throws IllegalArgumentException if {@code alleleList} is {@code null} or {@code likelihoods} is {@code null}
     *     or the alleleList size does not match the allele-count of this calculator, or there are missing allele vs
     *     read combinations in {@code likelihoods}.
     *
     * @return never {@code null}.
     */
    public <EVIDENCE, A extends Allele> GenotypeLikelihoods genotypeLikelihoods(final LikelihoodMatrix<EVIDENCE, A> likelihoods) {
        final double[] readLikelihoodsByGenotypeIndex = getReadRawReadLikelihoodsByGenotypeIndex(likelihoods);
        return GenotypeLikelihoods.fromLog10Likelihoods(readLikelihoodsByGenotypeIndex);
    }

    /**
     * A helper method that actually does the matrix operations but returns the raw values.
     *
     * @param likelihoods   log 10 likelihood matrix indexed by allele, then read
     * @return the raw array (in log10 likelihoods space) of the GL for each genotype
     */
    <EVIDENCE, A extends Allele> double[] getReadRawReadLikelihoodsByGenotypeIndex(final LikelihoodMatrix<EVIDENCE, A> likelihoods) {
        Utils.nonNull(likelihoods);
        Utils.validateArg(likelihoods.numberOfAlleles() == alleleCount, "mismatch between allele list and alleleCount");
        final int readCount = likelihoods.evidenceCount();
        ensureReadCapacity(readCount);

        final double[][] likelihoodsByAlleleAndRead = likelihoods.asRealMatrix().getData();

        final boolean triallelicGenotypesPossible = alleleCount > 2 && ploidy > 2;

        // non-log space likelihoods for multiallelic computation
        final Pair<double[][], Double> rescaledNonLogLikelihoodsAndCorrection = !triallelicGenotypesPossible ? null :
                rescaledNonLogLikelihoods(readCount, likelihoods);

        final double[] result = new double[genotypeCount];

        Utils.stream(iterator()).forEach(alleleCounts -> {
            final int componentCount = alleleCounts.distinctAlleleCount();
            final int genotypeIndex = alleleCounts.index();
            if (componentCount == 1) {
                // homozygous case: log P(reads|AAAAA. . .) = sum_{reads} log P(read|A)
                final int allele = alleleCounts.alleleIndexAt(0);
                result[genotypeIndex] = MathUtils.sum(likelihoodsByAlleleAndRead[allele]);
            } else if (componentCount == 2) {
                // biallelic het case: log P(reads | nA copies of A, nB copies of B) = sum_{reads} (log[(nA * P(read | A) + nB * P(read | B))] -log(ploidy))
                final double[] log10ReadLks1 = likelihoodsByAlleleAndRead[alleleCounts.alleleIndexAt(0)];
                final int count1 = alleleCounts.alleleCountAt(0);
                final double log10Count1 = MathUtils.log10(count1);
                final double[] log10ReadLks2  = likelihoodsByAlleleAndRead[alleleCounts.alleleIndexAt(1)];
                final double log10Count2 = MathUtils.log10(ploidy - count1);

                result[genotypeIndex] = new IndexRange(0, readCount).sum(r -> MathUtils.approximateLog10SumLog10(log10ReadLks1[r] + log10Count1, log10ReadLks2[r] + log10Count2))
                        - readCount * MathUtils.log10(ploidy);
            } else {
                // the multiallelic case is conceptually the same as the biallelic case but done in non-log space
                // We implement in a cache-friendly way by adding nA * P(read|A) to the per-read buffer for all reads (inner loop), and all alleles (outer loop)
                Arrays.fill(perReadBuffer,0, readCount, 0);
                final double[][] rescaledNonLogLikelihoods = rescaledNonLogLikelihoodsAndCorrection.getLeft();
                final double log10Rescaling = rescaledNonLogLikelihoodsAndCorrection.getRight();
                alleleCounts.forEachAlleleIndexAndCount((a, f) -> new IndexRange(0, readCount).forEach(r -> perReadBuffer[r] += f * rescaledNonLogLikelihoods[a][r]));
                result[genotypeIndex] = new IndexRange(0, readCount).sum(r -> FastMath.log10(perReadBuffer[r])) - readCount * MathUtils.log10(ploidy) + log10Rescaling;
            }
        });
        return result;
    }


    /**
     * Given an input log10 log10Likelihoods matrix, subtract off the maximum of each read column so that each column's maximum is zero for numerical
     * stability.  (This is akin to dividing each read column by its maximum in non-log space).  Then exponentiate to enter non-log space, mutating
     * the log10Likelihoods matrix in-place.  Finally, record the sum of all log-10 subtractions, which is the total amount in log10 space
     * that must later be added to the overall likelihood, which is a sum over all reads (product in npon-log space).
     * @param readCount
     * @param log10Likelihoods    and input log-10 likelihoods matrix
     */
    private <EVIDENCE, A extends Allele> Pair<double[][], Double> rescaledNonLogLikelihoods(final int readCount, final LikelihoodMatrix<EVIDENCE,A> log10Likelihoods) {
        final double[][] log10LikelihoodsByAlleleAndRead = log10Likelihoods.asRealMatrix().getData();
        // fill the per-read buffer with the maximum log-likelihood per read
        // note how we traverse by read for cache-friendliness
        Arrays.fill(perReadBuffer, 0, readCount, Double.NEGATIVE_INFINITY);
        for (int a = 0; a < alleleCount; a++) {
            for (int r = 0; r < readCount; r++) {
                perReadBuffer[r] = FastMath.max(perReadBuffer[r], log10LikelihoodsByAlleleAndRead[a][r]);
            }
        }

        // subtract these maxima
        for (int a = 0; a < alleleCount; a++) {
            for (int r = 0; r < readCount; r++) {
                log10LikelihoodsByAlleleAndRead[a][r] -= perReadBuffer[r];
            }
        }

        final double scaleFactor = MathUtils.sum(perReadBuffer, 0, readCount);

        // switch to non-log now that we have rescaled for numerical stability
        new IndexRange(0, alleleCount).forEach(a -> MathUtils.applyToArrayInPlace(log10LikelihoodsByAlleleAndRead[a], x -> Math.pow(10.0, x)));

        // note that the variable name is now wrong
        return ImmutablePair.of(log10LikelihoodsByAlleleAndRead, scaleFactor);


    }

    private void ensureReadCapacity(final int requestedCapacity) {
        if (readCapacity < requestedCapacity) {
            readCapacity = 2 * requestedCapacity;
            perReadBuffer = new double[readCapacity];
        }
    }


    // note that if the input has a high index that is not cached, it will be mutated in order to form the output
    private GenotypeAlleleCounts nextGenotypeAlleleCounts(final GenotypeAlleleCounts alleleCounts) {
        final int index = alleleCounts.index();
        if (index < (GenotypeLikelihoodCalculators.MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR - 1)) {
            return genotypeAlleleCounts[index + 1];
        } else if (index == GenotypeLikelihoodCalculators.MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR - 1) {
            return genotypeAlleleCounts[index].copy().increase();
        } else {
            return alleleCounts.increase();
        }
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
    public int[] newToOldGenotypeMap(final int[] newToOldAlleleMap, final GenotypeLikelihoodCalculators glCalcs) {
        Utils.nonNull(newToOldAlleleMap);
        final int newAlleleCount = newToOldAlleleMap.length;
        Utils.validateArg(newAlleleCount <= alleleCount,
                () -> String.format("New allele count %d exceeds old allele count %d.", newAlleleCount, alleleCount));
                ;
        final int[] result = new int[GenotypeIndexCalculator.genotypeCount(ploidy, newAlleleCount)];
        for (final GenotypeAlleleCounts newGAC : glCalcs.getInstance(ploidy, newAlleleCount)) {
            result[newGAC.index()] = GenotypeIndexCalculator.alleleCountsToIndex(newGAC, newToOldAlleleMap);
        }

        return result;
    }

    @Override
    public Iterator<GenotypeAlleleCounts> iterator() {
        return new Iterator<GenotypeAlleleCounts>() {
            private int index = 0;
            private GenotypeAlleleCounts alleleCounts = genotypeAlleleCounts[0];

            @Override
            public boolean hasNext() {
                return index < genotypeCount;
            }

            @Override
            public GenotypeAlleleCounts next() {
                alleleCounts = index++ == 0 ? genotypeAlleleCounts[0] : nextGenotypeAlleleCounts(alleleCounts);
                return alleleCounts;
            }
        };
    }
}
