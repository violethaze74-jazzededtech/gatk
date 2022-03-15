package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Iterator;

/**
 * This class holds caches of {@link GenotypeAlleleCounts} for multiple fixed ploidy, allele count pairs,
 * allowing for fast random access of genotypes.  Note that the increment method of GenotypeAlleleCounts is always fast,
 * so the caches here are only necessary when incremental traversal over genotypes in the canonical order is not possible.
 *
 *
 * This class is thread-safe since modifying the caches is synchronized.
 */
public final class GenotypesCache {

    private static final Logger logger = LogManager.getLogger(GenotypesCache.class);

    /**
     * The current maximum ploidy supported by the tables.  Initial value may be anything positive.
     */
    private static int maximumPloidy = 2;

    /**
     * The current number of cached genotypes per ploidy
     */
    private static int cacheSizePerPloidy = 10;

    /**
     * Maximum possible number of cached {@link GenotypeAlleleCounts} for each fixed ploidy and allele count.
     */
    public static final int MAX_CACHE_SIZE_PER_PLOIDY = 1000;

    /**
     * Cache of GenotypeAlleleCounts objects by ploidy.  Format is table[p][n] = nth genotype of ploidy p in canonical order,
     * with p up to the current maximum ploidy and n up to the maximum number of cached genotypes per table.
     */
    private static GenotypeAlleleCounts[][] caches = createCaches(maximumPloidy, cacheSizePerPloidy);

    public GenotypesCache(){

    }

    /**
     * Returns the GenotypeAlleleCounts associated to a particular ploidy and genotype index.
     *
     *  If the requested index is larger than {@link GenotypesCache#MAX_CACHE_SIZE_PER_PLOIDY},
     *  this method will construct the result iteratively from the largest cached object.  Thus if you are iterating
     *  through all genotype-allele-counts you should do sequentially using the iterator method to avoid a large efficiency drop.
     *
     * @param ploidy the ploidy
     * @param genotypeIndex  the genotype index in the canonical order
     * @return never {@code null}.
     */
    public static GenotypeAlleleCounts get(final int ploidy, final int genotypeIndex) {
        ensureCapacity(genotypeIndex, ploidy);
        Utils.validateArg(ploidy >= 0, "ploidy may not be negative");
        Utils.validateArg(genotypeIndex >= 0, "genotype index may not be negative");
        if (genotypeIndex < cacheSizePerPloidy) {
            return caches[ploidy][genotypeIndex];
        } else {
            final GenotypeAlleleCounts result = caches[ploidy][cacheSizePerPloidy - 1].copy();
            result.increase(genotypeIndex + 1 - cacheSizePerPloidy);
            return result;
        }
    }

    /**
     * Composes a table with the lists of all possible genotype allele counts given the the ploidy and maximum allele index.
     *
     * This method is synchronized since it modifies the shared cache.
     */
    private static synchronized GenotypeAlleleCounts[][] createCaches(final int maximumPloidy, final int cachedGenotypesPerPloidy) {
        Utils.validateArg(maximumPloidy >= 0, () -> "the ploidy provided cannot be negative: " + maximumPloidy);
        Utils.validateArg(cachedGenotypesPerPloidy >= 0, () -> "the cache size provided cannot be negative: " + cachedGenotypesPerPloidy);
        final GenotypeAlleleCounts[][] result = new GenotypeAlleleCounts[maximumPloidy + 1][];

        for (int ploidy = 0; ploidy <= maximumPloidy; ploidy++) {
            final int numberOfCachedGenotypes = Math.min(cachedGenotypesPerPloidy, MAX_CACHE_SIZE_PER_PLOIDY);
            final GenotypeAlleleCounts[] cache = new GenotypeAlleleCounts[numberOfCachedGenotypes];
            cache[0] = GenotypeAlleleCounts.first(ploidy);
            for (int genotypeIndex = 1; genotypeIndex < numberOfCachedGenotypes; genotypeIndex++) {
                cache[genotypeIndex] = cache[genotypeIndex - 1].next();
            }
            result[ploidy] = cache;
        }

        return result;
    }


    /**
     * Returns an instance given its ploidy and the number of alleles.
     *
     * @param alleleCount the required allele-count.
     * @param ploidy the required ploidy-count.
     *
     * @throws IllegalArgumentException if either {@code ploidy} or {@code alleleCount} is negative, or the resulting number of genotypes is too large.
     *
     * @return never {@code null}.
     */
    public synchronized GenotypeLikelihoodCalculator getInstance(final int ploidy, final int alleleCount) {
        //ensureCapacity(alleleCount, ploidy);
        return new GenotypeLikelihoodCalculator(ploidy, alleleCount, caches);
    }

    /**
     * Returns an instance of the DRAGEN genotypeLikelihoodCalculator given its ploidy and the number of alleles.
     *
     * @param alleleCount the required allele-count.
     * @param ploidy the required ploidy-count.
     *
     * @throws IllegalArgumentException if either {@code ploidy} or {@code alleleCount} is negative, or the resulting number of genotypes is too large.
     *
     * @return never {@code null}.
     */
    public synchronized GenotypeLikelihoodCalculatorDRAGEN getInstanceDRAGEN(final int ploidy, final int alleleCount) {
        Utils.validate(ploidy == 2, "DRAGEN genotyping mode currently only supports diploid samples");
        //ensureCapacity(alleleCount, ploidy);
        return new GenotypeLikelihoodCalculatorDRAGEN(ploidy, alleleCount, caches);
    }


    /**
     * Update cache if necessary
     */
    private static void ensureCapacity(final int genotypeIndex, final int requestedMaximumPloidy) {
        if ((genotypeIndex < cacheSizePerPloidy || cacheSizePerPloidy == MAX_CACHE_SIZE_PER_PLOIDY) && requestedMaximumPloidy <= maximumPloidy) {
            return;
        }

        maximumPloidy = Math.max(maximumPloidy, requestedMaximumPloidy);
        cacheSizePerPloidy = Math.min(Math.max(cacheSizePerPloidy, genotypeIndex + 1), MAX_CACHE_SIZE_PER_PLOIDY);
        caches = createCaches(maximumPloidy, cacheSizePerPloidy);
    }

}