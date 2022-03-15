package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

/**
 * This class holds caches of {@link GenotypeAlleleCounts} for multiple fixed ploidy, allele count pairs,
 * allowing for fast random access of genotypes.  Note that the increment method of GenotypeAlleleCounts is always fast,
 * so the caches here are only necessary when incremental traversal over genotypes in the canonical order is not possible.
 *
 *
 * This class is thread-safe since modifying the caches is synchronized.
 */
public final class GenotypeLikelihoodCalculators {

    private static final Logger logger = LogManager.getLogger(GenotypeLikelihoodCalculators.class);

    /**
     * The current maximum ploidy supported by the tables.  Initial value may be anything positive.
     */
    private int maximumPloidy = 2;

    /**
     * Maximum possible number of cached {@link GenotypeAlleleCounts} for each fixed ploidy and allele count.
     */
    public static final int MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR = 1000;
    

    /**
     * The current maximum allele index supported by the tables.
     */
    private int maximumAllele = 1;


    /**
     * Cache of GenotypeAlleleCounts objects by ploidy.  Format is table[p][n] = nth genotype of ploidy p in canonical order,
     * with p up to the current maximmum ploidy and n up to the maximum number of cached genotypes per table.
     */
    private GenotypeAlleleCounts[][] genotypeTableByPloidy = createCaches(maximumPloidy,maximumAllele);

    public GenotypeLikelihoodCalculators(){

    }

    /**
     * Composes a table with the lists of all possible genotype allele counts given the the ploidy and maximum allele index.
     *
     * This method is synchronized since it modifies the shared cache.
     */
    private static synchronized GenotypeAlleleCounts[][] createCaches(final int maximumPloidy, final int maximumAllele) {
        Utils.validateArg(maximumPloidy >= 0, () -> "the ploidy provided cannot be negative: " + maximumPloidy);
        Utils.validateArg(maximumAllele >= 0, () -> "the maximum allele index provided cannot be negative: " + maximumAllele);
        final GenotypeAlleleCounts[][] result = new GenotypeAlleleCounts[maximumPloidy + 1][]; // each row has a different number of columns.

        for (int ploidy = 0; ploidy <= maximumPloidy; ploidy++) {
            final int numberOfCachedGenotypes = Math.min(GenotypeIndexCalculator.genotypeCount(ploidy, maximumAllele), MAXIMUM_CACHED_GENOTYPES_PER_CALCULATOR);
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
        ensureCapacity(alleleCount, ploidy);
        return new GenotypeLikelihoodCalculator(ploidy, alleleCount, genotypeTableByPloidy);
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
        ensureCapacity(alleleCount, ploidy);
        return new GenotypeLikelihoodCalculatorDRAGEN(ploidy, alleleCount, genotypeTableByPloidy);
    }


    /**
     * Update cache if necessary
     */
    private void ensureCapacity(final int requestedMaximumAllele, final int requestedMaximumPloidy) {
        if (requestedMaximumAllele <= maximumAllele && requestedMaximumPloidy <= maximumPloidy) {
            return;
        }

        maximumPloidy = Math.max(maximumPloidy, requestedMaximumPloidy);
        maximumAllele = Math.max(maximumAllele, requestedMaximumAllele);
        logger.debug("Expanding capacity ploidy:" + maximumPloidy + "->" + maximumPloidy + " allele:" +  maximumAllele +"->" + maximumAllele );
        genotypeTableByPloidy = createCaches(maximumPloidy, maximumAllele);
    }

}