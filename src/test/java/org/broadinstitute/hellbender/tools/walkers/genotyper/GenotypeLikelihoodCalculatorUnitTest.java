package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoodsUnitTester;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.*;

/**
 * Tests {@link GenotypesCache}.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public final class GenotypeLikelihoodCalculatorUnitTest {
    
    @Test(dataProvider = "ploidyAndMaximumAlleleAndReadCountsData")
    public void testLikelihoodCalculation(final int ploidy, final int alleleCount, final int[] readCount) {
        final AlleleLikelihoods<GATKRead, Allele> readLikelihoods = ReadLikelihoodsUnitTester.readLikelihoods(alleleCount, readCount);
        final GenotypeLikelihoodCalculator calculator = new GenotypesCache().getInstance(ploidy, alleleCount);
        final int genotypeCount = calculator.genotypeCount();
        final int testGenotypeCount = Math.min(30000, genotypeCount);
        final int sampleCount = readCount.length;
        for (int s = 0; s < sampleCount ; s++) {
            final LikelihoodMatrix<GATKRead, Allele> sampleLikelihoods = readLikelihoods.sampleMatrix(s);
            final GenotypeLikelihoods genotypeLikelihoods = calculator.log10GenotypeLikelihoods(sampleLikelihoods);
            final double[] genotypeLikelihoodsDoubles = genotypeLikelihoods.getAsVector();
            Assert.assertEquals(genotypeLikelihoodsDoubles.length, genotypeCount);
            for (final GenotypeAlleleCounts gac : GenotypeAlleleCounts.iterable(ploidy, alleleCount)) {
                Assert.assertNotNull(genotypeLikelihoods);
                final double[] readGenotypeLikelihoods = new double[sampleLikelihoods.evidenceCount()];
                for (int r = 0; r < sampleLikelihoods.evidenceCount(); r++) {
                    final double[] compoments = new double[gac.distinctAlleleCount()];
                    for (int ar = 0; ar < gac.distinctAlleleCount(); ar++) {
                        final int a = gac.alleleIndexAt(ar);
                        final int aCount = gac.alleleCountAt(ar);
                        final double readLk = sampleLikelihoods.get(a, r);
                        compoments[ar] = readLk + Math.log10(aCount);
                    }
                    readGenotypeLikelihoods[r] = MathUtils.approximateLog10SumLog10(compoments) - Math.log10(ploidy);
                }
                final double genotypeLikelihood = MathUtils.sum(readGenotypeLikelihoods);
                Assert.assertEquals(genotypeLikelihoodsDoubles[gac.index()], genotypeLikelihood, 0.0001 * Math.abs(genotypeLikelihood));
            }
        }
    }

    @Test(dataProvider = "ploidyAndMaximumAlleleAndNewMaximumAlleleData")
    public void testGenotypeIndexMap(final int ploidy, final int oldAlleleCount, final int newAlleleCount) {

        // TODO: placeholder to cause failure because this method does not belong in this class
        Assert.assertTrue(false);


        Utils.resetRandomGenerator();
        final Random rnd = Utils.getRandomGenerator();
        final int maxAlleleCount = Math.max(oldAlleleCount, newAlleleCount);
        final int[] alleleMap = new int[newAlleleCount];
        final Map<Integer,Set<Integer>> reverseMap = new LinkedHashMap<>(oldAlleleCount);
        for (int i = 0; i < alleleMap.length; i++) {
            alleleMap[i] = rnd.nextInt(oldAlleleCount);
            if (reverseMap.get(alleleMap[i]) == null) reverseMap.put(alleleMap[i],new LinkedHashSet<>(6));
            reverseMap.get(alleleMap[i]).add(i);
        }
        final GenotypesCache calculators = new GenotypesCache();
        final GenotypeLikelihoodCalculator calculator = calculators.getInstance(ploidy, maxAlleleCount);

        final int[] genotypeIndexMap = calculator.newToOldGenotypeMap(alleleMap, calculators);
        Assert.assertNotNull(genotypeIndexMap);
        Assert.assertEquals(genotypeIndexMap.length, GenotypeIndexCalculator.genotypeCount(ploidy, newAlleleCount));

        final GenotypeLikelihoodCalculator oldCalculator = calculators.getInstance(ploidy, oldAlleleCount);
        final GenotypeLikelihoodCalculator newCalculator = calculators.getInstance(ploidy, newAlleleCount);

        for (int i = 0; i < genotypeIndexMap.length; i++) {
            final GenotypeAlleleCounts oldCounts = null;
                    //oldCalculator.genotypeAlleleCountsAt(genotypeIndexMap[i]);
            final GenotypeAlleleCounts newCounts = null;
                    //newCalculator.genotypeAlleleCountsAt(i);
            final int[] reverseCounts = new int[oldAlleleCount];
            for (int j = 0; j < newCounts.distinctAlleleCount(); j++) {
                final int newIndex = newCounts.alleleIndexAt(j);
                final int newRepeats = newCounts.alleleCountAt(j);
                final int expectedOldIndex = alleleMap[newIndex];
                final int oldIndexRank = oldCounts.alleleRankFor(expectedOldIndex);
                Assert.assertNotEquals(oldIndexRank, -1);
                final int oldIndex = oldCounts.alleleIndexAt(oldIndexRank);
                final int oldRepeats = oldCounts.alleleCountAt(oldIndexRank);
                Assert.assertEquals(oldIndex, expectedOldIndex);
                // not necessarily the same count if two or more new alleles map the same old allele.
                Assert.assertTrue(oldRepeats >= newRepeats);
                reverseCounts[oldIndex] += newRepeats;
            }
            for (int j = 0; j < oldAlleleCount; j++)
                Assert.assertEquals(oldCounts.alleleCountFor(j), reverseCounts[j]);
        }
    }


    private static final int[] MAXIMUM_ALLELE = { 1, 2, 5, 6};

    private static final int[] PLOIDY = { 1, 2, 3, 20 };

    private static final int[][] READ_COUNTS = {
            { 10 , 100, 50 },
            { 0, 100, 10, 1 , 50 },
            { 1, 2, 3, 4, 20 },
            { 10, 0 },
    };

    @DataProvider(name="ploidyAndMaximumAlleleAndReadCountsData")
    public Object[][] ploidyAndMaximumAlleleAndReadCountsData() {
        final Object[][] result = new Object[PLOIDY.length * MAXIMUM_ALLELE.length * READ_COUNTS.length][];
        int index = 0;
        for (final int i : PLOIDY)
            for (final int j : MAXIMUM_ALLELE)
                for (final int[] k : READ_COUNTS)
                result[index++] = new Object[] { i, j, k };
        return result;
    }

    @DataProvider(name="ploidyAndMaximumAlleleAndNewMaximumAlleleData")
    public Object[][] ploidyAndMaximumAlleleAndNewMaximumAlleleData() {
        final List<Object[]> result = new ArrayList<>(PLOIDY.length * MAXIMUM_ALLELE.length * 20);
        for (final int i : PLOIDY)
            for (final int j : MAXIMUM_ALLELE)
                for (int k = 0; k < (i < 10? j * 2 : j + 1); k++)
                    result.add(new Object[] { i, j, k });
        return result.toArray(new Object[result.size()][]);
    }
}
