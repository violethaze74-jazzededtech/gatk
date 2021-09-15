package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.apache.commons.collections4.CollectionUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.genotyper.AlleleSubsettingUtils;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAssignmentMethod;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.util.*;
import java.util.stream.Collectors;

/**
 * input variants should be tagged "actual" or "expected" -- single sample ONLY
 */
@BetaFeature
@CommandLineProgramProperties(
        summary = "Compare two VCFs, as for pipeline test updates",
        oneLineSummary = "Compare two VCFs",
        programGroup = VariantEvaluationProgramGroup.class
)
public class VCFComparator extends MultiVariantWalkerGroupedByOverlap {
    private boolean sampleCheckDone = false;
    private boolean isSingleSample = true;
    private boolean alleleNumberIsDifferent = false;
    private boolean inbreedingCoeffIsDifferent = false;
    private List<String> annotationsThatVaryWithNoCalls = Arrays.asList(VCFConstants.ALLELE_NUMBER_KEY, GATKVCFConstants.INBREEDING_COEFFICIENT_KEY,
            GATKVCFConstants.AS_INBREEDING_COEFFICIENT_KEY, GATKVCFConstants.EXCESS_HET_KEY,
            GATKVCFConstants.MLE_ALLELE_COUNT_KEY, GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY, VCFConstants.ALLELE_FREQUENCY_KEY,
            VCFConstants.DEPTH_KEY, GATKVCFConstants.QUAL_BY_DEPTH_KEY, GATKVCFConstants.AS_QUAL_BY_DEPTH_KEY,
            GATKVCFConstants.VQS_LOD_KEY, GATKVCFConstants.AS_VQS_LOD_KEY, GATKVCFConstants.AS_FILTER_STATUS_KEY,
            GATKVCFConstants.CULPRIT_KEY, GATKVCFConstants.AS_CULPRIT_KEY); //AC should be the same

    @Argument(fullName = "warn-on-errors",
            shortName = "warn-on-errors",
            doc = "just emit warnings on errors instead of terminating the run at the first instance",
            optional = true)
    Boolean WARN_ON_ERROR = false;

    @Argument(fullName = "ignore-quals", optional = true)
    Boolean IGNORE_QUALS = false;

    @Argument(fullName = "dp-change-allowed", doc = "Note that this is a signed change (actual - expected)")
    Integer DP_CHANGE = 0;

    @Argument(fullName = "ignore-annotations", doc = "Only match on position and alleles, ignoring annotations")
    Boolean IGNORE_ANNOTATIONS = false;

    @Argument(fullName = "ignore-genotype-annotations", doc = "Only match on genotype call")
    Boolean IGNORE_GENOTYPE_ANNOTATIONS = false;

    @Argument(fullName = "ignore-filters", doc = "Ignore filter status when comparing variants")
    Boolean IGNORE_FILTERS = false;

    @Argument(fullName = "ignore-attribute", doc = "Ignore INFO attributes with this key")
    List<String> IGNORE_ATTRIBUTES = new ArrayList<>(5);

    @Argument(fullName = "positions-only", doc = "Only match on position, ignoring alleles and annotations")
    Boolean POSITONS_ONLY = false;

    @Argument(fullName = "allow-new-stars", doc = "Allow additional * alleles in actual if there is a corresponding deletion")
    Boolean ALLOW_NEW_STARS = false;

    //TODO: this should have a mutex wrt the above
    @Argument(fullName = "allow-extra-alleles", doc = "Allow extra alleles in actual provided actual is a superset of expected")
    Boolean ALLOW_EXTRA_ALLELES = false;

    @Argument(fullName = "allow-missing-stars", doc = "Allow extra * in expected if actual has no corresponding deletions")
    Boolean ALLOW_MISSING_STARS = false;

    @Override
    public boolean doDictionaryCrossValidation() {
        return false;
    }

    @Override
    public void onTraversalStart() {
        if (getDrivingVariantsFeatureInputs().size() != 2) {
            throw new GATKException("VCFComparator expects exactly two inputs -- one actual and one expected.");
        }
    }

    @Override
    public void apply(final List<VariantContext> variantContexts, final ReferenceContext referenceContext, final List<ReadsContext> readsContexts) {
        if (!sampleCheckDone) {
            if (variantContexts.get(0).getGenotypes().size() > 1) {
                isSingleSample = false;
            }
            sampleCheckDone = true;
        }

        /*if (variantContexts.size() == 1) {
            final VariantContext vc = variantContexts.get(0);
            if (vc.getPhredScaledQual() > 60.0 && vc.isBiallelic() && vc.getAttributeAsInt(VCFConstants.ALLELE_COUNT_KEY, 0) > 1) {
                throwOrWarn(new UserException("Unmatched variant in " + vc.getSource() + " at position " + vc.getContig() + ":" + vc.getStart()));
            }
        }*/

        //ideally we have all variants overlapping a given position
        //if there's no actual for the expected (at the same position), see if the expected has an overlapping deletion
        for (final VariantContext vc : variantContexts) {
            if (vc.getSource().equals("actual")) {  //there may be more expected than actual variants, but only compare once
                continue;
            }
            final List<VariantContext> matches = variantContexts.stream().filter(v -> v.getStart() == vc.getStart()).collect(Collectors.toList());
            if (matches.size() != 2 && isHighQuality(vc)
            && !vc.getGenotype(0).getAlleles().contains(Allele.SPAN_DEL)) {
                throwOrWarn(new UserException("Apparent unmatched high quality variant in " + vc.getSource() + " at " + vc.getContig() + ":" + vc.getStart()));
            } else {
                if (!isHighQuality(vc) ||  (matches.size() < 2 && vc.getGenotype(0).getAlleles().contains(Allele.SPAN_DEL))) { //TODO: check * genotypes somehow
                    return;
                }
                final VariantContext match;
                if (!matches.get(0).getSource().equals(vc.getSource())) {
                    match = matches.get(0);
                } else {
                    match = matches.get(1);
                }

                final VariantContext trimmed1 = trimAlleles(vc);
                final VariantContext trimmed2 = trimAlleles(match);

                final List<VariantContext> overlappingDels = variantContexts.stream().filter(v -> v.getStart() < vc.getStart() && v.overlaps(vc))
                        .collect(Collectors.toList());

                //do single-sample GVCF checks, including deletion trimming and dropping
                if (isSingleSample) {
                    final List<Allele> gt1 = trimmed1.getGenotype(0).getAlleles();
                    final List<Allele> gt2 = trimmed2.getGenotype(0).getAlleles();

                    if (!gt1.get(0).equals(gt2.get(0)) || !gt1.get(1).equals(gt2.get(1))) {  //do use order here so we can check phasing
                        //this is okay if a star got corrected by dropping an upstream hom ref del
                        //could be dropped hom ref
                        if (overlappingDels.stream().anyMatch(v -> v.getGenotype(0).isHomRef())) {
                            if (!areAttributesEqual(vc.getAttributes(), match.getAttributes())) {
                                throwOrWarn(new UserException("INFO attributes do not match at " + vc.getContig() + ":" + vc.getStart()));
                            }
                            if (!areGenotypesEqual(trimmed1.getGenotype(0), trimmed2.getGenotype(0))) {
                                throwOrWarn(new UserException("Genotypes do not match at " + vc.getContig() + ":" + vc.getStart()));
                            }
                            return;
                        } else if (variantContexts.size() == 2 && gt1.contains(Allele.SPAN_DEL) && !gt2.contains(Allele.SPAN_DEL)) {
                            return;  //there never was an overlapping deletion in the expected
                        } else {
                            if (overlappingDels.size() != 1) {  //there should be one overlapping that got trimmed and the other shouldn't overlap anymore

                                throwOrWarn(new UserException("Genotype alleles do not match at " + vc.getContig() + ":" + vc.getStart()
                                        + ". " + vc.getSource() + " has " + trimmed1.getGenotype(0).toString() + " and " + match.getSource() + " has "
                                        + trimmed2.getGenotype(0).toString()));
                                return;
                            } else {
                                final VariantContext overlapper = overlappingDels.get(0);
                                final VariantContext trimmedOverlapper = trimAlleles(overlapper);
                                if (overlapper.overlaps(vc) && !trimmedOverlapper.overlaps(vc)) {
                                    return;
                                }
                            }

                        }
                    }
                }
                if (POSITONS_ONLY) {
                    return;
                }

                //more rigorous checks
                final VariantContext actual_trimmed, expected_trimmed;
                if (trimmed1.getSource().equals("actual")) {
                    actual_trimmed = trimmed1;
                    expected_trimmed = trimmed2;
                } else {
                    actual_trimmed = trimmed2;
                    expected_trimmed = trimmed1;
                }
                try {
                    if (!areVariantContextsEqualOrderIndependent(actual_trimmed, expected_trimmed, overlappingDels) && !alleleNumberIsDifferent
                    && !expected_trimmed.getGenotypes().stream().anyMatch(g -> g.getGQ() < 20)) {
                        throwOrWarn(new UserException("Variant contexts do not match: " + actual_trimmed.toStringDecodeGenotypes() + " versus " + expected_trimmed.toStringDecodeGenotypes()));
                    }
                } catch (UserException e) {
                    throwOrWarn(new UserException("Variant contexts do not match at " + actual_trimmed.getContig() + ":" + actual_trimmed.getStart(), e));
                }
                if (alleleNumberIsDifferent) {
                    logger.warn("Observed allele number differed at position " + vc.getContig() + ":" + vc.getStart());
                }
                if (inbreedingCoeffIsDifferent && expected_trimmed.getGenotypes().stream().anyMatch(g -> g.getGQ() < 20)) {
                    logger.warn("Low quality genotype may have caused inbreeing coeff differences at position " + vc.getContig() + ":" + vc.getStart());
                }
            }
        }
    }

    private boolean areVariantContextsEqualOrderIndependent(final VariantContext actual, final VariantContext expected,
                                                            final List<VariantContext> overlappingDels) {
        if (!actual.getContig().equals(expected.getContig())) {
            return false;
        }
        if (actual.getStart() != expected.getStart()) {
            return false;
        }

        //don't check end in case we're being lenient about alleles

        if (!ALLOW_EXTRA_ALLELES && !actual.getID().equals(expected.getID())) {  //more alleles might mean more dbSNP matches
            return false;
        }

        if (!expected.getAlternateAlleles().stream().allMatch(a -> GATKVariantContextUtils.isAlleleInList(
                expected.getReference(), a, actual.getReference(), actual.getAlternateAlleles()))) {
            if (!ALLOW_EXTRA_ALLELES && !actual.getAlternateAlleles().stream().allMatch(a -> GATKVariantContextUtils.isAlleleInList(
                    actual.getReference(), a, expected.getReference(), expected.getAlternateAlleles()))) {
                throwOrWarn(new UserException("Alleles are mismatched at " + actual.getContig() + ":" + actual.getStart() + ": actual has "
                        + actual.getAlternateAlleles() + " and expected has " + expected.getAlternateAlleles()));

            } else if (ALLOW_NEW_STARS && actual.getAlleles().contains(Allele.SPAN_DEL) && !expected.getAlleles().contains(Allele.SPAN_DEL)) {
                if (overlappingDels.size() == 0 || !overlappingDels.stream().anyMatch(vc -> vc.getSource().equals("actual"))) {
                    throwOrWarn(new UserException("Actual has new unmatched * allele. Alleles are mismatched at " + actual.getContig() + ":" + actual.getStart() + ": actual has "
                            + actual.getAlternateAlleles() + " and expected has " + expected.getAlternateAlleles()));
                }
            //this is a GenomicsDB/CombineGVCFs bug -- there is an overlapping deletion and * should be output, but those tools don't account for multiple haplotypes and upstream variant "ends" the deletion
            } else if (ALLOW_MISSING_STARS && expected.getAlleles().contains(Allele.SPAN_DEL) && !actual.getAlleles().contains(Allele.SPAN_DEL)) {
                final Set<Allele> remainder = new LinkedHashSet<>(expected.getAlleles());
                remainder.removeAll(actual.getAlleles());
                if (remainder.size() > 1 || !remainder.contains(Allele.SPAN_DEL)) {
                    throwOrWarn(new UserException("Actual missing * allele. Alleles are mismatched at " + actual.getContig() + ":" + actual.getStart() + ": actual has "
                            + actual.getAlternateAlleles() + " and expected has " + expected.getAlternateAlleles()));
                }
        } else {
                throwOrWarn(new UserException("Alleles are mismatched at " + actual.getContig() + ":" + actual.getStart() + ": actual has "
                        + actual.getAlternateAlleles() + " and expected has " + expected.getAlternateAlleles()));
            }
        }

        if (IGNORE_ANNOTATIONS) {
            return true;
        }

        if (Math.abs(actual.getPhredScaledQual() - expected.getPhredScaledQual()) > 0.001) {
            if (!IGNORE_QUALS) {
                throw new UserException("Qual scores differ by more than 0.001");
            }
        }

        if (!areAttributesEqual(actual.getAttributes(), expected.getAttributes())) {
            return false;
        }

        if (!alleleNumberIsDifferent) {
            if (!IGNORE_FILTERS) {
                if (actual.filtersWereApplied() != expected.filtersWereApplied()) {
                    return false;
                }
                if (actual.isFiltered() != expected.isFiltered()) {
                    return false;
                }
                if (!actual.getFilters().equals(expected.getFilters())) {
                    return false;
                }
            }


            if (!areGenotypesEqual(actual.getGenotype(0), expected.getGenotype(0))) {
                throw new UserException("Genotypes are not equal: " + actual.getGenotype(0) + " versus " + expected.getGenotype(0));
            }
        }
        return true;
    }

    @SuppressWarnings("unchecked")
    private boolean areAttributesEqual(final Map<String, Object> actual, final Map<String, Object> expected) {
        final Set<String> expectedKeys = new LinkedHashSet<>(expected.keySet());

        //do a precheck on AN because then we can't expect the rest of the annotations to match
        if (!isAttributeValueEqual(VCFConstants.ALLELE_NUMBER_KEY, actual.get(VCFConstants.ALLELE_NUMBER_KEY),
                expected.get(VCFConstants.ALLELE_NUMBER_KEY))) {
            alleleNumberIsDifferent = true;
        } else {
            alleleNumberIsDifferent = false;
        }

        if (actual.containsKey(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY)
                &&!isAttributeEqualDoubleSmart(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY, Double.parseDouble(actual.get(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY).toString()),
                Double.parseDouble(expected.get(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY).toString()), 0.001)) {
            inbreedingCoeffIsDifferent = true;
        } else {
            inbreedingCoeffIsDifferent = false;
        }

        for ( final Map.Entry<String, Object> exp : expected.entrySet() ) {
           final Object expectedValue = exp.getValue();
            final String key = exp.getKey();
            if ( actual.containsKey(key) && actual.get(key) != null ) {
                //we already checked these
                if (key.equals(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY) || key.equals(VCFConstants.ALLELE_NUMBER_KEY)
                || IGNORE_ATTRIBUTES.contains(key)) {
                    continue;
                }
                final Object actualValue = actual.get(key);
                if (expectedValue instanceof List && actualValue instanceof List) {
                    // both values are lists, compare element by element
                    List<Object> expectedList = (List<Object>) expectedValue;
                    List<Object> actualList = (List<Object>) actualValue;
                    if (actualList.size() != expectedList.size()) {
                        throw new UserException("Variant contexts have different attributes: actual has " + actual.keySet().toString()
                                + "\n expected has " + expected.keySet().toString());
                    }
                    for (int i = 0; i < expectedList.size(); i++) {
                        if (!isAttributeValueEqual(key, actualList.get(i), expectedList.get(i))) {
                            if (key.equals(VCFConstants.ALLELE_COUNT_KEY) || key.equals(GATKVCFConstants.MLE_ALLELE_COUNT_KEY)) {
                                return isACEqualEnough(key, )
                            } else {
                                return false;
                            }
                        }
                    }
                } else {
                    if (!isAttributeValueEqual(key, actualValue, expectedValue)) {
                        if (annotationsThatVaryWithNoCalls.contains(key) && alleleNumberIsDifferent) {
                            continue;
                        } else {
                            if (key.equals(GATKVCFConstants.QUAL_BY_DEPTH_KEY)) {
                                final double actualDouble = Double.parseDouble(actualValue.toString());
                                final double expectedDouble = Double.parseDouble(expectedValue.toString());
                                final double diff = Math.abs(expectedDouble - actualDouble);
                                final double relativeDiff = diff / (expectedDouble);
                                if (expectedDouble > 25.0 || relativeDiff < 0.01 || diff < 0.5) {  //25 is in the "jitter" zone
                                    continue;
                                } else {
                                    return false;
                                }
                            }
                        }
                    }
                }
            } else {
                // it's ok if there are new things in actual
            }
            expectedKeys.remove(key);
        }
        return true;
    }

    private boolean isAttributeValueEqual(final String key, final Object actual, final Object expected) {
           if (!actual.toString().equals(expected.toString())) {
               throwOrWarn(new UserException("Variant contexts have different attribute values for " + key + ": actual has " + actual.toString()
                       + " and expected has " + expected.toString()));
           }
           return actual.toString().equals(expected.toString());

    }

    private boolean isAttributeEqualDoubleSmart(final String key, final double actual, final double expected, final double tolerance) {
        final double diff = Math.abs(actual-expected);
        if (diff > tolerance) {
            throwOrWarn(new UserException("Attribute " + key + " has difference " + diff + ", which is larger difference than allowed delta " + tolerance));
            return false;
        }
        return true;
    }

    private boolean areGenotypesEqual(final Genotype actual, final Genotype expected) {
        if (!actual.getSampleName().equals(expected.getSampleName())) {
            return false;
        }
        if (!CollectionUtils.isEqualCollection(actual.getAlleles(), expected.getAlleles())) {
            return false;
        }
        if (!actual.getGenotypeString(false).equals(expected.getGenotypeString(false))){
            return false;
        }
        if (actual.isPhased() != expected.isPhased()) {
            return false;
        }
        if (IGNORE_GENOTYPE_ANNOTATIONS) {
            return true;
        }
        if (actual.hasDP() != expected.hasDP()) {
            return false;
        }
        //allow DP difference due to non-ref AD correction
        if (actual.hasDP() && expected.hasDP()) {
            if (actual.getDP() != expected.getDP() ) {
                if (DP_CHANGE == 0 && expected.getDP() - actual.getDP() != expected.getAD()[expected.getAD().length-1]) {
                    return false;
                }
                if (actual.getDP() - expected.getDP() > DP_CHANGE) {
                    return false;
                }
            }

        }
        if (actual.hasAD() != expected.hasAD()) {
            return false;
        }
        if (actual.hasAD() && (!Arrays.equals(actual.getAD(), expected.getAD()))) {
            if (DP_CHANGE == 0 && actual.getAD()[actual.getAD().length-1] == expected.getAD()[expected.getAD().length-1]) {  //non-ref differences are okay, but if there's an AD different that's not non-ref we're in trouble
                return false;
            }
            for (int i = 0; i < actual.getAD().length-1; i++) {
                if (actual.getAD()[i] - expected.getAD()[i] > DP_CHANGE) {
                    return false;
                }
            }
        }
        if (actual.hasGQ() != expected.hasGQ()) {
            return false;
        }
        if (actual.hasGQ() && (actual.getGQ() != expected.getGQ())) {
            return false;
        }
       if (actual.hasPL() != expected.hasPL()) {
            return false;
        }
        if (actual.hasPL() && (!Arrays.equals(actual.getPL(), (expected.getPL())))) {
            return false;
        }
        return true;
    }

    private void throwOrWarn(UserException e) {
        if (WARN_ON_ERROR) {
            logger.warn("***** " + e.getMessage() + " *****");
        } else {
            throw e;
        }
    }

    private VariantContext trimAlleles(final VariantContext variant) {
        final Allele ref = variant.getReference();
        final Set<Allele> relevantAlleles = new LinkedHashSet<>();
        for (final Genotype g : variant.getGenotypes()) {
            final List<Allele> gtAlleles = g.getAlleles();
            relevantAlleles.add(ref);
            if (!gtAlleles.contains(Allele.NO_CALL)) {  //in theory there could be half no-calls, but I've never seen one
                relevantAlleles.addAll(gtAlleles);
            }
        }
        final VariantContextBuilder vcBuilder = new VariantContextBuilder(variant);
        vcBuilder.alleles(relevantAlleles);
        final List<Allele> orderedRelevantAlleles = new ArrayList<>();
        orderedRelevantAlleles.addAll(relevantAlleles);
        if (vcBuilder.getGenotypes().size() == 1) {
            orderedRelevantAlleles.add(Allele.NON_REF_ALLELE);
        }
        //NOTE that we use BEST_MATCH_TO_ORIGINAL for post-reblocked VCFs with no hom ref PLs
        final GenotypesContext gc = AlleleSubsettingUtils.subsetAlleles(variant.getGenotypes(), variant.getGenotype(0).getPloidy(),  //mixed ploidy will be a problem, but we don't do that in practice
                variant.getAlleles(), orderedRelevantAlleles, null, GenotypeAssignmentMethod.BEST_MATCH_TO_ORIGINAL,
                variant.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0), false, false);
        vcBuilder.genotypes(gc);
        return GATKVariantContextUtils.reverseTrimAlleles(vcBuilder.make());
    }

    private boolean isHighQuality(final VariantContext vc) {
        if (vc.getGenotypes().size() == 1) {
            //genotyping engine returns 0.01 for hom-ref SNPs, probably because of prior
            final Genotype genotype = vc.getGenotype(0);
            if (!genotype.hasPL()) {
                return false;
            }
            if (genotype.isHomRef()) {
                return false;
            }
            if (genotype.getPL()[0] == 0) {
                return false;
            }
            if (genotype.getAlleles().contains(Allele.NON_REF_ALLELE)) {
                return false;
            }
            return vc.getPhredScaledQual() > 0.01 && (MathUtils.arrayMax(genotype.getPL()) > 0);
        } else {
            for (final Genotype genotype : vc.getGenotypes()) {
                if (genotype.isHomRef()) {
                    continue;
                }
                if (passesGnomadAdjCriteria(genotype)) {
                    return true;
                }
            }
            return false;
        }
    }

    //https://gnomad.broadinstitute.org/news/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/
    private boolean passesGnomadAdjCriteria(final Genotype g) {
        if (!g.hasGQ() || g.getGQ() < 20) {
            return false;
        }
        if (!g.hasDP() || g.getDP() < 10) {
            return false;
        }
        if (g.isHet()) {
            if (!g.hasAD()) {
                return false;
            }
            final int[] ad = g.getAD();
            final double alleleBalance;
            if (g.isHetNonRef()) {
                alleleBalance = ad[2] / ((double) ad[1]);

            } else {
                alleleBalance = ad[1] / ((double) ad[0]);
            }
            if (alleleBalance < 0.2 || alleleBalance > 0.8) {
                return false;
            }
        }
        return true;
    }
}
