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
    @Argument(fullName = "warn-on-errors",
            shortName = "warn-on-errors",
            doc = "just emit warnings on errors instead of terminating the run at the first instance",
            optional = true)
    Boolean WARN_ON_ERROR = false;

    @Argument(fullName = "ignore-quals", optional = true)
    Boolean IGNORE_QUALS = false;

    @Argument(fullName = "dp-change-allowed", doc = "Note that this is a signed change (actual - expected)")
    Integer DP_CHANGE = 0;

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
                final List<Allele> gt1 = trimmed1.getGenotype(0).getAlleles();

                final VariantContext trimmed2 = trimAlleles(match);
                final List<Allele> gt2 = trimmed2.getGenotype(0).getAlleles();

                if (!gt1.get(0).equals(gt2.get(0)) || !gt1.get(1).equals(gt2.get(1))) {  //do use order here so we can check phasing
                    //this is okay if a star got corrected by dropping an upstream hom ref del
                    final List<VariantContext> overlappingDels = variantContexts.stream().filter(v -> v.getStart() < vc.getStart() && v.overlaps(vc))
                            .collect(Collectors.toList());
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
                    }
                    else {
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
                final VariantContext actual_trimmed, expected_trimmed;
                if (trimmed1.getSource().equals("actual")) {
                    actual_trimmed = trimmed1;
                    expected_trimmed = trimmed2;
                } else {
                    actual_trimmed = trimmed2;
                    expected_trimmed = trimmed1;
                }
                try {
                    if (!areVariantContextsEqualOrderIndependent(actual_trimmed, expected_trimmed)) {
                        throwOrWarn(new UserException("Variant contexts do not match: " + actual_trimmed.toStringDecodeGenotypes() + " versus " + expected_trimmed.toStringDecodeGenotypes()));
                    }
                } catch (UserException e) {
                    throwOrWarn(new UserException("Variant contexts do not match at " + actual_trimmed.getContig() + ":" + actual_trimmed.getStart(), e));
                }

            }
        }
    }

    private boolean areVariantContextsEqualOrderIndependent(final VariantContext actual, final VariantContext expected) {
        if (!actual.getContig().equals(expected.getContig())) {
            return false;
        }
        if (actual.getStart() != expected.getStart()) {
            return false;
        }
        if (actual.getEnd() != expected.getEnd()) {
            return false;
        }
        if (!actual.getID().equals(expected.getID())) {
            return false;
        }
        if (!actual.getAlleles().equals(expected.getAlleles())) {
            return false;
        }

        if (!areAttributesEqual(actual.getAttributes(), expected.getAttributes())) {
            return false;
        }

        if (actual.filtersWereApplied() != expected.filtersWereApplied()) {
            return false;
        }
        if (actual.isFiltered() != expected.isFiltered()) {
            return false;
        }
        if (!actual.getFilters().equals(expected.getFilters())){
            return false;
        }
        if (Math.abs(actual.getPhredScaledQual() - expected.getPhredScaledQual()) > 0.001) {
            if (!IGNORE_QUALS) {
                throw new UserException("Qual scores differ by more than 0.001");
            }
        }

        if (!areGenotypesEqual(actual.getGenotype(0), expected.getGenotype(0))) {
            throw new UserException("Genotypes are not equal: " + actual.getGenotype(0) + " versus " + expected.getGenotype(0));
        }
        return true;
    }

    @SuppressWarnings("unchecked")
    private boolean areAttributesEqual(final Map<String, Object> actual, final Map<String, Object> expected) {
        final Set<String> expectedKeys = new LinkedHashSet<>(expected.keySet());

        for ( final Map.Entry<String, Object> exp : expected.entrySet() ) {
           final Object expectedValue = exp.getValue();
            final String key = exp.getKey();
            if ( actual.containsKey(key) && actual.get(key) != null ) {
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
                            return false;
                        }
                    }
                } else {
                    if (!isAttributeValueEqual(key, actualValue, expectedValue)) {
                        return false;
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
               throw new UserException("Variant contexts have different attribute values for + " + key + ": actual has " + actual.toString()
                       + "\n expected has " + expected.toString());
           } else {
               return actual.toString().equals(expected.toString());
           }
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
        if (actual.isPhased() != expected.isPhased()) {
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
