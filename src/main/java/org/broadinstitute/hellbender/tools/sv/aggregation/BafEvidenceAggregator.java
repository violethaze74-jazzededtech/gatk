package org.broadinstitute.hellbender.tools.sv.aggregation;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.tools.sv.BafEvidence;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;
import org.broadinstitute.hellbender.utils.SimpleInterval;

public class BafEvidenceAggregator extends SVEvidenceAggregator<BafEvidence> {

    private final int minSize;
    private final double paddingFraction;

    public BafEvidenceAggregator(final FeatureDataSource<BafEvidence> source,
                                 final SAMSequenceDictionary dictionary,
                                 final int minSize,
                                 final double paddingFraction) {
        super(source, dictionary);
        this.minSize = minSize;
        this.paddingFraction = paddingFraction;
    }

    @Override
    public SimpleInterval getEvidenceQueryInterval(final SVCallRecord call) {
        if (call.isSimpleCNV() && call.getLength() >= minSize) {
            final int padding = (int) Math.ceil(call.getLength() * paddingFraction);
            return new SimpleInterval(call.getContigA(), call.getPositionA(), call.getPositionB())
                    .expandWithinContig(padding, dictionary);
        } else {
            return null;
        }
    }

    @Override
    public boolean evidenceFilter(final SVCallRecord record, final BafEvidence evidence) {
        return true;
    }
}
