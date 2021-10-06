package org.broadinstitute.hellbender.tools.walkers.qc;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.ArrayList;
import java.util.List;

public class ReadPair {
    private GATKRead firstOfPair;
    private GATKRead secondOfPair;
    private List<GATKRead> secondaryAlignments = new ArrayList<>(10);
    private List<GATKRead> supplementaryAlignments = new ArrayList<>(10); // Finally understand the difference

    ReadPair() { }

    public GATKRead getFirstOfPair(){
        return firstOfPair;
    }

    public GATKRead getSecondOfPair(){
        return secondOfPair;
    }

    public void add(final GATKRead read) {
        if (read.isFirstOfPair()) {
            this.firstOfPair = read;
        } else if (read.isSecondOfPair()) {
            this.secondOfPair = read;
        } else if (read.isSupplementaryAlignment()) {
            this.secondaryAlignments.add(read);
        } else {
            int d = 3;
            throw new UserException("Unknown read type");
        }
    }

    public boolean isDuplicateMarked() {
        // Doing some investigation
        if (firstOfPair.isDuplicate()) {
            // Make sure the rest is duplicate-marked
            if (!secondOfPair.isDuplicate() || secondaryAlignments.stream().anyMatch(r -> !r.isDuplicate())) {
                throw new UserException("First of pair a duplicate but the rest is not" + secondOfPair.getName());
            }
        } else {
            // Make sure the rest is not duplicate-marked
            if (secondOfPair.isDuplicate() || secondaryAlignments.stream().anyMatch(r -> r.isDuplicate())) {
                throw new UserException("First of pair a not duplicate but the rest is " + secondOfPair.getName());
            }
        }
        return firstOfPair.isDuplicate();
    }
}
