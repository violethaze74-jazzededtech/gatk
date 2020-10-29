package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;


public class GenomeVetTsvCreator extends VetTsvCreator {
    static final Logger logger = LogManager.getLogger(GenomeVetTsvCreator.class);

    public GenomeVetTsvCreator(String sampleName, String sampleId, String tableNumberPrefix, final File outputDirectory) {
        super(sampleName, sampleId, tableNumberPrefix, outputDirectory, getHeaders());
    }

    @Override
    public List<String> createRow(final long location, final VariantContext variant, final String sampleId) {
        List<String> row = new ArrayList<>();
        for ( final GenomeFieldEnum fieldEnum : GenomeFieldEnum.values() ) {
            if (fieldEnum.equals(GenomeFieldEnum.location)) {
                row.add(String.valueOf(location));
            } else if (fieldEnum.equals(GenomeFieldEnum.sample)) {
                row.add(sampleId);
            } else {
                row.add(fieldEnum.getColumnValue(variant));
            }
        }
        return row;
    }

    public static List<String> getHeaders() {
        return Arrays.stream(GenomeFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }
}
