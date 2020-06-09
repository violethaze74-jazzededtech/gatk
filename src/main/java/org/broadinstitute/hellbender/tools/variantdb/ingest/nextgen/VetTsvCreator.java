package org.broadinstitute.hellbender.tools.variantdb.ingest.nextgen;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.tools.variantdb.ingest.IngestConstants;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

public abstract class VetTsvCreator {

    static final Logger logger = LogManager.getLogger(ExomeVetTsvCreator.class);

    private SimpleXSVWriter vetWriter = null;
    private final String sampleId;


    public VetTsvCreator(String sampleName, String sampleId, Path sampleDirectoryPath) {
        this.sampleId = sampleId;
        // If the vet directory inside it doesn't exist yet -- create it
        final String vetDirectoryName = "vet";
        final Path vetDirectoryPath = sampleDirectoryPath.resolve(vetDirectoryName);
        final File vetDirectory = new File(vetDirectoryPath.toString());
        if (!vetDirectory.exists()) {
            vetDirectory.mkdir();
        }

        // if the pet & vet tsvs don't exist yet -- create them
        try {
            // Create a vet file to go into the pet dir for _this_ sample
            final String vetOutputName = sampleName + vetDirectoryName + IngestConstants.FILETYPE;
            final Path vetOutputPath = vetDirectoryPath.resolve(vetOutputName);
            // write header to it
            List<String> vetHeader = ExomeVetTsvCreator.getHeaders();
            vetWriter = new SimpleXSVWriter(vetOutputPath, IngestConstants.SEPARATOR);
            vetWriter.setHeaderLine(vetHeader);
        } catch (final IOException e) {
            throw new UserException("Could not create vet outputs", e);
        }

    }

    public abstract List<String> createRow(long start, VariantContext variant, String sampleId);

    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        int start = variant.getStart();
        int end = variant.getEnd();
        // check to see if this is an array
        // else, it must be an exome or genome!
        List<String> row = createRow(
                SchemaUtils.encodeLocation(variant.getContig(), start),
                variant,
                sampleId
        );

        // write the variant to the XSV
        SimpleXSVWriter.LineBuilder vetLine = vetWriter.getNewLineBuilder();
        vetLine.setRow(row);
        vetLine.write();

    }

    public void closeTool() {
        if (vetWriter != null) {
            try {
                vetWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close VET writer", e);
            }
        }

    }
}