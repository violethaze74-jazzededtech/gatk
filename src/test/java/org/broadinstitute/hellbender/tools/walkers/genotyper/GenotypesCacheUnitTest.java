package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.apache.commons.math3.exception.MathArithmeticException;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public final class GenotypesCacheUnitTest extends GATKBaseTest {

    GenotypesCache calcs;

    @BeforeClass
    public void init(){
        calcs = new GenotypesCache();
    }

}
