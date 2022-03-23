package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.testutils.ArgumentsBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.function.Supplier;

/**
 * TODO
 */
public final class ExtractVariantAnnotationsIntegrationTest extends CommandLineProgramTest {

    private static final List<String> NON_ALLELE_SPECIFIC_ANNOTATIONS = Arrays.asList(
            "DP", "FS", "MQ", "MQRankSum", "QD", "ReadPosRankSum", "SOR");

    private static final List<String> ALLELE_SPECIFIC_ANNOTATIONS = Arrays.asList(
            "DP", "AS_FS", "AS_MQ", "AS_MQRankSum", "AS_QD", "AS_ReadPosRankSum", "AS_SOR");

    private static final File TEST_FILES_DIR = new File(largeFileTestDir,
            "org/broadinstitute/hellbender/tools/walkers/vqsr/scalable/");
    private static final File INPUT_VCF = new File(TEST_FILES_DIR, "stroke_vqsr_magic_as.chr1.1-10M.vcf.gz");
    private static final File SNP_RESOURCE_VCF = new File(TEST_FILES_DIR, "1000G_omni2.5.hg38.chr1.1-10M.vcf.gz");
    private static final File INDEL_RESOURCE_VCF = new File(TEST_FILES_DIR, "Mills_and_1000G_gold_standard.indels.hg38.chr1.1-10M.vcf.gz");

    @DataProvider(name = "dataValidInputs")
    public Object[][] dataValidInputs() {
        final Function<List<String>, ArgumentsBuilder> argumentsBuilderFromAnnotations = (annotations) -> {
            final ArgumentsBuilder argsBuilder = new ArgumentsBuilder();
            argsBuilder.add(StandardArgumentDefinitions.VARIANT_LONG_NAME, INPUT_VCF);
            annotations.forEach(a -> argsBuilder.add(StandardArgumentDefinitions.ANNOTATION_LONG_NAME, a));
            return argsBuilder;
        };

        return new Object[][]{
                {
                        argumentsBuilderFromAnnotations.apply(NON_ALLELE_SPECIFIC_ANNOTATIONS)
                                .add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, "SNP")
                                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":omni,training=true,truth=true", SNP_RESOURCE_VCF)
                },
                {
                        argumentsBuilderFromAnnotations.apply(NON_ALLELE_SPECIFIC_ANNOTATIONS)
                                .add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, "INDEL")
                                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":mills,training=true,truth=true", INDEL_RESOURCE_VCF)
                },
                {
                        argumentsBuilderFromAnnotations.apply(NON_ALLELE_SPECIFIC_ANNOTATIONS)
                                .add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, "SNP")
                                .add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, "INDEL")
                                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":omni,training=true,truth=true", SNP_RESOURCE_VCF)
                                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":mills,training=true,truth=true", INDEL_RESOURCE_VCF)
                },
                {
                        argumentsBuilderFromAnnotations.apply(ALLELE_SPECIFIC_ANNOTATIONS)
                                .addFlag(LabeledVariantAnnotationsWalker.USE_ALLELE_SPECIFIC_ANNOTATIONS_LONG_NAME)
                                .add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, "SNP")
                                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":omni,training=true,truth=true", SNP_RESOURCE_VCF)
                },
                {
                        argumentsBuilderFromAnnotations.apply(ALLELE_SPECIFIC_ANNOTATIONS)
                                .addFlag(LabeledVariantAnnotationsWalker.USE_ALLELE_SPECIFIC_ANNOTATIONS_LONG_NAME)
                                .add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, "INDEL")
                                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":mills,training=true,truth=true", INDEL_RESOURCE_VCF)
                },
                {
                        argumentsBuilderFromAnnotations.apply(ALLELE_SPECIFIC_ANNOTATIONS)
                                .addFlag(LabeledVariantAnnotationsWalker.USE_ALLELE_SPECIFIC_ANNOTATIONS_LONG_NAME)
                                .add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, "SNP")
                                .add(LabeledVariantAnnotationsWalker.MODE_LONG_NAME, "INDEL")
                                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":omni,training=true,truth=true", SNP_RESOURCE_VCF)
                                .add(StandardArgumentDefinitions.RESOURCE_LONG_NAME + ":mills,training=true,truth=true", INDEL_RESOURCE_VCF)
                }
        };
    }

    @Test(dataProvider = "dataValidInputs")
    public void testValidInputs(final ArgumentsBuilder argsBuilder) {
        final File outputDir = createTempDir("testDir");
        final String outputPrefix = outputDir + "/test";
        argsBuilder.add(StandardArgumentDefinitions.OUTPUT_LONG_NAME, outputPrefix);
        argsBuilder.add(StandardArgumentDefinitions.VERBOSITY_NAME, "INFO");
        runCommandLine(argsBuilder);
    }

    @Test
    public void test1kgp50ExomesAll() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.all",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesAllUnlabeled() throws IOException {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled",
                "--maximum-number-of-unlabeled-variants", "10000",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);

        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled.annot.hdf5 /home/slee/working/vqsr/scalable/extract-test/expected/test.all-unlabeled.annot.hdf5");
        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled.unlabeled.annot.hdf5 /home/slee/working/vqsr/scalable/extract-test/expected/test.all-unlabeled.unlabeled.annot.hdf5");
        runSystemCommand("diff /home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled.unlabeled.vcf.gz /home/slee/working/vqsr/scalable/extract-test/expected/test.all-unlabeled.unlabeled.vcf.gz");
        runSystemCommand("diff /home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled.unlabeled.vcf.gz.tbi /home/slee/working/vqsr/scalable/extract-test/expected/test.all-unlabeled.unlabeled.vcf.gz.tbi");
    }

    private static void runSystemCommand(final String command) {
        try {
            Process process = Runtime.getRuntime().exec(command);

            BufferedReader reader = new BufferedReader(
                    new InputStreamReader(process.getInputStream()));
            while (reader.readLine() != null) {
                Assert.fail(command);
            }

            reader.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void test1kgp50ExomesSNP() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.snp",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testSNPAS() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", largeFileTestDir + "VQSR/chr1snippet.doctoredMQ.doctoredAS.sites_only.vcf",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.snp.as",
                "--use-allele-specific-annotations",
                "-A", "AS_FS",
                "-A", "AS_ReadPosRankSum",
                "-A", "AS_MQRankSum",
                "-A", "AS_QD",
                "-A", "AS_SOR",
                "-A", "AS_MQ",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--resource:doctored,training=true,truth=true", largeFileTestDir + "VQSR/chr1snippet.doctoredMQ.doctoredAS.sites_only.vcf",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testSNPNonAS() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", largeFileTestDir + "VQSR/chr1snippet.doctoredMQ.doctoredAS.sites_only.vcf",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.snp.non-as",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--resource:doctored,training=true,truth=true", largeFileTestDir + "VQSR/chr1snippet.doctoredMQ.doctoredAS.sites_only.vcf",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void test1kgp50ExomesIndel() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/extract-test/test.indel",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "--trust-all-polymorphic",
                "--mode", "INDEL",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

//    @Test
//    public void test1kgp50ExomesSNPExactMatch() {
//        final File outputDir = createTempDir("extract-test");
//        final String[] arguments = {
//                "-L", "chr1",
//                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
//                "-O", new File(outputDir, "test.snp").getAbsolutePath(),
//                "-A", "FS",
//                "-A", "ReadPosRankSum",
//                "-A", "MQRankSum",
//                "-A", "QD",
//                "-A", "SOR",
//                "-A", "MQ",
//                "--trust-all-polymorphic",
//                "-mode", "SNP",
//                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
//                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
//                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//        assertAnnotationFilesEqual(new File(outputDir, "test.snp.annot.hdf5"), new File("/home/slee/working/vqsr/scalable/extract-exact-match", "test.snp.annot.hdf5"));
//    }
//
//    @Test
//    public void testSNPASExactMatch() {
//        final File outputDir = createTempDir("extract-test");
//        final String[] arguments = {
//                "-L", "chr1",
//                "-V", largeFileTestDir + "VQSR/chr1snippet.doctoredMQ.sites_only.vcf.gz",
//                "-O", new File(outputDir, "test.snp.as").getAbsolutePath(),
//                "--use-allele-specific-annotations",
//                "-A", "AS_FS",
//                "-A", "AS_ReadPosRankSum",
//                "-A", "AS_MQRankSum",
//                "-A", "AS_QD",
//                "-A", "AS_SOR",
//                "-A", "AS_MQ",
//                "--trust-all-polymorphic",
//                "--mode", "SNP",
//                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
//                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
//                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
//                "--verbosity", "INFO"
//        };
//        runCommandLine(arguments);
//        assertAnnotationFilesEqual(new File(outputDir, "test.snp.as.annot.hdf5"), new File("/home/slee/working/vqsr/scalable/extract-exact-match", "test.snp.as.annot.hdf5"));
//    }
//
//    private static void assertAnnotationFilesEqual(final File annotationsFile1,
//                                                   final File annotationsFile2) {
//        try (final HDF5File annotationsHDF5File1 = new HDF5File(annotationsFile1, HDF5File.OpenMode.READ_ONLY);
//             final HDF5File annotationsHDF5File2 = new HDF5File(annotationsFile2, HDF5File.OpenMode.READ_ONLY)) {
//            Assert.assertEquals(
//                    annotationsHDF5File1.readStringArray("/annotations/names"),
//                    annotationsHDF5File2.readStringArray("/annotations/names"));
//            for (final String label : Arrays.asList("training", "truth")) {
//                Assert.assertEquals(
//                        LabeledVariantAnnotationsData.readChunkedDoubleArray(annotationsHDF5File1, String.format("/labels/%s", label)),
//                        LabeledVariantAnnotationsData.readChunkedDoubleArray(annotationsHDF5File2, String.format("/labels/%s", label)));
//            }
//            Assert.assertEquals(
//                    HDF5Utils.readChunkedDoubleMatrix(annotationsHDF5File1, "/annotations"),
//                    HDF5Utils.readChunkedDoubleMatrix(annotationsHDF5File2, "/annotations"));
//        } catch (final HDF5LibException exception) {
//            Assert.fail("Exception encountered during reading of annotations:", exception);
//        }
//    }

    @Test
    public void testJbxAll() {
        final String[] arguments = {
                "-V", "/home/slee/working/vqsr/scalable/jbx/resources/Test50Callset.annoated_pids.sites-only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all.extract",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "-A", "COMBINED_TREE_SCORE",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxAllUnlabeled() {
        final String[] arguments = {
                "-V", "/home/slee/working/vqsr/scalable/jbx/resources/Test50Callset.annoated_pids.sites-only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all-unlabeled.extract",
                "--maximum-number-of-unlabeled-variants", "10000000",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "-A", "COMBINED_TREE_SCORE",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxSNP() {
        final String[] arguments = {
                "-V", "/home/slee/working/vqsr/scalable/jbx/resources/Test50Callset.annoated_pids.sites-only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.snp.extract",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "-A", "COMBINED_TREE_SCORE",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }
}