package org.broadinstitute.hellbender.tools.walkers.vqsr.scalable;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * See documentation for {@link ExtractVariantAnnotationsIntegrationTest} for information about how inputs and
 * expected outputs used there are related to those used here and in {@link TrainVariantAnnotationsModelIntegrationTest}.
 */
public final class ScoreVariantAnnotationsIntegrationTest extends CommandLineProgramTest {

    private static final String PYTHON_SCRIPT = packageMainResourcesDir + "tools/walkers/vqsr/scalable/isolation-forest.py";

    @Test
    public void test1kgp50ExomesAll() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/score-test/test.all",
                "--python-script", PYTHON_SCRIPT,
                "--model-prefix", "/home/slee/working/vqsr/scalable/train-test/test.all",
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
    public void test1kgp50ExomesAllUnlabeled() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/score-test/test.all-unlabeled",
                "--python-script", PYTHON_SCRIPT,
                "--model-prefix", "/home/slee/working/vqsr/scalable/train-test/test.all-unlabeled",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--mode", "INDEL",
                "--resource:extracted,extracted=true", "/home/slee/working/vqsr/scalable/extract-test/test.all-unlabeled.vcf.gz",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);

        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/score-test/test.all-unlabeled.annot.hdf5 /home/slee/working/vqsr/scalable/score-test/expected/test.all-unlabeled.annot.hdf5");
        runSystemCommand("h5diff /home/slee/working/vqsr/scalable/score-test/test.all-unlabeled.scores.hdf5 /home/slee/working/vqsr/scalable/score-test/expected/test.all-unlabeled.scores.hdf5");
        runSystemCommand("bcftools view -H /home/slee/working/vqsr/scalable/score-test/test.all-unlabeled.vcf.gz > /tmp/1.vcf; bcftools view -H /home/slee/working/vqsr/scalable/score-test/expected/test.all-unlabeled.vcf.gz > /tmp/2.vcf; diff /tmp/1.vcf /tmp/2.vcf; rm /tmp/1.vcf /tmp/2.vcf");
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
                "-O", "/home/slee/working/vqsr/scalable/score-test/test.snp",
                "--python-script", PYTHON_SCRIPT,
                "--model-prefix", "/home/slee/working/vqsr/scalable/train-test/test",
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
    public void test1kgp50ExomesIndel() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/score-test/test.indel",
                "--python-script", PYTHON_SCRIPT,
                "--model-prefix", "/home/slee/working/vqsr/scalable/train-test/test",
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

    @Test
    public void test1kgp50ExomesBGMMSNP() {
        final String[] arguments = {
                "-L", "chr1",
                "-V", "/home/slee/working/vqsr/1kgp-50-exomes/resources/1kgp-50-exomes.sites_only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/score-test/test.bgmm.snp",
                "--model-prefix", "/home/slee/working/vqsr/scalable/train-test/test.bgmm",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "--trust-all-polymorphic",
                "-mode", "SNP",
                "--resource:extracted-training,training=true,truth=false", "/home/slee/working/vqsr/scalable/extract-exact-match/test.snp.vcf",
                "--resource:hapmap,training=false,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=false,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxAll() {
        final String[] arguments = {
                "-L", "chr1",
                "-L", "chr2",
                "-V", "/home/slee/working/vqsr/scalable/jbx/resources/Test50Callset.annoated_pids.sites-only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all.score",
                "--python-script", PYTHON_SCRIPT,
                "--model-prefix", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all.train",
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
                "--resource:extracted,extracted=true", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all.extract.vcf",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxAllUnlabeledGC() {
        final String[] arguments = {
                "--omit-alleles-in-hdf5",
                "-L", "chr1",
                "-L", "chr2",
                "-V", "/home/slee/working/vqsr/scalable/jbx/resources/Test50Callset.annoated_pids.sites-only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/jbx/GC/Test50Callset.all-unlabeled.score",
                "--python-script", PYTHON_SCRIPT,
                "--model-prefix", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all-unlabeled.train",
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
                "--resource:extracted,extracted=true", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all-unlabeled.extract.vcf.gz",
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
        for (int i = 1; i <= 22; i++) {
            final String[] arguments = {
                    "--omit-alleles-in-hdf5",
                    "-L", "chr" + i,
                    "-V", "/home/slee/working/vqsr/scalable/jbx/resources/Test50Callset.annoated_pids.sites-only.vcf.gz",
                    "-O", String.format("/home/slee/working/vqsr/scalable/jbx/per-chr/Test50Callset.all-unlabeled.score.chr%02d", i),
                    "--python-script", PYTHON_SCRIPT,
                    "--model-prefix", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all-unlabeled.train",
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
                    "--resource:extracted,extracted=true", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all-unlabeled.extract.vcf.gz",
                    "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                    "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                    "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                    "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                    "--verbosity", "INFO"
            };
            runCommandLine(arguments);
        }
    }

    @Test
    public void testJbxSNP() {
        final String[] arguments = {
                "-L", "chr1:1-100000000",
                "-V", "/home/slee/working/vqsr/scalable/jbx/resources/Test50Callset.annoated_pids.sites-only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.snp.test",
                "--python-script", PYTHON_SCRIPT,
                "--model-prefix", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.snp.train",
                "-A", "FS",
                "-A", "ReadPosRankSum",
                "-A", "MQRankSum",
                "-A", "QD",
                "-A", "SOR",
                "-A", "MQ",
                "-A", "COMBINED_TREE_SCORE",
                "--trust-all-polymorphic",
                "--mode", "SNP",
                "--resource:extracted,extracted=true", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.snp.extract.vcf.gz",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--verbosity", "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testJbxBGMMAll() {
        final String[] arguments = {
                "-L", "chr1",
                "-L", "chr2",
                "-V", "/home/slee/working/vqsr/scalable/jbx/resources/Test50Callset.annoated_pids.sites-only.vcf.gz",
                "-O", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.bgmm.all.score",
                "--model-prefix", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.bgmm.all.train",
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
                "--resource:extracted,extracted=true", "/home/slee/working/vqsr/scalable/jbx/Test50Callset.all.extract.vcf",
                "--resource:hapmap,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/hapmap_3.3.hg38.vcf.gz",
                "--resource:omni,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/1000G_omni2.5.hg38.vcf.gz",
                "--resource:1000G,training=true,truth=false", "/mnt/4AB658D7B658C4DB/working/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
                "--resource:mills,training=true,truth=true", "/mnt/4AB658D7B658C4DB/working/ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                "--verbosity", "INFO"
        };
        runCommandLine(arguments);
    }
}