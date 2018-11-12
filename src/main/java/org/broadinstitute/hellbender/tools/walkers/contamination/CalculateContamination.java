package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.lang.mutable.MutableDouble;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import org.broadinstitute.hellbender.tools.walkers.mutect.FilterMutectCalls;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>
 *     Calculates the fraction of reads coming from cross-sample contamination, given results from {@link GetPileupSummaries}.
 *     The resulting contamination table is used with {@link FilterMutectCalls}.
 * </p>
 * <p>
 *     This tool is featured in the Somatic Short Mutation calling Best Practice Workflow.
 *     See <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11136">Tutorial#11136</a> for a
 *     step-by-step description of the workflow and <a href="https://software.broadinstitute.org/gatk/documentation/article?id=11127">Article#11127</a>
 *     for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the
 *     <a href="https://github.com/broadinstitute/gatk/tree/master/scripts/mutect2_wdl">Mutect2 WDL scripts directory</a>.
 * </p>
 *
 * <p>This tool borrows from <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3167057/">ContEst</a> by Cibulskis et al the idea of estimating contamination
 * from ref reads at hom alt sites.  However, ContEst uses a probabilistic model that assumes a diploid genotype with no copy number
 * variation and independent contaminating reads.  That is, ContEst assumes that each contaminating read is drawn randomly and
 * independently from a different human.  This tool uses a simpler estimate of contamination that relaxes these assumptions.  In particular,
 * it works in the presence of copy number variations and with an arbitrary number of contaminating samples.  In addition, this tool
 * is designed to work well with no matched normal data.  However, one can run {@link GetPileupSummaries} on a matched normal bam file
 * and input the result to this tool.</p>
 *
 * <h3>Usage examples</h3>
 *
 * <h4>Tumor-only mode</h4>
 *
 * <pre>
 * gatk CalculateContamination \
 *   -I pileups.table \
 *   -O contamination.table
 * </pre>
 *
 * <h4>Matched normal mode</h4>
 *
 * <pre>
 * gatk CalculateContamination \
 *   -I tumor-pileups.table \
 *   -matched normal-pileups.table \
 *   -O contamination.table
 * </pre>
 * <p>
 *     The resulting table provides the fraction contamination, one line per sample, e.g. SampleID--TAB--Contamination.
 *     The file has no header.
 * </p>
 *
 */
@CommandLineProgramProperties(
        summary = "Calculate the fraction of reads coming from cross-sample contamination",
        oneLineSummary = "Calculate the fraction of reads coming from cross-sample contamination",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class CalculateContamination extends CommandLineProgram {

    public static final Logger logger = LogManager.getLogger(CalculateContamination.class);

    // our analysis only cares about hom alt and het sites, so we throw away hom refs with a very conservative heuristic
    private static final double ALT_FRACTION_OF_DEFINITE_HOM_REF = 0.05;

    private static final double STRICT_LOH_MAF_THRESHOLD = 0.4;

    private static final double INITIAL_CONTAMINATION_GUESS = 0.05;
    private static final int MAX_ITERATIONS = 10;
    private static final double CONTAMINATION_CONVERGENCE_THRESHOLD = 0.001;

    private static final int MIN_COVERAGE = 10;
    private static final double DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD = 1.0/2;
    private static final double DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD = 3.0;
    public static final int DESIRED_MINIMUM_HOM_ALT_COUNT = 50;
    public static final double MINOR_ALLELE_FRACTION_STEP_SIZE = 0.05;

    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc="The input table")
    private File inputPileupSummariesTable;

    public static final String MATCHED_NORMAL_LONG_NAME = "matched-normal";
    public static final String MATCHED_NORMAL_SHORT_NAME = "matched";
    @Argument(fullName = MATCHED_NORMAL_LONG_NAME,
            shortName = MATCHED_NORMAL_SHORT_NAME,
            doc="The matched normal input table", optional = true)
    private File matchedPileupSummariesTable = null;

    @Argument(fullName= StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output table")
    private final File outputTable = null;

    public static final String TUMOR_SEGMENTATION_LONG_NAME = "tumor-segmentation";
    public static final String TUMOR_SEGMENTATION_SHORT_NAME = "segments";
    @Argument(fullName= TUMOR_SEGMENTATION_LONG_NAME,
            shortName= TUMOR_SEGMENTATION_SHORT_NAME,
            doc="The output table containing segmentation of the tumor by minor allele fraction", optional = true)
    private final File outputTumorSegmentation = null;

    public static final String LOW_COVERAGE_RATIO_THRESHOLD_NAME = "low-coverage-ratio-threshold";
    @Argument(fullName = LOW_COVERAGE_RATIO_THRESHOLD_NAME,
            doc="The minimum coverage relative to the median.", optional = true)
    private final double lowCoverageRatioThreshold = DEFAULT_LOW_COVERAGE_RATIO_THRESHOLD;

    public static final String HIGH_COVERAGE_RATIO_THRESHOLD_NAME = "high-coverage-ratio-threshold";
    @Argument(fullName= HIGH_COVERAGE_RATIO_THRESHOLD_NAME,
            doc="The maximum coverage relative to the mean.", optional = true)
    private final double highCoverageRatioThreshold = DEFAULT_HIGH_COVERAGE_RATIO_THRESHOLD;

    @Override
    public Object doWork() {
        final List<PileupSummary> sites = filterSitesByCoverage(PileupSummary.readFromFile(inputPileupSummariesTable));

        // used the matched normal to genotype (i.e. find hom alt sites) if available
        final List<PileupSummary> genotypingSites = matchedPileupSummariesTable == null ? sites :
                filterSitesByCoverage(PileupSummary.readFromFile(matchedPileupSummariesTable));

        // partition genome into minor allele fraction (MAF) segments to better distinguish hom alts from LoH hets.
        final List<List<PileupSummary>> genotypingSegments = ContaminationSegmenter.findSegments(genotypingSites);

        // TODO: 1. Estimate contamination and best-fit contamination model from least LoH segments
        // TODO:    A good way to estimate LoH is to take the ratio of to expected hets = sum_sites 2f(1-f)
        // TODO: 2. Estimate MAF of each segment


        List<PileupSummary> homAltGenotypingSites = new ArrayList<>();
        final MutableDouble genotypingContamination = new MutableDouble(INITIAL_CONTAMINATION_GUESS);

        for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
            List<List<PileupSummary>> homAltSitesBySegment = Arrays.asList(new ArrayList<>());
            final MutableDouble minorAlleleFractionThreshold = new MutableDouble(STRICT_LOH_MAF_THRESHOLD);
            while (homAltSitesBySegment.stream().mapToInt(List::size).sum() < DESIRED_MINIMUM_HOM_ALT_COUNT && minorAlleleFractionThreshold.doubleValue() > 0) {
                homAltSitesBySegment = genotypingSegments.stream()
                        .map(segment -> ContaminationEngine.segmentHomAlts(segment, genotypingContamination.doubleValue(), minorAlleleFractionThreshold.doubleValue()))
                        .collect(Collectors.toList());
                minorAlleleFractionThreshold.subtract(MINOR_ALLELE_FRACTION_STEP_SIZE);
            }
            homAltGenotypingSites = homAltSitesBySegment.stream().flatMap(List::stream).collect(Collectors.toList());
            final double newGenotypingContamination = ContaminationEngine.calculateContamination(homAltGenotypingSites, errorRate(genotypingSites)).getLeft();
            if (Math.abs(newGenotypingContamination - genotypingContamination.doubleValue()) < CONTAMINATION_CONVERGENCE_THRESHOLD) {
                break;
            }
            genotypingContamination.setValue(newGenotypingContamination);
        }

        if (outputTumorSegmentation != null) {
            final List<List<PileupSummary>> tumorSegments = matchedPileupSummariesTable == null ?
                    genotypingSegments : ContaminationSegmenter.findSegments(sites);
            List<MinorAlleleFractionRecord> tumorMinorAlleleFractions = tumorSegments.stream()
                    .map(this::makeMinorAlleleFractionRecord).collect(Collectors.toList());
            MinorAlleleFractionRecord.writeToFile(tumorMinorAlleleFractions, outputTumorSegmentation);

        }

        final List<PileupSummary> homAltSites = subsetSites(sites, homAltGenotypingSites);
        final Pair<Double, Double> contaminationAndError = ContaminationEngine.calculateContamination(homAltSites, errorRate(sites));
        final double contamination = contaminationAndError.getLeft();
        final double error = contaminationAndError.getRight();
        ContaminationRecord.writeToFile(Arrays.asList(new ContaminationRecord(ContaminationRecord.Level.WHOLE_BAM.toString(), contamination, error)), outputTable);

        return "SUCCESS";
    }

    // in a biallelic site, essentially every non-ref, non-primary alt base is an error, since there are 2 such possible
    // errors out of 3 total, we multiply by 3/2 to get the total base error rate
    private double errorRate(List<PileupSummary> sites) {
        final long totalBases = sites.stream().mapToInt(PileupSummary::getTotalCount).sum();
        final long otherAltBases = sites.stream().mapToInt(PileupSummary::getOtherAltCount).sum();
        return 1.5 * ((double) otherAltBases / totalBases);
    }

    private static List<PileupSummary> subsetSites(final List<PileupSummary> sites, final List<PileupSummary> subsetLoci) {
        final OverlapDetector<PileupSummary> od = OverlapDetector.create(subsetLoci);
        return sites.stream().filter(od::overlapsAny).collect(Collectors.toList());
    }

    private MinorAlleleFractionRecord makeMinorAlleleFractionRecord(final List<PileupSummary> segment) {
        final String contig = segment.get(0).getContig();
        final int start = segment.get(0).getStart();
        final int end = segment.get(segment.size() - 1).getEnd();
        final double minorAlleleFraction = ContaminationEngine.calculateMinorAlleleFraction(segment);
        return new MinorAlleleFractionRecord(new SimpleInterval(contig, start, end), minorAlleleFraction);
    }


    private List<PileupSummary> filterSitesByCoverage(final List<PileupSummary> allSites) {
        // Just in case the intervals given to GetPileupSummaries contained un-covered sites, we remove them
        // so that a bunch of zeroes don't throw off the median coverage
        final List<PileupSummary> coveredSites = allSites.stream().filter(s -> s.getTotalCount() > MIN_COVERAGE).collect(Collectors.toList());
        final double[] coverage = coveredSites.stream().mapToDouble(PileupSummary::getTotalCount).toArray();
        final double medianCoverage = new Median().evaluate(coverage);
        final double meanCoverage = new Mean().evaluate(coverage);
        final double lowCoverageThreshold = medianCoverage * lowCoverageRatioThreshold;
        final double highCoverageThreshold = meanCoverage * highCoverageRatioThreshold;
        return coveredSites.stream()
                .filter(ps -> ps.getTotalCount() > lowCoverageThreshold && ps.getTotalCount() < highCoverageThreshold)
                .collect(Collectors.toList());
    }

}
