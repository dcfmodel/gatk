package org.broadinstitute.hellbender.tools.walkers.contamination;

import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.lang3.Range;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ContaminationSegmenter {
    public static final Range<Double> ALT_FRACTIONS_FOR_SEGMENTATION = Range.between(0.1, 0.9);
    public static final double KERNEL_SEGMENTER_LINEAR_COST = 1.0;
    public static final double KERNEL_SEGMENTER_LOG_LINEAR_COST = 1.0;
    public static final int KERNEL_SEGMENTER_DIMENSION = 100;
    public static final int POINTS_PER_SEGMENTATION_WINDOW = 50;
    public static final int MAX_CHANGEPOINTS_PER_CHROMOSOME = 10;
    public static final int MIN_SITES_PER_SEGMENT = 5;
    private static final double SEGMENTATION_KERNEL_VARIANCE = 0.025;

    static final BiFunction<PileupSummary, PileupSummary, Double> SEGMENTATION_KERNEL = (ps1, ps2) -> {
        final double maf1 = FastMath.min(ps1.getAltFraction(), 1 - ps1.getAltFraction());
        final double maf2 = FastMath.min(ps2.getAltFraction(), 1 - ps2.getAltFraction());
        return FastMath.exp(-MathUtils.square(maf1 - maf2)/(2 * SEGMENTATION_KERNEL_VARIANCE));
    };

    private ContaminationSegmenter() {}

    public static List<PileupSummary> getLikelyHetsBasedOnAlleleFraction(List<PileupSummary> sites) {
        return sites.stream()
                    .filter(ps -> ALT_FRACTIONS_FOR_SEGMENTATION.contains(ps.getAltFraction()))
                    .collect(Collectors.toList());
    }

    // we want log(1/2 (likelihood of alt minor + likelihood of alt major))
    //         =  logSumLog(log likelihood of alt minor, log likelihood of alt major) - log(2)
    public static double logLikelihoodOfHetsInSegment(final List<PileupSummary> hets, final double minorAlleleFraction) {
        return hets.stream().mapToDouble(het -> {
            final int n = het.getTotalCount();
            final int a = het.getAltCount();
            final double altMinorLogLikelihood = new BinomialDistribution(null, n, minorAlleleFraction).logProbability(a);
            final double altMajorLogLikelihood = new BinomialDistribution(null, n, 1 - minorAlleleFraction).logProbability(a);

            return MathUtils.logSumLog(altMinorLogLikelihood, altMajorLogLikelihood) + MathUtils.LOG_ONE_HALF;
        }).sum();
    }

    public static double calculateMinorAlleleFraction(final List<PileupSummary> segment) {
        final List<PileupSummary> hets = getLikelyHetsBasedOnAlleleFraction(segment);
        final Function<Double, Double> objective = maf -> logLikelihoodOfHetsInSegment(hets, maf);
        return OptimizationUtils.argmax(objective, ALT_FRACTIONS_FOR_SEGMENTATION.getMinimum(), 0.5, 0.4, 0.01, 0.01, 20);
    }

    public static List<PileupSummary> segmentHomAlts(final List<PileupSummary> segment, final double contamination, double minimiumMinorAlleleFraction) {
        final double minorAlleleFraction = calculateMinorAlleleFraction(segment);
        return minorAlleleFraction < minimiumMinorAlleleFraction ? Collections.emptyList() :
                segment.stream().filter(site -> ContaminationEngine.homAltProbability(site, minorAlleleFraction, contamination) > 0.5).collect(Collectors.toList());
    }

    public static List<List<PileupSummary>> findSegments(final List<PileupSummary> sites) {
        final Map<String, List<PileupSummary>> sitesByContig = sites.stream().collect(Collectors.groupingBy(PileupSummary::getContig));

        return sitesByContig.values().stream()
                .flatMap(contig -> findContigSegments(contig).stream())
                .filter(segment -> segment.size() >= MIN_SITES_PER_SEGMENT)
                .collect(Collectors.toList());
    }

    public static List<List<PileupSummary>> findContigSegments(List<PileupSummary> sites) {
        // segment based on obvious hets
        final List<PileupSummary> hetSites = getLikelyHetsBasedOnAlleleFraction(sites);

        if (hetSites.isEmpty()) {
            return Collections.emptyList();
        }

        final List<Integer> changepoints = new ArrayList<>();
        // when the kernel segmenter finds a changepoint at index n, that means index n belongs to the left segment, which goes
        // against the usual end-exclusive intervals of IndexRange etc.  This explains adding in the first changepoint of -1
        // instead of 0 and all the "changepoint + 1" constructions below
        changepoints.add(-1);
        final KernelSegmenter<PileupSummary> segmenter = new KernelSegmenter<>(hetSites);
        changepoints.addAll(segmenter.findChangepoints(MAX_CHANGEPOINTS_PER_CHROMOSOME, SEGMENTATION_KERNEL, KERNEL_SEGMENTER_DIMENSION,
                Arrays.asList(POINTS_PER_SEGMENTATION_WINDOW), KERNEL_SEGMENTER_LINEAR_COST, KERNEL_SEGMENTER_LOG_LINEAR_COST, KernelSegmenter.ChangepointSortOrder.INDEX));
        changepoints.add(hetSites.size()-1);

        final List<SimpleInterval> segments = IntStream.range(0, changepoints.size() - 1)
                .mapToObj(n -> {
                    final PileupSummary firstSiteInSegment = hetSites.get(changepoints.get(n) + 1);
                    final PileupSummary lastSiteInSegment = hetSites.get(changepoints.get(n+1));
                    return new SimpleInterval(firstSiteInSegment.getContig(), firstSiteInSegment.getStart(), lastSiteInSegment.getEnd());
                }).collect(Collectors.toList());

        final OverlapDetector<PileupSummary> od = OverlapDetector.create(sites);

        // for each segment, find overlapping sites and sort by coordinate
        return segments.stream()
                .map(segment -> od.getOverlaps(segment).stream().sorted(Comparator.comparingInt(PileupSummary::getStart)).collect(Collectors.toList()))
                .collect(Collectors.toList());
    }
}
