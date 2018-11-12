package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;

import java.util.Collections;
import java.util.List;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Collectors;

public class ContaminationEngine {

    public static final Logger logger = LogManager.getLogger(ContaminationEngine.class);

    private ContaminationEngine() { }

    public static double calculateMinorAlleleFraction(final ContaminationModel model, final List<PileupSummary> segment) {
        final DoubleUnaryOperator objective = maf -> segment.stream().mapToDouble(site -> model.logLikelihood(site, maf)).sum();
        return OptimizationUtils.argmax(objective, 0.1, 0.5, 0.4, 0.01, 0.01, 20);
    }

    public static List<PileupSummary> segmentHomAlts(final List<PileupSummary> segment, final double contamination, double minimiumMinorAlleleFraction) {
        final double minorAlleleFraction = calculateMinorAlleleFraction(segment);
        return minorAlleleFraction < minimiumMinorAlleleFraction ? Collections.emptyList() :
                segment.stream().filter(site -> homAltProbability(site, minorAlleleFraction, contamination) > 0.5).collect(Collectors.toList());
    }

    //private static final List<SampleGenotype> SAMPLE_GENOTYPES = Arrays.asList(

    public static double homAltProbability(final PileupSummary site, final double minorAlleleFraction, final double contamination) {
        final double alleleFrequency = site.getAlleleFrequency();
        final double homAltPrior = MathUtils.square(alleleFrequency);
        final double hetPrior = 2 * alleleFrequency * (1 - alleleFrequency);

        final int altCount = site.getAltCount();
        final int totalCount = altCount + site.getRefCount();

        if (altCount < totalCount / 2) {
            return 0;
        }

        final double homAltLikelihood = new BinomialDistribution(null, totalCount, 1 - contamination).probability(altCount);
        final double hetLikelihood = new BinomialDistribution(null, totalCount, 1 - minorAlleleFraction).probability(altCount);

        final double unnormalizedHomAltProbability = homAltPrior * homAltLikelihood;
        final double unnormalizedHetProbability = hetPrior * hetLikelihood;

        final double result = unnormalizedHomAltProbability / (unnormalizedHetProbability + unnormalizedHomAltProbability);

        return result;

    }

    public static Pair<Double, Double> calculateContamination(List<PileupSummary> homAltSites, final double errorRate) {
        if (homAltSites.isEmpty()) {
            logger.warn("No hom alt sites found!  Perhaps GetPileupSummaries was run on too small of an interval, or perhaps the sample was extremely inbred or haploid.");
            return Pair.of(0.0, 1.0);
        }

        final long totalReadCount = homAltSites.stream().mapToLong(PileupSummary::getTotalCount).sum();
        final long totalRefCount = homAltSites.stream().mapToLong(PileupSummary::getRefCount).sum();

        // if eg ref is A, alt is C, then # of ref reads due to error is roughly (# of G read + # of T reads)/2
        final long errorRefCount = Math.round(totalReadCount * errorRate / 3);
        final long contaminationRefCount = Math.max(totalRefCount - errorRefCount, 0);
        final double totalDepthWeightedByRefFrequency = homAltSites.stream()
                .mapToDouble(ps -> ps.getTotalCount() * (1 - ps.getAlleleFrequency()))
                .sum();
        final double contamination = contaminationRefCount / totalDepthWeightedByRefFrequency;
        final double standardError = Math.sqrt(contamination / totalDepthWeightedByRefFrequency);

        logger.info(String.format("In %d homozygous variant sites we find %d reference reads due to contamination and %d" +
                    " due to to sequencing error out of a total %d reads.", homAltSites.size(), contaminationRefCount, errorRefCount, totalReadCount));
        logger.info(String.format("Based on population data, we would expect %d reference reads in a contaminant with equal depths at these sites.", (long) totalDepthWeightedByRefFrequency));
        logger.info(String.format("Therefore, we estimate a contamination of %.3f.", contamination));
        logger.info(String.format("The error bars on this estimate are %.5f.", standardError));

        return Pair.of(Math.min(contamination, 1.0), standardError);
    }

}
