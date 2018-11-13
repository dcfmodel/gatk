package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.OptimizationUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.DoubleUnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ContaminationEngine {

    private static final Logger logger = LogManager.getLogger(ContaminationEngine.class);
    private static final int NUM_ITERATIONS = 3;
    private static final List<Double> CONTAMINATIONS_FOR_COMPARISON =
            Arrays.asList(0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2);

    private ContaminationEngine() { }

    private static double calculateMinorAlleleFraction(final ContaminationModel model, final List<PileupSummary> segment) {
        final DoubleUnaryOperator objective = maf -> model.logLikelihood(segment, maf);
        return OptimizationUtils.argmax(objective, 0.1, 0.5, 0.4, 0.01, 0.01, 20);
    }

    private static ContaminationModel chooseBestModel(final List<ContaminationModel> models, final List<List<PileupSummary>> segments, final List<Double> mafs) {
        final double[] logLikelihoods = models.stream().mapToDouble(model -> model.logLikelihood(segments, mafs)).toArray();
        return models.get(MathUtils.maxElementIndex(logLikelihoods));
    }

    public static Pair<Double, Double> calculateContaminationFromHomAlts(List<PileupSummary> homAltSites, final double errorRate) {
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

    // in a biallelic site, essentially every non-ref, non-primary alt base is an error, since there are 2 such possible
    // errors out of 3 total, we multiply by 3/2 to get the total base error rate
    private static double errorRate(List<PileupSummary> sites) {
        final long totalBases = sites.stream().mapToInt(PileupSummary::getTotalCount).sum();
        final long otherAltBases = sites.stream().mapToInt(PileupSummary::getOtherAltCount).sum();
        return 1.5 * ((double) otherAltBases / totalBases);
    }

    public static Info getContaminationInfo(List<PileupSummary> genotypingSites) {
        final double genotypingErrorRate = errorRate(genotypingSites);

        // partition genome into minor allele fraction (MAF) segments to better distinguish hom alts from LoH hets.
        final List<List<PileupSummary>> genotypingSegments = ContaminationSegmenter.findSegments(genotypingSites);

        List<ContaminationModel> genotypingModelsToCompare = makeModelsToCompare(genotypingErrorRate);
        List<Double> genotypingMAFs = genotypingSegments.stream().map(segment -> 0.5).collect(Collectors.toList());
        ContaminationModel genotypingModel = ContaminationModel.createNoContaminantModel(genotypingErrorRate);

        for (int n = 0; n < NUM_ITERATIONS; n++) {
            final ContaminationModel genotypingModelForLambda = genotypingModel;
            genotypingMAFs = genotypingSegments.stream()
                    .map(segment -> calculateMinorAlleleFraction(genotypingModelForLambda, segment))
                    .collect(Collectors.toList());
            genotypingModel = chooseBestModel(genotypingModelsToCompare, genotypingSegments, genotypingMAFs);
            // TODO: logging
        }

        return new Info(genotypingErrorRate, genotypingSegments, genotypingModel, genotypingMAFs);
    }

    private static List<ContaminationModel> makeModelsToCompare(final double eps) {
        final List<ContaminationModel> result = new ArrayList<>();
        for (final double c : CONTAMINATIONS_FOR_COMPARISON) {
            result.add(ContaminationModel.createInfiniteContaminantModel(eps, c));
            result.add(ContaminationModel.createSingleContaminantModel(eps, c));
            result.add(ContaminationModel.createTwoContaminantModel(eps, 0.5 * c, 0.5 * c));
            result.add(ContaminationModel.createTwoContaminantModel(eps, 0.25 * c, 0.75 * c));
        }
        return result;
    }

    public static class Info {
        private final double errorRate;
        private final List<List<PileupSummary>> segments;
        private final ContaminationModel contaminationModel;
        private final List<Double> minorAlleleFractions;

        public Info(double errorRate, List<List<PileupSummary>> segments, ContaminationModel contaminationModel, List<Double> minorAlleleFractions) {
            this.errorRate = errorRate;
            this.segments = segments;
            this.contaminationModel = contaminationModel;
            this.minorAlleleFractions = minorAlleleFractions;
        }

        public double getErrorRate() { return errorRate; }

        public List<List<PileupSummary>> getSegments() { return segments; }

        public ContaminationModel getContaminationModel() { return contaminationModel; }

        private List<PileupSummary> getType(final ContaminationModel.SampleGenotype type) {
            return IntStream.range(0, segments.size())
                    .mapToObj(n -> segments.get(n).stream().filter(site -> contaminationModel.probability(site, minorAlleleFractions.get(n), type) > 0.95))
                    .flatMap(stream -> stream)
                    .collect(Collectors.toList());
        }

        public List<PileupSummary> homAlts() { return getType(ContaminationModel.SampleGenotype.HOM_ALT); }
        public List<PileupSummary> homRefs() { return getType(ContaminationModel.SampleGenotype.HOM_REF); }

        public List<MinorAlleleFractionRecord> segmentationRecords() {
            return IntStream.range(0, segments.size()).mapToObj(n -> {
                        final List<PileupSummary> segment = segments.get(n);
                        final String contig = segment.get(0).getContig();
                        final int start = segment.get(0).getStart();
                        final int end = segment.get(segment.size() - 1).getEnd();
                        final double maf = minorAlleleFractions.get(n);
                        return new MinorAlleleFractionRecord(new SimpleInterval(contig, start, end), maf);
                    }).collect(Collectors.toList());
        }
    }
}
