package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Arrays;
import java.util.List;
import java.util.function.DoubleUnaryOperator;
import java.util.function.ToDoubleFunction;

public class ContaminationEngine {

    public static final Logger logger = LogManager.getLogger(ContaminationEngine.class);

    private ContaminationEngine() { }

    private enum SampleGenotype {
        HOM_REF(maf -> 0,f -> (1 - f) * (1 - f)),
        ALT_MINOR(maf -> maf, f -> f * (1 - f)),
        ALT_MAJOR(maf -> 1 - maf, f -> f * (1 - f)),
        HOM_ALT(maf -> 1, f -> f * f);

        final DoubleUnaryOperator alleleFractionFunction;
        final DoubleUnaryOperator priorFunction;

        SampleGenotype(final DoubleUnaryOperator alleleFractionFunction, final DoubleUnaryOperator priorFunction) {
            this.alleleFractionFunction = alleleFractionFunction;
            this.priorFunction = priorFunction;
        }

        public double af(final double maf) {
            return alleleFractionFunction.applyAsDouble(maf);
        }

        public double prior(final double alleleFrequency) {
            return priorFunction.applyAsDouble(alleleFrequency);
        }
    }

    //private static final List<SampleGenotype> SAMPLE_GENOTYPES = Arrays.asList(

    private enum ContaminantGenotype {
        HOM_REF(0, f -> (1 - f) * ( 1 - f)), HET(0.5, f -> 2 * f * (1 - f)), HOM_ALT(1, f -> f * f);

        final double af;
        final DoubleUnaryOperator priorFunction;

        ContaminantGenotype(final double af, final DoubleUnaryOperator priorFunction) {
            this.af = af;
            this.priorFunction = priorFunction;
        }

        public double af() {
            return af;
        }

        public double prior(final double alleleFrequency) {
            return priorFunction.applyAsDouble(alleleFrequency);
        }
    }

    private static double contaminantSum(final ToDoubleFunction<ContaminantGenotype> func) {
        double result = 0;
        for (final ContaminantGenotype genotype : ContaminantGenotype.values()) {
            result += func.applyAsDouble(genotype);
        }
        return result;
    }

    public static double infiniteContaminantLikelihood(final PileupSummary site, final double maf, final double c) {
        final double f = site.getAlleleFrequency();
        final int k = site.getAltCount();
        final int n = k + site.getRefCount();

        return Arrays.stream(SampleGenotype.values())
                .mapToDouble(g -> g.prior(f) * binom(k, n, (1 - c) * g.af(maf) + c * f)).sum();
    }

    public static double singleContaminantLikelihood(final PileupSummary site, final double maf, final double c) {
        final double f = site.getAlleleFrequency();
        final int k = site.getAltCount();
        final int n = k + site.getRefCount();

        return Arrays.stream(SampleGenotype.values())
                .mapToDouble(g -> g.prior(f) * contaminantSum(h -> h.prior(f) * binom(k, n, (1 - c) * g.af(maf) + c * h.af()))).sum();
    }

    public static double twoContaminantLikelihood(final PileupSummary site, final double maf, final double c1, final double c2) {
        final double f = site.getAlleleFrequency();
        final int k = site.getAltCount();
        final int n = k + site.getRefCount();

        return Arrays.stream(SampleGenotype.values())
                .mapToDouble(g -> g.prior(f) *
                        contaminantSum(h -> h.prior(f) *
                                contaminantSum(i -> i.prior(f) * binom(k, n, (1 - c1 - c2) * g.af(maf) + c1 * h.af() + c1 * i.af())))).sum();
    }

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

    private static double binom(final int k, final int n, final double p) {
        return new BinomialDistribution(null, n, p).probability(k);
    }
}
