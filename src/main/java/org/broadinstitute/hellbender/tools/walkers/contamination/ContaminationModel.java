package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.functional.DoubleToDoubleBiFunction;

import java.util.Arrays;
import java.util.List;
import java.util.function.DoubleUnaryOperator;
import java.util.function.ToDoubleFunction;

public class ContaminationModel {
    private ContaminationModelType type;
    private double eps;
    private double c;
    private double c1;
    private double c2;

    public static double IRRELEVANT = Double.NEGATIVE_INFINITY;

    private ContaminationModel(final ContaminationModelType type, final double eps, final double c, final double c1, final double c2) {
        this.type = type;
        this.eps = eps;
        this.c = c;
        this.c1 = c1;
        this.c2 = c2;
    }

    public static ContaminationModel createNoContaminantModel(final double eps) {
        return new ContaminationModel(ContaminationModelType.NO_CONTAMINANT, eps, IRRELEVANT, IRRELEVANT, IRRELEVANT);
    }

    public static ContaminationModel createInfiniteContaminantModel(final double eps, final double c) {
        return new ContaminationModel(ContaminationModelType.INFINITE_CONTAMINANT, eps, c, IRRELEVANT, IRRELEVANT);
    }

    public static ContaminationModel createSingleContaminantModel(final double eps, final double c) {
        return new ContaminationModel(ContaminationModelType.SINGLE_CONTAMINANT, eps, c, c, 0);
    }

    public static ContaminationModel createTwoContaminantModel(final double eps, final double c1, final double c2) {
        return new ContaminationModel(ContaminationModelType.TWO_CONTAMINANT, eps, c1 + c2, c1, c2);
    }

    private static double contaminantSum(final ToDoubleFunction<ContaminantGenotype> func) {
        double result = 0;
        for (final ContaminantGenotype genotype : ContaminantGenotype.ALL_CONTAMINANT_GENOTYPES) {
            result += func.applyAsDouble(genotype);
        }
        return result;
    }

    private static double[] uncontaminatedLikelihoods(final PileupSummary site, final double maf, final double eps) {
        return infiniteContaminantLikelihoods(site, maf, eps, 0);
    }

    private static double[] infiniteContaminantLikelihoods(final PileupSummary site, final double maf, final double eps, final double c) {
        final double f = site.getAlleleFrequency();
        final int k = site.getAltCount();
        final int n = k + site.getRefCount();

        return applyOverSampleGenotypes(g -> g.prior(f) * binom(k, n, (1 - c) * g.af(maf, eps) + c * f));
    }

    private static double[] singleContaminantLikelihoods(final PileupSummary site, final double maf, final double eps, final double c) {
        final double f = site.getAlleleFrequency();
        final int k = site.getAltCount();
        final int n = k + site.getRefCount();

        return applyOverSampleGenotypes(g -> g.prior(f) * contaminantSum(h -> h.prior(f) * binom(k, n, (1 - c) * g.af(maf, eps) + c * h.af(eps))));
    }

    private static double[] twoContaminantLikelihoods(final PileupSummary site, final double maf, final double eps, final double c1, final double c2) {
        final double f = site.getAlleleFrequency();
        final int k = site.getAltCount();
        final int n = k + site.getRefCount();

        return applyOverSampleGenotypes(g -> g.prior(f) *
                        contaminantSum(h -> h.prior(f) *
                                contaminantSum(i -> i.prior(f) * binom(k, n, (1 - c1 - c2) * g.af(maf, eps) + c1 * h.af(eps) + c1 * i.af(eps)))));
    }

    private static double binom(final int k, final int n, final double p) {
        return new BinomialDistribution(null, n, p).probability(k);
    }



    private double[] sampleGenotypeLikelihoods(PileupSummary site, double maf) {
        final double[] likelihoods;
        if(type == ContaminationModelType.NO_CONTAMINANT) {
            likelihoods = uncontaminatedLikelihoods(site, maf, eps);
        } else if (type == ContaminationModelType.INFINITE_CONTAMINANT) {
            likelihoods = infiniteContaminantLikelihoods(site, maf, eps, c);
        } else if (type == ContaminationModelType.SINGLE_CONTAMINANT) {
            likelihoods = singleContaminantLikelihoods(site, maf, eps, c);
        } else if (type == ContaminationModelType.TWO_CONTAMINANT) {
            likelihoods = twoContaminantLikelihoods(site, maf, eps, c1, c2);
        } else {
            throw new GATKException.ShouldNeverReachHereException("Unexpected ContaminationModelType");
        }
        return likelihoods;
    }

    private double likelihood(final PileupSummary site, final double maf) {
        return MathUtils.sum(sampleGenotypeLikelihoods(site, maf));
    }

    public double probability(final PileupSummary site, final double maf, final SampleGenotype sampleGenotype) {
        final double[] likelihoods = sampleGenotypeLikelihoods(site, maf);
        return likelihoods[sampleGenotype.ordinal()] / MathUtils.sum(likelihoods);
    }

    private double logLikelihood(final PileupSummary site, final double maf) {
        return FastMath.log(likelihood(site, maf));
    }

    public double logLikelihood(final List<PileupSummary> segment, final double maf) {
        return segment.stream().mapToDouble(site -> logLikelihood(site, maf)).sum();
    }

    public double logLikelihood(final List<List<PileupSummary>> segments, final List<Double> mafs) {
        Utils.validate(segments.size() == mafs.size(), " Must have one MAF per segment");
        return new IndexRange(0, segments.size()).sum(n -> logLikelihood(segments.get(n), mafs.get(n)));
    }

    private static double[] applyOverSampleGenotypes(final ToDoubleFunction<SampleGenotype> func) {
        final double[] result = new double[SampleGenotype.ALL_SAMPLE_GENOTYPES.size()];
        for (final SampleGenotype sampleGenotype : SampleGenotype.ALL_SAMPLE_GENOTYPES) {
            result[sampleGenotype.ordinal()] = func.applyAsDouble(sampleGenotype);
        }
        return result;
    }

    public enum SampleGenotype {
        HOM_REF((maf, eps) -> eps,f -> (1 - f) * (1 - f)),
        ALT_MINOR((maf, eps) -> maf, f -> f * (1 - f)),
        ALT_MAJOR((maf, eps) -> 1 - maf, f -> f * (1 - f)),
        HOM_ALT((maf, eps) -> 1 - eps, f -> f * f);

        private final DoubleToDoubleBiFunction alleleFractionFunction;
        private final DoubleUnaryOperator priorFunction;
        public static final List<SampleGenotype> ALL_SAMPLE_GENOTYPES = Arrays.asList(values());

        SampleGenotype(final DoubleToDoubleBiFunction alleleFractionFunction, final DoubleUnaryOperator priorFunction) {
            this.alleleFractionFunction = alleleFractionFunction;
            this.priorFunction = priorFunction;
        }

        public double af(final double maf, final double eps) {
            return alleleFractionFunction.apply(maf, eps);
        }

        public double prior(final double alleleFrequency) {
            return priorFunction.applyAsDouble(alleleFrequency);
        }
    }

    public enum ContaminantGenotype {
        HOM_REF(eps -> eps, f -> (1 - f) * ( 1 - f)), HET(eps -> 0.5, f -> 2 * f * (1 - f)), HOM_ALT(eps -> 1 - eps, f -> f * f);

        private final DoubleUnaryOperator alleleFractionFunction;   // read af as function of error rate
        private final DoubleUnaryOperator priorFunction;
        public static final List<ContaminantGenotype> ALL_CONTAMINANT_GENOTYPES = Arrays.asList(values());

        ContaminantGenotype(final DoubleUnaryOperator alleleFractionFunction, final DoubleUnaryOperator priorFunction) {
            this.alleleFractionFunction = alleleFractionFunction;
            this.priorFunction = priorFunction;
        }

        public double af(final double eps) {
            return alleleFractionFunction.applyAsDouble(eps);
        }

        public double prior(final double alleleFrequency) {
            return priorFunction.applyAsDouble(alleleFrequency);
        }
    }
}
