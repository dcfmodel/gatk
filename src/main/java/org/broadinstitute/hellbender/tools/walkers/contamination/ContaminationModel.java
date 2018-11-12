package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.functional.DoubleToDoubleBiFunction;

import java.util.Arrays;
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
        for (final ContaminantGenotype genotype : ContaminantGenotype.values()) {
            result += func.applyAsDouble(genotype);
        }
        return result;
    }

    private static double uncontaminatedLikelihood(final PileupSummary site, final double maf, final double eps) {
        return infiniteContaminantLikelihood(site, maf, eps, 0);
    }

    private static double infiniteContaminantLikelihood(final PileupSummary site, final double maf, final double eps, final double c) {
        final double f = site.getAlleleFrequency();
        final int k = site.getAltCount();
        final int n = k + site.getRefCount();

        return Arrays.stream(SampleGenotype.values())
                .mapToDouble(g -> g.prior(f) * binom(k, n, (1 - c) * g.af(maf, eps) + c * f)).sum();
    }

    private static double singleContaminantLikelihood(final PileupSummary site, final double maf, final double eps, final double c) {
        final double f = site.getAlleleFrequency();
        final int k = site.getAltCount();
        final int n = k + site.getRefCount();

        return Arrays.stream(SampleGenotype.values())
                .mapToDouble(g -> g.prior(f) * contaminantSum(h -> h.prior(f) * binom(k, n, (1 - c) * g.af(maf, eps) + c * h.af(eps)))).sum();
    }

    private static double twoContaminantLikelihood(final PileupSummary site, final double maf, final double eps, final double c1, final double c2) {
        final double f = site.getAlleleFrequency();
        final int k = site.getAltCount();
        final int n = k + site.getRefCount();

        return Arrays.stream(SampleGenotype.values())
                .mapToDouble(g -> g.prior(f) *
                        contaminantSum(h -> h.prior(f) *
                                contaminantSum(i -> i.prior(f) * binom(k, n, (1 - c1 - c2) * g.af(maf, eps) + c1 * h.af(eps) + c1 * i.af(eps))))).sum();
    }

    private static double binom(final int k, final int n, final double p) {
        return new BinomialDistribution(null, n, p).probability(k);
    }

    public double likelihood(final PileupSummary site, final double maf) {
        if(type == ContaminationModelType.NO_CONTAMINANT) {
            return uncontaminatedLikelihood(site, maf, eps);
        } else if (type == ContaminationModelType.INFINITE_CONTAMINANT) {
            return infiniteContaminantLikelihood(site, maf, eps, c);
        } else if (type == ContaminationModelType.SINGLE_CONTAMINANT) {
            return singleContaminantLikelihood(site, maf, eps, c);
        } else if (type == ContaminationModelType.TWO_CONTAMINANT) {
            return twoContaminantLikelihood(site, maf, eps, c1, c2);
        } else {
            throw new GATKException.ShouldNeverReachHereException("Unexpected ContaminationModelType");
        }
    }

    public double logLikelihood(final PileupSummary site, final double maf) {
        return FastMath.log(likelihood(site, maf));
    }

    private enum SampleGenotype {
        HOM_REF((maf, eps) -> eps,f -> (1 - f) * (1 - f)),
        ALT_MINOR((maf, eps) -> maf, f -> f * (1 - f)),
        ALT_MAJOR((maf, eps) -> 1 - maf, f -> f * (1 - f)),
        HOM_ALT((maf, eps) -> 1 - eps, f -> f * f);

        final DoubleToDoubleBiFunction alleleFractionFunction;
        final DoubleUnaryOperator priorFunction;

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

    private enum ContaminantGenotype {
        HOM_REF(eps -> eps, f -> (1 - f) * ( 1 - f)), HET(eps -> 0.5, f -> 2 * f * (1 - f)), HOM_ALT(eps -> 1 - eps, f -> f * f);

        final DoubleUnaryOperator alleleFractionFunction;   // read af as function of error rate
        final DoubleUnaryOperator priorFunction;

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
