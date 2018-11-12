package org.broadinstitute.hellbender.utils.functional;

@FunctionalInterface
public interface DoubleToDoubleBiFunction {
    double apply(final double x, final double y);
}
