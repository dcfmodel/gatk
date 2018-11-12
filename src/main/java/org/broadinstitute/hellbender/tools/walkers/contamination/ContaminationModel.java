package org.broadinstitute.hellbender.tools.walkers.contamination;

import org.broadinstitute.hellbender.exceptions.GATKException;

public class ContaminationModel {
    private ContaminationModelType type;
    private double c;
    private double c1;
    private double c2;

    public static double IRRELEVANT = Double.NEGATIVE_INFINITY;

    private ContaminationModel(final ContaminationModelType type, final double c, final double c1, final double c2) {
        this.type = type;
        this.c = c;
        this.c1 = c1;
        this.c2 = c2;
    }

    public static ContaminationModel createNoContaminantModel() {
        return new ContaminationModel(ContaminationModelType.NO_CONTAMINANT, IRRELEVANT, IRRELEVANT, IRRELEVANT);
    }

    public static ContaminationModel createInfiniteContaminantModel(final double c) {
        return new ContaminationModel(ContaminationModelType.INFINITE_CONTAMINANT, c, IRRELEVANT, IRRELEVANT);
    }

    public static ContaminationModel createSingleContaminantModel(final double c) {
        return new ContaminationModel(ContaminationModelType.SINGLE_CONTAMINANT, c, c, 0);
    }

    public static ContaminationModel createTwoContaminantModel(final double c1, final double c2) {
        return new ContaminationModel(ContaminationModelType.TWO_CONTAMINANT, c1 + c2, c1, c2);
    }

    public double likelihood(final PileupSummary site, final double maf) {
        if(type == ContaminationModelType.NO_CONTAMINANT) {
            return ContaminationEngine.uncontaminatedLikelihood(site, maf);
        } else if (type == ContaminationModelType.INFINITE_CONTAMINANT) {
            return ContaminationEngine.infiniteContaminantLikelihood(site, maf, c);
        } else if (type == ContaminationModelType.SINGLE_CONTAMINANT) {
            return ContaminationEngine.singleContaminantLikelihood(site, maf, c);
        } else if (type == ContaminationModelType.TWO_CONTAMINANT) {
            return ContaminationEngine.twoContaminantLikelihood(site, maf, c1, c2);
        } else {
            throw new GATKException.ShouldNeverReachHereException("Unexpected ContaminationModelType");
        }
    }
}
