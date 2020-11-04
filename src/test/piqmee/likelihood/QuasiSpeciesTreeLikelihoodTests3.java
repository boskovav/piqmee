package test.piqmee.likelihood;


import beast.evolution.likelihood.GenericTreeLikelihood;
import piqmee.likelihood.QuasiSpeciesTreeLikelihood3;

public class QuasiSpeciesTreeLikelihoodTests3 extends QuasiSpeciesTreeLikelihoodTests2 {

    protected GenericTreeLikelihood newQSTreeLikelihood() {
        System.setProperty("java.only","true");
        return new QuasiSpeciesTreeLikelihood3();
    }

}
