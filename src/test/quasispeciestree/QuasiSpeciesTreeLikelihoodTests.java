package test.quasispeciestree;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Test;
import quasispeciestree.likelihood.QuasiSpeciesTreeLikelihood;
import quasispeciestree.tree.QuasiSpeciesTree;
import test.beast.BEASTTestCase;
import test.quasispeciestree.QuasiSpeciesTestCase;

import static org.junit.Assert.assertEquals;

/**
 * @author Veronika Boskova created on 09/02/2018 finished on 0?/0?/201?
 */
public class QuasiSpeciesTreeLikelihoodTests {

    /**
     *
     * Likelihood calculation P(data|tree) testing
     *
     */


    protected QuasiSpeciesTreeLikelihood newQSTreeLikelihood() {
        System.setProperty("java.only","true");
        return new QuasiSpeciesTreeLikelihood();
    }

    protected TreeLikelihood newTreeLikelihood() {
        System.setProperty("java.only","true");
        return new TreeLikelihood();
    }

    public double[] testJC69Likelihood(Alignment data, Tree tree) throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "0.6", "substModel", JC);
        // NB The rate in the JC model used here is actually alpha * 3 in the usual sense, because
        // it's divided by 3 before multiplying in the exponent (not sure why)

        System.out.println("Without tip likelihoods:");
        QuasiSpeciesTreeLikelihood likelihood = newQSTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "scaling", TreeLikelihood.Scaling.none);
        double[] logP = new double[2];
        logP[0] = likelihood.calculateLogP();
        System.out.println(logP[0]);

        System.out.println("With tip likelihoods:");
        likelihood.initByName("useTipLikelihoods", true, "data", data, "tree", tree, "siteModel", siteModel, "scaling", TreeLikelihood.Scaling.none);
        logP[1]= likelihood.calculateLogP();
        System.out.println(logP[1]);

        return logP;
    }

    @Test
    public void testJC69Likelihood() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("((t2 : 0.5, t3 : 0.5) : 0.5, (t0 : 0.25, t1 : 0.25) : 0.25);", new String[] {"A", "A", "C", "C"});
        Tree treeNormal = new TreeParser("( t1 : 0.25, t0 : 0.5 ) : 1;",false);

        Alignment data = QuasiSpeciesTestCase.getAlignment(new String[] {"A", "A", "C", "C"});
        Alignment dataNormal = QuasiSpeciesTestCase.getAlignment(new String[] {"A", "C"});

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);

        // QS likelihood
        QuasiSpeciesTreeLikelihood likelihood = newQSTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "origin", new RealParameter("1.000000000000000000000000000000000001"));
        double logQSP = 0;
        logQSP = likelihood.calculateLogP();

        // normal likelihood
        TreeLikelihood likelihoodNormal = newTreeLikelihood();
        likelihoodNormal.initByName("data", dataNormal, "tree", treeNormal, "siteModel", siteModel);
        double logP = 0;
        logP = likelihoodNormal.calculateLogP();

        double[] rates = new double[4];
        likelihood.getNoChangeRates(rates);

        // compare the two
        assertEquals(logP+Math.log(Math.exp(rates[3] * 0.5))+Math.log(Math.exp(rates[1] * 1)), logQSP, BEASTTestCase.PRECISION);

        // with ambiguities
        likelihood.initByName("useAmbiguities", true, "data", data, "tree", tree, "siteModel", siteModel, "origin", new RealParameter("1.000000000000000000000000000000000001"));
        logQSP = likelihood.calculateLogP();
        likelihoodNormal.initByName("useAmbiguities", true, "data", data, "tree", tree, "siteModel", siteModel);
        logP = likelihoodNormal.calculateLogP();
        assertEquals(logP+Math.log(Math.exp(rates[3] * 0.5))+Math.log(Math.exp(rates[1] * 1)), logQSP, BEASTTestCase.PRECISION);
    }

//    @Test
//    public void testJC69LikelihoodWithUncertainCharacters() throws Exception {
//
//        Alignment data = UncertainAlignmentTest.getAlignment();
//        Alignment data2 = UncertainAlignmentTest.getUncertainAlignment();
//        double[] logL, logL_uncertain;
//
//        System.out.println("\nTree A:");
//        Tree tree = UncertainAlignmentTest.getTreeA(data2);
//        logL = testJC69Likelihood(data,tree);
//        logL_uncertain = testJC69Likelihood(data2,tree);
//        double x1 = -11.853202336328778;
//        double x2 = -12.069603116476458;
//        assertEquals(logL[0], x1, BEASTTestCase.PRECISION);
//        assertEquals(logL[1], x1, BEASTTestCase.PRECISION);
//        assertEquals(logL_uncertain[0], x1, BEASTTestCase.PRECISION);
//        assertEquals(logL_uncertain[1], x2, BEASTTestCase.PRECISION);
//
//        System.out.println("\nTree B:");
//        tree = UncertainAlignmentTest.getTreeB(data2);
//        logL = testJC69Likelihood(data,tree);
//        logL_uncertain = testJC69Likelihood(data2,tree);
//        double x3 = -12.421114302827698;
//        double x4 = -11.62105662310513;
//        assertEquals(logL[0], x3, BEASTTestCase.PRECISION);
//        assertEquals(logL[1], x3, BEASTTestCase.PRECISION);
//        assertEquals(logL_uncertain[0], x3, BEASTTestCase.PRECISION);
//        assertEquals(logL_uncertain[1], x4, BEASTTestCase.PRECISION);
//
//        System.out.println("\nTesting alignment doubling:");
//        Alignment data3 = UncertainAlignmentTest.getUncertainAlignmentDoubled();
//        logL_uncertain = testJC69Likelihood(data3,tree);
//        assertEquals(logL_uncertain[0], 2 * x3, BEASTTestCase.PRECISION);
//        assertEquals(logL_uncertain[1], 2 * x4, BEASTTestCase.PRECISION);
//
//    }

}
