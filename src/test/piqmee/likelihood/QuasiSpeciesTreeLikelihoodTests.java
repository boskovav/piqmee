package test.piqmee.likelihood;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import piqmee.evolution.branchratemodel.QuasiSpeciesUCRelaxedClockModel;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.Test;
import piqmee.likelihood.QuasiSpeciesTreeLikelihood;
import piqmee.tree.QuasiSpeciesNode;
import piqmee.tree.QuasiSpeciesTree;
import test.beast.BEASTTestCase;
import test.piqmee.QuasiSpeciesTestCase;
import beast.base.inference.distribution.LogNormalDistributionModel;
import test.piqmee.branchratemodel.UCRelaxedClockModel_setCategories;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

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
        siteModel.initByName("mutationRate", "13.0", "gammaCategoryCount", 1, "substModel", JC);

        // QS likelihood
        QuasiSpeciesTreeLikelihood likelihood = newQSTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);
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
        assertEquals(logP + (rates[3] * 0.5 * 13) + (rates[1] * 1 * 13), logQSP, BEASTTestCase.PRECISION);

        // with ambiguities
        likelihood.initByName("useAmbiguities", true, "data", data, "tree", tree, "siteModel", siteModel);
        logQSP = likelihood.calculateLogP();
        likelihoodNormal.initByName("useAmbiguities", true, "data", dataNormal, "tree", treeNormal, "siteModel", siteModel);
        logP = likelihoodNormal.calculateLogP();
        assertEquals(logP+ (rates[3] * 0.5 * 13)+ (rates[1] * 1 * 13), logQSP, BEASTTestCase.PRECISION);
    }

// todo how to test uncertain characters in our case? Use orig beast tree? or make a new one, and then what would be the likelihood? calc by hand?
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


    // todo make a test where one qs is above another === and use ambiguities (todo from 26.02)


    /**
     * To make the following 6 test work we needed to copy the UCRelaxedClockModel class and
     * add the following lines to it:
     *
     *
     // this is for testing purposes only
     // set categories
     public void setCategories(int position, int category) {
     categories.setValue(position,category);
     }

     // set categories
     public int getCategories(int position) {
     return categories.getValue(position);
     }

     */

    /**
     * test that the two likelihoods are different even if all rates are the same
     */
    @Test
    public void testJC69AndRelaxedClockLikelihood() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((t3:1.,t4:1.):2.0,(t5:1.5,t6:0.5):0.5):1.0,((t0:1.5,t1:0.5):1.,t2:0.5):1.5);", new String[] {"C", "G", "T", "A", "A", "A", "A"});
        Tree treeNormal = new TreeParser("(((t3:1.,t4:1.):2.0,(t5:1.5,t6:0.5):0.5):1.0,((t0:1.5,t1:0.5):1.,t2:0.5):1.5);",false);

        Alignment data = QuasiSpeciesTestCase.getAlignment(new String[] {"C", "G", "T", "A", "A", "A", "A"});

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "13.0", "gammaCategoryCount", "1", "substModel", JC);

        LogNormalDistributionModel uldistr = new LogNormalDistributionModel();
        uldistr.initByName("M","1.0","S","0.1","meanInRealSpace",true);

        UCRelaxedClockModel_setCategories branchModelNormal = new UCRelaxedClockModel_setCategories();
        branchModelNormal.initByName("distr",uldistr,"rateCategories","13","numberOfDiscreteRates","3","tree",treeNormal);
        // the 13th branch (root-orig branch) has rate 0 so do not set to 1
        for (int i=0; i<12; i++) {
            branchModelNormal.setCategories(i,1);
        }

        QuasiSpeciesUCRelaxedClockModel branchModel = new QuasiSpeciesUCRelaxedClockModel();
        branchModel.initByName("distr",uldistr,"rateCategories","11","numberOfDiscreteRates","3","tree",tree);
        // the 11th branch (root-orig branch) has rate 0 so do not set to 1
        for (int i=0; i<10; i++) {
            branchModel.setCategories(i,1);
        }

        // QS likelihood
        QuasiSpeciesTreeLikelihood likelihood = newQSTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "branchRateModel", branchModel);
        double logQSP = 0;
        logQSP = likelihood.calculateLogP();

        // normal likelihood
        TreeLikelihood likelihoodNormal = newTreeLikelihood();
        likelihoodNormal.initByName("data", data, "tree", treeNormal, "siteModel", siteModel, "branchRateModel", branchModelNormal);
        double logP = 0;
        logP = likelihoodNormal.calculateLogP();

        // compare the two
        assertTrue(logP != logQSP);
    }

    /**
     * test how much different the PIQMEE is from the CLASSIC model, given the qs tree structure
     */
    @Test
    public void testJC69AndRelaxedClockLikelihood1() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((t3:1.,t4:1.):2.0,(t5:1.5,t6:0.5):0.5):1.0,((t0:1.5,t1:0.5):1.,t2:0.5):1.5);", new String[] {"C", "G", "T", "A", "A", "A", "A"});
        Tree treeNormal = new TreeParser("(t3:1.0,((t0:1.5,t1:0.5):1.,t2:0.5):1.5);",false);

        Alignment data = QuasiSpeciesTestCase.getAlignment(new String[] {"C", "G", "T", "A", "A", "A", "A"});
        Alignment dataNormal = QuasiSpeciesTestCase.getAlignment(new String[] {"C", "G", "T", "A"});

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "13.0", "gammaCategoryCount", "1", "substModel", JC);

        LogNormalDistributionModel uldistr = new LogNormalDistributionModel();
        uldistr.initByName("M","1.0","S","0.1","meanInRealSpace",true);

        UCRelaxedClockModel_setCategories branchModelNormal = new UCRelaxedClockModel_setCategories();
        branchModelNormal.initByName("distr",uldistr,"rateCategories","7","numberOfDiscreteRates","3","tree",treeNormal);
        // the 7th branch (root-orig branch) has rate 0 so do not set to 1
        for (int i=0; i<6; i++) {
            branchModelNormal.setCategories(i,1);
        }

        QuasiSpeciesUCRelaxedClockModel branchModel = new QuasiSpeciesUCRelaxedClockModel();
        branchModel.initByName("distr",uldistr,"rateCategories","11","numberOfDiscreteRates","3","tree",tree);
        // the 10th branch (root-orig branch) has rate 0 so do not set to 1
        for (int i=0; i<10; i++) {
            branchModel.setCategories(i,1);
        }

        // QS likelihood
        QuasiSpeciesTreeLikelihood likelihood = newQSTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "branchRateModel", branchModel);
        double logQSP = 0;
        logQSP = likelihood.calculateLogP();

        // normal likelihood
        TreeLikelihood likelihoodNormal = newTreeLikelihood();
        likelihoodNormal.initByName("data", dataNormal, "tree", treeNormal, "siteModel", siteModel, "branchRateModel", branchModelNormal);
        double logP = 0;
        logP = likelihoodNormal.calculateLogP();

        double[] rates = new double[4];
        likelihood.getNoChangeRates(rates);

        // compare the two
        QuasiSpeciesNode newNode = new QuasiSpeciesNode();
        newNode.setNr(0);
        assertEquals(logP+ (rates[1] * 6.5 * 13 * branchModel.getRateForBranch(newNode)), logQSP, BEASTTestCase.PRECISION);
    }

    /**
     * test at the ML value for the CLASSIC model
     */
    @Test
    public void testJC69AndRelaxedClockLikelihoodMLforCLASSIC() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((t3:1.,t4:1.):2.0,(t5:1.5,t6:0.5):0.5):1.0,((t0:1.5,t1:0.5):1.,t2:0.5):1.5);", new String[] {"C", "G", "T", "A", "A", "A", "A"});
        Tree treeNormal = new TreeParser("(((t3:1.,t4:1.):2.0,(t5:1.5,t6:0.5):0.5):1.0,((t0:1.5,t1:0.5):1.,t2:0.5):1.5);",false);

        Alignment data = QuasiSpeciesTestCase.getAlignment(new String[] {"C", "G", "T", "A", "A", "A", "A"});
        Alignment dataNormal = QuasiSpeciesTestCase.getAlignment(new String[] {"C", "G", "T", "A", "A", "A", "A"});

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", "4", "shape", "3.4870306588330533", "substModel", JC);

        LogNormalDistributionModel uldistr = new LogNormalDistributionModel();
        uldistr.initByName("M","1.0","S","3.0427110046894157","meanInRealSpace",true);

        UCRelaxedClockModel_setCategories branchModelNormal = new UCRelaxedClockModel_setCategories();
        branchModelNormal.initByName("distr",uldistr,"clock.rate","9.6561609453625","rateCategories","13","numberOfDiscreteRates","3","tree",treeNormal);
        //from log file categories for branches are 2,2,0,0,0,0,0,2,2,0,0,2,0
        // branches are t0,t1,t2,t3,t4,t5,t6, 7(parent t0,t1),8(parent 7,t2), 9(parent t3,t4), 10(parent t5,t6), 11(parent 9,10), 12 root
        //
        // order of branches in our read in tree is different, so we need to shuffle the rates
        // order of branches here is:
        //              t0,t1,t2,t3,t4,t5,t6, 7(parent t3,t4),8(parent t5,t6), 9(parent 7,8), 10(parent t0,t1), 11(parent 10,t2), 12 root
        //
        // thus new order of categories is          2,2,0,0,0,0,0,0,0,2,2,2,0
        int[] categoriesNormal = {2,2,0,0,0,0,0,0,0,2,2,2,0};
        for (int i=0; i<12; i++) {
            branchModelNormal.setCategories(i,categoriesNormal[i]);
        }

        QuasiSpeciesUCRelaxedClockModel branchModel = new QuasiSpeciesUCRelaxedClockModel();
        branchModel.initByName("distr",uldistr,"clock.rate","9.6561609453625","rateCategories","11","numberOfDiscreteRates","3","tree",tree);
        // from the CLASSIC model above we have
        // tips: t0=2,t1=2,t2=0,t3=0,t4=0,t5=0,t6=0
        //
        // internal nodes: 7(parent of t3,t4)=0; 8(parent of t5,t6)=0, 9(parent of 7,8 === partial branch above A)=2,
        //                 10(parent of t0,t1)=2, 11(parent of 10,t2)=2, 12 (root of the tree)=0
        // average for subtree of haplotype A == 2.5*0+4*0=0 === so first category of 3 --- so category 0
        //                                       here 2.5 is the branch length of branch 7+8 and 4 is the branch length of 3+4+5+6
        int[] categories = {0,0,0,0,2,2,0,2,2,0,2};
        for (int i=0; i<11; i++) {
            branchModel.setCategories(i,categories[i]);
        }

        // QS likelihood
        QuasiSpeciesTreeLikelihood likelihood = newQSTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "branchRateModel", branchModel);
        double logQSP = 0;
        logQSP = likelihood.calculateLogP();

        // normal likelihood
        TreeLikelihood likelihoodNormal = newTreeLikelihood();
        likelihoodNormal.initByName("data", dataNormal, "tree", treeNormal, "siteModel", siteModel, "branchRateModel", branchModelNormal);
        double logP = 0;
        logP = likelihoodNormal.calculateLogP();

        // compare the two
        assertEquals(logP, -5.704853191310973, BEASTTestCase.PRECISION);
        assertEquals(logQSP, -5.704910906518559, BEASTTestCase.PRECISION);

        assertTrue(logP > logQSP);
    }

    /**
     * test at the ML value for the CLASSIC model
     * FOR THIS TEST WE NEEDED TO MAKE calculateIntegratePartials IN LikelihoodCore AND IN BeerLikelihoodCore CLASS
     *  AND likelihoodCore IN TreeLikelihood CLASS PUBLIC
     *  or copy these classes locally to the test folder --- which is what was done here
     */
    @Test
    public void testJC69AndRelaxedClockLikelihoodMLforCLASSIC1() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((t3:1.,t4:1.):2.0,(t5:1.5,t6:0.5):0.5):1.0,((t0:1.5,t1:0.5):1.,t2:0.5):1.5);", new String[] {"C", "G", "T", "A", "A", "A", "A"});
        Tree treeNormal = new TreeParser("(t3:1.0,((t0:1.5,t1:0.5):1.,t2:0.5):1.5);",false);

        Alignment data = QuasiSpeciesTestCase.getAlignment(new String[] {"C", "G", "T", "A", "A", "A", "A"});
        Alignment dataNormal = QuasiSpeciesTestCase.getAlignment(new String[] {"C", "G", "T", "A"});

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", "4", "shape", "3.4870306588330533", "substModel", JC);

        LogNormalDistributionModel uldistr = new LogNormalDistributionModel();
        uldistr.initByName("M","1.0","S","3.0427110046894157","meanInRealSpace",true);

        UCRelaxedClockModel_setCategories branchModelNormal = new UCRelaxedClockModel_setCategories();
        branchModelNormal.initByName("distr",uldistr,"clock.rate","9.6561609453625","rateCategories","7","numberOfDiscreteRates","3","tree",treeNormal);
        //from log file categories for branches are 2,2,0,0,0,0,0,2,2,0,0,2,0
        // branches are t0,t1,t2,t3,t4,t5,t6, 7(parent t0,t1),8(parent 7,t2), 9(parent t3,t4), 10(parent t5,t6), 11(parent 9,10), 12 root
        //
        // order of branches in our read in tree is different, so we need to shuffle the rates
        // order of branches here is:
        //              t0,t1,t2,t3,t4,t5,t6, 7(parent t3,t4),8(parent t5,t6), 9(parent 7,8), 10(parent t0,t1), 11(parent 10,t2), 12 root
        //
        // thus new order of categories is          2,2,0,0,0,0,0,0,0,2,2,2,0
        int[] categoriesNormal = {2,2,0,2,2,2,0};
        for (int i=0; i<6; i++) {
            branchModelNormal.setCategories(i,categoriesNormal[i]);
        }

        QuasiSpeciesUCRelaxedClockModel branchModel = new QuasiSpeciesUCRelaxedClockModel();
        branchModel.initByName("distr",uldistr,"clock.rate","9.6561609453625","rateCategories","11","numberOfDiscreteRates","3","tree",tree);
        // from the CLASSIC model above we have
        // tips: t0=2,t1=2,t2=0,t3=0,t4=0,t5=0,t6=0
        //
        // internal nodes: 7(parent of t3,t4)=0; 8(parent of t5,t6)=0, 9(parent of 7,8 === partial branch above A)=2,
        //                 10(parent of t0,t1)=2, 11(parent of 10,t2)=2, 12 (root of the tree)=0
        // average for subtree of haplotype A == 2.5*0+4*0=0 === so first catefory of 3 --- so category 0
        //                                       here 2.5 is the branch length of branch 7+8 and 4 is the branch length of 3+4+5+6
        int[] categories = {0,0,0,0,2,2,0,2,2,0,2};
        for (int i=0; i<11; i++) {
            branchModel.setCategories(i,categories[i]);
        }

        // QS likelihood
        QuasiSpeciesTreeLikelihood likelihood = newQSTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "branchRateModel", branchModel);
        double logQSP = 0;
        logQSP = likelihood.calculateLogP();

        // normal likelihood
        test.piqmee.likelihood.TreeLikelihood likelihoodNormal = new test.piqmee.likelihood.TreeLikelihood();
        likelihoodNormal.initByName("data", dataNormal, "tree", treeNormal, "siteModel", siteModel, "branchRateModel", branchModelNormal);
        double logP = 0;
        logP = likelihoodNormal.calculateLogP();

        double[] rates = new double[4];
        likelihood.getNoChangeRates(rates);

        // compare the two
        QuasiSpeciesNode newNode = new QuasiSpeciesNode();
        newNode.setNr(3);
        double[] finalPartials = new double[4];
        double[] rootPartials = new double[16];
        double[] partials = new double[16];
        likelihoodNormal.likelihoodCore.getNodePartials(treeNormal.getRoot().getNr(),partials);
        int v=0;
        for (int i=0; i<4; i++) {
            for (int j=0; j < siteModel.getCategoryCount(); j++) {
                rootPartials[v] = partials[v] * Math.exp(rates[1] * 6.5 * branchModel.getRateForBranch(newNode) * siteModel.getCategoryRates(newNode)[i]);
                v++;
            }
        }

        final double[] proportions = siteModel.getCategoryProportions(treeNormal.getRoot());
        likelihoodNormal.likelihoodCore.calculateIntegratePartials(rootPartials, proportions, finalPartials);
        double[] frequencies = siteModel.substModelInput.get().getFrequencies();
        double[] logPnew = new double[1];
        likelihoodNormal.likelihoodCore.calculateLogLikelihoods(finalPartials, frequencies, logPnew);

        assertEquals(logPnew[0], logQSP, BEASTTestCase.PRECISION);
    }

    /**
     * test at the ML value for the CLASSIC model
     */
    @Test
    public void testJC69AndRelaxedClockLikelihoodMLforPIQMEE() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((t3:1.,t4:1.):2.0,(t5:1.5,t6:0.5):0.5):1.0,((t0:1.5,t1:0.5):1.,t2:0.5):1.5);", new String[] {"C", "G", "T", "A", "A", "A", "A"});
        Tree treeNormal = new TreeParser("(((t3:1.,t4:1.):2.0,(t5:1.5,t6:0.5):0.5):1.0,((t0:1.5,t1:0.5):1.,t2:0.5):1.5);",false);

        Alignment data = QuasiSpeciesTestCase.getAlignment(new String[] {"C", "G", "T", "A", "A", "A", "A"});
        Alignment dataNormal = QuasiSpeciesTestCase.getAlignment(new String[] {"C", "G", "T", "A", "A", "A", "A"});

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", "4", "shape", "4.063052723241771", "substModel", JC);

        LogNormalDistributionModel uldistr = new LogNormalDistributionModel();
        uldistr.initByName("M","1.0","S","2.81780094730661","meanInRealSpace",true);

        QuasiSpeciesUCRelaxedClockModel branchModel = new QuasiSpeciesUCRelaxedClockModel();
        branchModel.initByName("distr",uldistr,"clock.rate","7.462630258553472","rateCategories","11","numberOfDiscreteRates","3","tree",tree);
        //from log file categories for branches are 1,0,1,0,2,2,2,2,2,2,2
        // branches are t0,t1,t2,t3, 4(parent t0,t1),5(parent 4,t2), 6 root, 7/8/9/10 branch above haplo 0/1/2/3
        //
        // order of branches in our read in tree the same
        int[] categories = {1,0,1,0,2,2,2,2,2,2,2};
        for (int i=0; i<11; i++) {
            branchModel.setCategories(i,categories[i]);
        }

        UCRelaxedClockModel_setCategories branchModelNormal = new UCRelaxedClockModel_setCategories();
        branchModelNormal.initByName("distr",uldistr,"clock.rate","7.462630258553472","rateCategories","13","numberOfDiscreteRates","3","tree",treeNormal);
        // from the PIQMEE model above we have
        // tips: t0=1,t1=0,t2=1, are irrelevant, because they are not used
        //   true rate for tips is the rate above these tips so tip t0(node nr7)=2,t1(node nr8)=2,t2(node nr9)=2
        //   for tip representing the subtree of haplotype A we have t3=0 and rate of the branch above t3 (node nr10)=2
        //
        // internal nodes: 4(parent of t0,t1)=2; 5(parent of 2,t2)=2, 6 (root of the tree)=2 but is irrelevant
        int[] categoriesNormal = {2,2,2,0,0,0,0,0,0,2,2,2,0};
        for (int i=0; i<12; i++) {
            branchModelNormal.setCategories(i,categoriesNormal[i]);
        }

        // QS likelihood
        QuasiSpeciesTreeLikelihood likelihood = newQSTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel, "branchRateModel", branchModel);
        double logQSP = 0;
        logQSP = likelihood.calculateLogP();

        // normal likelihood
        TreeLikelihood likelihoodNormal = newTreeLikelihood();
        likelihoodNormal.initByName("data", dataNormal, "tree", treeNormal, "siteModel", siteModel, "branchRateModel", branchModelNormal);
        double logP = 0;
        logP = likelihoodNormal.calculateLogP();

        // compare the two
        assertEquals(logP, -5.646908588976165, BEASTTestCase.PRECISION);
        assertEquals(logQSP, -5.647090960029157, BEASTTestCase.PRECISION);
        // note that logP > logQSP by very little but still, the main difference is that at the ML value of CLASSIC model (test just above)
        // PIQMEE has a way lower likelihood value, while CLASSIC has just a slightly lower value
    }

    /**
     * test at the median likelihood value for the CLASSIC model
     */
    @Test
    public void testJC69AndRelaxedClockLikelihoodMedianforCLASSIC() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((t3:1.,t4:1.):2.0,(t5:1.5,t6:0.5):0.5):1.0,((t0:1.5,t1:0.5):1.,t2:0.5):1.5);", new String[] {"C", "G", "T", "A", "A", "A", "A"});
        Tree treeNormal = new TreeParser("(((t3:1.,t4:1.):2.0,(t5:1.5,t6:0.5):0.5):1.0,((t0:1.5,t1:0.5):1.,t2:0.5):1.5);",false);

        Alignment data = QuasiSpeciesTestCase.getAlignment(new String[] {"C", "G", "T", "A", "A", "A", "A"});
        Alignment dataNormal = QuasiSpeciesTestCase.getAlignment(new String[] {"C", "G", "T", "A", "A", "A", "A"});

        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel1 = new SiteModel();
        siteModel1.initByName("mutationRate", "1.0", "gammaCategoryCount", "4", "shape", "0.5612565047930533", "substModel", JC);

        SiteModel siteModel2 = new SiteModel();
        siteModel2.initByName("mutationRate", "1.0", "gammaCategoryCount", "4", "shape", "0.7755492449693009", "substModel", JC);


        LogNormalDistributionModel uldistr1 = new LogNormalDistributionModel();
        uldistr1.initByName("M","1.0","S","0.8806322765682721","meanInRealSpace",true);

        LogNormalDistributionModel uldistr2 = new LogNormalDistributionModel();
        uldistr2.initByName("M","1.0","S","2.9324292919674435","meanInRealSpace",true);

        UCRelaxedClockModel_setCategories branchModelNormal1 = new UCRelaxedClockModel_setCategories();
        branchModelNormal1.initByName("distr",uldistr1,"clock.rate","7.835061389720933","rateCategories","13","numberOfDiscreteRates","3","tree",treeNormal);
        UCRelaxedClockModel_setCategories branchModelNormal2 = new UCRelaxedClockModel_setCategories();
        branchModelNormal2.initByName("distr",uldistr2,"clock.rate","9.349074968454529","rateCategories","13","numberOfDiscreteRates","3","tree",treeNormal);
        // likelihood value -9.1312
        // Sample	posterior	likelihood	prior	BDTreelikelihood	TreeHeight	ucldMean	ucldStdev	rate.mean	rate.variance	rate.coefficientOfVariation	categoryForBranchNr1	categoryForBranchNr2	categoryForBranchNr3	categoryForBranchNr4	categoryForBranchNr5	categoryForBranchNr6	categoryForBranchNr7	categoryForBranchNr8	categoryForBranchNr9	categoryForBranchNr10	categoryForBranchNr11	categoryForBranchNr12	categoryForBranchNr13	gammaShape	origin	reproductiveNumber	becomeUninfectiousRate	samplingProportion	rho	birth	growth-rate	birthRho	death	sampling
        // 117380000	-1090.188086168279	-9.131287623950689	-1081.0567985443283	-1045.2636184640642	4.0	7.835061389720933	0.8806322765682721	5.320670415838548	11.620493928657105	0.5932991089843568	1	0	2	1	0	1	1	2	1	0	1	1	0	0.5612565047930533	1000.0	2.0	1.0	1.0E-4	1.0	2.0	1.0	2.0	0.9999	1.0E-4
        // 431800000	-1090.4023084330631	-9.13121714855862	-1081.2710912845046	-1045.2636184640642	4.0	9.349074968454529	2.9324292919674435	0.7553291837896942	1.054574266563593	1.3224439761724611	2	2	0	2	0	1	1	2	1	1	1	0	0	0.7755492449693009	1000.0	2.0	1.0	1.0E-4	1.0	2.0	1.0	2.0	0.9999	1.0E-4
        //from log file categories for branches are 1	0	2	1	0	1	1	2	1	0	1	1	0 (step 117380000)
        //                                      and 2	2	0	2	0	1	1	2	1	1	1	0	0 (step 431800000)
        // branches are t0,t1,t2,t3,t4,t5,t6, 7(parent t0,t1),8(parent 7,t2), 9(parent t3,t4), 10(parent t5,t6), 11(parent 9,10), 12 root
        //
        // order of branches in our read in tree is different, so we need to shuffle the rates
        // order of branches here is:
        //              t0,t1,t2,t3,t4,t5,t6, 7(parent t3,t4),8(parent t5,t6), 9(parent 7,8), 10(parent t0,t1), 11(parent 10,t2), 12 root
        //
        // thus new order of categories is          1,0,2,1,0,1,1,0,1,1,2,1,0  &   2,2,0,2,0,1,1,1,1,0,2,1,0
        int[] categoriesNormal1   = {1,0,2,1,0,1,1,0,1,1,2,1,0};
        int[] categoriesNormal2   = {2,2,0,2,0,1,1,1,1,0,2,1,0};
        for (int i=0; i<12; i++) {
            branchModelNormal1.setCategories(i,categoriesNormal1[i]);
            branchModelNormal2.setCategories(i,categoriesNormal2[i]);
        }

        QuasiSpeciesUCRelaxedClockModel branchModel1 = new QuasiSpeciesUCRelaxedClockModel();
        branchModel1.initByName("distr",uldistr1,"clock.rate","9.6561609453625","rateCategories","11","numberOfDiscreteRates","3","tree",tree);
        QuasiSpeciesUCRelaxedClockModel branchModel2 = new QuasiSpeciesUCRelaxedClockModel();
        branchModel2.initByName("distr",uldistr2,"clock.rate","9.6561609453625","rateCategories","11","numberOfDiscreteRates","3","tree",tree);
        // search through PIQMEE lig file for likelihood values -9.1312 gives zero hits
        //      for value -9.131 gives 4 hits
        //  2	2	2	0	0	0	0	0	1	0	1
        //  2	0	1	1	1	2	2	1	2	1	0
        //  1	1	0	1	1	2	0	0	2	1	2
        //  1	2	0	0	2	1	2	0	1	2	1
        //
        // from the CLASSIC model above we have for the first set: 1,0,2,1,0,1,1,0,1,1,2,1,0
        // tips: t0=1,t1=0,t2=2,t3=1,t4=0,t5=1,t6=1  - so for piqmee tip t3=1*1 t4=1*0 t5=1.5*1 t6=0.5*1
        //
        // internal nodes: 7(parent of t3,t4)=0; 8(parent of t5,t6)=1, 9(parent of 7,8 === partial branch above A)=1,
        //                 10(parent of t0,t1)=2, 11(parent of 10,t2)=1, 12 (root of the tree)=0 but is irrelevant
        // average for subtree of haplotype A == 2.5*0 + 1*1 + 1*0 + 1.5*1 + 0.5*1 = 0.46 === so first category of 3 --- so category 0
        //                                       here 2.5 is the branch length of branch 7+8 and 4 is the branch length of 3+4+5+6
        int[] categories1 = {0,0,0,0,2,1,0,1,0,2,1};
        // from the CLASSIC model above we have for the first set: 2,2,0,2,0,1,1,1,1,0,2,1,0
        // tips: t0=2,t1=2,t2=0,t3=2,t4=0,t5=1,t6=1  - so for piqmee tip t3=1*2 t4=1*0 t5=1.5*1 t6=0.5*1
        //
        // internal nodes: 7(parent of t3,t4)=1; 8(parent of t5,t6)=1, 9(parent of 7,8 === partial branch above A)=0,
        //                 10(parent of t0,t1)=2, 11(parent of 10,t2)=1, 12 (root of the tree)=0 but is irrelevant
        // average for subtree of haplotype A == 2.5*1 + 1*2 + 1*0 + 1.5*1 + 0.5*1 = 1 === so second category of 3 --- so category 1
        //                                       here 2.5 is the branch length of branch 7+8 and 4 is the branch length of 3+4+5+6
        int[] categories2 = {0,0,0,1,2,1,0,2,2,0,0};
        for (int i=0; i<11; i++) {
            branchModel1.setCategories(i,categories1[i]);
            branchModel2.setCategories(i,categories2[i]);
        }

        // QS likelihood
        QuasiSpeciesTreeLikelihood likelihood1 = newQSTreeLikelihood();
        likelihood1.initByName("data", data, "tree", tree, "siteModel", siteModel1, "branchRateModel", branchModel1);
        double logQSP1 = 0;
        logQSP1 = likelihood1.calculateLogP();
        QuasiSpeciesTreeLikelihood likelihood2 = newQSTreeLikelihood();
        likelihood2.initByName("data", data, "tree", tree, "siteModel", siteModel2, "branchRateModel", branchModel2);
        double logQSP2 = 0;
        logQSP2 = likelihood2.calculateLogP();

        // normal likelihood
        TreeLikelihood likelihoodNormal1 = newTreeLikelihood();
        likelihoodNormal1.initByName("data", dataNormal, "tree", treeNormal, "siteModel", siteModel1, "branchRateModel", branchModelNormal1);
        double logP1 = 0;
        logP1 = likelihoodNormal1.calculateLogP();
        TreeLikelihood likelihoodNormal2 = newTreeLikelihood();
        likelihoodNormal2.initByName("data", dataNormal, "tree", treeNormal, "siteModel", siteModel2, "branchRateModel", branchModelNormal2);
        double logP2 = 0;
        logP2 = likelihoodNormal2.calculateLogP();

        // compare the two
        assertEquals(logP1, -9.131287623950689, BEASTTestCase.PRECISION);
        assertEquals(logP2, -9.13121714855862, BEASTTestCase.PRECISION);
        assertEquals(logQSP1, -9.071017893578265, BEASTTestCase.PRECISION);
        assertEquals(logQSP2, -8.409375927260333, BEASTTestCase.PRECISION);
        assertTrue(logP1 < logQSP1);
        assertTrue(logP2 < logQSP2);
    }
}
