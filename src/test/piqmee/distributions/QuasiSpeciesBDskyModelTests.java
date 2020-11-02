package piqmee.distributions;

import static org.junit.Assert.assertEquals;
import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Test;

import piqmee.distributions.BirthDeathSkylineModel;
import piqmee.distributions.QuasiSpeciesBirthDeathSkylineModel;
import piqmee.tree.QuasiSpeciesTree;
import test.piqmee.QuasiSpeciesTestCase;

/**
 * @author Veronika Boskova created on 15/07/16
 */

@Description("Test the piqmee birth-death model with a small tree example")
public class QuasiSpeciesBDskyModelTests {


    /**
     *
     * Tree prior calculation P(tree|tree parameters) testing - 4 haplo with duplicates
     *
     */

    private QuasiSpeciesBirthDeathSkylineModel getQSBDSKYmodel(QuasiSpeciesTree tree, RealParameter origin, boolean conditionOnSurvival,
                                                               RealParameter birth, RealParameter death, RealParameter sampling) {

        QuasiSpeciesBirthDeathSkylineModel bdsqs =  new QuasiSpeciesBirthDeathSkylineModel();

        bdsqs.setInputValue("tree", tree);
        bdsqs.setInputValue("origin", origin);
        bdsqs.setInputValue("conditionOnSurvival", conditionOnSurvival);

        // test without rate change
        bdsqs.setInputValue("birthRate", birth);
        bdsqs.setInputValue("deathRate", death);
        bdsqs.setInputValue("samplingRate", sampling);

        bdsqs.initAndValidate();
        bdsqs.printTempResults = true;

        return bdsqs;
    }

    private BirthDeathSkylineModel getBDSKYmodel(Tree tree, RealParameter origin, boolean conditionOnSurvival,
                                                 RealParameter birth, RealParameter death, RealParameter sampling) {
        BirthDeathSkylineModel bdssm = new BirthDeathSkylineModel();

        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", origin);
        bdssm.setInputValue("conditionOnSurvival", conditionOnSurvival);

        // test without rate change
        bdssm.setInputValue("birthRate", birth);
        bdssm.setInputValue("deathRate", death);
        bdssm.setInputValue("samplingRate", sampling);

        bdssm.initAndValidate();
        bdssm.printTempResults = true;

        return bdssm;
    }

    // test without rate change
    @Test
    public void testLikelihoodCalculation1() throws Exception {

        // Assemble BEASTObjects needed by QuasiSpeciesTree
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromNewick("((t3 : 1.5, t0 : 0.5) : 1 , (t1 : 2, t2 : 1) : 3);", new String[] {"A","C","G","T"});

        // calculate the tree likelihood with QS algorithm
        QuasiSpeciesBirthDeathSkylineModel bdsqs =  this.getQSBDSKYmodel(tree, new RealParameter("6.0"),
                false, new RealParameter("2.0"), new RealParameter("1.0"),
                new RealParameter("0.5"));

        // calculate the basis of the tree likelihood with bdsky
        BirthDeathSkylineModel bdssm =  this.getBDSKYmodel(tree, new RealParameter("6.0"),
                false, new RealParameter("2.0"), new RealParameter("1.0"),
                new RealParameter("0.5"));

        double bdsky=bdssm.calculateTreeLogLikelihood(tree);
        double t1_1 = 0.5663474;
        double t2_1 = -0.1769383;
        double t2_2 = -0.60098;
        double t3_1 = -3.147534;
        double t3_2 = -3.761309;
        double t3_3 = -4.377242;
        double t1_tips = 3*2.5377465380121595;
        double t2_tips = 1*-0.6931471805599453;
        double t3_tips = 2*-0.19427342018656657;
        // Since all tips have haplo duplicates sampled only at a single time point, no factor for number of full trees = constant
        // double nrFullTrees = Math.log(4*3*3*2*2*1 * 3*2*2*1 * 2*1);
        assertEquals(bdsky+bdsqs.logNumberOfQSTrees(tree)+t1_1+t2_1+t2_2+t3_1+t3_2+t3_3+t1_tips+t2_tips+t3_tips, bdsqs.calculateTreeLogLikelihood(tree), 1e-5);
    }

    @Test
    public void testLikelihoodCalculation2() throws Exception {

        // Assemble BEASTObjects needed by QuasiSpeciesTree
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromNewick("((t3 : 1.5, t0 : 0.5) : 1 , (t1 : 2, t2 : 1) : 3);", new String[] {"A","C","G","T"});

        // calculate the tree likelihood with QS algorithm
        QuasiSpeciesBirthDeathSkylineModel bdsqs =  this.getQSBDSKYmodel(tree, new RealParameter("6.0"),
                false, new RealParameter("2.0"), new RealParameter("1.0"),
                new RealParameter("0.5"));

        // calculate the basis of the tree likelihood with bdsky
        BirthDeathSkylineModel bdssm =  this.getBDSKYmodel(tree, new RealParameter("6.0"),
                false, new RealParameter("2.0"), new RealParameter("1.0"),
                new RealParameter("0.5"));

        double bdsky=bdssm.calculateTreeLogLikelihood(tree);
        double t1_1 = 0.5663474;
        double t2_1 = -0.1769383;
        double t2_2 = -0.60098;
        double t3_1 = -3.147534;
        double t3_2 = -3.761309;
        double t3_3 = -4.377242;
        double t1_tips = 3*2.5377465380121595;
        double t2_tips = 1*-0.6931471805599453;
        double t3_tips = 2*-0.19427342018656657;
        assertEquals(bdsky+bdsqs.logNumberOfQSTrees(tree)+t1_1+t2_1+t2_2+t3_1+t3_2+t3_3+t1_tips+t2_tips+t3_tips, bdsqs.calculateTreeLogLikelihood(tree), 1e-5);
    }

    /**
     *
     * Tree prior calculation P(tree|tree parameters) testing - lineageCountAtTime
     *
     */

    // test without rate change
    @Test
    public void testLineageCountAtTimeFunction() throws Exception {

        // Assemble BEASTObjects needed by QuasiSpeciesTree
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((((((((((t0 : 1, t1 : 1) : 1, t2 : 2) : 0.5, t3 : 1) : 1, t4 : 2) : 1, t5 : 3) : 1, t6 : 4) : 1, t7 : 5) : 0.5, t8 : 2) : 1, t9 : 3) : 1, t10 : 4 ) : 1, t11 : 10 ) : 1;",
                new String[] {"A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "G"});

        // calculate the tree likelihood with QS algorithm
        QuasiSpeciesBirthDeathSkylineModel bdsqs = this.getQSBDSKYmodel(tree, new RealParameter("11.0"),
                false, new RealParameter("2.0"), new RealParameter("1.0"),
                new RealParameter("0.5"));

        assertEquals(bdsqs.lineageCountAtTime(4.25,tree),5);
        assertEquals(bdsqs.lineageCountAtTime(1.25,tree),3);

    }

    /**
     *
     * Tree prior calculation P(tree|tree parameters) testing - number of full trees from qs tree
     *
     */

    @Test
    public void testCounterOfFullTreesRepresentedByQsTree1() throws Exception {

        // Assemble BEASTObjects needed by QuasiSpeciesTree
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((((((((((t0 : 1, t1 : 1) : 1, t2 : 2) : 0.5, t3 : 1) : 1, t4 : 2) : 1, t5 : 3) : 1, t6 : 4) : 1, t7 : 5) : 0.5, t8 : 2) : 1, t9 : 3) : 1, t10 : 4 ) : 1, t11 : 10 ) : 1;",
                new String[] {"A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "G"});
        tree.setEverythingDirty(true);

        // calculate the tree likelihood with QS algorithm
        QuasiSpeciesBirthDeathSkylineModel bdsqs = this.getQSBDSKYmodel(tree, new RealParameter("11.0"),
                false, new RealParameter("2.0"), new RealParameter("1.0"),
                new RealParameter("0.5"));

        // calculate the basis of the tree likelihood with bdsky
        Tree treeNormal = new TreeParser("(((((((((((t0 : 1, t1 : 1) : 1, t2 : 2) : 0.5, t3 : 1) : 1, t4 : 2) : 1, t5 : 3) : 1, t6 : 4) : 1, t7 : 5) : 0.5, t8 : 2) : 1, t9 : 3) : 1, t10 : 4 ) : 1, t11 : 10 ) : 1;",false);
        BirthDeathSkylineModel bdssm = this.getBDSKYmodel(treeNormal, new RealParameter("11.0"),
                false, new RealParameter("2.0"), new RealParameter("1.0"),
                new RealParameter("0.5"));

        double qsbdsky = bdsqs.calculateTreeLogLikelihood(tree);
        double bdsky = bdssm.calculateTreeLogLikelihood(treeNormal)
                // no qs to be attached to for G (node t11) so no gamma contribution here
                // before height 5 there are new 5 attachment points, and one lineage existing from height A
                + Math.log(6*5*5*4*4*3*3*2*2*1)
                // at the same time, at height 5 there is in total 6 + 1 (G) lineages
//                    - Math.log(7*6*6*5*5*4*4*3*3*2)
//                - Math.log(4*3*3*2)
                // before height 1.5 there are 3 surviving A lineages from previous interval, plus
                // 4 new attachment points
                + Math.log(7*6*6*5*5*4*4*3)
                // at the same time, at height 1.5 there is in total 7 + 1 (G) lineages
//                    - Math.log(8*7*7*6*6*5*5*4)
//                - Math.log(6*5*5*4*4*3*3*2)
                // before height 0 there are 2 surviving A lineages from previous interval, plus
                // 1 new attachment points
                + Math.log(3*2)
                // at the same time, at height 0 there is in total 3 + 1 (G) lineages
//                    - Math.log(4*3)
//                - Math.log(3*2)
                ;

        assertEquals(qsbdsky,bdsky,1e-10);

    }

    /**
     *
     * Tree prior calculation P(tree|tree parameters) testing - number of full trees from qs tree
     *
     */

    @Test
    public void testCounterOfFullTreesRepresentedByQsTree2() throws Exception {

        // Assemble BEASTObjects needed by QuasiSpeciesTree
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((((((((((t0 : 1, t1 : 1) : 1, t2 : 2) : 0.5, t3 : 1) : 1, t4 : 2) : 1, t5 : 3) : 1, t6 : 4) : 1, t7 : 5) : 0.5, t8 : 2) : 1, t9 : 3) : 1, t10 : 4 ) : 1, t11 : 10 ) : 1;",
                new String[] {"A", "A", "G", "A", "A", "A", "A", "A", "A", "A", "A", "A"});
        tree.setEverythingDirty(true);

        // calculate the tree likelihood with QS algorithm
        QuasiSpeciesBirthDeathSkylineModel bdsqs = this.getQSBDSKYmodel(tree, new RealParameter("11.0"),
                false, new RealParameter("2.0"), new RealParameter("1.0"),
                new RealParameter("0.5"));

        // calculate the basis of the tree likelihood with bdsky
        Tree treeNormal = new TreeParser("(((((((((((t0 : 1, t1 : 1) : 1, t2 : 2) : 0.5, t3 : 1) : 1, t4 : 2) : 1, t5 : 3) : 1, t6 : 4) : 1, t7 : 5) : 0.5, t8 : 2) : 1, t9 : 3) : 1, t10 : 4 ) : 1, t11 : 10 ) : 1;",false);
        BirthDeathSkylineModel bdssm = this.getBDSKYmodel(treeNormal, new RealParameter("11.0"),
                false, new RealParameter("2.0"), new RealParameter("1.0"),
                new RealParameter("0.5"));

        double qsbdsky = bdsqs.calculateTreeLogLikelihood(tree);
        double bdsky = bdssm.calculateTreeLogLikelihood(treeNormal)
                // before height 5 there are new 6 attachment points, and one lineage existing from height A
                + Math.log(7*6*6*5*5*4*4*3*3*2*2*1)
                // at the same time, at height 5 there is in total 7 (A) lineages
//                    - Math.log(7*6*6*5*5*4*4*3*3*2*2*1)
                // before height 2 there are 4 surviving A lineages from previous interval, plus
                // 3 new attachment points
                + Math.log(7*6*6*5*5*4)
                // at the same time, at height 2 there is in total 7 (A) lineages
//                    - Math.log(7*6*6*5*5*4)
                // at height 2, G lineage is born through true internal node
                // there are thus 7 qs A lineages to be attached to for G (node t2) so there is gamma contribution
                + Math.log(7)
                // before height 1.5 there are 2 surviving A lineages from previous interval, plus
                // 0 new attachment points
                //+ Math.log(0)
                // at the same time, at height 1.5 there is in total 2 + 1 (G) lineages
                    //- Math.log(0)
                // before height 0 there are 2 surviving A lineages from previous interval, plus
                // 1 new attachment points
                + Math.log(3*2)
                // at the same time, at height 0 there is in total 3 + 1 (G) lineages
//                    - Math.log(4*3)
                ;

        assertEquals(qsbdsky,bdsky,1e-10);

    }


    /**
     *
     * Tree prior calculation P(tree|tree parameters) testing - number of full trees from qs tree
     *
     */

    @Test
    public void testCounterOfFullTreesRepresentedByQsTree3() throws Exception {

        // Assemble BEASTObjects needed by QuasiSpeciesTree
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((((((((((t0 : 1, t1 : 1) : 1, t2 : 2) : 0.5, t3 : 1) : 1, t4 : 2) : 1, t5 : 3) : 1, t6 : 4) : 1, t7 : 5) : 0.5, t8 : 2) : 1, t9 : 3) : 1, t10 : 4 ) : 1, t11 : 10 ) : 1;",
                new String[] {"A", "A", "G", "A", "T", "A", "A", "A", "A", "A", "A", "A"});
        tree.setEverythingDirty(true);

        // calculate the tree likelihood with QS algorithm
        QuasiSpeciesBirthDeathSkylineModel bdsqs = this.getQSBDSKYmodel(tree, new RealParameter("11.0"),
                false, new RealParameter("2.0"), new RealParameter("1.0"),
                new RealParameter("0.5"));

        // calculate the basis of the tree likelihood with bdsky
        Tree treeNormal = new TreeParser("(((((((((((t0 : 1, t1 : 1) : 1, t2 : 2) : 0.5, t3 : 1) : 1, t4 : 2) : 1, t5 : 3) : 1, t6 : 4) : 1, t7 : 5) : 0.5, t8 : 2) : 1, t9 : 3) : 1, t10 : 4 ) : 1, t11 : 10 ) : 1;",false);
        BirthDeathSkylineModel bdssm = this.getBDSKYmodel(treeNormal, new RealParameter("11.0"),
                false, new RealParameter("2.0"), new RealParameter("1.0"),
                new RealParameter("0.5"));

        double qsbdsky = bdsqs.calculateTreeLogLikelihood(tree);
        double bdsky = bdssm.calculateTreeLogLikelihood(treeNormal)
                // before height 5 there are new 6 attachment points, and one lineage existing from height A
                + Math.log(7*6*6*5*5*4*4*3*3*2*2*1)
                // at the same time, at height 5 there is in total 7 (A) lineages
//                    - Math.log(7*6*6*5*5*4*4*3*3*2*2*1)
                // before height 3.5 there are 4 surviving A lineages from previous interval, plus
                // 1 new attachment points
                + Math.log(5*4)
                // at the same time, at height 3.5 there are in total 5 (A) lineages
//                    - Math.log(5*4)
                // at height 3.5, T lineage is born through true internal node
                // there are thus 5 qs A lineages to be attached to for T (node t4) so there is gamma contribution
                + Math.log(5)
                // before height 2 there are 5 surviving A lineages from previous interval, plus
                // 1 new attachment points
                + Math.log(6*5)
                // at the same time, at height 2 there is in total 6 (A) lineages and one T lineage
//                    - Math.log(7*6)
                // at height 2, G lineage is born through true internal node
                // there are thus 6 qs A lineages to be attached to for G (node t2) so there is gamma contribution
                + Math.log(6)
                // before height 1.5 there are 2 surviving A lineages from previous interval, plus
                // 0 new attachment points
                //+ Math.log(0)
                // at the same time, at height 1.5 there is in total 2 + 1 (G) lineages
                //- Math.log(0)
                // before height 0 there are 2 surviving A lineages from previous interval, plus
                // 1 new attachment points
                + Math.log(3*2)
                // at the same time, at height 0 there is in total 3 + 1 (G) lineages
//                    - Math.log(4*3)
                ;

        assertEquals(qsbdsky,bdsky,1e-10);

    }

    /**
     *
     * Tree prior calculation P(tree|tree parameters) testing - number of full trees from qs tree
     *
     */

    @Test
    public void testCounterOfFullTreesRepresentedByQsTreeSAP() throws Exception {

        // Assemble BEASTObjects needed by QuasiSpeciesTree
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((((((((((t0 : 1, t1 : 1) : 1, t2 : 2) : 0.5, t3 : 2.5) : 1, t4 : 3.5) : 1, t5 : 4.5) : 1, t6 : 5.5) : 1, t7 : 6.5) : 0.5, t8 : 7) : 1, t9 : 8) : 1, t10 : 9 ) : 1, t11 : 10 ) : 1;",
                new String[] {"A", "A", "G", "A", "T", "A", "A", "A", "A", "A", "A", "A"});
        tree.setEverythingDirty(true);

        // calculate the tree likelihood with QS algorithm
        QuasiSpeciesBirthDeathSkylineModel bdsqs = this.getQSBDSKYmodel(tree, new RealParameter("11.0"),
                false, new RealParameter("2.0"), new RealParameter("1.0"),
                new RealParameter("0.5"));

        // calculate the basis of the tree likelihood with bdsky
        Tree treeNormal = new TreeParser("(((((((((((t0 : 1, t1 : 1) : 1, t2 : 2) : 0.5, t3 : 2.5) : 1, t4 : 3.5) : 1, t5 : 4.5) : 1, t6 : 5.5) : 1, t7 : 6.5) : 0.5, t8 : 7) : 1, t9 : 8) : 1, t10 : 9 ) : 1, t11 : 10 ) : 1;",false);
        BirthDeathSkylineModel bdssm = this.getBDSKYmodel(treeNormal, new RealParameter("11.0"),
                false, new RealParameter("2.0"), new RealParameter("1.0"),
                new RealParameter("0.5"));

        double qsbdsky = bdsqs.calculateTreeLogLikelihood(tree);
        double bdsky = bdssm.calculateTreeLogLikelihood(treeNormal)
                // before height 3.5 there are new 7 attachment points, and one lineage existing from height A
//                + Math.log(8*7*7*6*6*5*5*4*4*3*3*2*2*1)
                // at the same time, at height 3.5 there are in total 8 (A) lineages
//                    - Math.log(8*7*7*6*6*5*5*4*4*3*3*2*2*1)
                // at height 3.5, T lineage is born through true internal node
                // there are thus 8 qs A lineages to be attached to for T (node t4) so there is gamma contribution
                + Math.log(8)
                // before height 2 there are 8 surviving A lineages from previous interval, plus
                // 1 new attachment points
//                + Math.log(9*8)
                // at the same time, at height 2 there is in total 9 (A) lineages and one T lineage
//                    - Math.log(10*9)
                // at height 2, G lineage is born through true internal node
                // there are thus 9 qs A lineages to be attached to for G (node t2) so there is gamma contribution
                + Math.log(9)
                // before height 0 there are 9 surviving A lineages from previous interval, plus
                // 1 new attachment points
//                + Math.log(10*9)
                // at the same time, at height 0 there is in total 10 + 1 (G) + 1 (T) lineages
//                    - Math.log(12*11)
                ;

        assertEquals(qsbdsky,bdsky,1e-10);

    }
}
