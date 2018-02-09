package test.quasispeciestree;

import static org.junit.Assert.assertEquals;
import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.speciation.BirthDeathSkylineModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.HKY;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import org.junit.Test;
import quasispeciestree.distributions.BirthDeathSkylineQuasiSpeciesModel;
import quasispeciestree.tree.QuasiSpeciesTree;
import quasispeciestree.tree.QuasiSpeciesTreeFromNewick;
import quasispeciestree.tree.QuasiSpeciesTreeFromFullNewick;
import test.beast.BEASTTestCase;
import test.beast.evolution.alignment.UncertainAlignmentTest;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Veronika Boskova created on 15/07/16 finished on 0?/0?/2016
 */

@Description("Test the quasispeciestree birth-death model with a small tree example")
public class QuasiSpeciesBDModelTests {


    /**
     *
     * Tree prior calculation P(tree|tree parameters) testing - 4 haplo with duplicates
     *
     */

    @Test
    public void testLikelihoodCalculation1() throws Exception {

        // calculate the tree likelihood with QS algorithm
        BirthDeathSkylineQuasiSpeciesModel bdsqs =  new BirthDeathSkylineQuasiSpeciesModel();

        // Assemble BEASTObjects needed by QuasiSpeciesTree
        int Nleaves = 4;
        StringBuilder traitSB = new StringBuilder();
        List<Sequence> seqList = new ArrayList<Sequence>();

        String[] data = {"A","C","G","T"};
        for (int i=0; i<Nleaves; i++) {
            String taxonID = "t" + i;
            seqList.add(new Sequence(taxonID, data[i]));

            if (i>0)
                traitSB.append(",");
            traitSB.append(taxonID).append("=").append(i+1);
        }

        Alignment alignment = new Alignment(seqList, "nucleotide");
        TaxonSet taxonSet = new TaxonSet(alignment);
        TraitSet haploCounts = new TraitSet();

        haploCounts.initByName(
                "traitname", "qscounts",
                "taxa", taxonSet,
                "value", traitSB.toString());

        QuasiSpeciesTree tree = new QuasiSpeciesTreeFromNewick();
        // for QuasiSpeciesTreeFromNewick() class
        tree.setInputValue("newick","((t3 : 1.5, t0 : 0.5) : 1 , (t1 : 2, t2 : 1) : 3);");
        tree.setInputValue("adjustTipHeights","false");
        // for QuasiSpeciesTree() class
        //tree.setInputValue("origin",new RealParameter("6."));
        tree.setInputValue("haplotypeCounts",haploCounts);
        tree.setInputValue("data",alignment);
        tree.initAndValidate();

        // for BirthDeathSkylineQuasiSpeciesModel() class
        bdsqs.setInputValue("tree", tree);
        bdsqs.setInputValue("origin", new RealParameter("6."));
        bdsqs.setInputValue("conditionOnSurvival", false);

        // test without rate change
        bdsqs.setInputValue("birthRate", new RealParameter("2."));
        bdsqs.setInputValue("deathRate", new RealParameter("1."));
        bdsqs.setInputValue("samplingRate", new RealParameter("0.5"));

        bdsqs.initAndValidate();
        bdsqs.printTempResults = true;


        // calculate the basis of the tree likelihood with bdsky
        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", false);

        // test without rate change
        bdssm.setInputValue("birthRate", new RealParameter("2."));
        bdssm.setInputValue("deathRate", new RealParameter("1."));
        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));

        bdssm.initAndValidate();
        bdssm.printTempResults = true;

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
        double nrFullTrees = bdsqs.logNumberOfTrees(4,3);
        nrFullTrees += bdsqs.logNumberOfTrees(3,2);
        nrFullTrees += bdsqs.logNumberOfTrees(2,1);
        nrFullTrees += bdsqs.logNumberOfTrees(1,0);
        assertEquals(bdsky+t1_1+t2_1+t2_2+t3_1+t3_2+t3_3+t1_tips+t2_tips+t3_tips+nrFullTrees, bdsqs.calculateTreeLogLikelihood(tree), 1e-5);
    }

    /**
     *
     * Tree prior calculation P(tree|tree parameters) testing - 1 haplo 100 duplicates -- TODO
     *
     */

    @Test
    public void testLikelihoodCalculation2() throws Exception {

        // calculate the tree likelihood with QS algorithm
        BirthDeathSkylineQuasiSpeciesModel bdsqs =  new BirthDeathSkylineQuasiSpeciesModel();

        // Assemble BEASTObjects needed by QuasiSpeciesTree
        int Nleaves = 4;
        StringBuilder traitSB = new StringBuilder();
        List<Sequence> seqList = new ArrayList<Sequence>();

        String[] data = {"A","C","G","T"};
        for (int i=0; i<Nleaves; i++) {
            String taxonID = "t" + i;
            // this was originally for not data[i], but all tips with "?"
            seqList.add(new Sequence(taxonID, data[i]));

            if (i>0)
                traitSB.append(",");
            traitSB.append(taxonID).append("=").append(i+1);
        }

        Alignment alignment = new Alignment(seqList, "nucleotide");
        TaxonSet taxonSet = new TaxonSet(alignment);
        TraitSet haploCounts = new TraitSet();

        haploCounts.initByName(
                "traitname", "qscounts",
                "taxa", taxonSet,
                "value", traitSB.toString());
//        Tree treeT = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
//        QuasiSpeciesTree tree = (QuasiSpeciesTree) treeT;

        QuasiSpeciesTree tree = new QuasiSpeciesTreeFromNewick();
        // for QuasiSpeciesTreeFromNewick() class
        tree.setInputValue("newick","((t3 : 1.5, t0 : 0.5) : 1 , (t1 : 2, t2 : 1) : 3);");
        tree.setInputValue("adjustTipHeights","false");
        // for QuasiSpeciesTree() class
        tree.setInputValue("haplotypeCounts",haploCounts);
        tree.setInputValue("data",alignment);
        tree.initAndValidate();

        // for BirthDeathSkylineQuasiSpeciesModel() class
        bdsqs.setInputValue("tree", tree);
        bdsqs.setInputValue("origin", new RealParameter("6."));
        bdsqs.setInputValue("conditionOnSurvival", false);

        // test without rate change
        bdsqs.setInputValue("birthRate", new RealParameter("2."));
        bdsqs.setInputValue("deathRate", new RealParameter("1."));
        bdsqs.setInputValue("samplingRate", new RealParameter("0.5"));

        bdsqs.initAndValidate();
        bdsqs.printTempResults = true;


        // calculate the basis of the tree likelihood with bdsky
        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        bdssm.setInputValue("tree", tree);
        bdssm.setInputValue("origin", new RealParameter("6."));
        bdssm.setInputValue("conditionOnSurvival", false);

        // test without rate change
        bdssm.setInputValue("birthRate", new RealParameter("2."));
        bdssm.setInputValue("deathRate", new RealParameter("1."));
        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));

        bdssm.initAndValidate();
        bdssm.printTempResults = true;

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
        double nrFullTrees = bdsqs.logNumberOfTrees(4,3);
        nrFullTrees += bdsqs.logNumberOfTrees(3,2);
        nrFullTrees += bdsqs.logNumberOfTrees(2,1);
        nrFullTrees += bdsqs.logNumberOfTrees(1,0);
        assertEquals(bdsky+t1_1+t2_1+t2_2+t3_1+t3_2+t3_3+t1_tips+t2_tips+t3_tips+nrFullTrees, bdsqs.calculateTreeLogLikelihood(tree), 1e-5);
    }

    /**
     *
     * Tree prior calculation P(tree|tree parameters) testing - lineageCountAtTime
     *
     */

    @Test
    public void testLineageCountAtTimeFunction() throws Exception {
        // calculate the tree likelihood with QS algorithm
        BirthDeathSkylineQuasiSpeciesModel bdsqs = new BirthDeathSkylineQuasiSpeciesModel();

        // Assemble BEASTObjects needed by QuasiSpeciesTree
        int Nleaves = 12;
        StringBuilder traitSB = new StringBuilder();
        List<Sequence> seqList = new ArrayList<Sequence>();

        String[] data = {"A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "G"};

        for (int i = 0; i < Nleaves; i++) {
            String taxonID = "t" + i;
            seqList.add(new Sequence(taxonID, data[i]));

            if (i > 0)
                traitSB.append(",");
            traitSB.append(taxonID).append("=").append(1);
        }

        Alignment alignment = new Alignment(seqList, "nucleotide");
        TaxonSet taxonSet = new TaxonSet(alignment);
        TraitSet haploCounts = new TraitSet();

        haploCounts.initByName("traitname", "qscounts", "taxa", taxonSet, "value", traitSB.toString());
        //        Tree treeT = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        //        QuasiSpeciesTree tree = (QuasiSpeciesTree) treeT;

        QuasiSpeciesTree tree = new QuasiSpeciesTreeFromFullNewick();
        // for QuasiSpeciesTreeFromNewick() class
        tree.setInputValue("newick", "(((((((((((t0 : 1, t1 : 1) : 1, t2 : 2) : 0.5, t3 : 1) : 1, t4 : 2) : 1, t5 : 3) : 1, t6 : 4) : 1, t7 : 5) : 0.5, t8 : 2) : 1, t9 : 3) : 1, t10 : 4 ) : 1, t11 : 10 ) : 1;");
        tree.setInputValue("adjustTipHeights", "false");
        // for QuasiSpeciesTree() class
        //tree.setInputValue("haplotypeCounts", haploCounts);
        tree.setInputValue("data", alignment);
        tree.initAndValidate();

        // for BirthDeathSkylineQuasiSpeciesModel() class
        bdsqs.setInputValue("tree", tree);
        bdsqs.setInputValue("origin", new RealParameter("11."));
        bdsqs.setInputValue("conditionOnSurvival", false);

        // test without rate change
        bdsqs.setInputValue("birthRate", new RealParameter("2."));
        bdsqs.setInputValue("deathRate", new RealParameter("1."));
        bdsqs.setInputValue("samplingRate", new RealParameter("0.5"));

        bdsqs.initAndValidate();

        assertEquals(bdsqs.lineageCountAtTime(4.25,tree),5);
        assertEquals(bdsqs.lineageCountAtTime(1.25,tree),3);

    }

    /**
     *
     * Tree prior calculation P(tree|tree parameters) testing - number of full trees from qs tree
     *
     */

    @Test
    public void testCounterOfFullTreesRepresentedByQsTree() throws Exception {
        // calculate the tree likelihood with QS algorithm
        BirthDeathSkylineQuasiSpeciesModel bdsqs = new BirthDeathSkylineQuasiSpeciesModel();

        // Assemble BEASTObjects needed by QuasiSpeciesTree
        int Nleaves = 12;
        StringBuilder traitSB = new StringBuilder();
        List<Sequence> seqList = new ArrayList<Sequence>();

        String[] data = {"A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "G"};

        for (int i = 0; i < Nleaves; i++) {
            String taxonID = "t" + i;
            seqList.add(new Sequence(taxonID, data[i]));

            if (i > 0)
                traitSB.append(",");
            traitSB.append(taxonID).append("=").append(1);
        }

        Alignment alignment = new Alignment(seqList, "nucleotide");
        TaxonSet taxonSet = new TaxonSet(alignment);
        TraitSet haploCounts = new TraitSet();

        haploCounts.initByName("traitname", "qscounts", "taxa", taxonSet, "value", traitSB.toString());
        //        Tree treeT = new TreeParser("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);",false);
        //        QuasiSpeciesTree tree = (QuasiSpeciesTree) treeT;

        QuasiSpeciesTree tree = new QuasiSpeciesTreeFromFullNewick();
        // for QuasiSpeciesTreeFromNewick() class
        tree.setInputValue("newick", "(((((((((((t0 : 1, t1 : 1) : 1, t2 : 2) : 0.5, t3 : 1) : 1, t4 : 2) : 1, t5 : 3) : 1, t6 : 4) : 1, t7 : 5) : 0.5, t8 : 2) : 1, t9 : 3) : 1, t10 : 4 ) : 1, t11 : 10 ) : 1;");
        tree.setInputValue("adjustTipHeights", "false");
        // for QuasiSpeciesTree() class
        //tree.setInputValue("haplotypeCounts", haploCounts);
        tree.setInputValue("data", alignment);
        tree.initAndValidate();

        // for BirthDeathSkylineQuasiSpeciesModel() class
        bdsqs.setInputValue("tree", tree);
        bdsqs.setInputValue("origin", new RealParameter("11."));
        bdsqs.setInputValue("conditionOnSurvival", false);

        // test without rate change
        bdsqs.setInputValue("birthRate", new RealParameter("2."));
        bdsqs.setInputValue("deathRate", new RealParameter("1."));
        bdsqs.setInputValue("samplingRate", new RealParameter("0.5"));

        bdsqs.initAndValidate();

        Tree treeNormal = new TreeParser("(((((((((((t0 : 1, t1 : 1) : 1, t2 : 2) : 0.5, t3 : 1) : 1, t4 : 2) : 1, t5 : 3) : 1, t6 : 4) : 1, t7 : 5) : 0.5, t8 : 2) : 1, t9 : 3) : 1, t10 : 4 ) : 1, t11 : 10 ) : 1;",false);

        // calculate the basis of the tree likelihood with bdsky
        BirthDeathSkylineModel bdssm =  new BirthDeathSkylineModel();

        bdssm.setInputValue("tree", treeNormal);
        bdssm.setInputValue("origin", new RealParameter("11."));
        bdssm.setInputValue("conditionOnSurvival", false);

        // test without rate change
        bdssm.setInputValue("birthRate", new RealParameter("2."));
        bdssm.setInputValue("deathRate", new RealParameter("1."));
        bdssm.setInputValue("samplingRate", new RealParameter("0.5"));

        bdssm.initAndValidate();

        bdsqs.printTempResults = true;
        double qsbdsky = bdsqs.calculateTreeLogLikelihood(tree);
        bdssm.printTempResults = true;
        double bdsky = bdssm.calculateTreeLogLikelihood(treeNormal)
                +bdsqs.logNumberOfTrees(3,1)
                +bdsqs.logNumberOfTrees(7,4)
                +bdsqs.logNumberOfTrees(6,5);

        assertEquals(qsbdsky,bdsky,1e-10);

    }
}
