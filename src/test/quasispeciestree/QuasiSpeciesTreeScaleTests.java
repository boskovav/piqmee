package test.quasispeciestree;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import org.junit.Test;
import beast.core.Description;
import quasispeciestree.operators.QuasiSpeciesTreeScale;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;
import quasispeciestree.tree.QuasiSpeciesTreeFromFullNewick;

import java.util.ArrayList;
import java.util.List;


/**
 * @author Veronika Boskova created on 23/10/17
 */

@Description("Test the quasispeciestree tree scale operator class functions")
public class QuasiSpeciesTreeScaleTests {

    /**
     *
     * Tree scale operator: root scale testing - 2 haplo with duplicates
     *
     */

    private QuasiSpeciesTree setTree(String inputTree) {
        int Nleaves = 4;
        StringBuilder traitSB = new StringBuilder();
        List<Sequence> seqList = new ArrayList<Sequence>();

        String[] data = {"A", "A", "C", "C"};
        for (int i = 0; i < Nleaves; i++) {
            String taxonID = "t" + (i);
            seqList.add(new Sequence(taxonID, data[i]));

            if (i > 0)
                traitSB.append(",");
            traitSB.append(taxonID).append("=").append(i+1);
        }

        Alignment alignment = new Alignment(seqList, "nucleotide");
        TaxonSet taxonSet = new TaxonSet(alignment);
        TraitSet haploCounts = new TraitSet();

        QuasiSpeciesTree tree = new QuasiSpeciesTreeFromFullNewick();

        tree.setInputValue("newick", inputTree);
        tree.setInputValue("adjustTipHeights", "false");
        // for QuasiSpeciesTree() class
        tree.setInputValue("origin", new Double("6."));
        tree.setInputValue("data", alignment);
        tree.initAndValidate();

        return tree;
    }

    private double getTree(QuasiSpeciesTree tree, double scale) {
        QuasiSpeciesTreeScale treeScale = new QuasiSpeciesTreeScale(new DeterministicRandomGenerator(1));
        treeScale.setInputValue( "quasiSpeciesTree", tree);
        treeScale.setInputValue( "scaleFactor", scale);
        treeScale.setInputValue( "rootOnly", true);
        treeScale.initAndValidate();
        double result = treeScale.proposal();
        return result;
    }

    @Test
    public void testQuasiSpeciesTreeScale0() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 0.75 and 0.5
        QuasiSpeciesTree tree = setTree("((t2 : 0.5, t3 : 0.5) : 0.5, (t0 : 0.25, t1 : 0.25) : 0.25);");

        assertTrue(getTree(tree,(1-0.25)/1.0*1.0) ==Double.NEGATIVE_INFINITY);

    }



    // testing scaling of root down



    @Test
    public void testQuasiSpeciesTreeScaleDown1() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 0.75 and 0.5
        QuasiSpeciesTree tree = setTree("((t2 : 0.5, t3 : 0.5) : 0.5, (t0 : 0.25, t1 : 0.25) : 0.25);");
        double result = getTree(tree,(1-0.25)/1.0*1.1);

        assertTrue(((QuasiSpeciesNode)tree.getRoot()).getHaploAboveName() == -1);
    }

    @Test
    public void testQuasiSpeciesTreeScaleDown2() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 0.75 and 0.5
        QuasiSpeciesTree tree = setTree("((t2 : 0.5, t3 : 0.5) : 0.5, (t0 : 0.25, t1 : 0.25) : 0.25);");

        double result = getTree(tree,(1-0.25)/1.0*0.9);

        assertTrue(((QuasiSpeciesNode)tree.getRoot()).getHaploAboveName() != -1);

    }

    @Test
    public void testQuasiSpeciesTreeScaleDown3() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 0.75 and 0.5
        QuasiSpeciesTree tree = setTree("((t2 : 0.5, t3 : 0.5) : 0.5, (t0 : 0.25, t1 : 0.25) : 0.25);");

        assertTrue(getTree(tree,(1-0.25)/1.0*0.5) == Double.NEGATIVE_INFINITY);


    }

    @Test
    public void testQuasiSpeciesTreeScaleDown4() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 1.25 and 0.5
        QuasiSpeciesTree tree = setTree("(((t2 : 0.5, t3 : 0.5) : 0.5, t0 : 0.5) : 0.25, t1 : 0.75);");
        int haploAbove = ((QuasiSpeciesNode) tree.getRoot()).getHaploAboveName();
        double result = getTree(tree,(1-0.5)/1.0*1.0001);

        assertTrue(((QuasiSpeciesNode)tree.getRoot()).getHaploAboveName() == haploAbove);

    }

    @Test
    public void testQuasiSpeciesTreeScaleDown5() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 1.25 and 0.5
        QuasiSpeciesTree tree = setTree("(((t2 : 0.5, t3 : 0.5) : 0.5, t0 : 0.5) : 0.25, t1 : 0.75);");
        double result = getTree(tree,(1-0.5)/1.0*0.9);

        assertTrue(result == Double.NEGATIVE_INFINITY);

    }



    // testing scaling of root up



    @Test
    public void testQuasiSpeciesTreeScaleUp1() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 0.75 and 0.5
        QuasiSpeciesTree tree = setTree("((t2 : 0.5, t3 : 0.5) : 0.5, (t0 : 0.25, t1 : 0.25) : 0.25);");
        double result = getTree(tree,(6.0-1.0)/0.5);

        assertTrue(((QuasiSpeciesNode)tree.getRoot()).getHaploAboveName() == -1);
    }

    @Test
    public void testQuasiSpeciesTreeScaleUp2() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 1.25 and 0.5
        QuasiSpeciesTree tree = setTree("(((t2 : 0.5, t3 : 0.5) : 0.5, t0 : 0.5) : 0.25, t1 : 0.75);");
        int haploAbove = ((QuasiSpeciesNode) tree.getRoot()).getHaploAboveName();
        double result = getTree(tree,(1.25-1.0)/1.0*0.9);

        assertTrue(((QuasiSpeciesNode)tree.getRoot()).getHaploAboveName() == haploAbove);

    }

    @Test
    public void testQuasiSpeciesTreeScaleUp3() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 1.25 and 0.5
        QuasiSpeciesTree tree = setTree("(((t2 : 0.5, t3 : 0.5) : 0.5, t0 : 0.5) : 0.25, t1 : 0.75);");
        int haploAbove = ((QuasiSpeciesNode) tree.getRoot()).getHaploAboveName();
        double result = getTree(tree,1.25/1.0*1.1);

        assertTrue(((QuasiSpeciesNode)tree.getRoot()).getHaploAboveName() == -1);

    }

}




//    int Nleaves = 4;
//    StringBuilder traitSB = new StringBuilder();
//    List<Sequence> seqList = new ArrayList<Sequence>();
//
//    String[] data = {"A", "C", "G", "T"};
//        for (int i = 0; i < Nleaves; i++) {
//        String taxonID = "t" + (i);
//        seqList.add(new Sequence(taxonID, data[i]));
//
//        if (i > 0)
//        traitSB.append(",");
//        traitSB.append(taxonID).append("=").append(i+1);
//        }
//
//        Alignment alignment = new Alignment(seqList, "nucleotide");
//        TaxonSet taxonSet = new TaxonSet(alignment);
//        TraitSet haploCounts = new TraitSet();
//
//        haploCounts.initByName("traitname", "qscounts", "taxa", taxonSet, "value", traitSB.toString());
//
//        QuasiSpeciesTree tree = new QuasiSpeciesTreeFromNewick();
//        // for QuasiSpeciesTreeFromNewick() class
//
//        //        QuasiSpeciesBuilder().setNewick(7).build()
//
//
//
//        tree.setInputValue("newick", "((t3 : 1.5, t0 : 0.5) : 1 , (t1 : 2, t2 : 1) : 3);");
//        tree.setInputValue("adjustTipHeights", "false");
//        // for QuasiSpeciesTree() class
//        tree.setInputValue("origin", new Double("6."));
//        tree.setInputValue("haplotypeCounts", haploCounts);
//        tree.setInputValue("data", alignment);
//        tree.initAndValidate();
//
//        QuasiSpeciesTreeScale treeScale = new QuasiSpeciesTreeScale();
//        treeScale.setInputValue( "scaleFactor", 1.5);
//        treeScale.initAndValidate();
//        treeScale.proposal();
//
