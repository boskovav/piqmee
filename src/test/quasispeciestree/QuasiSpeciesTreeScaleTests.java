package test.quasispeciestree;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import org.junit.Test;
import beast.core.Description;
import quasispeciestree.operators.QuasiSpeciesTreeScale;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;
import test.quasispeciestree.DeterministicRandomGenerator;
import test.quasispeciestree.QuasiSpeciesTestCase;


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

    private double getTree(QuasiSpeciesTree tree, double scale) {
        QuasiSpeciesTreeScale treeScale = new QuasiSpeciesTreeScale(new DeterministicRandomGenerator(1));
        treeScale.setInputValue( "quasiSpeciesTree", tree);
        treeScale.setInputValue( "scaleFactor", scale);
        treeScale.setInputValue( "rootOnly", true);
        treeScale.initAndValidate();
        double result = treeScale.proposal();
        return result;
    }



    // testing scaling of root down



    @Test
    public void testQuasiSpeciesTreeScale0() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 0.75 and 0.5
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("((t2 : 0.5, t3 : 0.5) : 0.5, (t0 : 0.25, t1 : 0.25) : 0.25);", new String[] {"A", "A", "C", "C"});

        assertTrue(getTree(tree,(1-0.25)/1.0*1.0) ==Double.NEGATIVE_INFINITY);

    }

    @Test
    public void testQuasiSpeciesTreeScaleDown1() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 0.75 and 0.5
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("((t2 : 0.5, t3 : 0.5) : 0.5, (t0 : 0.25, t1 : 0.25) : 0.25);", new String[] {"A", "A", "C", "C"});
        double result = getTree(tree,(1-0.25)/1.0*1.1);

        assertTrue(((QuasiSpeciesNode)tree.getRoot()).getHaploAboveName() == -1);
    }

    @Test
    public void testQuasiSpeciesTreeScaleDown2() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 0.75 and 0.5
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("((t2 : 0.5, t3 : 0.5) : 0.5, (t0 : 0.25, t1 : 0.25) : 0.25);", new String[] {"A", "A", "C", "C"});

        double result = getTree(tree,(1-0.25)/1.0*0.9);

        assertTrue(((QuasiSpeciesNode)tree.getRoot()).getHaploAboveName() != -1);

    }

    @Test
    public void testQuasiSpeciesTreeScaleDown3() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 0.75 and 0.5
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("((t2 : 0.5, t3 : 0.5) : 0.5, (t0 : 0.25, t1 : 0.25) : 0.25);", new String[] {"A", "A", "C", "C"});

        assertTrue(getTree(tree,(1-0.25)/1.0*0.5) == Double.NEGATIVE_INFINITY);


    }

    @Test
    public void testQuasiSpeciesTreeScaleDown4() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 1.25 and 0.5
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((t2 : 0.5, t3 : 0.5) : 0.5, t0 : 0.5) : 0.25, t1 : 0.75);", new String[] {"A", "A", "C", "C"});
        int haploAbove = ((QuasiSpeciesNode) tree.getRoot()).getHaploAboveName();
        double result = getTree(tree,(1-0.5)/1.0*1.0001);

        assertTrue(((QuasiSpeciesNode)tree.getRoot()).getHaploAboveName() == haploAbove);

    }

    @Test
    public void testQuasiSpeciesTreeScaleDown5() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 1.25 and 0.5
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((t2 : 0.5, t3 : 0.5) : 0.5, t0 : 0.5) : 0.25, t1 : 0.75);", new String[] {"A", "A", "C", "C"});
        double result = getTree(tree,(1-0.5)/1.0*0.9);

        assertTrue(result == Double.NEGATIVE_INFINITY);

    }



    // testing scaling of root up



    @Test
    public void testQuasiSpeciesTreeScaleUp1() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 0.75 and 0.5
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("((t2 : 0.5, t3 : 0.5) : 0.5, (t0 : 0.25, t1 : 0.25) : 0.25);", new String[] {"A", "A", "C", "C"});
        double result = getTree(tree,(6.0-1.0)/0.5);

        assertTrue(((QuasiSpeciesNode)tree.getRoot()).getHaploAboveName() == -1);
    }

    @Test
    public void testQuasiSpeciesTreeScaleUp2() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 1.25 and 0.5
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((t2 : 0.5, t3 : 0.5) : 0.5, t0 : 0.5) : 0.25, t1 : 0.75);", new String[] {"A", "A", "C", "C"});
        int haploAbove = ((QuasiSpeciesNode) tree.getRoot()).getHaploAboveName();
        double result = getTree(tree,(1.25-1.0)/1.0*0.9);

        assertTrue(((QuasiSpeciesNode)tree.getRoot()).getHaploAboveName() == haploAbove);

    }

    @Test
    public void testQuasiSpeciesTreeScaleUp3() throws Exception {
        // Assemble BEASTObjects needed by QuasiSpeciesTree
        // root at 1.0, haplo attach at 1.25 and 0.5
        QuasiSpeciesTree tree = QuasiSpeciesTestCase.setTreeFromFullNewick("(((t2 : 0.5, t3 : 0.5) : 0.5, t0 : 0.5) : 0.25, t1 : 0.75);", new String[] {"A", "A", "C", "C"});
        int haploAbove = ((QuasiSpeciesNode) tree.getRoot()).getHaploAboveName();
        double result = getTree(tree,1.25/1.0*1.1);

        assertTrue(((QuasiSpeciesNode)tree.getRoot()).getHaploAboveName() == -1);

    }

}