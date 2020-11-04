package test.piqmee.tree;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import org.junit.Test;
import beast.core.Description;
import piqmee.tree.QuasiSpeciesNode;


/**
 * @author Veronika Boskova created on 23/10/17
 */

@Description("Test the piqmee node class functions")
public class QuasiSpeciesNodeTests {

    /**
     *
     * Node scale
     *
     */

    private QuasiSpeciesNode createQuasiSpeciesNode() {
        QuasiSpeciesNode node = new QuasiSpeciesNode();
        node.setHeight(1.0);
        node.setAttachmentTimesList(new double[] {4.0, 4.0, 2.0, 1.5});
        node.setTipTimesList(new double[] {node.getHeight()});
        node.setTipTimesCountList(new int[] {node.getAttachmentTimesList().length});

        return node;
    }

    private QuasiSpeciesNode createQuasiSpeciesNodeSampledThroughTime() {
        QuasiSpeciesNode node = new QuasiSpeciesNode();
        node.setHeight(1.0);
        node.setAttachmentTimesList(new double[] {4.0, 4.0, 2.0, 1.5});
        node.setTipTimesList(new double[] {node.getHeight(),node.getHeight()*2});
        node.setTipTimesCountList(new int[] {node.getAttachmentTimesList().length});

        return node;
    }

    @Test
    public void testScaleOneNode() throws Exception {

        // define an example node
        QuasiSpeciesNode node =  createQuasiSpeciesNode();

        // scale the node and attachment time of it's duplicates
        node.scale(2);

        // check if values changed correctly
        assertArrayEquals(new double[] {8, 8, 4, 3}, node.getAttachmentTimesList(), 1e-100);

    }

    @Test(expected = IllegalArgumentException.class)
    public void testScaleOneNodeWithNegativeScaleFactor() throws Exception {

        // define an example node
        QuasiSpeciesNode node = createQuasiSpeciesNode();

        // scale the node and attachment time of it's duplicates
        node.scale(-2);

    }

    @Test(expected = IllegalArgumentException.class)
    public void testScaleOneNodeWithScaleFactorZero() throws Exception {

        // define an example node
        QuasiSpeciesNode node = createQuasiSpeciesNode();

        // scale the node and attachment time of it's duplicates
        node.scale(0);

    }


    @Test(expected = IllegalArgumentException.class)
    public void testScaleOneNodeAttachHaploAboveNode() throws Exception {

        // define an example node
        QuasiSpeciesNode node = createQuasiSpeciesNode();

        // scale the node and attachment time of it's duplicates
        node.scale(0.1);

    }

    @Test(expected = IllegalArgumentException.class)
    public void testScaleOneNodeAttachHaploAtNode() throws Exception {

        // define an example node
        QuasiSpeciesNode node = createQuasiSpeciesNode();

        // scale the node and attachment time of it's duplicates
        node.scale(1 / 1.5);

    }





//    // check if values changed correctly
//    assertTrue(1.0==node.getHeight());
//    assertTrue(-1==node.getHaploAboveName());
//    assertTrue(0==node.getContinuingHaploName());
//    assertTrue(1==node.getStartBranchCounts());
//    assertArrayEquals(new double[] {8, 8, 4, 3}, node.getAttachmentTimesList(), 1e-100);
//    assertArrayEquals(new double[] {1.0}, node.getTipTimesList(), 1e-100);
//    assertArrayEquals(new int[] {4}, node.getTipTimesCountList());
//    assertTrue(-1 == node.getParentHaplo());

}
