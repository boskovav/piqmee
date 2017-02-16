package quasispeciestree.operators;

import beast.core.Description;
import beast.util.Randomizer;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;

/**
 * @author Veronika Boskova created on 18/01/2017
 */
@Description("Randomly selects true internal tree node or a root) and move node height uniformly in interval " +
             "restricted by the nodes parent or an attachment time of a haplo going through it " +
             "and children or attachment time of a haplo going through.")
public class QuasiSpeciesUniform extends QuasiSpeciesTreeOperator {

    // empty constructor to facilitate construction by XML + initAndValidate
    public QuasiSpeciesUniform() {
    }

    public QuasiSpeciesUniform(QuasiSpeciesTree qsTree) {
        try {
            initByName(quasiSpeciesTreeInput.getName(), qsTree);
        } catch (Exception e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            throw new RuntimeException("Failed to construct Uniform Tree Operator.");
        }
    }

    @Override
    public void initAndValidate() {
    }

    /**
     * change the parameter and return the hastings ratio.
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        final QuasiSpeciesTree qsTree = quasiSpeciesTreeInput.get(this);

        // randomly select internal node
        final int nodeCount = qsTree.getNodeCount();

        // Abort if no non-root internal nodes
//        if (qsTree.getInternalNodeCount()==1)
//            return Double.NEGATIVE_INFINITY;

        QuasiSpeciesNode node;
        do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = (QuasiSpeciesNode) qsTree.getNode(nodeNr);
        //} while (node.isRoot() || node.isLeaf());
        } while (node.isLeaf());
        double oldValue = node.getHeight();
        QuasiSpeciesNode lowerBoundNode;
        double upperQS = Double.POSITIVE_INFINITY;
        double lowerQS = Double.NEGATIVE_INFINITY;
        double upperHard = Double.POSITIVE_INFINITY;
        double lowerHard = Double.NEGATIVE_INFINITY;
        double upper;
        double lower;

        // find the parent node (Hard upper)
        if (node.isRoot()){
            upperHard = originInput.get().getValue();
        }
        else {
            upperHard = node.getParent().getHeight();
        }
        // find the higher child node (Hard lower)
        if (node.getLeft().getHeight()>node.getRight().getHeight()) {
            lowerBoundNode = (QuasiSpeciesNode) node.getLeft();
        }
        else {
            lowerBoundNode = (QuasiSpeciesNode) node.getRight();
        }
        lowerHard = lowerBoundNode.getHeight();

        // find if there is QS passing through the node
        int haplo = node.getContinuingHaploName();
        if (haplo != -1 && qsTree.getAttachmentTimesList(haplo).length > 1){
            Double[] tempqstimes = qsTree.getAttachmentTimesList(haplo).clone();

            for (int i=1; i<=tempqstimes.length-1; i++){
                if (tempqstimes[i]>oldValue)
                    upperQS = tempqstimes[i];
                else {
                    lowerQS = tempqstimes[i];
                    break;
                }
            }

            // check which is the used lower/upper bounds
            if(lowerQS<lowerHard)
                lower = lowerHard;
            else lower = lowerQS;

            if (upperQS>upperHard)
                upper = upperHard;
            else upper = upperQS;
        }
        else{
            lower = lowerHard;
            upper = upperHard;
        }
        // there can still be (another) QS on the left or right from the node -- important for the LOWER bound
        int haploleft = ((QuasiSpeciesNode) node.getLeft()).getHaploAboveName();
        int haploright = ((QuasiSpeciesNode) node.getRight()).getHaploAboveName();
        double lowerQSleft = Double.NEGATIVE_INFINITY;
        double lowerQSright = Double.NEGATIVE_INFINITY;
        if (haploleft != -1 && qsTree.getAttachmentTimesList(haploleft).length > 1)
            lowerQSleft = qsTree.getAttachmentTimesList(haploleft).clone()[1];
        if (haploright != -1 && qsTree.getAttachmentTimesList(haploright).length > 1)
            lowerQSright = qsTree.getAttachmentTimesList(haploright).clone()[1];
        if (lowerQSleft<lowerQSright){
            if (lower<lowerQSright)
                lower = lowerQSright;
        }
        else {
            if (lower < lowerQSleft)
                lower = lowerQSleft;
        }

        // assign a new value to the node
        final double newValue = (Randomizer.nextDouble() * (upper - lower)) + lower;
        node.setHeight(newValue);

        return 0.0;
    }

}
