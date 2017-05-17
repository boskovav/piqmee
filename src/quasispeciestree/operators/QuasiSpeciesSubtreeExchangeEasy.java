package quasispeciestree.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import quasispeciestree.tree.QuasiSpeciesNode;

import java.util.ArrayList;

/**
 *  @author Veronika Boskova created on 26/04/17
 */
@Description("Subtree/branch exchange operations for quasispecies tree")
public class QuasiSpeciesSubtreeExchangeEasy extends QuasiSpeciesTreeOperator{

    final public Input<Boolean> isNarrowInput = new Input<>("isNarrow",
            "if true (default) a narrow exchange is performed, otherwise a wide exchange", true);

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        if (qsTree.getLeafNodeCount() < 4)
            throw new IllegalArgumentException("QuasiSpeciesSubtreeExchange operator cannot be " +
                    "used since there are only 4 tips in the tree! Remove this " +
                    "operator from your XML file.");
    }

    @Override
    public double proposal() {

        double logHastingsRatio = 0.0;

        // Select source and destination nodes:

        Node srcNode, srcNodeParent, destNode, destNodeParent;

        // Narrow exchange selection:

        if (isNarrowInput.get()) {
            // check that this operator can actually perform a move
            int cousinhaplobelowparentnode = 0;
            for (int i = 0; i < qsTree.getNodeCount(); i++){
                srcNode = qsTree.getNode(i);
                if (isaValidNodeForNarrow(srcNode)){
                    cousinhaplobelowparentnode = 1;
                    break;
                }
            }
            if (cousinhaplobelowparentnode == 0)
                return Double.NEGATIVE_INFINITY;
            // select a node for narrow exchange
            do {
                srcNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
            } while ( !isaValidNodeForNarrow(srcNode) );
            srcNodeParent = srcNode.getParent();
            destNodeParent = srcNodeParent.getParent();
            destNode = getOtherChild(destNodeParent, srcNodeParent);
        }

        // Wide exchange selection:

        else {
            // check that this operator can actually perform a move
            int haplobelowparentnode = 0;
            for (int i = 0; i < qsTree.getNodeCount(); i++){
                for (int j = i + 1; j < qsTree.getNodeCount(); j++) {
                    srcNode = qsTree.getNode(i);
                    destNode = qsTree.getNode(j);
                    if (isaValidNodePairForWide(srcNode, destNode)) {
                        haplobelowparentnode = 1;
                        break;
                    }
                }
                if (haplobelowparentnode==1)
                    break;
            }
            if (haplobelowparentnode == 0)
                return Double.NEGATIVE_INFINITY;
            // select a node for wide exchange
            do {
                srcNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
            } while ( !isaValidSrcNodeForWide(srcNode) );
            srcNodeParent = srcNode.getParent();

            do {
                destNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
            } while( !isaValidDestNodeForWide(srcNode, srcNodeParent, destNode) );
            destNodeParent = destNode.getParent();

            // Reject if substitution would result in negative branch lengths:
            int srcHaplo = ((QuasiSpeciesNode) srcNode).getHaploAboveName();
            int destHaplo = ((QuasiSpeciesNode) destNode).getHaploAboveName();
            if ( ( srcHaplo != -1 && ((QuasiSpeciesNode) qsTree.getNode(srcHaplo)).getAttachmentTimesList()[0] > destNodeParent.getHeight())
                    ||
                    (destHaplo != -1 && ((QuasiSpeciesNode) qsTree.getNode(destHaplo)).getAttachmentTimesList()[0] > srcNodeParent.getHeight()))
                return Double.NEGATIVE_INFINITY;
        }


        // Reject if substitution would result in negative branch lengths:
        if (destNode.getHeight() >= srcNodeParent.getHeight() || srcNode.getHeight() >= destNodeParent.getHeight())
            return Double.NEGATIVE_INFINITY;

        // Make changes to tree topology:
        replace(srcNodeParent, srcNode, destNode);
        replace(destNodeParent, destNode, srcNode);

        // Recalculate continuingHaplo and HaploAbove arrays
        recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) qsTree.getRoot());

        // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
        qsTree.countAndSetPossibleStartBranches();

        // RETURN log(HASTINGS RATIO)
        return logHastingsRatio;
    }

    /**
     * Function to decide if the source node is a valid node for narrow easy exchange
     *
     * @param srcNode
     */
    private boolean isaValidNodeForNarrow(Node srcNode) {

        Node srcNodeParent = srcNode.getParent();

        if (srcNode.isRoot() || srcNodeParent.isRoot())
            return false;

        int srcParentHaplo = ((QuasiSpeciesNode) srcNodeParent).getHaploAboveName();
        Node destNode = getOtherChild(srcNodeParent.getParent(), srcNodeParent);
        int destHaplo = ((QuasiSpeciesNode) destNode).getHaploAboveName();

        if ( (srcParentHaplo == -1
                && (destHaplo == -1 || ((QuasiSpeciesNode) qsTree.getNode(destHaplo)).getAttachmentTimesList()[0] < srcNodeParent.getHeight()) )
            ||
             (srcParentHaplo == ((QuasiSpeciesNode)getOtherChild(srcNodeParent, srcNode)).getContinuingHaploName()
                && (destHaplo == ((QuasiSpeciesNode) destNode).getContinuingHaploName()
                && (destHaplo == -1 || ((QuasiSpeciesNode) qsTree.getNode(destHaplo)).getAttachmentTimesList()[0] < srcNodeParent.getHeight())) )
            )
            return true;
        else
            return false;
    }

    /**
     * Function to decide if the source node is a valid node for wide easy exchange
     *
     * @param srcNode
     */
    private boolean isaValidSrcNodeForWide(Node srcNode) {
        if (srcNode.isRoot())
            return false;

        int srcHaplo = ((QuasiSpeciesNode) srcNode).getHaploAboveName();
        if (srcHaplo != ((QuasiSpeciesNode) srcNode).getContinuingHaploName())
            return false;

        return true;
    }

    /**
     * Function to decide if the source node is a valid node for wide easy exchange
     *
     * @param srcNode
     * @param srcNodeParent
     * @param destNode
     */
    private boolean isaValidDestNodeForWide(Node srcNode, Node srcNodeParent, Node destNode) {
        if ( destNode == srcNode || destNode.isRoot() || destNode.getParent() == srcNodeParent
                || srcNodeParent == destNode || destNode.getParent() == srcNode)
            return false;

        int destHaplo = ((QuasiSpeciesNode) destNode).getHaploAboveName();
        if (destHaplo != ((QuasiSpeciesNode) destNode).getContinuingHaploName())
            return false;
        return true;
    }

    /**
     * Function to decide if the source and destination nodes are a valid node pair for wide easy exchange
     *
     * @param srcNode
     * @param destNode
     */
    private boolean isaValidNodePairForWide(Node srcNode, Node destNode) {
        boolean isValid = true;
        Node srcNodeParent = srcNode.getParent();
        Node destNodeParent = destNode.getParent();
        isValid = isaValidSrcNodeForWide(srcNode);
        if (isValid)
            isValid = isaValidDestNodeForWide(srcNode,srcNodeParent,destNode);
        else
            return isValid;

        if (isValid){
            int srcHaplo = ((QuasiSpeciesNode) srcNode).getHaploAboveName();
            int destHaplo = ((QuasiSpeciesNode) destNode).getHaploAboveName();
            if ( srcHaplo != -1 && ((QuasiSpeciesNode) qsTree.getNode(srcHaplo)).getAttachmentTimesList()[0] > destNodeParent.getHeight())
                return false;
            if ( destHaplo != -1 && ((QuasiSpeciesNode) qsTree.getNode(destHaplo)).getAttachmentTimesList()[0] > srcNodeParent.getHeight())
                return false;
        }
        return isValid;
    }
}
