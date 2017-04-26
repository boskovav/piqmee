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
            int parenthaplo = -1;
            int desthaplo = -1;
            for (int i = 0; i < qsTree.getNodeCount(); i++){
                srcNode = qsTree.getNode(i);
                srcNodeParent = srcNode.getParent();
                destNodeParent = srcNodeParent.getParent();
                destNode = getOtherChild(destNodeParent, srcNodeParent);
                parenthaplo = ((QuasiSpeciesNode) srcNodeParent).getHaploAboveName();
                desthaplo = ((QuasiSpeciesNode) destNode).getHaploAboveName();
                if ( !srcNode.isRoot() && !srcNodeParent.isRoot()
                        && (   ( parenthaplo == -1
                                 &&
                                 (desthaplo == -1 || ((QuasiSpeciesNode) qsTree.getNode(desthaplo)).getAttachmentTimesList()[1] < srcNodeParent.getHeight()) )
                            || ( parenthaplo == ((QuasiSpeciesNode)getOtherChild(srcNodeParent, srcNode)).getContinuingHaploName()
                                 &&
                                 ((QuasiSpeciesNode) qsTree.getNode(desthaplo)).getAttachmentTimesList()[1] < srcNodeParent.getHeight() )
                            )
                   ){
                    cousinhaplobelowparentnode = 1;
                    break;
                }
            }
            if (cousinhaplobelowparentnode == 0)
                return Double.NEGATIVE_INFINITY;
            // select a node for narrow exchange
            do {
                srcNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
                srcNodeParent = srcNode.getParent();
                destNodeParent = srcNodeParent.getParent();
                destNode = getOtherChild(destNodeParent, srcNodeParent);
                parenthaplo = ((QuasiSpeciesNode) srcNodeParent).getHaploAboveName();
                desthaplo = ((QuasiSpeciesNode) destNode).getHaploAboveName();
            } while (srcNode.isRoot() || srcNode.getParent().isRoot()
                    || (   ( parenthaplo != -1
                             ||
                             (desthaplo != -1 && ((QuasiSpeciesNode) qsTree.getNode(desthaplo)).getAttachmentTimesList()[1] >= srcNodeParent.getHeight()) )
                        && ( parenthaplo != ((QuasiSpeciesNode)getOtherChild(srcNodeParent, srcNode)).getContinuingHaploName()
                             ||
                             ((QuasiSpeciesNode) qsTree.getNode(desthaplo)).getAttachmentTimesList()[1] >= srcNodeParent.getHeight() )
                        )
            );
        }

        // Wide exchange selection:

        else {
            // check that this operator can actually perform a move
            int haplobelowparentnode = 0;
            int srchaplo = -1;
            int desthaplo = -1;
            for (int i = 0; i < qsTree.getNodeCount(); i++){
                for (int j = i + 1; j < qsTree.getNodeCount(); j++) {
                    srcNode = qsTree.getNode(i);
                    srcNodeParent = srcNode.getParent();
                    destNode = qsTree.getNode(j);
                    destNodeParent = destNode.getParent();
                    srchaplo = ((QuasiSpeciesNode) srcNode).getHaploAboveName();
                    desthaplo = ((QuasiSpeciesNode) destNode).getHaploAboveName();
                    if (!srcNode.isRoot() && destNode != srcNode && !destNode.isRoot() && destNode.getParent() != srcNodeParent
                            && srcNodeParent != destNode && destNode.getParent() != srcNode
                            && ( (srchaplo == -1 || (srchaplo == ((QuasiSpeciesNode) srcNode).getContinuingHaploName() && ((QuasiSpeciesNode) qsTree.getNode(srchaplo)).getAttachmentTimesList()[1] < destNodeParent.getHeight()) )
                                 &&
                                 (desthaplo == -1 || (desthaplo == ((QuasiSpeciesNode) destNode).getContinuingHaploName() && ((QuasiSpeciesNode) qsTree.getNode(desthaplo)).getAttachmentTimesList()[1] < srcNodeParent.getHeight()) )
                                )
                            ) {
                        haplobelowparentnode = 1;
                        break;
                    }
                }
            }
            if (haplobelowparentnode == 0)
                return Double.NEGATIVE_INFINITY;
            // select a node for wide exchange
            do {
                srcNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
            } while (srcNode.isRoot()   || (srchaplo != -1 && srchaplo != ((QuasiSpeciesNode) srcNode).getContinuingHaploName()));
            srcNodeParent = srcNode.getParent();

            do {
                destNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
            } while(destNode == srcNode || destNode.isRoot() || destNode.getParent() == srcNodeParent
                    || srcNodeParent == destNode || destNode.getParent() == srcNode
                                        || (desthaplo != -1 && desthaplo != ((QuasiSpeciesNode) destNode).getContinuingHaploName()));
            destNodeParent = destNode.getParent();

            // Reject if substitution would result in negative branch lengths:
            if (((QuasiSpeciesNode) qsTree.getNode(srchaplo)).getAttachmentTimesList()[1] > destNodeParent.getHeight()
                    || ((QuasiSpeciesNode) qsTree.getNode(desthaplo)).getAttachmentTimesList()[1] > srcNodeParent.getHeight())
                return Double.NEGATIVE_INFINITY;
        }


        // Reject if substitution would result in negative branch lengths:
        if (destNode.getHeight() > srcNodeParent.getHeight() || srcNode.getHeight() > destNodeParent.getHeight())
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
}
