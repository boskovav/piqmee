package quasispeciestree.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;

import java.util.ArrayList;


/**
 * @author Veronika Boskova created on 02/02/2016 finished on 20/01/2017
 *      largely based on TypedWilsonBalding
 */
@Description("Implements the unweighted Wilson-Balding branch"
        +" swapping move.  This move is similar to one proposed by WILSON"
        +" and BALDING 1998 and involves removing a subtree and"
        +" re-attaching it on a new parent branch. "
        +" See <a href='http://www.genetics.org/cgi/content/full/161/3/1307/F1'>picture</a>."
        +" This version only allows moves if haplotypes are not 'interrupted' ")
public class QuasiSpeciesWilsonBaldingEasy extends QuasiSpeciesTreeOperator{

    public Input<Double> alphaInput = new Input<>("alpha",
            "Root height proposal parameter", Input.Validate.REQUIRED);
    private double alpha;

    @Override
    public void initAndValidate(){
        super.initAndValidate();

        alpha = alphaInput.get();

        // Check that operator can be applied to tree:
        if (qsTree.getLeafNodeCount() < 3)
            throw new IllegalArgumentException("Tree is too small for"
                    +" QuasiSpeciesWilsonBaldingEasy operator.");
    }

    /**
     * Change the start time and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {

        // count number of pairs valid for WB (may not be the same back and forth without root moves)
        int numberofpairs = countValidPairsForWBEasy();

        // Select source node:
        Node srcNode;
        do {
            srcNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
        } while (invalidSrcNode(srcNode) || (((QuasiSpeciesNode) srcNode).getContinuingHaploName() != ((QuasiSpeciesNode) srcNode).getHaploAboveName()) || srcNode.getParent().isRoot());

        Node srcNodeP = srcNode.getParent();
        Node srcNodeS = getOtherChild(srcNodeP, srcNode);
        double t_srcNode = srcNode.getHeight();
        double t_srcNodeP = srcNodeP.getHeight();
        double t_srcNodeS = srcNodeS.getHeight();

        int srcHaplo = ((QuasiSpeciesNode) srcNode).getContinuingHaploName();

        // Select destination branch node:
        Node destNode;
        do {
            destNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
        } while (invalidDestNode(srcNode, destNode) || destNode.isRoot());

        Node destNodeP = destNode.getParent();
        double t_destNode = destNode.getHeight();

        assert destNodeP != null;

        double scalefactor=0;
        double toldQSstart=0;
        double tnewQSstart=0;

//        // Handle special cases involving root:
//        if (destNode.isRoot()) {
//            // FORWARD ROOT MOVE
//
//            double logHastingsRatio = 0.0;
//
//            // Record srcNode grandmother height:
//            double t_srcNodeG = srcNodeP.getParent().getHeight();
//
//            // Choose new root height:
//            double newTime = t_destNode+Randomizer.nextExponential(1.0 / (alpha * t_destNode));
//
//            if (srcHaplo != -1){
//                QuasiSpeciesNode node = (QuasiSpeciesNode) qsTree.getNode(srcHaplo);
//
//                // scale the haplotype
//                if (node.getAttachmentTimesList().length > 1) {
//
//                    // Find current boundaries for haplotype start times
//                    double haploStartMin = t_srcNode;
//                    double haploStartMax = t_srcNodeP;
//
//                    double logHastingsRatioContribution = scaleThisHaplo(node, newTime, haploStartMin, node.getAttachmentTimesList()[1],
//                                                                               haploStartMax, haploStartMin, -1);
//                    if (logHastingsRatioContribution == Double.NEGATIVE_INFINITY)
//                        return Double.NEGATIVE_INFINITY;
//                    else logHastingsRatio += logHastingsRatioContribution;
//                }
//
//                tnewQSstart = node.getAttachmentTimesList()[0];
//
//                // Ensure BEAST knows to recalculate affected likelihood:
//                node.makeDirty(QuasiSpeciesTree.IS_FILTHY);
//            }
//
//            // Implement tree changes:
//            QuasiSpeciesNode correctTreeFromThisNode1 = disconnectBranch((QuasiSpeciesNode) srcNode, srcHaplo);
//            QuasiSpeciesNode correctTreeFromThisNode2 = connectBranchToRoot((QuasiSpeciesNode) srcNode,
//                    (QuasiSpeciesNode) destNode, newTime, srcHaplo, tnewQSstart);
//            qsTree.setRoot((QuasiSpeciesNode) srcNodeP);
//
//            // Recalculate continuingHaplo and HaploAbove arrays
//            recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) srcNodeP);
//            // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
//            qsTree.countAndSetPossibleStartBranches();
//
//            // HR contribution of topology and node height changes:
//            // P(moving back)=1/(time span moving back)
//            // P(moving forth)=(1/(alpha*t_destNode)*exp(-1/(alpha*t_destNode)*(newTime-t_destNode))) // exponentially distr
//            // HR=P(back)/P(forth)=(alpha*t_destNode)*exp(1/(alpha*t_destNode)*(newTime-t_destNode))/(time span moving back)
//            // log(HR)=log(alpha*t_destNode)+(1/(alpha*t_destNode)*(newTime-t_destNode))-log(time span moving back)
//            logHastingsRatio += Math.log(alpha * t_destNode) + (1.0 / alpha) * (newTime / t_destNode - 1.0)
//                              - Math.log(t_srcNodeG - Math.max(t_srcNode, t_srcNodeS));
//
//            if ((origin.getValue() - Math.max(t_srcNode, t_destNode)) == 0 || (t_srcNodeG - Math.max(t_srcNode, t_srcNodeS)) == 0) {
//                // This happens when some branch lengths are zero.
//                // If oldRange = 0, hastingsRatio == Double.POSITIVE_INFINITY and
//                // node i can be catapulted anywhere in the tree, resulting in
//                // very bad trees that are always accepted.
//                // For symmetry, newRange = 0 should therefore be ruled out as well
//                //                return Double.NEGATIVE_INFINITY;
//                throw new IllegalStateException("problem in hereeeeee: QuasiSpeciesWilsonBaldingEasy - some branch lengths are 0?");
//            }
//
//            // RETURN log(HASTINGS RATIO)
//            return logHastingsRatio;
//        }
//
//        if (srcNodeP.isRoot()) {
//            // BACKWARD ROOT MOVE
//
//            double logHastingsRatio = 0.0;
//
//            // Record old srcNode parent height
//            double oldTime = t_srcNodeP;
//
//            // Choose height of new attachment point:
//            double min_newTime = Math.max(t_srcNode, t_destNode);
//            double t_destNodeP = destNodeP.getHeight();
//            double span = t_destNodeP - min_newTime;
//            double newTime = min_newTime + span * Randomizer.nextDouble();
//
//            if (srcHaplo != -1){
//                QuasiSpeciesNode node = (QuasiSpeciesNode) qsTree.getNode(srcHaplo);
//
//                // scale the haplotype
//                if (node.getAttachmentTimesList().length > 1) {
//
//                    // Find current boundaries for haplotype start times
//                    double haploStartMin = t_srcNode;
//                    double haploStartMax = t_srcNodeP;
//
//                    double logHastingsRatioContribution = scaleThisHaplo(node, newTime, haploStartMin, node.getAttachmentTimesList()[1],
//                                                                               haploStartMax, haploStartMin, -1);
//                    if (logHastingsRatioContribution == Double.NEGATIVE_INFINITY)
//                        return Double.NEGATIVE_INFINITY;
//                    else logHastingsRatio += logHastingsRatioContribution;
//                }
//
//                tnewQSstart = node.getAttachmentTimesList()[0];
//
//                // Ensure BEAST knows to recalculate affected likelihood:
//                node.makeDirty(QuasiSpeciesTree.IS_FILTHY);
//            }
//
//            // Implement tree changes:
//            QuasiSpeciesNode correctTreeFromThisNode1 = disconnectBranchFromRoot((QuasiSpeciesNode) srcNode, srcHaplo);
//            QuasiSpeciesNode correctTreeFromThisNode2 = connectBranch((QuasiSpeciesNode) srcNode,
//                    (QuasiSpeciesNode) destNode, newTime, srcHaplo, tnewQSstart);
//            qsTree.setRoot((QuasiSpeciesNode) srcNodeS);
//
//            // Recalculate continuingHaplo and HaploAbove arrays
//            recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) qsTree.getRoot());
//            // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
//            qsTree.countAndSetPossibleStartBranches();
//
//            // Return HR:
//            // P(moving back)=(1/(alpha*t_srcNodeS)*exp(-1/(alpha*t_srcNodeS)*(oldTime-srcNodeS))) // exponentially distr
//            // P(moving forth)=1/(time span moving forth)
//            // HR=P(back)/P(forth)=(time span moving forth)/(alpha*t_srcNodeS)*exp(1/(alpha*t_srcNodeS)*(oldTime-t_srcNodeS))
//            // log(HR)=log(time span moving forth)-log(alpha*t_srcNodeS)-(1/(alpha*t_srcNodeS)*(oldTime-t_srcNodeS))
//            logHastingsRatio += Math.log(t_destNodeP - Math.max(t_srcNode, t_destNode))
//                              - Math.log(alpha * t_srcNodeS) - (1.0 / alpha) * (oldTime / t_srcNodeS - 1.0);
//            if ((t_destNodeP - Math.max(t_srcNode, t_destNode)) == 0 || (origin.getValue() - t_srcNodeS) == 0) {
//                // This happens when some branch lengths are zero.
//                // If oldRange = 0, hastingsRatio == Double.POSITIVE_INFINITY and
//                // node i can be catapulted anywhere in the tree, resulting in
//                // very bad trees that are always accepted.
//                // For symmetry, newRange = 0 should therefore be ruled out as well
//                //                return Double.NEGATIVE_INFINITY;
//                throw new IllegalStateException("problem in hereeeeee: QuasiSpeciesWilsonBaldingEasy - some branch lengths are 0?");
//
//            }
//
//            // RETURN log(HASTINGS RATIO)
//            return logHastingsRatio;
//        }

        // NON-ROOT MOVE

        double logHastingsRatio = 0.0;

        // Record srcNode grandmother height:
        double t_srcNodeG = srcNodeP.getParent().getHeight();

        // Choose height of new attachment point:
        double min_newTime = Math.max(t_destNode, t_srcNode);
        double t_destNodeP = destNodeP.getHeight();
        double span = t_destNodeP - min_newTime;
        double newTime = min_newTime + span * Randomizer.nextDouble();

        if (srcHaplo != -1){
            QuasiSpeciesNode node = (QuasiSpeciesNode) qsTree.getNode(srcHaplo);

            // scale the haplotype
            if (node.getAttachmentTimesList().length > 1) {

                // Find current boundaries for haplotype start times
                double haploStartMin = t_srcNode;
                double haploStartMax = t_srcNodeP;

                double logHastingsRatioContribution = scaleThisHaplo(node, newTime, haploStartMin, node.getAttachmentTimesList()[1],
                                                                           haploStartMax, haploStartMin, -1);
                if (logHastingsRatioContribution == Double.NEGATIVE_INFINITY)
                    return Double.NEGATIVE_INFINITY;
                else logHastingsRatio += logHastingsRatioContribution;
            }

            tnewQSstart = node.getAttachmentTimesList()[0];

            // Ensure BEAST knows to recalculate affected likelihood:
            node.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        }

        // Implement tree changes:
        QuasiSpeciesNode correctTreeFromThisNode1 = disconnectBranch((QuasiSpeciesNode) srcNode, srcHaplo);
        QuasiSpeciesNode correctTreeFromThisNode2 = connectBranch((QuasiSpeciesNode) srcNode,
                (QuasiSpeciesNode) destNode, newTime, srcHaplo, tnewQSstart);

        // Recalculate continuingHaplo and HaploAbove arrays
        recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) qsTree.getRoot());

        // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
        qsTree.countAndSetPossibleStartBranches();

        // HR contribution of topology and node height changes:
        logHastingsRatio += Math.log(t_destNodeP - Math.max(t_srcNode, t_destNode))
                           -Math.log(t_srcNodeG - Math.max(t_srcNode, t_srcNodeS));

        if ((t_destNodeP - Math.max(t_srcNode, t_destNode)) == 0 || (t_srcNodeG - Math.max(t_srcNode, t_srcNodeS)) == 0) {
            // This happens when some branch lengths are zero.
            // If oldRange = 0, hastingsRatio == Double.POSITIVE_INFINITY and
            // node i can be catapulted anywhere in the tree, resulting in
            // very bad trees that are always accepted.
            // For symmetry, newRange = 0 should therefore be ruled out as well
            throw new IllegalStateException("problem in hereeeeee: QuasiSpeciesWilsonBaldingEasy - some branch lengths are 0?");
        }

        // count number of pairs valid for WB (may not be the same back and forth without root moves)
        int numberofpairsback = countValidPairsForWBEasy();

        if (numberofpairs != numberofpairsback){
            logHastingsRatio += Math.log(numberofpairs);
            logHastingsRatio -= Math.log(numberofpairsback);
        }

        // RETURN log(HASTINGS RATIO)
        return logHastingsRatio;
    }

    /**
     * Returns true if srcNode CANNOT be used for the Wilson-Balding Easy move.
     *
     * @param srcNode
     * @return True if srcNode invalid.
     */
    private boolean invalidSrcNode(Node srcNode) {

        if (srcNode.isRoot())
            return true;

        Node parent = srcNode.getParent();

        // This check is important for avoiding situations where it is
        // impossible to choose a valid destNode:
        if (parent.isRoot()) {

            Node sister = getOtherChild(parent, srcNode);

            if (sister.isLeaf())
                return true;
            // make sure that the target branch is above the subtree being moved
            // if the parent is root then the only possibility is the sister branch tip/node,
            // and if this one is below the source node, we cannot attach to it
            // and attaching to root would just swap left/right = not much use
            if (srcNode.getHeight() >= sister.getHeight())
                return true;
        }

        // avoid root moves
        if (srcNode.getParent().isRoot())
            return true;

        return false;
    }

    /**
     * Returns true if destNode CANNOT be used for the Wilson-Balding Easy move in conjunction
     * with srcNode.
     *
     * @param srcNode
     * @param destNode
     * @return True if destNode invalid.
     */
    private boolean invalidDestNode(Node srcNode, Node destNode) {

        if (    // cannot attach to the same place
                destNode == srcNode
                        || destNode == srcNode.getParent()
                        // changing height of attachment, does not change much
                        || destNode.getParent() == srcNode.getParent())
            return true;

        Node destNodeP = destNode.getParent();
        // make sure that the target branch is above the subtree being moved
        if (destNodeP != null && (destNodeP.getHeight() <= srcNode.getHeight()))
            return true;

        // avoid root moves
        if (destNode.isRoot())
            return true;

        return false;
    }

    /**
     * Function to count how many valid source and destination node pairs there are for wilson-balding operator
     *
     */
    private int countValidPairsForWBEasy() {
        int count = 0;
        for (int i = 0; i < qsTree.getNodeCount(); i++) {
            for (int j = i + 1; j < qsTree.getNodeCount(); j++) {
                Node srcNode = qsTree.getNode(i);
                Node destNode = qsTree.getNode(j);
                if (!invalidSrcNode(srcNode) && !invalidDestNode(srcNode, destNode)
//                        && !(((QuasiSpeciesNode) srcNode).getContinuingHaploName() != ((QuasiSpeciesNode) srcNode).getHaploAboveName())
                        )
                    count += 1;
                srcNode = qsTree.getNode(j);
                destNode = qsTree.getNode(i);
                if (!invalidSrcNode(srcNode) && !invalidDestNode(srcNode, destNode)
//                        && !(((QuasiSpeciesNode) srcNode).getContinuingHaploName() != ((QuasiSpeciesNode) srcNode).getHaploAboveName())
                        )
                    count += 1;
            }
        }
        return count;
    }
}