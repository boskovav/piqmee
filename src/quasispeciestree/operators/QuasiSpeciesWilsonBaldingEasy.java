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
    }

    /**
     * Change the start time and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {

        // Check that operator can be applied to tree:
        if (qsTree.getLeafNodeCount()<3)
            throw new IllegalStateException("Tree too small for"
                    +" QuasiSpeciesWilsonBalding operator.");

        // Select source node:
        Node srcNode;
        do {
            srcNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
        } while (invalidSrcNode(srcNode) || (((QuasiSpeciesNode) srcNode).getContinuingHaploName()!=((QuasiSpeciesNode) srcNode).getHaploAboveName()));
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
        } while (invalidDestNode(srcNode, destNode));
        Node destNodeP = destNode.getParent();
        double t_destNode = destNode.getHeight();

        // TODO THIS IS TO AVOID ROOT MOVES FOR NOW
        // disallow moves that change the root.
        if (destNode.isRoot() || srcNodeP.isRoot()) {
            return Double.NEGATIVE_INFINITY;
        }
        assert destNodeP != null;

        double scalefactor=0;
        double toldQSstart=0;
        double tnewQSstart=0;

        // NON-ROOT MOVE

        double logHastingsRatio = 0.0;

        // Find current boundaries for haplotype start times
//        double haploStartMin = qsTree.getNode(srcHaplo).getHeight();
        double haploStartMin = t_srcNode;
        double haploStartMax = getMaxPossibleHaploAttachTime(srcHaplo,(QuasiSpeciesNode) srcNode);

        // Record srcNode grandmother height:
        double t_srcNodeG = srcNodeP.getParent().getHeight();

        // Choose height of new attachment point:
        double min_newTime = Math.max(t_destNode, t_srcNode);
        double t_destNodeP = destNodeP.getHeight();
        double span = t_destNodeP-min_newTime;
        double newTime = min_newTime+span*Randomizer.nextDouble();

        if (srcHaplo != -1){
            // Choose attachment point for the moved haplotype
            double haploStartMaxNew = newTime;
            // choose new time to attach
            // get a random number deciding where the current haplo will be moved
            double u = Randomizer.nextDouble();
//        double tnew = u*haploStartMin + (1-u)*haploStartMaxNew;
            double tnew = 0;
            // reposition attachment times: attach ((time - haploStartMin) * (tHaploNew/told)) + haploStartMin
            double[] tempqstimes=qsTree.getAttachmentTimesList(srcHaplo).clone();
            // get the haplotype's starting time
            toldQSstart = tempqstimes[0];
            double told=0;
            if (tempqstimes.length > 1){
                told = tempqstimes[1];
                // Scale the haplotype strains
                // scale all the other positions in the array but the 0 position (haplo start time)
//            scalefactor = (tnew - haploStartMin)/(told - haploStartMin);
                scalefactor = u*(haploStartMin/told) + (1.0-u)*(haploStartMaxNew/told);
                tnew = tempqstimes[1] * scalefactor;
                for (int i=1; i<tempqstimes.length; i++) {
//                tempqstimes[i] = ((tempqstimes[i] - haploStartMin) * scalefactor) + haploStartMin;
                    tempqstimes[i] = tempqstimes[i] * scalefactor;
                }
                // Reject invalid haplotype scalings:
                if (scalefactor<1.0 && tempqstimes[tempqstimes.length-1]<qsTree.getNode(srcHaplo).getHeight())
                    return Double.NEGATIVE_INFINITY;

                // set the haplotype's starting time to the new time
//            double x = Randomizer.nextDouble();
//            tnewQSstart = x*tempqstimes[1] + (1-x)*haploStartMaxNew;
                tnewQSstart = tempqstimes[1];
                tempqstimes[0] = tnewQSstart;
                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//            logHastingsRatio -= Math.log(haploStartMax - told);
//            logHastingsRatio += Math.log(haploStartMaxNew - tempqstimes[1]);
            }
            else {
                // set the haplotype's starting time to the new time
//            double x = Randomizer.nextDouble();
//            tnewQSstart = x*haploStartMin + (1-x)*haploStartMaxNew;
//            tnew = u*haploStartMin + (1-u)*haploStartMaxNew;
//            tnewQSstart = tnew;
                tnewQSstart = haploStartMin;
                tempqstimes[0] = tnewQSstart;
                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//            logHastingsRatio -= Math.log(haploStartMax - haploStartMin);
//            logHastingsRatio += Math.log(haploStartMaxNew - haploStartMin);
            }
            // rewrite the attachment times array
            qsTree.setAttachmentTimesList(srcHaplo, tempqstimes);

            // Incorporate probability of current haplotype to move
            if (tempqstimes.length > 1){
                // scaling the attachment times (see Scaling Operator)
                // assign contribution to the Hastings ratio for having different possible scales for told
                logHastingsRatio += Math.log(haploStartMaxNew/told - haploStartMin/told);
                logHastingsRatio -= Math.log(haploStartMax/tnew - haploStartMin/tnew);
                // assign contribution of each scaled attachment time
//            logHastingsRatio -= 2 * (Math.log(scalefactor * (haploStartMaxNew - haploStartMin) + haploStartMin) - Math.log(haploStartMaxNew));
//            logHastingsRatio += (tempqstimes.length-1) * Math.log(scalefactor);
                logHastingsRatio += (tempqstimes.length-3) * Math.log(scalefactor);
            }
        }


        // Implement tree changes:
        // TODO disconnect and connect branches need the correction for the fact that QS starts on the passed nodes could move up/down
        QuasiSpeciesNode correctTreeFromThisNode1 = disconnectBranch((QuasiSpeciesNode) srcNode, srcHaplo);
        QuasiSpeciesNode correctTreeFromThisNode2 = connectBranch((QuasiSpeciesNode) srcNode,
                (QuasiSpeciesNode) destNode, newTime, srcHaplo, tnewQSstart);

        // Recalculate continuingHaplo and HaploAbove arrays
        recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) qsTree.getRoot());

        // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
        int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
        qsTree.setStartBranchCounts(startBranchCountsArray);

        // HR contribution of topology and node height changes:
        logHastingsRatio += Math.log(t_destNodeP-Math.max(t_srcNode, t_destNode))
                           -Math.log(t_srcNodeG-Math.max(t_srcNode, t_srcNodeS));

        if ((t_destNodeP-Math.max(t_srcNode, t_destNode)) == 0 || (t_srcNodeG-Math.max(t_srcNode, t_srcNodeS)) == 0) {
            // This happens when some branch lengths are zero.
            // If oldRange = 0, hastingsRatio == Double.POSITIVE_INFINITY and
            // node i can be catapulted anywhere in the tree, resulting in
            // very bad trees that are always accepted.
            // For symmetry, newRange = 0 should therefore be ruled out as well
//                return Double.NEGATIVE_INFINITY;
            System.out.println("problem in hereeeeee 3");
            System.exit(0);
        }

        return logHastingsRatio;
    }

    /**
     * Returns true if srcNode CANNOT be used for the CWBR ? Wilson-Balding ? move.
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
            // and if this one is below the source node, we cannot attach it to it
            // and attaching to root would just swap left/right = not much use
            if (srcNode.getHeight()>=sister.getHeight())
                return true;
        }

        return false;
    }

    /**
     * Returns true if destNode CANNOT be used for the CWBR move in conjunction
     * with srcNode.
     *
     * @param srcNode
     * @param destNode
     * @return True if destNode invalid.
     */
    private boolean invalidDestNode(Node srcNode, Node destNode) {

        if (    // cannot attach to the same place
                destNode==srcNode
                        ||destNode==srcNode.getParent()
                        // changing height of attachment, does not change much
                        ||destNode.getParent()==srcNode.getParent())
            return true;

        Node destNodeP = destNode.getParent();
        // make sure that the target branch is above the subtree being moved
        if (destNodeP!=null&&(destNodeP.getHeight()<=srcNode.getHeight()))
            return true;

        return false;
    }

}