package quasispeciestree.operators;


import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import quasispeciestree.tree.QuasiSpeciesNode;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * @author Veronika Boskova created on 02/02/2016.
 */
@Description("Implements the unweighted Wilson-Balding branch"
        +" swapping move.  This move is similar to one proposed by WILSON"
        +" and BALDING 1998 and involves removing a subtree and"
        +" re-attaching it on a new parent branch. "
        +" See <a href='http://www.genetics.org/cgi/content/full/161/3/1307/F1'>picture</a>."
        +" This version selects a new starting point for the 'interrupted' haplotype"
        + "rescales all the attachment times of the corresponding haplotype on a newly "
        + "generated branch")

// TODO 1 for the normal Wilson-balding need to implement pushing down of all QS in tips of chosen subtree
// TODO    and then moving and reseeding QS on that sub-tree.  -->> actually, I only need to choose a new attachemtn time
// TODO    which is delimited with "new" parent haplo and the current haplotype's tip height!

// TODO 2 make a function in QS tree to re-seed start points of QS and their copy sequences --- can be used by above WB
// TODO    -->> re-seeding not necessary, keep everything as it was, just randomly choose haplo start point, if possible and scale all the copies

// TODO 3 Do I need to re-calculate the array holding the number of possible attachment branches of each internal node?
// TODO    -->> I think YES


/*
* TODO: need to implement check for whether the proposed branching times for each duplicate of a haplotype are
       in accordance with "appear-only-once-in-tree-history" assumption --- where to do it?  -- in operator? to be efficient
       -- randomly choose a haplotype, propose its attachment times from root to tip, if another haplotype attaches to the path (or below),
       where the QS of the first haplotype chose to attach, then allow these to only propose within attachmet node to first haplotype
       to the tip of the second haplotype, or to the tip of the haplotypes below.

       Also check in tree operators whether the QS starting/attachment times are still compatible, if not propose new QS times...
 *
 */


public class QuasiSpeciesWilsonBalding extends QuasiSpeciesTreeOperator{
////public class TypedWilsonBalding extends UniformizationRetypeOperator {
    /**
     * Change the start time and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {

        // mark the tree as dirty (startEditing)
        qsTree.startEditing(null);

        // Check that operator can be applied to tree:
        if (qsTree.getLeafNodeCount()<3)
            throw new IllegalStateException("Tree too small for"
                    +" TypedWilsonBalding operator.");

        // Select source node:
        Node srcNode;
        do {
            srcNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
        } while (invalidSrcNode(srcNode));
        Node srcNodeP = srcNode.getParent();
        Node srcNodeS = getOtherChild(srcNodeP, srcNode);
        int srcHaplo = ((QuasiSpeciesNode) srcNode).getContinuingHaploName();
        // if there is no haplotype passing the current chosen node,
        //  get all the possible haplotypes that can be moved up
        ArrayList<Integer> possibleHaplo = new ArrayList<>();
        int[] oldParentHaplo = qsTree.getParentHaplo();
        if (srcHaplo == -1){
            checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, oldParentHaplo, possibleHaplo);
            // choose randomly from the array of haplotypes one to move up
            srcHaplo = possibleHaplo.get(Randomizer.nextInt(possibleHaplo.size()));
        }
        else {
            possibleHaplo.add(srcHaplo);
        }
        double t_srcNode = srcNode.getHeight();
        double t_srcNodeP = srcNodeP.getHeight();
        double t_srcNodeS = srcNodeS.getHeight();

        // Select destination branch node:
        Node destNode;
        do {
            destNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
        } while (invalidDestNode(srcNode, destNode));
        Node destNodeP = destNode.getParent();
        double t_destNode = destNode.getHeight();


        // Handle special cases involving root:
//        if (destNode.isRoot()) {
//            // FORWARD ROOT MOVE
//
//            double logHastingsRatio = 0.0;
//
//
//            // get probability that the current haplotype attaches where it does now ..
//            // get branch length to MRCA
//            // Record probability of current colouring:
//            logHastingsRatio += getBranchTypeProb(srcNode);
//
//            // Record srcNode grandmother height:
//            double t_srcNodeG = srcNodeP.getParent().getHeight();
//
//            // Choose new root height:
//            double newTime = qs_destNode+Randomizer.nextExponential(1.0/(alpha*qs_destNode));
//
//            // Implement tree changes:
//            disconnectBranch(srcNode);
//            // TODO return which haplotype is moved --- we have to re-attach and rescale it!
//            // TODO and re-assign continuing haplotypes/parenthaploarray from parent node onwards
//            connectBranchToRoot(srcNode, destNode, newTime);
//            qsTree.setRoot(srcNodeP);
//
//            // Recolour root branches:
//            try {
//                logHastingsRatio -= retypeRootBranches(srcNode);
//            } catch (NoValidPathException e) {
//                return Double.NEGATIVE_INFINITY;
//            }
//
//            // Return HR:
//            logHastingsRatio += Math.log(alpha*qs_destNode)
//                    +(1.0/alpha)*(newTime/qs_destNode-1.0)
//                    -Math.log(t_srcNodeG-Math.max(qs_srcNode, qs_srcNodeS));
//
//            return logHastingsRatio;
//        }






//
//        if (srcNodeP.isRoot()) {
//            // BACKWARD ROOT MOVE
//
//            double logHastingsRatio = 0.0;
//
//            // get probability that the current haplotype attaches where it does now ..
//            // get branch length to MRCA
//            logHastingsRatio += getRootBranchTypeProb(srcNode);
//
//            // Record old srcNode parent height
//            double oldTime = t_srcNodeP;
//
//            // Choose height of new attachement point:
//            double min_newTime = Math.max(t_srcNode, t_destNode);
//            double t_destNodeP = destNodeP.getHeight();
//            double span = t_destNodeP-min_newTime;
//            double newTime = min_newTime+span*Randomizer.nextDouble();
//
//            // Implement tree changes:
//            disconnectBranchFromRoot(srcNode);
//            connectBranch(srcNode, destNode, newTime);
//            srcNodeS.setParent(null);
//            mtTree.setRoot(srcNodeS);
//
//            // Recolour new branch:
//            try {
//                logHastingsRatio -= retypeBranch(srcNode);
//            } catch (NoValidPathException e) {
//                return Double.NEGATIVE_INFINITY;
//            }
//
//            // Return HR:
//            logHastingsRatio += Math.log(t_destNodeP-Math.max(t_srcNode, t_destNode))
//                    -Math.log(alpha*t_srcNodeS)
//                    -(1.0/alpha)*(oldTime/t_srcNodeS-1.0);
//
//            return logHastingsRatio;
//        }










        // NON-ROOT MOVE

        double logHastingsRatio = 0.0;
        // Incorporate probability of choosing current haplotype to move
        logHastingsRatio += Math.log(possibleHaplo.size());

        // Incorporate probability of current positioning of the start time of the haplotype.
        double haploStartMin = qsTree.getNode(srcHaplo).getHeight();
        double haploStartMax = getMaxPossibleHaploAttachTime(srcHaplo,(QuasiSpeciesNode) srcNode);

        logHastingsRatio -= Math.log(haploStartMax-haploStartMin);

        // Record srcNode grandmother height:
        double t_srcNodeG = srcNodeP.getParent().getHeight();

        // Choose height of new attachment point:
        double min_newTime = Math.max(t_destNode, t_srcNode);
        double t_destNodeP = destNodeP.getHeight();
        double span = t_destNodeP-min_newTime;
        double newTime = min_newTime+span*Randomizer.nextDouble();

        // Choose attachment point for the moved haplotype
        QuasiSpeciesNode haploStartMaxNode = getMaxPossibleHaploAttachNode((QuasiSpeciesNode) destNode);
        double haploStartMaxNew = getMaxPossibleHaploAttachTime((QuasiSpeciesNode) destNode);
        // choose new time to attach
        // get a random number deciding where the current haplo will be moved
        double u = Randomizer.nextDouble();
        double tHaploNew = u*haploStartMin + (1-u)*haploStartMin;
        // reposition attachment times: attach ((time - haploStartMin) * (tHaploNew/told)) + haploStartMin
        Double[] tempqstimes=qsTree.getAttachmentTimesList(srcHaplo).clone();
        // get the haplotype's starting time
        double told = tempqstimes[0];
        // Scale the haplotype strains
        // scale all the other positions in the array but the 0 position (haplo start time)
        double scalefactor = (tHaploNew - haploStartMin)/(told - haploStartMaxNew);
        for (int i=1; i<tempqstimes.length; i++) {
            tempqstimes[i] = ((tempqstimes[i] - haploStartMin) * scalefactor) + haploStartMin;
        }
        // set the haplotype's starting time to the new time
        tempqstimes[0] = tHaploNew;
        // rewrite the attachment times array
        qsTree.setAttachmentTimesList(srcHaplo, tempqstimes);
        // Incorporate probability of current positioning of the start time of the haplotype.
        logHastingsRatio += Math.log(haploStartMaxNew-haploStartMin);


        // Implement tree changes:
        QuasiSpeciesNode correctTreeFromThisNode1 = disconnectBranch((QuasiSpeciesNode) srcNode, srcHaplo);
        QuasiSpeciesNode correctTreeFromThisNode2 = connectBranch((QuasiSpeciesNode) srcNode,
                (QuasiSpeciesNode) destNode, newTime, srcHaplo, tHaploNew);

        // Recalculate continuingHaplo and HaploAbove arrays
        if(correctTreeFromThisNode1 == qsTree.getRoot()){
            recalculateParentHaploAndCorrectContinuingHaploName(-1, correctTreeFromThisNode1);
        }
        else if (correctTreeFromThisNode2 == qsTree.getRoot()){
            recalculateParentHaploAndCorrectContinuingHaploName(-1, correctTreeFromThisNode2);
        }
        else {
            // for source node, need to find what the parent haplo was and start at the changed node
            recalculateParentHaploAndCorrectContinuingHaploName(oldParentHaplo[srcHaplo], correctTreeFromThisNode1);
            // for the destination node, need to find out new parent haplo and start at the changed node
            recalculateParentHaploAndCorrectContinuingHaploName(oldParentHaplo[haploStartMaxNode.getContinuingHaploName()], correctTreeFromThisNode2);
        }

        // Incorporate probability of choosing current haplotype to move back
        if (((QuasiSpeciesNode) srcNode).getContinuingHaploName() != srcHaplo){
            // check how many haplotypes can be pulled up
            ArrayList<Integer> backPossibleHaplo = new ArrayList<>();
            int[] newParentHaplo = qsTree.getParentHaplo();
            checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode,newParentHaplo, backPossibleHaplo);
            // account for this in the hastings ratio
            logHastingsRatio -= Math.log(backPossibleHaplo.size());
        }
        // no need to be adding 0.0, so just commented out
        // else{
        //  logHastingsRatio -= Math.log(1);
        // }

        // HR contribution of topology and node height changes:
        logHastingsRatio += Math.log(t_destNodeP-Math.max(t_srcNode, t_destNode))
                           -Math.log(t_srcNodeG-Math.max(t_srcNode, t_srcNodeS));

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
                // swapping left right, does not change anything
                ||destNode.getParent()==srcNode.getParent())
            return true;

        Node destNodeP = destNode.getParent();
        // make sure that the target branch is above the subtree being moved
        if (destNodeP!=null&&(destNodeP.getHeight()<=srcNode.getHeight()))
            return true;

        return false;
    }

}