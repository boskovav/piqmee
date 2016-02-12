package quasispeciestree.operators;


import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import quasispeciestree.tree.QuasiSpeciesNode;
import java.util.ArrayList;


/**
 * @author Veronika Boskova created on 02/02/2016 finished on 11.02.2016.
 *      largely based on TypedWilsonBalding
 */
@Description("Implements the unweighted Wilson-Balding branch"
        +" swapping move.  This move is similar to one proposed by WILSON"
        +" and BALDING 1998 and involves removing a subtree and"
        +" re-attaching it on a new parent branch. "
        +" See <a href='http://www.genetics.org/cgi/content/full/161/3/1307/F1'>picture</a>."
        +" This version selects a new starting point for the 'interrupted' haplotype"
        + "rescales all the attachment times of the corresponding haplotype on a newly "
        + "generated branch")
public class QuasiSpeciesWilsonBalding extends QuasiSpeciesTreeOperator{


    public Input<Double> alphaInput = new Input<>("alpha",
        "Root height proposal parameter", Validate.REQUIRED);
    private double alpha;

    @Override
    public void initAndValidate() throws Exception {
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
        int[] oldParentHaplo = qsTree.getParentHaplo().clone();
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
        if (destNode.isRoot()) {
            // FORWARD ROOT MOVE

            double logHastingsRatio = 0.0;

            // Incorporate probability of choosing current haplotype to move
            logHastingsRatio += Math.log(possibleHaplo.size());

            // Incorporate probability of current positioning of the start time of the haplotype.
            double haploStartMin = qsTree.getNode(srcHaplo).getHeight();
            double haploStartMax = getMaxPossibleHaploAttachTime(srcHaplo,(QuasiSpeciesNode) srcNode);

            logHastingsRatio -= Math.log(haploStartMax-haploStartMin);

            // Record srcNode grandmother height:
            double t_srcNodeG = srcNodeP.getParent().getHeight();

            // Choose new root height:
            double newTime = t_destNode+Randomizer.nextExponential(1.0/(alpha*t_destNode));

            // Choose attachment point for the moved haplotype
            double haploStartMaxNew = qsTree.originInput.get().getValue();
            // choose new time to attach
            // get a random number deciding where the current haplo will be moved
            double u = Randomizer.nextDouble();
            double tHaploNew = u*haploStartMin + (1-u)*haploStartMaxNew;
            // reposition attachment times: attach ((time - haploStartMin) * (tHaploNew/told)) + haploStartMin
            Double[] tempqstimes=qsTree.getAttachmentTimesList(srcHaplo).clone();
            // get the haplotype's starting time
            double told = tempqstimes[0];
            // Scale the haplotype strains
            // scale all the other positions in the array but the 0 position (haplo start time)
            double scalefactor = (tHaploNew - haploStartMin)/(told - haploStartMin);
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
            QuasiSpeciesNode correctTreeFromThisNode2 = connectBranchToRoot((QuasiSpeciesNode) srcNode,
                    (QuasiSpeciesNode) destNode, newTime, srcHaplo, tHaploNew);
            qsTree.setRoot(srcNodeP);

            // Recalculate continuingHaplo and HaploAbove arrays
            recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) srcNodeP);

            // Incorporate probability of choosing current haplotype to move back
            if (((QuasiSpeciesNode) srcNode).getContinuingHaploName() != srcHaplo){
                // check how many haplotypes can be pulled up
                ArrayList<Integer> backPossibleHaplo = new ArrayList<>();
                int[] newParentHaplo = qsTree.getParentHaplo();
                checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, newParentHaplo, backPossibleHaplo);
                // account for this in the hastings ratio
                logHastingsRatio -= Math.log(backPossibleHaplo.size());
            }
            // no need to be adding 0.0, so just commented out
            // else{
            //  logHastingsRatio -= Math.log(1);
            // }

            // Return HR:
            // P(moving back)=1/(time span moving back)
            // P(moving forth)=(1/(alpha*t_destNode)*exp(-1/(alpha*t_destNode)*(newTime-t_destNode))) // exponentially distr
            // HR=P(back)/P(forth)=(alpha*t_destNode)*exp(1/(alpha*t_destNode)*(newTime-t_destNode))/(time span moving back)
            // log(HR)=log(alpha*t_destNode)+(1/(alpha*t_destNode)*(newTime-t_destNode))-log(time span moving back)
            logHastingsRatio += Math.log(alpha*t_destNode)+(1.0/alpha)*(newTime/t_destNode-1.0)
                               -Math.log(t_srcNodeG-Math.max(t_srcNode, t_srcNodeS));

            return logHastingsRatio;
        }

        if (srcNodeP.isRoot()) {
            // BACKWARD ROOT MOVE

            double logHastingsRatio = 0.0;

            // Incorporate probability of choosing current haplotype to move
            logHastingsRatio += Math.log(possibleHaplo.size());

            // Incorporate probability of current positioning of the start time of the haplotype.
            double haploStartMin = qsTree.getNode(srcHaplo).getHeight();
            double haploStartMax = getMaxPossibleHaploAttachTime(srcHaplo,(QuasiSpeciesNode) srcNode);

            logHastingsRatio -= Math.log(haploStartMax-haploStartMin);

            // Record old srcNode parent height
            double oldTime = t_srcNodeP;

            // Choose height of new attachement point:
            double min_newTime = Math.max(t_srcNode, t_destNode);
            double t_destNodeP = destNodeP.getHeight();
            double span = t_destNodeP-min_newTime;
            double newTime = min_newTime+span*Randomizer.nextDouble();

            // Choose attachment point for the moved haplotype
            QuasiSpeciesNode haploStartMaxNode = getMaxPossibleHaploAttachNode((QuasiSpeciesNode) srcNode, (QuasiSpeciesNode) destNode, newTime);
            double haploStartMaxNew = getMaxPossibleHaploAttachTime((QuasiSpeciesNode) destNode, newTime);
            // choose new time to attach
            // get a random number deciding where the current haplo will be moved
            double u = Randomizer.nextDouble();
            double tHaploNew = u*haploStartMin + (1-u)*haploStartMaxNew;
            // reposition attachment times: attach ((time - haploStartMin) * (tHaploNew/told)) + haploStartMin
            Double[] tempqstimes=qsTree.getAttachmentTimesList(srcHaplo).clone();
            // get the haplotype's starting time
            double told = tempqstimes[0];
            // Scale the haplotype strains
            // scale all the other positions in the array but the 0 position (haplo start time)
            double scalefactor = (tHaploNew - haploStartMin)/(told - haploStartMin);
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
            QuasiSpeciesNode correctTreeFromThisNode1 = disconnectBranchFromRoot((QuasiSpeciesNode) srcNode, srcHaplo);
            QuasiSpeciesNode correctTreeFromThisNode2 = connectBranch((QuasiSpeciesNode) srcNode,
                    (QuasiSpeciesNode) destNode, newTime, srcHaplo, tHaploNew);
            srcNodeS.setParent(null);
            qsTree.setRoot(srcNodeS);

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
                checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, newParentHaplo, backPossibleHaplo);
                // account for this in the hastings ratio
                logHastingsRatio -= Math.log(backPossibleHaplo.size());
            }
            // no need to be adding 0.0, so just commented out
            // else{
            //  logHastingsRatio -= Math.log(1);
            // }

            // Return HR:
            // P(moving back)=(1/(alpha*t_srcNodeS)*exp(-1/(alpha*t_srcNodeS)*(oldTime-srcNodeS))) // exponentially distr
            // P(moving forth)=1/(time span moving forth)
            // HR=P(back)/P(forth)=(time span moving forth)/(alpha*t_srcNodeS)*exp(1/(alpha*t_srcNodeS)*(oldTime-t_srcNodeS))
            // log(HR)=log(time span moving forth)-log(alpha*t_srcNodeS)-(1/(alpha*t_srcNodeS)*(oldTime-t_srcNodeS))
            logHastingsRatio += Math.log(t_destNodeP-Math.max(t_srcNode, t_destNode))
                    -Math.log(alpha*t_srcNodeS)
                    -(1.0/alpha)*(oldTime/t_srcNodeS-1.0);

            return logHastingsRatio;
        }

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
        QuasiSpeciesNode haploStartMaxNode = getMaxPossibleHaploAttachNode((QuasiSpeciesNode) srcNode, (QuasiSpeciesNode) destNode, newTime);
        double haploStartMaxNew = getMaxPossibleHaploAttachTime((QuasiSpeciesNode) destNode, newTime);
        // choose new time to attach
        // get a random number deciding where the current haplo will be moved
        double u = Randomizer.nextDouble();
        double tHaploNew = u*haploStartMin + (1-u)*haploStartMaxNew;
        // reposition attachment times: attach ((time - haploStartMin) * (tHaploNew/told)) + haploStartMin
        Double[] tempqstimes=qsTree.getAttachmentTimesList(srcHaplo).clone();
        // get the haplotype's starting time
        double told = tempqstimes[0];
        // Scale the haplotype strains
        // scale all the other positions in the array but the 0 position (haplo start time)
        double scalefactor = (tHaploNew - haploStartMin)/(told - haploStartMin);
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
            checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, newParentHaplo, backPossibleHaplo);
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