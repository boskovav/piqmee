package quasispeciestree.operators;


import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import beast.util.*;
import quasispeciestree.tree.QuasiSpeciesNode;
import java.util.ArrayList;


/**
 * @author Veronika Boskova created on 02/02/2016 finished on 27/05/2016
 *      largely based on TypedWilsonBalding
 */
@Description("Implements the unweighted Wilson-Balding branch"
        +" swapping move.  This move is similar to one proposed by WILSON"
        +" and BALDING 1998 and involves removing a subtree and"
        +" re-attaching it on a new parent branch. "
        +" See <a href='http://www.genetics.org/cgi/content/full/161/3/1307/F1'>picture</a>."
        +" This version selects a new starting point for the 'interrupted' or randomly chosen haplotype"
        + "from the moved subtree, rescales all the attachment times of the corresponding haplotype on a newly "
        + "generated branch with a scale factor chosen for the first attachment time and finds a new haplotype"
        + "starting point between the first attachment time and the next MRCA with another haplotype")
public class QuasiSpeciesWilsonBalding extends QuasiSpeciesTreeOperator{


    public Input<Double> alphaInput = new Input<>("alpha",
        "Root height proposal parameter", Validate.REQUIRED);
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
        //testing
        //for (int i=0; i<qsTree.getLeafNodeCount(); i++){
        //    if (qsTree.getParentHaplo()[i]==i){
        //        System.out.println("problem in hereeeeee");
        //        System.exit(0);
        //    }
        //    if (qsTree.getParentHaplo()[i]!=-1 && i==qsTree.getParentHaplo()[qsTree.getParentHaplo()[i]]){
        //        System.out.println("problem in hereeeeee");
        //        System.exit(0);
        //    }
        //}

        // Check that operator can be applied to tree:
        if (qsTree.getLeafNodeCount()<3)
            throw new IllegalStateException("Tree too small for"
                    +" TypedWilsonBalding operator.");

        // Select source node:
        Node srcNode;
        do {
            srcNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
        } while (srcNode.isRoot());
        // TODO THIS IS TO AVOID ROOT MOVES FOR NOW
//        } while (invalidSrcNode(srcNode));
        Node srcNodeP = srcNode.getParent();
        Node srcNodeS = getOtherChild(srcNodeP, srcNode);
        double t_srcNode = srcNode.getHeight();
        double t_srcNodeP = srcNodeP.getHeight();
        double t_srcNodeS = srcNodeS.getHeight();

        int srcHaplo = ((QuasiSpeciesNode) srcNode).getContinuingHaploName();
        // if there is no haplotype passing the current chosen node,
        //  get all the possible haplotypes that can be moved up
        ArrayList<Integer> possibleHaplo = new ArrayList<>();
//        int[] oldParentHaplo = qsTree.getParentHaplo().clone();
//        int[] parentHaplo = qsTree.getParentHaplo();
        if (srcHaplo == -1){
//            checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, oldParentHaplo, possibleHaplo);
            checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, qsTree.getParentHaplo(), possibleHaplo);
            // choose randomly from the array of haplotypes one to move up
            srcHaplo = possibleHaplo.get(Randomizer.nextInt(possibleHaplo.size()));
        }
        else {
            possibleHaplo.add(srcHaplo);
        }

        // Select destination branch node:
        Node destNode;
        do {
            destNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
        } while ((destNode.getParent() != null && destNode.getParent().getHeight() <= srcNode.getHeight()) || (srcNode.getNr() == destNode.getNr()));
        // TODO THIS IS TO AVOID ROOT MOVES FOR NOW
//        } while (invalidDestNode(srcNode, destNode));
        Node destNodeP = destNode.getParent();
        double t_destNode = destNode.getHeight();

        // TODO THIS IS TO AVOID ROOT MOVES FOR NOW
        // disallow moves that change the root.
        if (destNode.isRoot() || srcNodeP.isRoot()) {
            return Double.NEGATIVE_INFINITY;
        }
        assert destNodeP != null;
        if ( destNodeP.getNr() == srcNodeP.getNr() || destNode.getNr() == srcNodeP.getNr() || destNodeP.getNr() == srcNode.getNr())
            return Double.NEGATIVE_INFINITY;

        double scalefactor=0;
        double toldQSstart=0;
        double tnewQSstart=0;

        // Handle special cases involving root:
//        if (destNode.isRoot()) {
//            // FORWARD ROOT MOVE
//
//            double logHastingsRatio = 0.0;
//
//            // Find current boundaries for haplotype start times
//            double haploStartMin = qsTree.getNode(srcHaplo).getHeight();
//            double haploStartMax = getMaxPossibleHaploAttachTime(srcHaplo,(QuasiSpeciesNode) srcNode);
//
//            // Record srcNode grandmother height:
//            double t_srcNodeG = srcNodeP.getParent().getHeight();
//
//            // Choose new root height:
////            double newTime = t_destNode+Randomizer.nextExponential(1.0/(alpha*t_destNode));
////            double newTime = t_destNode + span*Randomizer.nextDouble();
//            double newTime = Randomizer.uniform(t_destNode,origin.getValue());
//
//            if (newTime>origin.getValue()){
//                        System.out.println("problem in hereeeeee");
//                        System.exit(0);
//            }
//            // Choose attachment point for the moved haplotype
////            double haploStartMaxNew = origin.getValue();
//            ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTime((QuasiSpeciesNode) destNode, srcHaplo, newTime);
//            double haploStartMaxNew = (double) haploStartMaxNewArray.get(1);
//            // choose new time to attach
//            // get a random number deciding where the current haplo will be moved
//            double u = Randomizer.nextDouble();
//            double tnew = u*haploStartMin + (1-u)*haploStartMaxNew;
//            // reposition attachment times: attach ((time - haploStartMin) * (tHaploNew/told)) + haploStartMin
//            Double[] tempqstimes=qsTree.getAttachmentTimesList(srcHaplo).clone();
//            // get the haplotype's starting time
//            toldQSstart = tempqstimes[0];
//            double told=0;
//            if (tempqstimes.length > 1){
//                told = tempqstimes[1];
//                // Scale the haplotype strains
//                // scale all the other positions in the array but the 0 position (haplo start time)
//                scalefactor = (tnew - haploStartMin)/(told - haploStartMin);
//                for (int i=1; i<tempqstimes.length; i++) {
//                    tempqstimes[i] = ((tempqstimes[i] - haploStartMin) * scalefactor) + haploStartMin;
//                }
//                // set the haplotype's starting time to the new time
//                double x = Randomizer.nextDouble();
//                tnewQSstart = x*tempqstimes[1] + (1-x)*haploStartMaxNew;
//                tempqstimes[0] = tnewQSstart;
//                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
////                logHastingsRatio -= Math.log(haploStartMax - told);
////                logHastingsRatio += Math.log(haploStartMaxNew - tempqstimes[1]);
//            }
//            else {
//                // set the haplotype's starting time to the new time
//                tnewQSstart = tnew;
//                tempqstimes[0] = tnewQSstart;
//                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
////                logHastingsRatio -= Math.log(haploStartMax - haploStartMin);
////                logHastingsRatio += Math.log(haploStartMaxNew - haploStartMin);
//            }
//            // rewrite the attachment times array
//            qsTree.setAttachmentTimesList(srcHaplo, tempqstimes);
//
//            // Implement tree changes:
//            // TODO disconnect and connect branches need the correction for the fact that QS starts on the passed nodes could move up/down
//            QuasiSpeciesNode correctTreeFromThisNode1 = disconnectBranch((QuasiSpeciesNode) srcNode, srcHaplo);
//            QuasiSpeciesNode correctTreeFromThisNode2 = connectBranchToRoot((QuasiSpeciesNode) srcNode,
//                    (QuasiSpeciesNode) destNode, newTime, srcHaplo, tnewQSstart);
//            qsTree.setRoot(srcNodeP);
//
//            // Recalculate continuingHaplo and HaploAbove arrays
//            recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) srcNodeP);
////            in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
//            int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
//            qsTree.setStartBranchCounts(startBranchCountsArray);
//
//            // Incorporate probability of choosing current haplotype to move --- only with Felsenstein
////            logHastingsRatio += Math.log(possibleHaplo.size());
//            if (tempqstimes.length>1){
//                // Incorporate probability of current positioning of the start time of the haplotype.
////                logHastingsRatio -= Math.log(haploStartMax-haploStartMin);
//                // Incorporate probability of selected positioning of the start time of the haplotype.
////                logHastingsRatio += Math.log(haploStartMaxNew-haploStartMin);
//                // scaling the attachment times (see Scaling Operator)
//                // assign contribution to the Hastings ratio for having different possible scales for told
//                logHastingsRatio += Math.log(haploStartMaxNew/told - haploStartMin/told);
//                logHastingsRatio -= Math.log(haploStartMax/tnew - haploStartMin/tnew);
//                // assign contribution of each scaled attachment time
//                logHastingsRatio -= 2 * (Math.log(scalefactor * (haploStartMaxNew - haploStartMin) + haploStartMin) - Math.log(haploStartMaxNew));
//                logHastingsRatio += (tempqstimes.length-1) * Math.log(scalefactor);
//            }
//            // Incorporate probability of choosing current haplotype to move back --- only with Felsenstein
//            if (((QuasiSpeciesNode) srcNode).getContinuingHaploName() != srcHaplo){
//                // check how many haplotypes can be pulled up
//                ArrayList<Integer> backPossibleHaplo = new ArrayList<>();
////                int[] newParentHaplo = qsTree.getParentHaplo();
////                checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, newParentHaplo, backPossibleHaplo);
//                checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, qsTree.getParentHaplo(), backPossibleHaplo);
//                // account for this in the hastings ratio
////                logHastingsRatio -= Math.log(backPossibleHaplo.size());
//            }
//            // no need to be adding 0.0, so just commented out
//            // else{
//            //  logHastingsRatio -= Math.log(1);
//            // }
//
//            // Return HR:
//            // P(moving back)=1/(time span moving back)
//            // P(moving forth)=(1/(alpha*t_destNode)*exp(-1/(alpha*t_destNode)*(newTime-t_destNode))) // exponentially distr
//            // HR=P(back)/P(forth)=(alpha*t_destNode)*exp(1/(alpha*t_destNode)*(newTime-t_destNode))/(time span moving back)
//            // log(HR)=log(alpha*t_destNode)+(1/(alpha*t_destNode)*(newTime-t_destNode))-log(time span moving back)
////            logHastingsRatio += Math.log(alpha*t_destNode)+(1.0/alpha)*(newTime/t_destNode-1.0)
//            logHastingsRatio += Math.log(origin.getValue()-t_destNode)
//                               -Math.log(t_srcNodeG-Math.max(t_srcNode, t_srcNodeS));
//
//            if ((t_srcNodeG-Math.max(t_srcNode, t_srcNodeS)) == 0 || (origin.getValue()-t_destNode) == 0) {
//                // This happens when some branch lengths are zero.
//                // If oldRange = 0, hastingsRatio == Double.POSITIVE_INFINITY and
//                // node i can be catapulted anywhere in the tree, resulting in
//                // very bad trees that are always accepted.
//                // For symmetry, newRange = 0 should therefore be ruled out as well
////                return Double.NEGATIVE_INFINITY;
//                System.out.println("problem in hereeeeee 1");
//                System.exit(0);
//            }
//
//            return logHastingsRatio;
//        }
//
//        if (srcNodeP.isRoot()) {
//            // BACKWARD ROOT MOVE
//
//            double logHastingsRatio = 0.0;
//
//            // Find current boundaries for haplotype start times
//            double haploStartMin = qsTree.getNode(srcHaplo).getHeight();
//            double haploStartMax = getMaxPossibleHaploAttachTime(srcHaplo,(QuasiSpeciesNode) srcNode);
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
//            // Choose attachment point for the moved haplotype
//            //QuasiSpeciesNode haploStartMaxNode = getMaxPossibleHaploAttachNode((QuasiSpeciesNode) srcNode, (QuasiSpeciesNode) destNode, newTime);
//            ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTime((QuasiSpeciesNode) destNode, srcHaplo, newTime);
//            double haploStartMaxNew = (double) haploStartMaxNewArray.get(1);
//            // choose new time to attach
//            // get a random number deciding where the current haplo will be moved
//            double u = Randomizer.nextDouble();
//            double tnew = u*haploStartMin + (1-u)*haploStartMaxNew;
//            // reposition attachment times: attach ((time - haploStartMin) * (tHaploNew/told)) + haploStartMin
//            Double[] tempqstimes=qsTree.getAttachmentTimesList(srcHaplo).clone();
//            // get the haplotype's starting time
//            toldQSstart = tempqstimes[0];
//            double told=0;
//            if (tempqstimes.length > 1){
//                told = tempqstimes[1];
//                // Scale the haplotype strains
//                // scale all the other positions in the array but the 0 position (haplo start time)
//                scalefactor = (tnew - haploStartMin)/(told - haploStartMin);
//                for (int i=1; i<tempqstimes.length; i++) {
//                    tempqstimes[i] = ((tempqstimes[i] - haploStartMin) * scalefactor) + haploStartMin;
//                }
//                // set the haplotype's starting time to the new time
//                double x = Randomizer.nextDouble();
//                tnewQSstart = x*tempqstimes[1] + (1-x)*haploStartMaxNew;
//                tempqstimes[0] = tnewQSstart;
//                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
////                logHastingsRatio -= Math.log(haploStartMax - told);
////                logHastingsRatio += Math.log(haploStartMaxNew - tempqstimes[1]);
//            }
//            else {
//                // set the haplotype's starting time to the new time
////                double x = Randomizer.nextDouble();
////                tnewQSstart = x*haploStartMin + (1-x)*haploStartMaxNew;
//                tnewQSstart = tnew;
//                tempqstimes[0] = tnewQSstart;
//                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
////                logHastingsRatio -= Math.log(haploStartMax - haploStartMin);
////                logHastingsRatio += Math.log(haploStartMaxNew - haploStartMin);
//            }
//            // rewrite the attachment times array
//            qsTree.setAttachmentTimesList(srcHaplo, tempqstimes);
//
//            // Implement tree changes:
//            // TODO disconnect and connect branches need the correction for the fact that QS starts on the passed nodes could move up/down
//            QuasiSpeciesNode correctTreeFromThisNode1 = disconnectBranchFromRoot((QuasiSpeciesNode) srcNode, srcHaplo);
//            QuasiSpeciesNode correctTreeFromThisNode2 = connectBranch((QuasiSpeciesNode) srcNode,
//                    (QuasiSpeciesNode) destNode, newTime, srcHaplo, tnewQSstart);
//            qsTree.setRoot(srcNodeS);
//
//            // Recalculate continuingHaplo and HaploAbove arrays
////            recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) qsTree.getRoot());
//            if(correctTreeFromThisNode1 == qsTree.getRoot()){
//                recalculateParentHaploAndCorrectContinuingHaploName(-1, correctTreeFromThisNode1);
//            }
//            else if (correctTreeFromThisNode2 == qsTree.getRoot()){
//                recalculateParentHaploAndCorrectContinuingHaploName(-1, correctTreeFromThisNode2);
//            }
//            else {
//                // for source node, need to find what the parent haplo was and start at the changed node
//                if (correctTreeFromThisNode1!=null)
////                    recalculateParentHaploAndCorrectContinuingHaploName(oldParentHaplo[srcHaplo], correctTreeFromThisNode1);
//                    recalculateParentHaploAndCorrectContinuingHaploName(qsTree.getParentHaplo()[srcHaplo], correctTreeFromThisNode1);
//                // for the attached haplotype, need to find out new parent haplo and start at the changed node
//                // changed node can be newly attached parent or any node above
//                // getMaxPossibleHaploAttachTime function returned array of Node, time, haplo (above in the tree)
//                if (correctTreeFromThisNode2.getHaploAboveName()!=-1){
//                    if ((int) haploStartMaxNewArray.get(2) != correctTreeFromThisNode2.getHaploAboveName())
//                        recalculateParentHaploAndCorrectContinuingHaploName((int) haploStartMaxNewArray.get(2), correctTreeFromThisNode2);
//                    else
////                        recalculateParentHaploAndCorrectContinuingHaploName(oldParentHaplo[(int) haploStartMaxNewArray.get(2)], correctTreeFromThisNode2);
//                        recalculateParentHaploAndCorrectContinuingHaploName(qsTree.getParentHaplo()[(int) haploStartMaxNewArray.get(2)], correctTreeFromThisNode2);
//                }
//                // if the continuingHaplo is -1, we must have attached the haplo below the moved parent
//                // there can however still be some other haplo up in the tree (and not coming to this part)
//                else{
//                    // we have to start correcting from the srcNode
//                    int contHaplo = (int) haploStartMaxNewArray.get(2);
//                    boolean contHaploHere = false;
//                    for (Node node : correctTreeFromThisNode2.getAllLeafNodes()){
//                        if (node.getNr() == contHaplo){
//                            contHaploHere = true;
//                            break;
//                        }
//                    }
//                    if (contHaploHere)
//                        correctTreeFromThisNode2.setContinuingHaploName(contHaplo);
//                    else
//                        correctTreeFromThisNode2.setContinuingHaploName(-1);
//                    recalculateParentHaploAndCorrectContinuingHaploName(contHaplo, (QuasiSpeciesNode) srcNode);
//                }
//            }
//            // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
//            int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
//            qsTree.setStartBranchCounts(startBranchCountsArray);
//
//            // Incorporate probability of choosing current haplotype to move --- only with Felsenstein
////            logHastingsRatio += Math.log(possibleHaplo.size());
//            if (tempqstimes.length>1){
//                // Incorporate probability of current positioning of the start time of the haplotype.
////                logHastingsRatio -= Math.log(haploStartMax-haploStartMin);
//                // Incorporate probability of selected positioning of the start time of the haplotype.
////                logHastingsRatio += Math.log(haploStartMaxNew-haploStartMin);
//                // scaling the attachment times (see Scaling Operator)
//                // assign contribution to the Hastings ratio for having different possible scales for told
//                logHastingsRatio += Math.log(haploStartMaxNew/told - haploStartMin/told);
//                logHastingsRatio -= Math.log(haploStartMax/tnew - haploStartMin/tnew);
//                // assign contribution of each scaled attachment time
//                logHastingsRatio -= 2 * (Math.log(scalefactor * (haploStartMaxNew - haploStartMin) + haploStartMin) - Math.log(haploStartMaxNew));
//                logHastingsRatio += (tempqstimes.length-1) * Math.log(scalefactor);
//            }
//            // Incorporate probability of choosing current haplotype to move back --- only with Felsenstein
//            if (((QuasiSpeciesNode) srcNode).getContinuingHaploName() != srcHaplo){
//                // check how many haplotypes can be pulled up
//                ArrayList<Integer> backPossibleHaplo = new ArrayList<>();
////                int[] newParentHaplo = qsTree.getParentHaplo();
////                checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, newParentHaplo, backPossibleHaplo);
//                checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, qsTree.getParentHaplo(), backPossibleHaplo);
//                // account for this in the hastings ratio
////                logHastingsRatio -= Math.log(backPossibleHaplo.size());
//            }
//            // no need to be adding 0.0, so just commented out
//            // else{
//            //  logHastingsRatio -= Math.log(1);
//            // }
//
//            // Return HR:
//            // P(moving back)=(1/(alpha*t_srcNodeS)*exp(-1/(alpha*t_srcNodeS)*(oldTime-srcNodeS))) // exponentially distr
//            // P(moving forth)=1/(time span moving forth)
//            // HR=P(back)/P(forth)=(time span moving forth)/(alpha*t_srcNodeS)*exp(1/(alpha*t_srcNodeS)*(oldTime-t_srcNodeS))
//            // log(HR)=log(time span moving forth)-log(alpha*t_srcNodeS)-(1/(alpha*t_srcNodeS)*(oldTime-t_srcNodeS))
//            logHastingsRatio += Math.log(t_destNodeP-Math.max(t_srcNode, t_destNode))
//                               -Math.log(origin.getValue() - t_srcNodeS);
////                    -Math.log(alpha*t_srcNodeS)
////                    -(1.0/alpha)*(oldTime/t_srcNodeS-1.0);
//            if ((t_destNodeP-Math.max(t_srcNode, t_destNode)) == 0 || (origin.getValue()-t_srcNodeS) == 0) {
//                // This happens when some branch lengths are zero.
//                // If oldRange = 0, hastingsRatio == Double.POSITIVE_INFINITY and
//                // node i can be catapulted anywhere in the tree, resulting in
//                // very bad trees that are always accepted.
//                // For symmetry, newRange = 0 should therefore be ruled out as well
////                return Double.NEGATIVE_INFINITY;
//                System.out.println("problem in hereeeeee 2");
//                        System.exit(0);
//            }
//
//            return logHastingsRatio;
//        }

        // NON-ROOT MOVE

        double logHastingsRatio = 0.0;

        // Find current boundaries for haplotype start times
        double haploStartMin = qsTree.getNode(srcHaplo).getHeight();
        double haploStartMax = getMaxPossibleHaploAttachTime(srcHaplo,(QuasiSpeciesNode) srcNode);

        // Record srcNode grandmother height:
        double t_srcNodeG = srcNodeP.getParent().getHeight();

        // Choose height of new attachment point:
        double min_newTime = Math.max(t_destNode, t_srcNode);
        double t_destNodeP = destNodeP.getHeight();
        double span = t_destNodeP-min_newTime;
        double newTime = min_newTime+span*Randomizer.nextDouble();

        // Choose attachment point for the moved haplotype
        //QuasiSpeciesNode haploStartMaxNode = getMaxPossibleHaploAttachNode((QuasiSpeciesNode) srcNode, (QuasiSpeciesNode) destNode, newTime);
        ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTimeForQSStart((QuasiSpeciesNode) destNode, srcHaplo, newTime);
        double haploStartMaxNew = (double) haploStartMaxNewArray.get(1);
        // choose new time to attach
        // get a random number deciding where the current haplo will be moved
        double u = Randomizer.nextDouble();
//        double tnew = u*haploStartMin + (1-u)*haploStartMaxNew;
        double tnew = 0;
        // reposition attachment times: attach ((time - haploStartMin) * (tHaploNew/told)) + haploStartMin
        Double[] tempqstimes=qsTree.getAttachmentTimesList(srcHaplo).clone();
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
            double x = Randomizer.nextDouble();
            tnewQSstart = x*tempqstimes[1] + (1-x)*haploStartMaxNew;
            tempqstimes[0] = tnewQSstart;
            // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//            logHastingsRatio -= Math.log(haploStartMax - told);
//            logHastingsRatio += Math.log(haploStartMaxNew - tempqstimes[1]);
        }
        else {
            // set the haplotype's starting time to the new time
//            double x = Randomizer.nextDouble();
//            tnewQSstart = x*haploStartMin + (1-x)*haploStartMaxNew;
            tnew = u*haploStartMin + (1-u)*haploStartMaxNew;
            tnewQSstart = tnew;
            tempqstimes[0] = tnewQSstart;
            // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//            logHastingsRatio -= Math.log(haploStartMax - haploStartMin);
//            logHastingsRatio += Math.log(haploStartMaxNew - haploStartMin);
        }
        // rewrite the attachment times array
        qsTree.setAttachmentTimesList(srcHaplo, tempqstimes);

        // Implement tree changes:
        // TODO disconnect and connect branches need the correction for the fact that QS starts on the passed nodes could move up/down
        QuasiSpeciesNode correctTreeFromThisNode1 = disconnectBranch((QuasiSpeciesNode) srcNode, srcHaplo);
        QuasiSpeciesNode correctTreeFromThisNode2 = connectBranch((QuasiSpeciesNode) srcNode,
                (QuasiSpeciesNode) destNode, newTime, srcHaplo, tnewQSstart);

        // Recalculate continuingHaplo and HaploAbove arrays
        recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) qsTree.getRoot());
//        if(correctTreeFromThisNode1 == qsTree.getRoot()){
//            recalculateParentHaploAndCorrectContinuingHaploName(-1, correctTreeFromThisNode1);
//        }
//        else if (correctTreeFromThisNode2 == qsTree.getRoot()){
//            recalculateParentHaploAndCorrectContinuingHaploName(-1, correctTreeFromThisNode2);
//        }
//        else {
//            // for source node, need to find what the parent haplo was and start at the changed node
//            if (correctTreeFromThisNode1!=null)
////                recalculateParentHaploAndCorrectContinuingHaploName(oldParentHaplo[srcHaplo], correctTreeFromThisNode1);
//                recalculateParentHaploAndCorrectContinuingHaploName(qsTree.getParentHaplo()[srcHaplo], correctTreeFromThisNode1);
//            // for the attached haplotype, need to find out new parent haplo and start at the changed node
//            // changed node can be newly attached parent or any node above
//            // getMaxPossibleHaploAttachTime function returned array of Node, time, haplo (above in the tree)
//            if (correctTreeFromThisNode2.getHaploAboveName()!=-1){
//                if ((int) haploStartMaxNewArray.get(2) != correctTreeFromThisNode2.getHaploAboveName())
//                    recalculateParentHaploAndCorrectContinuingHaploName((int) haploStartMaxNewArray.get(2), correctTreeFromThisNode2);
//                else
////                    recalculateParentHaploAndCorrectContinuingHaploName(oldParentHaplo[(int) haploStartMaxNewArray.get(2)], correctTreeFromThisNode2);
//                    recalculateParentHaploAndCorrectContinuingHaploName(qsTree.getParentHaplo()[(int) haploStartMaxNewArray.get(2)], correctTreeFromThisNode2);
//            }
//            // if the continuingHaplo is -1, we must have attached the haplo below the moved parent
//            // there can however still be some other haplo up in the tree (and not coming to this part)
//            else{
//                // we have to start correcting from the srcNode
//                int contHaplo = (int) haploStartMaxNewArray.get(2);
//                boolean contHaploHere = false;
//                for (Node node : correctTreeFromThisNode2.getAllLeafNodes()){
//                    if (node.getNr() == contHaplo){
//                        contHaploHere = true;
//                        break;
//                    }
//                }
//                if (contHaploHere)
//                    correctTreeFromThisNode2.setContinuingHaploName(contHaplo);
//                else
//                    correctTreeFromThisNode2.setContinuingHaploName(-1);
//                recalculateParentHaploAndCorrectContinuingHaploName(contHaplo, (QuasiSpeciesNode) srcNode);
//            }
//        }
        // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
        int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
        qsTree.setStartBranchCounts(startBranchCountsArray);

        // Incorporate probability of choosing current haplotype to move --- only with Felsenstein
//        logHastingsRatio += Math.log(possibleHaplo.size());
        if (tempqstimes.length > 1){
            // Incorporate probability of current positioning of the start time of the haplotype.
//            logHastingsRatio -= Math.log(haploStartMax-haploStartMin);
            // Incorporate probability of current positioning of the start time of the haplotype.
//            logHastingsRatio += Math.log(haploStartMaxNew-haploStartMin);
            // scaling the attachment times (see Scaling Operator)
            // assign contribution to the Hastings ratio for having different possible scales for told
            logHastingsRatio += Math.log(haploStartMaxNew/told - haploStartMin/told);
            logHastingsRatio -= Math.log(haploStartMax/tnew - haploStartMin/tnew);
            // assign contribution of each scaled attachment time
//            logHastingsRatio -= 2 * (Math.log(scalefactor * (haploStartMaxNew - haploStartMin) + haploStartMin) - Math.log(haploStartMaxNew));
//            logHastingsRatio += (tempqstimes.length-1) * Math.log(scalefactor);
            logHastingsRatio += (tempqstimes.length-3) * Math.log(scalefactor);
        }
        // Incorporate probability of choosing current haplotype to move back --- only with Felsenstein
        if (((QuasiSpeciesNode) srcNode).getContinuingHaploName() != srcHaplo){
            // check how many haplotypes can be pulled up
            ArrayList<Integer> backPossibleHaplo = new ArrayList<>();
//            int[] newParentHaplo = qsTree.getParentHaplo();
//            checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, newParentHaplo, backPossibleHaplo);
            checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, qsTree.getParentHaplo(), backPossibleHaplo);
            // account for this in the hastings ratio
//            logHastingsRatio -= Math.log(backPossibleHaplo.size());
        }
        // no need to be adding 0.0, so just commented out
        // else{
        //  logHastingsRatio -= Math.log(1);
        // }

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