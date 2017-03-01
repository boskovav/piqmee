package quasispeciestree.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesNode;
import beast.util.Randomizer;
import quasispeciestree.tree.QuasiSpeciesTree;

import java.util.ArrayList;


/**
 *  @author Veronika Boskova created on 09/06/2016 finished on 09/06/2016
 */
@Description("Chooses a haplotype at random and moves"
        + "its first attachment time uniformly in interval "
        + "restricted by the parent haplotype (start)"
        +" and tip time of the moved haplotype,"
        +" and scales all the attachment times accordingly.")
public class QuasiSpeciesHaplotypeScale extends QuasiSpeciesTreeOperator{

    /**
     * Change the start time and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {

//        qsTree.startEditing(this);

        if (qsTree.getTotalAttachmentCounts()==0){
            System.out.println("In QuasiSpeciesHaploScale operator --- "
                              +"there are no QS duplicates. The QuasiSpeciesHaplotypeScale "
                              +"operator cannot be used. Please remove it from your xml file.");
            System.exit(0);
        }

        // keep track of the Hastings ratio
        double logHastingsRatio = 0.0;
        // Randomly select event on tree:
        // weighted by the number of events (i.e. count of each haplotype??)
        // for now do random uniform from just the haplotype count (disregarding the counts)
        int haplo=-1;
        do {
            haplo = Randomizer.nextInt(qsTree.getLeafNodeCount());
        } while (qsTree.getAttachmentTimesList(haplo).length<2);

        QuasiSpeciesNode node = (QuasiSpeciesNode) qsTree.getNode(haplo);

        // if we did not assign the node at this stage, throw exception
        if (node == null)
            throw new IllegalStateException("Event selection loop fell through!");

        // get the attachment times array to be changed
        double[] tempqstimes=qsTree.getAttachmentTimesList(node).clone();

        // get a node above which the current haplotype arises
        QuasiSpeciesNode oldNodeBelowHaploMoved = findNodeBelowThisHaplo(node,haplo);

        // reposition the event (i.e. haplotype's first attachment time) uniformly at random between tmin and tmax
        // what is the parent first attachment time?
        ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTimeForQSStart(oldNodeBelowHaploMoved, haplo, 0);
        int haplotypesParentHaplo = (int) haploStartMaxNewArray.get(2);

        // get a random number deciding where the start point of the current haplo will be moved
        double u = Randomizer.nextDouble();

        // note down what needs to be found out to propose a new start time
        // and to reposition the event (i.e. haplotype's first attachment time) uniformly at random between tmin and tmax
        double scalefactor=0;
        double toldQSstart=tempqstimes[0];
        double tnewQSstart=0;

        double tmin = node.getHeight();
        double tmax;
        double told=0;
        //if (tempqstimes.length>1)
            told=tempqstimes[1];
        double tnew=0;

        // find out tmax
        tmax=(double)haploStartMaxNewArray.get(1);
        // reposition attachment times: attach ((time - tmin) * (tnew/told)) + tmin
        // scale all the other positions in the array but the 0 position (haplo start time)
        scalefactor = u*(tmin/told) + (1.0-u)*(tmax/told);
        // get new time to attach of first attachment time
        tnew = tempqstimes[1] * scalefactor;
        for (int i=1; i<tempqstimes.length; i++) {
            tempqstimes[i] = tempqstimes[i] * scalefactor;
        }

        if (tempqstimes[1]>tmax)
            return Double.NEGATIVE_INFINITY;

        // set the haplotype's starting time to the new time
//        double x = Randomizer.nextDouble();
//        tnewQSstart = x*tempqstimes[1] + (1.0-x)*tmax;
//        tempqstimes[0] = tnewQSstart;
        tnewQSstart = tempqstimes[1];
        tempqstimes[0] = tnewQSstart;

        // Reject invalid haplotype scalings:
        if (scalefactor<1.0 && tempqstimes.length>1){
            if (tempqstimes[0]<node.getHeight()){
                System.out.println("problem in hereeeeee you did not really scale the QS start apparently");
                System.exit(0);
            }
            else if (tempqstimes[tempqstimes.length-1]<node.getHeight())
                return Double.NEGATIVE_INFINITY;
        }

        // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                logHastingsRatio -= Math.log(tmax - told);
//                logHastingsRatio += Math.log(tmax - tnew);
        // assign contribution to the Hastings ratio for having different possible scales for told
        logHastingsRatio += Math.log(tmax/told - tmin/told);
        logHastingsRatio -= Math.log(tmax/tnew - tmin/tnew);
        // Incorporate the probability of scaling all the attachment times
        // assign contribution of each scaled attachment time
        logHastingsRatio += (tempqstimes.length-3) * Math.log(scalefactor);

        qsTree.setAttachmentTimesList(node, tempqstimes);

        // scale all QS start points for all haplotypes on a way from toldQSstart to tnewQSstart
        // find the new node above which the haplo arises
        QuasiSpeciesNode nodeBelowHaploMoved;
        if (tnewQSstart > toldQSstart) {
            nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(oldNodeBelowHaploMoved,tnewQSstart);
        }
        else {
            nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(node,tnewQSstart);
        }

        if (oldNodeBelowHaploMoved.getNr()!=nodeBelowHaploMoved.getNr()){
            // change internal node haploName - both remove and add if necessary
            // change the internal nodes' haploName, where necessary
            oldNodeBelowHaploMoved.setHaploAboveName(-1);
            nodeBelowHaploMoved.setHaploAboveName(node.getNr());
//            recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode)qsTree.getRoot());

            if (oldNodeBelowHaploMoved.getHeight() > nodeBelowHaploMoved.getHeight()){
                if (haplotypesParentHaplo == -1)
                    recalculateParentHaploAndCorrectContinuingHaploName(haplotypesParentHaplo, oldNodeBelowHaploMoved);
                else {
                    int parenthaplo = qsTree.getParentHaplo(haplotypesParentHaplo);
                    QuasiSpeciesNode oldNode = findNodeBelowThisHaplo(oldNodeBelowHaploMoved, haplotypesParentHaplo);
                    recalculateParentHaploAndCorrectContinuingHaploName(parenthaplo, oldNode);
                }
            }
            else
                recalculateParentHaploAndCorrectContinuingHaploName(haplotypesParentHaplo, nodeBelowHaploMoved);
        }

        // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
        int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
        qsTree.setStartBranchCounts(startBranchCountsArray);

        // Ensure BEAST knows to recalculate affected likelihood:
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);

    // RETURN log(HASTINGS RATIO)
    return logHastingsRatio; // proper hastings ratio!!!
    }

}
