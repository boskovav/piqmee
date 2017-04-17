package quasispeciestree.operators;

import beast.core.Description;
import beast.util.Randomizer;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;

import java.util.ArrayList;


/**
 *  @author Veronika Boskova created on 09/06/2016 finished on 27/03/2017
 */
@Description("Chooses a haplotype at random and moves"
        + "its first attachment time uniformly in interval "
        + "restricted by the parent haplotype (start)"
        +" and tip time of the moved haplotype,"
        +" and scales all the attachment times accordingly.")
public class QuasiSpeciesHaplotypeScale extends QuasiSpeciesTreeOperator{

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        if (qsTree.getTotalAttachmentCounts()==0){
            throw new IllegalArgumentException("In QuasiSpeciesHaploScale operator --- "
                    +"there are no QS duplicates. The QuasiSpeciesHaplotypeScale "
                    +"operator cannot be used. Please remove it from your xml file.");
        }
    }

    /**
     * Change the start time and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {

        // keep track of the Hastings ratio
        double logHastingsRatio = 0.0;
        // Randomly select event on tree:
        // weighted by the number of events (i.e. count of each haplotype??)
        // for now do random uniform from just the haplotype count (disregarding the counts)
        QuasiSpeciesNode node = null;
        do {
            node = (QuasiSpeciesNode) qsTree.getNode(Randomizer.nextInt(qsTree.getLeafNodeCount()));

        } while (node.getAttachmentTimesList().length < 2);

        // if we did not assign the node at this stage, throw exception
        if (node == null)
            throw new IllegalStateException("Event selection loop fell through!");

        int haplo = node.getNr();

        // get the attachment times array to be changed
        double[] tempqstimes = node.getAttachmentTimesList().clone();
        // get also tip times to help define max/min scalings
        double[] temptiptimes = node.getTipTimesList();
        int[] temptiptimescount = node.getTipTimesCountList();

        // get a node above which the current haplotype arises
        QuasiSpeciesNode oldNodeBelowHaploMoved = findNodeBelowThisHaplo(node,haplo);

        // reposition the event (i.e. haplotype's first attachment time) uniformly at random between tmin and tmax
        // what is the parent first attachment time?
        ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTime(oldNodeBelowHaploMoved, haplo, 0);
        int parentHaplo = (int) haploStartMaxNewArray.get(2);

        // get a random number deciding where the start point of the current haplo will be moved
        double u = Randomizer.nextDouble();

        // note down what needs to be found out to propose a new start time
        // and to reposition the event (i.e. haplotype's first attachment time) uniformly at random between tmin and tmax
        double scalefactor = 0;
        double toldQSstart = tempqstimes[0];
        double tnewQSstart = 0;

        // find out tmin/tmax/told
        double tmax = (double)haploStartMaxNewArray.get(1);
        double tminold = node.getHeight();
        double tminnew = node.getHeight();
        double toldtop = tempqstimes[1];
        double tnewtop = 0;
        double toldbottom = tempqstimes[tempqstimes.length-1];
        double tnewbottom = 0;
        // to get the tmin check for each sampling time of the haplo that the last possible attach time is above the
        //  the corresponding sampling time
        int currentPosition = tempqstimes.length-1-(temptiptimescount[0]-1);
        for (int i = 1; i < temptiptimes.length; i++){
            if (tminold/toldbottom < temptiptimes[i]/tempqstimes[currentPosition]) {
                tminold = temptiptimes[i];
                toldbottom = tempqstimes[currentPosition];
            }
            currentPosition -= temptiptimescount[i];
        }
        // reposition attachment times: attach ((time - tmin) * (tnew/told)) + tmin
        // scale all the other positions in the array but the 0 position (haplo start time)
        scalefactor = u*(tminold/toldbottom) + (1.0-u)*(tmax/toldtop);
        for (int i = 1; i < tempqstimes.length; i++) {
            tempqstimes[i] = tempqstimes[i] * scalefactor;
        }
        // get new time to attach of first attachment time
        tnewtop = tempqstimes[1];
        // set the haplotype's starting time to the new time
        tnewQSstart = tempqstimes[1];
        tempqstimes[0] = tnewQSstart;

        // check that the scaling was ok
        if (tmax/toldtop < tminold/toldbottom){
            throw new IllegalStateException("problem in hereeeeee: HaplotypeScale: The scaled values are not calculated properly?");
        }
        // check that the newly scaled attachment times do not go over the boundaries
        if (tnewQSstart > tmax)
            return Double.NEGATIVE_INFINITY;

        //get the possible back scale interval
        tnewbottom = tempqstimes[tempqstimes.length-1];
        currentPosition = tempqstimes.length-1-(temptiptimescount[0]-1);
        for (int i = 1; i < temptiptimes.length; i++){
            if(temptiptimes[i] > tempqstimes[currentPosition])
                return Double.NEGATIVE_INFINITY;
            if (tminnew/tnewbottom < temptiptimes[i]/tempqstimes[currentPosition]) {
                tminnew = temptiptimes[i];
                tnewbottom = tempqstimes[currentPosition];
            }
            currentPosition -= temptiptimescount[i];
        }

        // check that the scaling was ok
        if (tmax/tnewtop < tminnew/tnewbottom){
            throw new IllegalStateException("problem in hereeeeee: HaplotypeScale: The scaled values are not calculated properly?");
        }

        // Reject invalid haplotype scalings:
        if (scalefactor<1.0){
            if (tnewQSstart < node.getHeight()){
                throw new IllegalStateException("problem in hereeeeee you did not really scale the QS start apparently");
            }
            else if (tempqstimes[tempqstimes.length-1] < node.getHeight())
                return Double.NEGATIVE_INFINITY;
        }

        // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                logHastingsRatio -= Math.log(tmax - toldtop);
//                logHastingsRatio += Math.log(tmax - tnewtop);
        // assign contribution to the Hastings ratio for having different possible scales for toldtop
        logHastingsRatio += Math.log(tmax/toldtop - tminold/toldbottom);
        logHastingsRatio -= Math.log(tmax/tnewtop - tminnew/tnewbottom);
        // Incorporate the probability of scaling all the attachment times
        // assign contribution of each scaled attachment time
        logHastingsRatio += (tempqstimes.length - 3) * Math.log(scalefactor);

        node.setAttachmentTimesList(tempqstimes);

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
            //recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode)qsTree.getRoot());

            if (oldNodeBelowHaploMoved.getHeight() > nodeBelowHaploMoved.getHeight()){
                if (parentHaplo == -1)
                    recalculateParentHaploAndCorrectContinuingHaploName(parentHaplo, oldNodeBelowHaploMoved);
                else {
                    int grandParentHaplo = ((QuasiSpeciesNode) qsTree.getNode(parentHaplo)).getParentHaplo();
                    QuasiSpeciesNode oldNode = findNodeBelowThisHaplo(oldNodeBelowHaploMoved, parentHaplo);
                    recalculateParentHaploAndCorrectContinuingHaploName(grandParentHaplo, oldNode);
                }
            }
            else
                recalculateParentHaploAndCorrectContinuingHaploName(parentHaplo, nodeBelowHaploMoved);
        }

        // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
        qsTree.countAndSetPossibleStartBranches();

        // Ensure BEAST knows to recalculate affected likelihood:
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);

    // RETURN log(HASTINGS RATIO)
    return logHastingsRatio; // proper hastings ratio!!!
    }

}
