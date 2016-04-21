package quasispeciestree.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;
import beast.util.Randomizer;

import java.util.ArrayList;

/**
 *  @author Veronika Boskova created on 07/08/2015 finished on 20/04/2016
 */

@Description("Within given haplotype, randomly selects one sequence "
        + "attachment time, selects a new interval, restricted by the "
        + "clone start and clone tip times, where this attachment "
        + "time will be added to and attaches it uniformly in that interval.")
public class QuasiSpeciesSequenceAttachmentRandom extends QuasiSpeciesTreeOperator{

    /**
     * Change the attachment time and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {
        // Randomly select event on tree:
        // weighted by the number of events (i.e. count of each haplotype)
        int event = Randomizer.nextInt(qsTree.getTotalAttachmentCounts());


        QuasiSpeciesNode node = null;

        // index for the attachment time chosen to change
        int changeIdx = -1;

        // find the haplotype and the sequence position corresponding to the event number
        // backward search...
        for (Node thisNode : qsTree.getExternalNodes()) {
            double tempHaploCount = qsTree.getHaplotypeCounts((QuasiSpeciesNode) thisNode);
            if (event<tempHaploCount) {
                node = (QuasiSpeciesNode)thisNode;
                // change index is event +1 as our arrays are 1-#haplotype repetition instances
                // and position 0 in the array is the haplotype starting point
                // TODO haplotype starting time changes by another operator...
                changeIdx = event+1;
                break;
            }
            event -= tempHaploCount;
        }


        // if we did not assign the node at this stage, throw exception
        if (node == null)
            throw new IllegalStateException("Event selection loop fell through!");

        // reposition the event (i.e. haplotype sequence changeIdx attachment time)
        double tmin, tmax;
        int tminIdx, tmaxIdx;
        Double[] tempqstimes=qsTree.getAttachmentTimesList(node).clone();

        // choose new max index between 0 - #sequences of this haplotype
        // here we do not allow the event to be repositioned to the interval delimited by itself
            do { tmaxIdx = Randomizer.nextInt(tempqstimes.length); }
        while (tmaxIdx == changeIdx);

        ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTime(node, node.getNr(), 0);

        if (tmaxIdx==0)
            tmax = (double) haploStartMaxNewArray.get(1);
        else
            tmax = tempqstimes[tmaxIdx];

        if (tmaxIdx == changeIdx-1)
            tminIdx = tmaxIdx + 2;
        else
            tminIdx = tmaxIdx + 1;

        if (tminIdx == tempqstimes.length)
            tmin = node.getHeight();
        else
            tmin = tempqstimes[tminIdx];

        double newRange = Math.log(tmax - tmin);
        double oldRange;
        if (changeIdx+1 == tempqstimes.length)
            oldRange = Math.log(tempqstimes[changeIdx-1] - node.getHeight());
        else if (changeIdx == 1)
            oldRange = Math.log((double) haploStartMaxNewArray.get(1) - tempqstimes[changeIdx+1]);
        else
            oldRange = Math.log(tempqstimes[changeIdx-1] - tempqstimes[changeIdx+1]);

        double u = Randomizer.nextDouble();
        double tnew = u*tmin + (1-u)*tmax; // = u*(tmin-tmax)+tmax =u*(tmin-(tmin+tmax-tmin))+tmin+tmax-tmin
        // invert u and 1-u and have the same as u(tmax-tmin)+tmin
        // (1-u)*tmin + u*tmax-u*(tmax-tmin)-tmin=0 ???
        // indeed tmin-u*tmin+u*tmax-u*tmax+u*tmin-tmin=0

        // if we moved the first attachment time, then move also the QS start appropriately
        // generic move of the 1st attachment time
        double tnewQSstart = -1;
        if (changeIdx==1){
            double v = Randomizer.nextDouble();
            if (tmaxIdx == 0){
                tnewQSstart = (v * tnew) + ((1 - v) * (double) haploStartMaxNewArray.get(1));
                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                oldRange += Math.log((double) haploStartMaxNewArray.get(1) - tempqstimes[1]);
//                newRange += Math.log((double) haploStartMaxNewArray.get(1) - tnew);
            }
            else{
                tnewQSstart = (v * tempqstimes[2]) + ((1 - v) * (double) haploStartMaxNewArray.get(1));
                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                oldRange += Math.log((double) haploStartMaxNewArray.get(1) - tempqstimes[1]);
//                newRange += Math.log((double) haploStartMaxNewArray.get(1) - tempqstimes[2]);
            }
            tempqstimes[0] = tnewQSstart;
        }
        // creation of new 1st attachment time
        else if (tmaxIdx==0){
            double v = Randomizer.nextDouble();
            tnewQSstart = v*tnew + (1-v)*tmax;
            tempqstimes[0] = tnewQSstart;
            // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//            oldRange += Math.log(tmax - tempqstimes[1]);
//            newRange += Math.log(tmax - tnew);
        }

        tempqstimes[changeIdx]=tnew;
        if (changeIdx > tmaxIdx){
            for (int i=changeIdx; i>tmaxIdx+1; i--){
                double swaptime = tempqstimes[i];
                tempqstimes[i] = tempqstimes[i-1];
                tempqstimes[i-1] = swaptime;
            }
        } else {
            for (int i=changeIdx; i<tmaxIdx; i++){
                double swaptime = tempqstimes[i];
                tempqstimes[i] = tempqstimes[i+1];
                tempqstimes[i+1] = swaptime;
            }
        }
        qsTree.setAttachmentTimesList(node, tempqstimes);

        if (tnewQSstart != -1){
            QuasiSpeciesNode nodeBelowHaploMoved = node;
            while (nodeBelowHaploMoved != qsTree.getRoot() && tnewQSstart > nodeBelowHaploMoved.getParent().getHeight()){
                nodeBelowHaploMoved = (QuasiSpeciesNode) nodeBelowHaploMoved.getParent();
            }
            QuasiSpeciesNode oldNodeBelowHaploMoved = node;
            while (oldNodeBelowHaploMoved.getHaploAboveName()!=node.getNr()){
                oldNodeBelowHaploMoved = (QuasiSpeciesNode) oldNodeBelowHaploMoved.getParent();
            }
            oldNodeBelowHaploMoved.setHaploAboveName(-1);
            nodeBelowHaploMoved.setHaploAboveName(node.getNr());

            if (oldNodeBelowHaploMoved.getHeight() > nodeBelowHaploMoved.getHeight())
                recalculateParentHaploAndCorrectContinuingHaploName((int) haploStartMaxNewArray.get(2), oldNodeBelowHaploMoved);
            else
                recalculateParentHaploAndCorrectContinuingHaploName((int) haploStartMaxNewArray.get(2), nodeBelowHaploMoved);
        }

        // recalculate countPossibleStartBranches
        int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
        qsTree.setStartBranchCounts(startBranchCountsArray);

        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);

        return newRange-oldRange;

    }

}
