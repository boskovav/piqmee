package quasispeciestree.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;
import beast.util.Randomizer;

import java.util.ArrayList;

/**
 *  @author Veronika Boskova created on 07/08/2015 finished on 25/05/2016
 */

@Description("Within given haplotype, randomly selects one sequence "
        + "attachment time, selects a new interval, restricted by the "
        + "first branching event of haplotype above and haplotype tip times, where this attachment "
        + "time will be added to and attaches it uniformly in that interval.")
public class QuasiSpeciesSequenceAttachmentRandom extends QuasiSpeciesTreeOperator{

    /**
     * Change the attachment time and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {

        if (qsTree.getTotalAttachmentCounts()==0){
            System.out.println("In QuasiSpeciesSequenceAttachmentRandom operator --- "
                    +"there are no QS duplicates. The QuasiSpeciesSequenceAttachmentRandom "
                    +"operator cannot be used. Please remove it from your xml file.");
            System.exit(0);
        }

        double logHastingsRatio = 0.0;

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
                // NOTE: haplotype starting time changes by another operator QuasiSpeciesHaplotypeStartRandom
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
        double toldQSstart = tempqstimes[0];
        double tnewQSstart = -1;

        // choose new max index between 0 - #sequences of this haplotype
        // here we do not allow the event to be repositioned to the interval delimited by itself
            do { tmaxIdx = Randomizer.nextInt(tempqstimes.length); }
        while (tmaxIdx == changeIdx);

//        ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTime(node, node.getNr(), 0);
        ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTimeForQSStart(node, node.getNr(), 0);

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

        logHastingsRatio += Math.log(tmax - tmin);

        if (changeIdx+1 == tempqstimes.length && changeIdx == 1)
            logHastingsRatio -= Math.log((double) haploStartMaxNewArray.get(1) - node.getHeight());
        else if (changeIdx+1 == tempqstimes.length)
            logHastingsRatio -= Math.log(tempqstimes[changeIdx-1] - node.getHeight());
        else if (changeIdx == 1)
            logHastingsRatio -= Math.log((double) haploStartMaxNewArray.get(1) - tempqstimes[changeIdx+1]);
        else
            logHastingsRatio -= Math.log(tempqstimes[changeIdx-1] - tempqstimes[changeIdx+1]);

        double u = Randomizer.nextDouble();
        double tnew = u*tmin + (1-u)*tmax; // = u*(tmin-tmax)+tmax =u*(tmin-(tmin+tmax-tmin))+tmin+tmax-tmin
        // invert u and 1-u and have the same as u(tmax-tmin)+tmin
        // (1-u)*tmin + u*tmax-u*(tmax-tmin)-tmin=0 ???
        // indeed tmin-u*tmin+u*tmax-u*tmax+u*tmin-tmin=0

        // if we moved the first attachment time, then move also the QS start appropriately
        // generic move of the 1st attachment time
        if (changeIdx==1){
            double v = Randomizer.nextDouble();
            if (tmaxIdx == 0){
                tnewQSstart = (v * tnew) + ((1 - v) * (double) haploStartMaxNewArray.get(1));
                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                logHastingsRatio -= Math.log((double) haploStartMaxNewArray.get(1) - tempqstimes[1]);
//                logHastingsRatio += Math.log((double) haploStartMaxNewArray.get(1) - tnew);
            }
            else{
                tnewQSstart = (v * tempqstimes[2]) + ((1 - v) * (double) haploStartMaxNewArray.get(1));
                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                logHastingsRatio -= Math.log((double) haploStartMaxNewArray.get(1) - tempqstimes[1]);
//                logHastingsRatio += Math.log((double) haploStartMaxNewArray.get(1) - tempqstimes[2]);
            }
            tempqstimes[0] = tnewQSstart;
        }
        // creation of new 1st attachment time
        else if (tmaxIdx==0){
            double v = Randomizer.nextDouble();
            tnewQSstart = v*tnew + (1-v)*tmax;
            tempqstimes[0] = tnewQSstart;
            // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//            logHastingsRatio -= Math.log(tmax - tempqstimes[1]);
//            logHastingsRatio += Math.log(tmax - tnew);
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

        // account for the fact that we are changing the QS start
        if (tnewQSstart != -1){
            // get a node above which the current haplotype arises
            QuasiSpeciesNode oldNodeBelowHaploMoved = findNodeBelowThisHaplo(node,node.getNr());
            // scale all QS start points for all haplotypes on a way from toldQSstart to tnewQSstart
            // find out what nodes have been passed by changing the haplo height
//            ArrayList<QuasiSpeciesNode> nodesPassed = new ArrayList<>(qsTree.getInternalNodeCount());
            // find the new node above which the haplo arises
            QuasiSpeciesNode nodeBelowHaploMoved;
            if (tnewQSstart > toldQSstart) {
//                nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(oldNodeBelowHaploMoved,tnewQSstart,nodesPassed);
                nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(oldNodeBelowHaploMoved,tnewQSstart);
            }
            else {
                // first find the new node below
                nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(node,tnewQSstart);
                // then check out which nodes were passed
//                QuasiSpeciesNode intnodehelper=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(nodeBelowHaploMoved,toldQSstart,nodesPassed);
            }
//            nodesPassed.trimToSize();

            if (oldNodeBelowHaploMoved.getNr()!=nodeBelowHaploMoved.getNr()){
//                // reposition QS start points of the QS on the way and add contribution to the HR for QS start moves
//                int[] parentHaploArray = qsTree.getParentHaplo();
//                // calculate up to where the haplotypes from the passed nodes can max attach to
//                // also get the node up to which they can go to
//                ArrayList upToWhere;
//
//                if (oldNodeBelowHaploMoved.getHeight() > nodeBelowHaploMoved.getHeight())
//                    upToWhere = getMaxPossibleHaploAttachTimeForQSStart(oldNodeBelowHaploMoved, node.getNr(), 0);
//                else if (nodeBelowHaploMoved.getContinuingHaploName() != -1)
//                    upToWhere = getMaxPossibleHaploAttachTimeForQSStart(nodeBelowHaploMoved, nodeBelowHaploMoved.getContinuingHaploName(), 0);
//                else
//                // here the haplo is a dummy thing... we just need to find out up to where the QS start from the passed haplotypes can move
//                    upToWhere = getMaxPossibleHaploAttachTimeForQSStart(nodeBelowHaploMoved, node.getNr(), 0);
//
//                logHastingsRatio += moveQSInWayAfterMovingHaplotype((double)upToWhere.get(1),tmin,tnewQSstart,toldQSstart,nodesPassed,parentHaploArray);
//                if (oldNodeBelowHaploMoved.getHaploAboveName()==node.getNr())
                    oldNodeBelowHaploMoved.setHaploAboveName(-1);
                nodeBelowHaploMoved.setHaploAboveName(node.getNr());

                if (oldNodeBelowHaploMoved.getHeight() > nodeBelowHaploMoved.getHeight()){
                    if ((int) haploStartMaxNewArray.get(2) == -1)
                        recalculateParentHaploAndCorrectContinuingHaploName((int) haploStartMaxNewArray.get(2), oldNodeBelowHaploMoved);
                    else {
                        int parenthaplo = qsTree.getParentHaplo((int) haploStartMaxNewArray.get(2));
                        QuasiSpeciesNode oldNode = findNodeBelowThisHaplo(oldNodeBelowHaploMoved,(int) haploStartMaxNewArray.get(2));
                        recalculateParentHaploAndCorrectContinuingHaploName(parenthaplo, oldNode);
                    }
                }
                else
                    recalculateParentHaploAndCorrectContinuingHaploName((int) haploStartMaxNewArray.get(2), nodeBelowHaploMoved);

//                if ((double) upToWhere.get(1)==origin.getValue()){
//                    recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) qsTree.getRoot());
//                }
//                else if (((QuasiSpeciesNode) upToWhere.get(0)).getHaploAboveName() != -1)
//                    recalculateParentHaploAndCorrectContinuingHaploName(parentHaploArray[(int) upToWhere.get(2)], (QuasiSpeciesNode) upToWhere.get(0));
//                else {
//                    QuasiSpeciesNode nodetemp = findNodeBelowThisHaplo((QuasiSpeciesNode) upToWhere.get(0),(int) upToWhere.get(2));
//                    recalculateParentHaploAndCorrectContinuingHaploName((int) upToWhere.get(2), nodetemp);
//                }

                // mark passed nodes as dirty
//                if (tnewQSstart>toldQSstart){
//                    for (Node nodePassed : nodesPassed)
//                        nodePassed.makeDirty(QuasiSpeciesTree.IS_FILTHY);
//                }
            }
        }

        // recalculate countPossibleStartBranches
        int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
        qsTree.setStartBranchCounts(startBranchCountsArray);

        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);

        return logHastingsRatio;
    }
}