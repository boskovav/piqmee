package quasispeciestree.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;
import beast.util.Randomizer;

import java.util.ArrayList;

/**
 *  @author Veronika Boskova created on 27/07/2015
 */
@Description("Within given haplotype, randomly selects one sequence "
        + "attachment time and moves it uniformly in interval "
        + "restricted by the closest previous and next attachment times "
        + "of another sequence from the same haplotype.")
public class QuasiSpeciesSequenceAttachmentUniform extends QuasiSpeciesTreeOperator{

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        if (qsTree.getTotalAttachmentCounts() == 0){
            throw new IllegalArgumentException("In QuasiSpeciesSequenceAttachmentRandom operator --- "
                    + "there are no QS duplicates. The QuasiSpeciesSequenceAttachmentRandom "
                    + "operator cannot be used. Please remove it from your xml file.");
        }
    }

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
            int tempHaploCount = qsTree.getHaplotypeCounts(thisNode) - 1;
            if (event < tempHaploCount) {
                node = (QuasiSpeciesNode) thisNode;
                // change index is event +1 as our arrays are 1-#haplotype repetition instances
                // and position 0 in the array is the haplotype starting point
                // NOTE: haplotype starting time changes by another operator QuasiSpeciesHaplotypeStartRandom
                changeIdx = event + 1;
                break;
            }
            event -= tempHaploCount;
        }

        // if we did not assign the node at this stage, throw exception
        if (node == null)
            throw new IllegalStateException("Event selection loop fell through!");

        // reposition the event (i.e. haplotype sequence changeIdx attachment time)
        double tmin, tmax;
        double[] tempqstimes = node.getAttachmentTimesList().clone();
        double[] temptiptimes = node.getTipTimesList();
        int[] temptiptimescount = node.getTipTimesCountList();

        // as we choose index changeIdx between 1 - #sequences of this haplotype then changeIdx-1=0 at minimum
        // tmax can got up to the next real node or next attach time
        ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTime(node, node.getNr(), 0);
        if (changeIdx == 1){
            tmax = (double) haploStartMaxNewArray.get(1);
        }
        else
            tmax = tempqstimes[changeIdx-1];

        // tmin can go down to the previous attach time or is delimited by the tip sampling times
        if (changeIdx + 1 == tempqstimes.length)
            tmin = node.getHeight();
        else
            tmin = tempqstimes[changeIdx + 1];
        // check if changeIdx points to the attach time that needs to be above some tip time
        int currentPosition = tempqstimes.length-1-(temptiptimescount[0]-1);
        int tminIdxHard = tempqstimes.length;
        double tminHard = temptiptimes[0];
        for (int i = 1; i < temptiptimes.length; i++){
            if (changeIdx < currentPosition)
                break;
            else if (tempqstimes[currentPosition] < temptiptimes[i])
                throw new IllegalStateException("QuasiSpeciesSequenceAttachmentRandom: seems like the attachment time is below its tip.");
            else if (tempqstimes[currentPosition + 1] < temptiptimes[i]) {
                tminIdxHard = currentPosition;
                tminHard = temptiptimes[i];
            }
            currentPosition -= temptiptimescount[i];
        }
        if (changeIdx == tminIdxHard) {
            if (tmin < tminHard)
                tmin = tminHard;
        }

        double u = Randomizer.nextDouble();
        double tnew = u*tmin + (1-u)*tmax; // = u*(tmin-tmax)+tmax
                                            // invert u and 1-u and have the same as u(tmax-tmin)+tmin
                                            // (1-u)*tmin + u*tmax-u*(tmax-tmin)-tmin=0 ???
                                            // indeed tmin-u*tmin+u*tmax-u*tmax+u*tmin-tmin=0

        double told = tempqstimes[changeIdx];
        tempqstimes[changeIdx] = tnew;
        if (changeIdx == 1){
            tempqstimes[0] = tnew;
        }
        node.setAttachmentTimesList(tempqstimes);

        // account for the fact that we are changing the QS start
        if (changeIdx == 1){
            // get a node above which the current haplotype arises
            QuasiSpeciesNode oldNodeBelowHaploMoved = findNodeBelowThisHaplo(node,node.getNr());
            // find the new node above which the haplo arises
            QuasiSpeciesNode nodeBelowHaploMoved;
            if ( tnew > told) {
                nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(oldNodeBelowHaploMoved,tnew);
            }
            else {
                // first find the new node below
                nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(node,tnew);
            }

            if (oldNodeBelowHaploMoved.getNr()!=nodeBelowHaploMoved.getNr()){
                oldNodeBelowHaploMoved.setHaploAboveName(-1);
                nodeBelowHaploMoved.setHaploAboveName(node.getNr());

                if (oldNodeBelowHaploMoved.getHeight() > nodeBelowHaploMoved.getHeight()){
                    if ((int) haploStartMaxNewArray.get(2) == -1)
                        recalculateParentHaploAndCorrectContinuingHaploName((int) haploStartMaxNewArray.get(2), oldNodeBelowHaploMoved);
                    else {
                        int parenthaplo = ((QuasiSpeciesNode) qsTree.getNode((int) haploStartMaxNewArray.get(2))).getParentHaplo();
                        QuasiSpeciesNode oldNode = findNodeBelowThisHaplo(oldNodeBelowHaploMoved,(int) haploStartMaxNewArray.get(2));
                        recalculateParentHaploAndCorrectContinuingHaploName(parenthaplo, oldNode);
                    }
                }
                else
                    recalculateParentHaploAndCorrectContinuingHaploName((int) haploStartMaxNewArray.get(2), nodeBelowHaploMoved);
            }
        }

        // recalculate countPossibleStartBranches
        qsTree.countAndSetPossibleStartBranches();

        // Ensure BEAST knows to recalculate affected likelihood:
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);

        return 0.0;
    }
}
