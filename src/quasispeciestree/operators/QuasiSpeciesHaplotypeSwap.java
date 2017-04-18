package quasispeciestree.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;

import java.util.ArrayList;


/**
 *  @author Veronika Boskova created on 06/08/2015 finished on 27/03/2017
 */
@Description("Chooses a haplotype at random and moves"
        +" its first attachment time uniformly in interval"
        +" restricted by the grand-parent haplotype (start)"
        +" and a time of a common ancestor with the parent haplotype,"
        +" and scales all the current/parent haplo attachment"
        +" times accordingly."
        +" This move is performed ONLY if there is another haplotype"
        +" above the currently chosen one that we can SWAP with.")
public class QuasiSpeciesHaplotypeSwap extends QuasiSpeciesTreeOperator{

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        // check if this operator can be used at all?
        int canbeused = 0;
        for (Node node : qsTree.getExternalNodes()){
            if (qsTree.getHaplotypeCounts(node) > 0)
                canbeused++;
        }
        if (canbeused<1){
            throw new IllegalArgumentException("QuasiSpeciesHaplotypeSwap operator cannot be " +
                    "used since there are no or just one haplotype with duplicates! Remove this " +
                    "operator from your XML file.");
        }
    }

    /**
     * Change the start time and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {

        // check that this operator can actually perform a move
        int haplowithparents = 0;
        for (int i = 0; i < qsTree.getLeafNodeCount(); i++){
            if (((QuasiSpeciesNode)qsTree.getNode(i)).getParentHaplo() != -1)
                haplowithparents += 1;
        }
        if (haplowithparents == 0)
            return Double.NEGATIVE_INFINITY;

        // keep track of the Hastings ratio
        double logHastingsRatio = 0.0;
        // Randomly select event on tree:
        // for now choose uniformly at random from just the unique haplotype count (disregarding the counts)
        QuasiSpeciesNode node = null;
        do {
            node = (QuasiSpeciesNode) qsTree.getNode(Randomizer.nextInt(qsTree.getLeafNodeCount()));
        } while (node.getAttachmentTimesList().length < 2 || node.getParentHaplo() == -1);

        // if we did not assign the node at this stage, throw exception
        if (node == null)
            throw new IllegalStateException("Event selection loop fell through!");

        int haplo = node.getNr();

        // get the parentHaploArray and the attachment times array to be changed
        QuasiSpeciesNode parentHaplo = (QuasiSpeciesNode) qsTree.getNode(node.getParentHaplo());
        QuasiSpeciesNode grandParentHaplo = null;
        if (parentHaplo.getParentHaplo() != -1)
            grandParentHaplo = (QuasiSpeciesNode) qsTree.getNode(parentHaplo.getParentHaplo());;

        // get the attachment times for tip
        double[] tempqstimes = node.getAttachmentTimesList().clone();
        // get also tip times to help define max/min scalings
        double[] temptiptimes = node.getTipTimesList();
        int[] temptiptimescount = node.getTipTimesCountList();

        // get the attachment times for tip's parent
        double[] tempqstimesparent = parentHaplo.getAttachmentTimesList().clone();
        // get also tip times to help define maxparent/minparent scalings
        double[] temptiptimesparent = parentHaplo.getTipTimesList();
        int[] temptiptimescountparent = parentHaplo.getTipTimesCountList();

        // get a node above which the current/parent/grand parent haplotype arises
        QuasiSpeciesNode oldNodeBelowCurrentHaplo = findNodeBelowThisHaplo(node,haplo);
        QuasiSpeciesNode oldNodeBelowParentHaplo = findNodeBelowThisHaplo(node,parentHaplo.getNr());

        // reposition the event (i.e. haplotype's first attachment time) uniformly at random between tmin and tmax
        // what is the parent haplo start time?
        ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTime(oldNodeBelowCurrentHaplo, haplo, 0);


        // CURRENT
        // note down what needs to be found out to propose a new start time
        // and to reposition the event (i.e. haplotype's first attachment time) uniformly at random between tmin and tmax
        double scalefactor = 0;
        double toldQSstart = tempqstimes[0];
        double tnewQSstart = 0;
        double tmax;
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
            currentPosition = currentPosition - temptiptimescount[i];
        }
        // PARENT
        // get the haplotype's starting time
        double parentscalefactor = 0;
        double toldQSstartparent = tempqstimesparent[0];
        double tnewQSstartparent = 0;
        double tmaxparent = (double) haploStartMaxNewArray.get(1);
        double tminparentold = parentHaplo.getHeight();
        double tminparentnew = parentHaplo.getHeight();
        double toldparenttop = tempqstimesparent[1];
        double tnewparenttop = 0;
        double toldparentbottom = tempqstimesparent[tempqstimesparent.length-1];
        double tnewparentbottom = 0;
        // to get the tminparent check for each sampling time of the haplo that the last possible attach time is above the
        //  the corresponding sampling time
        currentPosition = tempqstimesparent.length-1-(temptiptimescountparent[0]-1);
        for (int i = 1; i < temptiptimesparent.length; i++){
            if (tminparentold/toldparentbottom < temptiptimesparent[i]/tempqstimesparent[currentPosition]) {
                tminparentold = temptiptimesparent[i];
                toldparentbottom = tempqstimesparent[currentPosition];
            }
            currentPosition -= temptiptimescountparent[i];
        }
        // check that the scaling was ok
        if (tmaxparent/toldparenttop < tminparentold/toldparentbottom){
            // in this case, it is possible that there is no acceptable scaling for the parent haplo
            // but then this move cannot be performed at the moment
            throw new IllegalStateException("problem in hereeeeee: HaplotypeSwap: The scaled values are not calculated properly?");
        }

// NO GRANDPARENT
        // find out tmax: grandparent haplotype is null, then tmax=origin
        if (grandParentHaplo == null){
            tmax = origin.getValue();
        }
// GRANDPARENT
        // otherwise tmax is the common ancestor of grandParentHaplo and parentHaplo
        else{
            // get a node above which the grandparent haplotype arises
            tmax = (double) getMaxPossibleHaploAttachTime(
                    (QuasiSpeciesNode)haploStartMaxNewArray.get(0),parentHaplo.getNr(),0).get(1);
        }
// PROPOSE NEW TIMES FOR CURRENT HAPLO
        // get a random number deciding where the start point of the current haplo will be moved
        double u = Randomizer.nextDouble();
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
            throw new IllegalStateException("problem in hereeeeee: HaplotypeSwap: The scaled values are not calculated properly?");

        }
        // check that the newly scaled attachment times do not go over the boundaries
        if (tnewQSstart > tmax || tnewQSstart < tmaxparent)
            return Double.NEGATIVE_INFINITY;
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
            throw new IllegalStateException("problem in hereeeeee: HaplotypeSwap: The scaled values are not calculated properly?");

        }
        // Reject invalid haplotype scalings:
        if (scalefactor < 1.0){
            if (tnewQSstart < node.getHeight()){
                throw new IllegalStateException("problem in hereeeeee you did not really scale the QS start apparently");
            }
            else if (tempqstimes[tempqstimes.length-1] < node.getHeight())
                return Double.NEGATIVE_INFINITY;
        }

// PROPOSE NEW TIMES FOR PARENT HAPLO
        // choose new time to start for the parent node uniformly at random between tminparent and tmaxparent
        //  i.e. below the common ancestor node and scale all attachment times of this haplo
        double v = Randomizer.nextDouble();
        // reposition parent haplo attachment times: attach ((time - tminparent) * (tnewparent/toldparent)) + tminparent
        // scale all the other positions in the array but the 0 position (haplo start time)
        parentscalefactor = v*(tminparentold/toldparentbottom) + (1.0-v)*(tmaxparent/toldparenttop);
        for (int i = 1; i < tempqstimesparent.length; i++) {
            tempqstimesparent[i] = tempqstimesparent[i] * parentscalefactor;
        }
        // get new time to attach of first attachment time
        tnewparenttop = tempqstimesparent[1];
        // set the haplotype's starting time to the new time
        tnewQSstartparent = tempqstimesparent[1];
        tempqstimesparent[0] = tnewQSstartparent;

        // check that the scaling was ok
        if (tmax/toldparenttop < tminparentold/toldparentbottom){
            throw new IllegalStateException("problem in hereeeeee: HaplotypeSwap: The scaled values are not calculated properly?");
        }
        // check that the newly scaled attachment times do not go over the boundaries
        if (tnewQSstartparent > tmaxparent || tnewQSstartparent > tmax)
            return Double.NEGATIVE_INFINITY;
        tnewparentbottom = tempqstimesparent[tempqstimesparent.length-1];
        currentPosition = tempqstimesparent.length-1-(temptiptimescountparent[0]-1);
        for (int i = 1; i < temptiptimesparent.length; i++){
            if(temptiptimesparent[i] > tempqstimesparent[currentPosition])
                return Double.NEGATIVE_INFINITY;
            if (tminparentnew/tnewparentbottom < temptiptimesparent[i]/tempqstimesparent[currentPosition]) {
                tminparentnew = temptiptimesparent[i];
                tnewparentbottom = tempqstimesparent[currentPosition];
            }
            currentPosition -= temptiptimescountparent[i];
        }
        // check that the scaling was ok
        if (tmaxparent/tnewparenttop < tminparentnew/tnewparentbottom){
            throw new IllegalStateException("problem in hereeeeee: HaplotypeSwap: The scaled values are not calculated properly?");
        }
        // Reject invalid haplotype scalings:
        if (parentscalefactor < 1.0){
            if (tnewQSstartparent < parentHaplo.getHeight()){
                throw new IllegalStateException("problem in hereeeeee you did not really scale the QS start apparently");
            }
            else if (tempqstimesparent[tempqstimesparent.length-1] < parentHaplo.getHeight())
                return Double.NEGATIVE_INFINITY;
        }

// HR AND CORRECTIONS TO THE TREE

        // HR contribution from CURRENT haplo
        // Incorporate the probability of scaling all the attachment times
        // assign contribution of each scaled attachment time
        logHastingsRatio += (tempqstimes.length - 3) * Math.log(scalefactor);

        // assign contribution to the Hastings ratio for having different possible scales for toldtop
        logHastingsRatio += Math.log(tmax/toldtop - tminold/toldbottom);
        logHastingsRatio -= Math.log(tmaxparent/tnewtop - tminold/tnewbottom);

        // HR contribution from PARENT haplo
        // Incorporate the probability of scaling all the attachment times
        // assign contribution of each scaled attachment time
        logHastingsRatio += (tempqstimesparent.length - 3) * Math.log(parentscalefactor);

        // assign contribution to the Hastings ratio for having different possible scales for toldparenttop
        logHastingsRatio += Math.log(tmaxparent/toldparenttop - tminparentold/toldparentbottom);
        logHastingsRatio -= Math.log(tmax/tnewparenttop - tminparentnew/tnewparentbottom);

        // rewrite the attachment times array
        node.setAttachmentTimesList(tempqstimes);

         // rewrite the attachment times array
        parentHaplo.setAttachmentTimesList(tempqstimesparent);

        // find out what nodes have been passed by changing the haplo height
        QuasiSpeciesNode nodeBelowHaploMoved = (QuasiSpeciesNode) findNodeBelowAfterRepositioningHaploStart(oldNodeBelowCurrentHaplo,tnewQSstart);

        // find the new node above which the haplo arises
        QuasiSpeciesNode nodeBelowHaploParent = (QuasiSpeciesNode) findNodeBelowAfterRepositioningHaploStart(parentHaplo,tnewQSstartparent);


        // calculate up to where the haplotypes from the passed nodes can max attach to
        // also get the node up to which they can go to
//        ArrayList upToWhere;
//
//        if (nodeBelowHaploMoved.getContinuingHaploName() != -1)
//            upToWhere = getMaxPossibleHaploAttachTimeForQSStart(nodeBelowHaploMoved, nodeBelowHaploMoved.getContinuingHaploName(), 0);
//        else
//        // here the haplo is a dummy thing... we just need to find out up to where the QS start from the passed haplotypes can move
//            upToWhere = getMaxPossibleHaploAttachTimeForQSStart(nodeBelowHaploMoved, haplo, 0);

        oldNodeBelowCurrentHaplo.setHaploAboveName(-1);
        nodeBelowHaploMoved.setHaploAboveName(haplo);

        if (oldNodeBelowParentHaplo.getHaploAboveName() == parentHaplo.getNr())
            oldNodeBelowParentHaplo.setHaploAboveName(-1);
        nodeBelowHaploParent.setHaploAboveName(parentHaplo.getNr());

        // recalculate parentHaplo array, re-assign continuingHaploName
//        if ((double) upToWhere.get(1)==origin.getValue()){
//            recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode)qsTree.getRoot());
//        }
//        else if (((QuasiSpeciesNode) upToWhere.get(0)).getHaploAboveName() != -1)
//            recalculateParentHaploAndCorrectContinuingHaploName(parentHaploArray[(int) upToWhere.get(2)], (QuasiSpeciesNode) upToWhere.get(0));
//        else {
//            QuasiSpeciesNode nodetemp = findNodeBelowThisHaplo((QuasiSpeciesNode) upToWhere.get(0),(int) upToWhere.get(2));
//            recalculateParentHaploAndCorrectContinuingHaploName((int) upToWhere.get(2), nodetemp);
//        }
        recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode)qsTree.getRoot());

// NODE ITSELF
// FOR ALL CASES: reposition the HAPLOTYPE START ITSELF

        // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
        qsTree.countAndSetPossibleStartBranches();

        // Ensure BEAST knows to recalculate affected likelihood:
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        parentHaplo.makeDirty(QuasiSpeciesTree.IS_FILTHY);

    // RETURN log(HASTINGS RATIO)
    return logHastingsRatio; // proper hastings ratio!!!
    }

}
