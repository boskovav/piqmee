package quasispeciestree.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesNode;
import beast.util.Randomizer;
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
            if (qsTree.getHaplotypeCounts((QuasiSpeciesNode) node)>0)
                canbeused++;
        }
        if (canbeused<1){
            System.out.println("QuasiSpeciesHaplotypeSwap operator cannot be used since there are no or just one "
                    + "haplotype with duplicates! Remove this operator from your XML file.");
            System.exit(0);
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
        for (int i = 0; i < qsTree.getParentHaplo().length; i++){
            if (qsTree.getParentHaplo(i)!=-1)
                haplowithparents++;
        }
        if (haplowithparents==0)
            return Double.NEGATIVE_INFINITY;

        // keep track of the Hastings ratio
        double logHastingsRatio = 0.0;
        // Randomly select event on tree:
        // for now choose uniformly at random from just the unique haplotype count (disregarding the counts)
        int haplo = -1;
        do {
            haplo = Randomizer.nextInt(qsTree.getLeafNodeCount());
        } while (qsTree.getAttachmentTimesList(haplo).length<2 || qsTree.getParentHaplo(haplo)==-1);

        QuasiSpeciesNode node = (QuasiSpeciesNode) qsTree.getNode(haplo);

        // if we did not assign the node at this stage, throw exception
        if (node == null)
            throw new IllegalStateException("Event selection loop fell through!");

        // get the parentHaploArray and the attachment times array to be changed
        int[] parentHaploArray = qsTree.getParentHaplo();
        int parentHaplo = parentHaploArray[haplo];
//        if (parentHaplo != qsTree.getParentHaplo(node.getNr())){
//            System.out.println("problem in hereeeeee: HaplotypeSwap: Why the found parent haplo is not the same as stored one?");
//            System.exit(0);
//        }
        int grandParentHaplo = parentHaploArray[parentHaplo];
//        if (grandParentHaplo != qsTree.getParentHaplo(parentHaplo)) {
//            System.out.println("problem in hereeeeee: HaplotypeSwap: Why the found grand parent haplo is not the same as stored one?");
//            System.exit(0);
//        }

        // get the attachment times
        double[] tempqstimes=qsTree.getAttachmentTimesList(haplo).clone();
        // get also tip times to help define max/min scalings
        double[] temptiptimes=qsTree.getTipTimesList(haplo);
        int[] temptiptimescount=qsTree.getTipTimesCountList(haplo);
        double[] tempqstimesparent=qsTree.getAttachmentTimesList(parentHaplo).clone();
        // get also tip times to help define maxparent/minparent scalings
        double[] temptiptimesparent=qsTree.getTipTimesList(parentHaplo);
        int[] temptiptimescountparent=qsTree.getTipTimesCountList(parentHaplo);

        // get a node above which the current/parent/grand parent haplotype arises
        QuasiSpeciesNode oldNodeBelowCurrentHaplo = findNodeBelowThisHaplo(node,haplo);
        QuasiSpeciesNode oldNodeBelowParentHaplo = findNodeBelowThisHaplo(node,parentHaplo);

        // reposition the event (i.e. haplotype's first attachment time) uniformly at random between tmin and tmax
        // what is the parent haplo start time?
        ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTimeForQSStart(oldNodeBelowCurrentHaplo, haplo, 0);


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
        int currentPosition=tempqstimes.length-1-(temptiptimescount[0]-1);
        for (int i = 1; i < temptiptimes.length; i++){
            if (tminold/toldbottom<temptiptimes[i]/tempqstimes[currentPosition]) {
                tminold=temptiptimes[i];
                toldbottom=tempqstimes[currentPosition];
            }
            currentPosition=currentPosition-temptiptimescount[i];
        }
        // PARENT
        // get the haplotype's starting time
        double parentscalefactor = 0;
        double toldQSstartparent = tempqstimesparent[0];
        double tnewQSstartparent = 0;
        double tmaxparent = (double) haploStartMaxNewArray.get(1);
        double tminparentold = qsTree.getNode(parentHaplo).getHeight();
        double tminparentnew = qsTree.getNode(parentHaplo).getHeight();
        double toldparenttop = tempqstimesparent[1];
        double tnewparenttop = 0;
        double toldparentbottom = tempqstimesparent[tempqstimesparent.length-1];
        double tnewparentbottom = 0;
        // to get the tminparent check for each sampling time of the haplo that the last possible attach time is above the
        //  the corresponding sampling time
        currentPosition=tempqstimesparent.length-1-(temptiptimescountparent[0]-1);
        for (int i = 1; i < temptiptimesparent.length; i++){
            if (tminparentold/toldparentbottom<temptiptimesparent[i]/tempqstimesparent[currentPosition]) {
                tminparentold=temptiptimesparent[i];
                toldparentbottom=tempqstimesparent[currentPosition];
            }
            currentPosition=currentPosition-temptiptimescountparent[i];
        }
        // check that the scaling was ok
        if (tmaxparent/toldparenttop < tminparentold/toldparentbottom){
            // in this case, it is possible that there is no acceptable scaling for the parent haplo
            // but then this move cannot be performed at the moment
            return Double.NEGATIVE_INFINITY;
        }

// NO GRANDPARENT
        // find out tmax: grandparent haplotype is null, then tmax=origin
        if (grandParentHaplo == -1){
            tmax = origin.getValue();
        }
// GRANDPARENT
        // otherwise tmax is the common ancestor of grandParentHaplo and parentHaplo
        else{
            // get a node above which the grandparent haplotype arises
            tmax = (double) getMaxPossibleHaploAttachTimeForQSStart(
                    (QuasiSpeciesNode)haploStartMaxNewArray.get(0),parentHaplo,0).get(1);
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
//        double x = Randomizer.nextDouble();
//        tnewQSstart = x*tnewtop + (1.0-x)*tmax;
        tempqstimes[0] = tnewQSstart;

        // check that the scaling was ok
        if (tmax/toldtop < tminold/toldbottom){
            System.out.println("problem in hereeeeee: HaplotypeSwap: The scalings are not calculated properly?");
            System.exit(0);
        }
        // check that the newly scaled attachment times do not go over the boundaries
        if (tnewQSstart > tmax || tnewQSstart < tmaxparent)
            return Double.NEGATIVE_INFINITY;
        tnewbottom=tempqstimes[tempqstimes.length-1];
        currentPosition=tempqstimes.length-1-(temptiptimescount[0]-1);
        for (int i = 1; i < temptiptimes.length; i++){
            if(temptiptimes[i]>tempqstimes[currentPosition])
                return Double.NEGATIVE_INFINITY;
            if (tminnew/tnewbottom<temptiptimes[i]/tempqstimes[currentPosition]) {
                tminnew=temptiptimes[i];
                tnewbottom=tempqstimes[currentPosition];
            }
            currentPosition=currentPosition-temptiptimescount[i];
        }
        // check that the scaling was ok
        if (tmax/tnewtop < tminnew/tnewbottom){
            System.out.println("problem in hereeeeee: HaplotypeSwap: The scalings are not calculated properly?");
            System.exit(0);
        }
        // Reject invalid haplotype scalings:
        if (scalefactor<1.0){
            if (tnewQSstart<node.getHeight()){
                System.out.println("problem in hereeeeee you did not really scale the QS start apparently");
                System.exit(0);
            }
            else if (tempqstimes[tempqstimes.length-1]<node.getHeight())
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
//        double y = Randomizer.nextDouble();
//        tnewQSstartparent = y*tnewparenttop + (1-y)*tmaxparent;
        tempqstimesparent[0] = tnewQSstartparent;

        // check that the scaling was ok
        if (tmax/toldparenttop < tminparentold/toldparentbottom){
            System.out.println("problem in hereeeeee: HaplotypeSwap: The scalings are not calculated properly?");
            System.exit(0);
        }
        // check that the newly scaled attachment times do not go over the boundaries
        if (tnewQSstartparent > tmaxparent || tnewQSstartparent > tmax)
            return Double.NEGATIVE_INFINITY;
        tnewparentbottom=tempqstimesparent[tempqstimesparent.length-1];
        currentPosition=tempqstimesparent.length-1-(temptiptimescountparent[0]-1);
        for (int i = 1; i < temptiptimesparent.length; i++){
            if(temptiptimesparent[i]>tempqstimesparent[currentPosition])
                return Double.NEGATIVE_INFINITY;
            if (tminparentnew/tnewparentbottom<temptiptimesparent[i]/tempqstimesparent[currentPosition]) {
                tminparentnew=temptiptimesparent[i];
                tnewparentbottom=tempqstimesparent[currentPosition];
            }
            currentPosition=currentPosition-temptiptimescountparent[i];
        }
        // check that the scaling was ok
        if (tmaxparent/tnewparenttop < tminparentnew/tnewparentbottom){
            System.out.println("problem in hereeeeee: HaplotypeScale: The scalings are not calculated properly?");
            System.exit(0);
        }
        // Reject invalid haplotype scalings:
        if (parentscalefactor<1.0){
            if (tnewQSstartparent<qsTree.getNode(parentHaplo).getHeight()){
                System.out.println("problem in hereeeeee you did not really scale the QS start apparently");
                System.exit(0);
            }
            else if (tempqstimesparent[tempqstimesparent.length-1]<qsTree.getNode(parentHaplo).getHeight())
                return Double.NEGATIVE_INFINITY;
        }

// HR AND CORRECTIONS TO THE TREE

        // HR contribution from CURRENT haplo
        // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//        logHastingsRatio -= Math.log(tmaxparent - toldtop);
//        logHastingsRatio += Math.log(tmax - tnewtop);

        // Incorporate the probability of scaling all the attachment times
        // assign contribution of each scaled attachment time
        logHastingsRatio += (tempqstimes.length-3) * Math.log(scalefactor);

        // assign contribution to the Hastings ratio for having different possible scales for toldtop
        logHastingsRatio += Math.log(tmax/toldtop - tminold/toldbottom);
        logHastingsRatio -= Math.log(tmaxparent/tnewtop - tminold/tnewbottom);

        // HR contribution from PARENT haplo
        // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//        logHastingsRatio -= Math.log(tmax - toldparenttop);
//        logHastingsRatio += Math.log(tmaxparent - tnewparenttop);
        // Incorporate the probability of scaling all the attachment times
        // assign contribution of each scaled attachment time
        logHastingsRatio += (tempqstimesparent.length-3) * Math.log(parentscalefactor);

        // assign contribution to the Hastings ratio for having different possible scales for toldparenttop
        logHastingsRatio += Math.log(tmaxparent/toldparenttop - tminparentold/toldparentbottom);
        logHastingsRatio -= Math.log(tmax/tnewparenttop - tminparentnew/tnewparentbottom);

        // rewrite the attachment times array
        qsTree.setAttachmentTimesList(node, tempqstimes);

         // rewrite the attachment times array
        qsTree.setAttachmentTimesList(parentHaplo, tempqstimesparent);

        // find out what nodes have been passed by changing the haplo height
        QuasiSpeciesNode nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(oldNodeBelowCurrentHaplo,tnewQSstart);

        // find the new node above which the haplo arises
        QuasiSpeciesNode nodeBelowHaploParent=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(qsTree.getNode(parentHaplo),tnewQSstartparent);


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

        if (oldNodeBelowParentHaplo.getHaploAboveName()==parentHaplo)
            oldNodeBelowParentHaplo.setHaploAboveName(-1);
        nodeBelowHaploParent.setHaploAboveName(parentHaplo);

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
        int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
        qsTree.setStartBranchCounts(startBranchCountsArray);

        // Ensure BEAST knows to recalculate affected likelihood:
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        qsTree.getNode(parentHaplo).makeDirty(QuasiSpeciesTree.IS_FILTHY);


    // RETURN log(HASTINGS RATIO)
    return logHastingsRatio; // proper hastings ratio!!!
    }

}
