package quasispeciestree.operators;

import beast.core.Description;
import beast.core.Input;
import quasispeciestree.tree.QuasiSpeciesNode;
import beast.util.Randomizer;
import java.util.ArrayList;


/**
 *  @author Veronika Boskova created on 06/08/2015 finished on 09/06/2016
 */
@Description("Chooses a haplotype at random and moves"
        + "its first attachment time uniformly in interval "
        + "restricted by the grand-parent haplotype (start)"
        +" and tip time of the moved haplotype,"
        +" and scales all the attachment times accordingly."
        +" If there is another haplotype in the way, this will"
        +" be scaled down below the MRCA of the two haplotypes.")
//  NOTE: we are only shifting the passed parent haplotype's first attachment time, if the first attachment
//        time of the moved haplotype is above their common ancestor. This is to make sure all is reversible.
//  QS start time swap is coded in the operator QuasiSpeciesHaplotypeStartSwap
public class QuasiSpeciesHaplotypeSwap extends QuasiSpeciesTreeOperator{

//***    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
//***            "Scaling is restricted to the range [1/scaleFactor, scaleFactor]");

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
        int haplo = Randomizer.nextInt(qsTree.getLeafNodeCount());

        QuasiSpeciesNode node = (QuasiSpeciesNode) qsTree.getNode(haplo);

        // if we did not assign the node at this stage, throw exception
        if (node == null)
            throw new IllegalStateException("Event selection loop fell through!");

        // get the parentHaploArray and the attachment times array to be changed
        int[] parentHaploArray = qsTree.getParentHaplo();
        Double[] tempqstimes=qsTree.getAttachmentTimesList(node).clone();

        // get a node above which the current haplotype arises
        QuasiSpeciesNode oldNodeBelowHaploMoved = findNodeBelowThisHaplo(node,haplo);
        QuasiSpeciesNode oldNodeBelowHaploParent=null;
        QuasiSpeciesNode oldNodeBelowHaploGrandParent=null;

        // reposition the event (i.e. haplotype's first attachment time) uniformly at random between tmin and tmax
        // what is the parent haplo start time?
        ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTimeForQSStart(oldNodeBelowHaploMoved, haplo, 0);
//        int haplotypesParentHaplo = (int) haploStartMaxNewArray.get(2);
        int haplotypesParentHaplo = parentHaploArray[haplo];
        // what is the grandparent haplo start time?
        ArrayList haploStartMaxNewArrayParent = new ArrayList(3);
        int haplotypesGrandParentHaplo = -1;

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
        if (tempqstimes.length>1)
            told=tempqstimes[1];
        double tnew=0;
        // get origin height
        double originheight = origin.getValue();

// NO PARENT
        if (haplotypesParentHaplo == -1){

            // find out tmax
            tmax = originheight;

            // reposition attachment times: attach ((time - tmin) * (tnew/told)) + tmin
            if (tempqstimes.length > 1){
                // scale all the other positions in the array but the 0 position (haplo start time)
                scalefactor = u*(tmin/told) + (1.0-u)*(tmax/told);
//***            scalefactor = u*scaleFactorInput.get()+(1.0-u)/scaleFactorInput.get();
                for (int i=1; i<tempqstimes.length; i++) {
                    tempqstimes[i] = tempqstimes[i] * scalefactor;
                }
                // get new time to attach of first attachment time
                tnew = tempqstimes[1];

//***                if (tnew>tmax)
//***                    return Double.NEGATIVE_INFINITY;

                // set the haplotype's starting time to the new time
                double x = Randomizer.nextDouble();
                tnewQSstart = x*tnew + (1.0-x)*tmax;
                tempqstimes[0] = tnewQSstart;

                // Reject invalid haplotype scalings:
                if (scalefactor<1.0 && tempqstimes.length>1){
                    if (tnewQSstart<node.getHeight()){
                        System.out.println("problem in hereeeeee you did not really scale the QS start apparently");
                        System.exit(0);
                    }
                    else if (tempqstimes[tempqstimes.length-1]<node.getHeight())
                        return Double.NEGATIVE_INFINITY;
                }

                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                logHastingsRatio -= Math.log(tmax - told);
//                logHastingsRatio += Math.log(tmax - tempqstimes[1]);
                // assign contribution to the Hastings ratio for having different possible scales for told
                logHastingsRatio += Math.log(tmax/told - tmin/told);
                logHastingsRatio -= Math.log(tmax/tnew - tmin/tnew);
                // Incorporate the probability of scaling all the attachment times
                // assign contribution of each scaled attachment time
                logHastingsRatio += (tempqstimes.length-3) * Math.log(scalefactor);
                }
                else {
                    // set the haplotype's starting time
                    tnew = u*tmin + (1-u)*tmax;
                    tnewQSstart = tnew;
                    tempqstimes[0] = tnewQSstart;
                    // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
                    //              logHastingsRatio -= Math.log(tmax - tmin);
                    //              logHastingsRatio += Math.log(tmax - tmin);
                    // since they are identical, no need to change Hastings ratio
                }
            // rewrite the attachment times array
            qsTree.setAttachmentTimesList(node, tempqstimes);

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

            // only shuffle haplotypes start times, if the current haplo moved such that there is a new node below it
            if (oldNodeBelowHaploMoved.getNr()!=nodeBelowHaploMoved.getNr()){
//                // reposition QS start points of the QS on the way and add contribution to the HR for QS start moves
//                // calculate up to where the haplotypes from the passed nodes can max attach to
//                // also get the node up to which they can go to
//                ArrayList upToWhere;
//
//                if (oldNodeBelowHaploMoved.getHeight() > nodeBelowHaploMoved.getHeight())
//                    upToWhere = getMaxPossibleHaploAttachTimeForQSStart(oldNodeBelowHaploMoved, haplo, 0);
//                else if (nodeBelowHaploMoved.getContinuingHaploName() != -1)
//                    upToWhere = getMaxPossibleHaploAttachTimeForQSStart(nodeBelowHaploMoved, nodeBelowHaploMoved.getContinuingHaploName(), 0);
//                else
//                    // here the haplo is a dummy thing... we just need to find out up to where the QS start from the passed haplotypes can move
//                    upToWhere = getMaxPossibleHaploAttachTimeForQSStart(nodeBelowHaploMoved, haplo, 0);
//
//                logHastingsRatio += moveQSInWayAfterMovingHaplotype((double)upToWhere.get(1),tmin,tnewQSstart,toldQSstart,nodesPassed,parentHaploArray);
                // at the same time change internal node haploName - both remove and add if necessary
                // change the internal nodes' haploName, where necessary
//                if (oldNodeBelowHaploMoved.getHaploAboveName()==node.getNr())
                    oldNodeBelowHaploMoved.setHaploAboveName(-1);
                nodeBelowHaploMoved.setHaploAboveName(haplo);

//                recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode)qsTree.getRoot());
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

//                // recalculate parentHaplo array, re-assign continuingHaploName
//                if ((double) upToWhere.get(1)==origin.getValue()){
//                    recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode)qsTree.getRoot());
//                }
//                else if (((QuasiSpeciesNode) upToWhere.get(0)).getHaploAboveName() != -1)
//                    recalculateParentHaploAndCorrectContinuingHaploName(parentHaploArray[(int) upToWhere.get(2)], (QuasiSpeciesNode) upToWhere.get(0));
//                else {
//                    QuasiSpeciesNode nodetemp = findNodeBelowThisHaplo((QuasiSpeciesNode) upToWhere.get(0),(int) upToWhere.get(2));
//                    recalculateParentHaploAndCorrectContinuingHaploName((int) upToWhere.get(2), nodetemp);
//                }
//
//                // mark passed nodes as dirty
//                if (tnewQSstart>toldQSstart){
//                    for (Node nodePassed : nodesPassed)
//                        nodePassed.makeDirty(QuasiSpeciesTree.IS_FILTHY);
//                }
            }
        }
// PARENT
        // if there is a parent haplotype, we need to change more stuff (also for the parent haplotype)
        else{
            // check what is the grand parent haplotype
            haploStartMaxNewArrayParent = getMaxPossibleHaploAttachTimeForQSStart(
                    (QuasiSpeciesNode)haploStartMaxNewArray.get(0),haplotypesParentHaplo ,0);
//            haplotypesGrandParentHaplo = (int) haploStartMaxNewArrayParent.get(2);
            haplotypesGrandParentHaplo = parentHaploArray[haplotypesParentHaplo];
// NO GRANDPARENT
            // if none, change only the arrays for the parent
            if (haplotypesGrandParentHaplo == -1){

                // get a node above which the parent haplotype arises
                oldNodeBelowHaploParent = findNodeBelowThisHaplo(node,haplotypesParentHaplo);

                // find out tmax: grandparent haplotype is null, then tmax=origin
                tmax = originheight;

            }
// GRANDPARENT
            // otherwise change stuff also for the grandparent
            else{
                oldNodeBelowHaploParent = findNodeBelowThisHaplo(node,haplotypesParentHaplo);
                // get a node above which the parent haplotype arises
//                oldNodeBelowHaploGrandParent = findNodeBelowThisHaplo(node,(int) haploStartMaxNewArrayParent.get(2));
                oldNodeBelowHaploGrandParent = findNodeBelowThisHaplo(node,haplotypesGrandParentHaplo);

                // find out tmax: restricted by the last common ancestor node of grandparent haplotype and chosen haplotype
                // the grandparent haplotype can be from the same subtree or from the other subtree
                // in either case (same or other subtree) the last common ancestor node is where the current haplo has tmax
                tmax = (double) haploStartMaxNewArrayParent.get(1);
            }
// PARENT AND GRANDPARENT
            //choose new time to attach
//            tnew = u*tmin + (1-u)*tmax;

            // propose a new start time for the parent node uniformly at random between tminparent and tmaxparent
            double tminparent = qsTree.getNode(haplotypesParentHaplo).getHeight();
            double tmaxparent = (double) haploStartMaxNewArray.get(1);
            double tnewparent = 0;
            double toldparent=0;

            // reposition attachment times: attach ((time - tmin) * (tnew/told)) + tmin
            if (tempqstimes.length > 1){
                // scale all the other positions in the array but the 0 position (haplo start time)
//                scalefactor = (tnew - tmin)/(told - tmin);
                scalefactor = u*(tmin/told) + (1.0-u)*(tmax/told);
//                tnew = tempqstimes[1] * scalefactor;
                for (int i=1; i<tempqstimes.length; i++) {
//                    tempqstimes[i] = ((tempqstimes[i] - tmin) * scalefactor) + tmin;
                    tempqstimes[i] = tempqstimes[i] * scalefactor;
                }
                // get new time to attach of first attachment time
                tnew = tempqstimes[1];

                // set the haplotype's starting time to the new time
                double x = Randomizer.nextDouble();
//                // if the first attachment time did no go above tmaxparent, then QS start time also cannot go beyond it
//                if (tempqstimes[1]<tmaxparent){
//                    tnewQSstart = x*tempqstimes[1] + (1.0-x)*tmaxparent;
//                    // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                    logHastingsRatio -= Math.log(tmaxparent - told);
//                    logHastingsRatio += Math.log(tmaxparent - tempqstimes[1]);
//                    // assign contribution to the Hastings ratio for having different possible scales for told
//                    logHastingsRatio += Math.log(tmax/told - tmin/told);
//                    logHastingsRatio -= Math.log(tmax/tnew - tmin/tnew);
//                }
//                else{
                    tnewQSstart = x*tnew + (1.0-x)*tmax;
                    // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                    logHastingsRatio -= Math.log(tmaxparent - told);
//                    logHastingsRatio += Math.log(tmax - tempqstimes[1]);
//                    //logHastingsRatio += Math.log(tmax-tmin);
//                    //logHastingsRatio -= Math.log(tmaxparent-tmin);

                // assign contribution to the Hastings ratio for having different possible scales for told
                logHastingsRatio += Math.log(tmax/told - tmin/told);
                logHastingsRatio -= Math.log(tmaxparent/tnew - tmin/tnew);
//                }
                // Incorporate the probability of scaling all the attachment times
                // assign contribution of each scaled attachment time
//                if (tempqstimes.length>1){
//                logHastingsRatio -= 2 * (Math.log (scalefactor * (told - tmin) + tmin) - Math.log(told));
//                logHastingsRatio -= 2 * Math.log(scalefactor);
//                logHastingsRatio += (tempqstimes.length-1) * Math.log(scalefactor);
                logHastingsRatio += (tempqstimes.length-3) * Math.log(scalefactor);
//                }

                tempqstimes[0] = tnewQSstart;

                //Reject invalid haplotype scalings:
                if (scalefactor<1.0 && tempqstimes.length>1){
                    if (tnewQSstart<node.getHeight()){
                        System.out.println("problem in hereeeeee you did not really scale the QS start apparently");
                        System.exit(0);
                    }
                    else if (tempqstimes[tempqstimes.length-1]<node.getHeight())
                        return Double.NEGATIVE_INFINITY;
                }
            }
            else {
                // set the haplotype's starting time to the new time
//                tnewQSstart = tnew;
                tnew = u*tmin + (1-u)*tmax;
                tnewQSstart = tnew;
                tempqstimes[0] = tnewQSstart;
                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                logHastingsRatio -= Math.log(tmaxparent - tmin);
//                logHastingsRatio += Math.log(tmax - tmin);
            }
            // rewrite the attachment times array
            qsTree.setAttachmentTimesList(node, tempqstimes);

            // move the parent haplotype at random below the common ancestor node
            // ! only if chosen haplo moves above last common ancestor of parent and current haplo

            if (tnewQSstart > tmaxparent) {
                Double[] tempqstimesparent=qsTree.getAttachmentTimesList(haplotypesParentHaplo).clone();


                // choose new time to attach for the moved parent node and scale all attachment times of this haplo
                double v = Randomizer.nextDouble();
                // choose new time to attach
//                double tnewparent = v*tminparent + (1-v)*tmaxparent;

                // get the haplotype's starting time
                double toldQSstartparent = tempqstimesparent[0];
                double parentscalefactor=0;
                double tnewQSstartparent;

                // reposition parent haplo attachment times: attach ((time - tminparent) * (tnew/told)) + tminparent
                if (tempqstimesparent.length > 1){
                    toldparent = tempqstimesparent[1];
                    // scale all the other positions in the array but the 0 position (haplo start time)
//                    parentscalefactor = (tnewparent - tminparent)/(toldparent - tminparent);
                    parentscalefactor = v*(tminparent/toldparent) + (1.0-v)*(tmaxparent/toldparent);

                    for (int i=1; i<tempqstimesparent.length; i++) {
//                        tempqstimesparent[i] = ((tempqstimesparent[i] - tminparent) * parentscalefactor) + tminparent;
                        tempqstimesparent[i] = tempqstimesparent[i] * parentscalefactor;
                    }
                    // get new time to attach of first attachment time
                    tnewparent = tempqstimesparent[1];

                    // set the haplotype's starting time to the new time
                    double y = Randomizer.nextDouble();
                    tnewQSstartparent = y*tnewparent + (1-y)*tmaxparent;
                    tempqstimesparent[0] = tnewQSstartparent;

                    // Reject invalid haplotype scalings:
                    if (parentscalefactor<1.0 && tempqstimesparent.length>1){
                        if (tnewQSstartparent<qsTree.getNode(haplotypesParentHaplo).getHeight()){
                            System.out.println("problem in hereeeeee you did not really scale the QS start apparently");
                            System.exit(0);
                        }
                        else if (tempqstimesparent[tempqstimesparent.length-1]<qsTree.getNode(haplotypesParentHaplo).getHeight())
                            return Double.NEGATIVE_INFINITY;
                    }
                    // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                    logHastingsRatio -= Math.log(tmax - toldparent);
//                    logHastingsRatio += Math.log(tmaxparent - tempqstimesparent[1]);
                    // assign contribution to the Hastings ratio for having different possible scales for told
//                    logHastingsRatio += Math.log((tmaxparent-tminparent)/(toldparent-tminparent) - (tminparent-tminparent)/(toldparent-tminparent));
                    logHastingsRatio += Math.log(tmaxparent/toldparent - tminparent/toldparent);
//                    logHastingsRatio -= Math.log((tmax-tminparent)/(tnewparent-tminparent) - (tminparent-tminparent)/(tnewparent-tminparent));
                    logHastingsRatio -= Math.log(tmax/tnewparent - tminparent/tnewparent);
                    // Incorporate the probability of scaling all the attachment times
//                    logHastingsRatio += Math.log(tmaxparent-tminparent);
//                    logHastingsRatio -= Math.log(tmax-tminparent);
                    // assign contribution of each scaled attachment time
//                    if (tempqstimesparent.length>1){
//                    logHastingsRatio -= 2 * (Math.log(parentscalefactor * (toldparent - tminparent) + tminparent) - Math.log(toldparent));
//                    logHastingsRatio -= 2 * Math.log(parentscalefactor);
//                    logHastingsRatio += (tempqstimesparent.length-1) * Math.log(parentscalefactor);
                    logHastingsRatio += (tempqstimesparent.length-3) * Math.log(parentscalefactor);
//                    }
                }
                else {
                    // set the haplotype's starting time to the new time
//                    tnewQSstartparent = tnewparent;
                    tnewparent = u*tminparent + (1-u)*tmaxparent;
                    tnewQSstartparent = tnewparent;
                    tempqstimesparent[0] = tnewQSstartparent;
                    // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                    logHastingsRatio -= Math.log(tmax - tminparent);
//                    logHastingsRatio += Math.log(tmaxparent - tminparent);
                }

                // rewrite the attachment times array
                qsTree.setAttachmentTimesList(haplotypesParentHaplo, tempqstimesparent);

                // add probability that tempqstimes[1] (and then also tnewQSstart) > tmaxparent
                logHastingsRatio += Math.log(tmax-tmaxparent) - Math.log(tmax-tminparent);
//                logHastingsRatio += Math.log(tmax/tnewparent-tmaxparent/tnewparent) - Math.log(tmax/tnewparent-tminparent/tnewparent);
                logHastingsRatio -= (Math.log(tmax-tmaxparent) - Math.log(tmax-tmin));
//                logHastingsRatio -= (Math.log(tmax/told-tmaxparent/told) - Math.log(tmax/told-tmin/told));

                // scale all QS start points for all haplotypes on a way from toldQSstart to tnewQSstart
                // find out what nodes have been passed by changing the haplo height
//                ArrayList<QuasiSpeciesNode> nodesPassed = new ArrayList<>(qsTree.getInternalNodeCount());
                // find the new node above which the haplo arises
//                QuasiSpeciesNode nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(
//                                                        oldNodeBelowHaploMoved,tnewQSstart,nodesPassed);
                QuasiSpeciesNode nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(
                                                        oldNodeBelowHaploMoved,tnewQSstart);
//                nodesPassed.trimToSize();

                // find out what nodes have been passed by changing the parenthaplo height
                //ArrayList<QuasiSpeciesNode> nodesPassedParent = new ArrayList<>(qsTree.getInternalNodeCount());
                // find the new node above which the haplo arises
                QuasiSpeciesNode nodeBelowHaploParent=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(
                                                        (QuasiSpeciesNode)qsTree.getNode(haplotypesParentHaplo),tnewQSstartparent);
                //NOT NECESSARY
                //// then check out which nodes were passed
                //QuasiSpeciesNode intnodehelper=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(
                //                                    nodeBelowHaploParent,toldQSstartparent,nodesPassedParent);
                //nodesPassedParent.trimToSize();

//                // only shuffle haplotypes start times, if the current haplo moved such that there is a new node below it
//                // IN THIS CASE WE KNOW THERE WERE SOME NODES PASSED
//                //if (oldNodeBelowHaploMoved.getNr()!=nodeBelowHaploMoved.getNr()){
//                    // reposition QS start points of the QS on the way and add contribution to the HR for QS start moves
//                    // FIRST contribution comes from nodes passed from oldNodeBelowHaploMoved to haploStartMaxNewArray.get(0)
//                    // find out what nodes have been passed by changing the haplo height
//                    ArrayList<QuasiSpeciesNode> nodesPassed1 = new ArrayList<>(qsTree.getInternalNodeCount());
//                    findNodesPassed(oldNodeBelowHaploMoved,(Node)haploStartMaxNewArray.get(0),nodesPassed1,false);
//                    nodesPassed1.trimToSize();
//                    // correct QS start times at these nodes
//                    logHastingsRatio += moveQSInWayAfterMovingHaplotype((double)haploStartMaxNewArray.get(1),tmin,tnewQSstart,toldQSstart,nodesPassed1,parentHaploArray);
//
//                    // SECOND contribution comes from nodes passed from oldNodeBelowHaploParent to haploStartMaxNewArrayParent.get(0)
//                    // calculate up to where the haplotypes from the passed nodes can max attach to
//                    // also get the node up to which they can go to
//                    ArrayList upToWhere;
//
//                    if (nodeBelowHaploMoved.getContinuingHaploName() != -1)
//                        upToWhere = getMaxPossibleHaploAttachTimeForQSStart(nodeBelowHaploMoved, nodeBelowHaploMoved.getContinuingHaploName(), 0);
//                    else
//                        // here the haplo is a dummy thing... we just need to find out up to where the QS start from the passed haplotypes can move
//                        upToWhere = getMaxPossibleHaploAttachTimeForQSStart(nodeBelowHaploMoved, haplo, 0);
//
//                    // find out what nodes have been passed by changing the haplo height
//                    ArrayList<QuasiSpeciesNode> nodesPassed2 = new ArrayList<>(qsTree.getInternalNodeCount());
//                    if ((double)upToWhere.get(1)==origin.getValue())
//                        findNodesPassed(oldNodeBelowHaploParent,qsTree.getRoot(),nodesPassed2,true);
//                    else
//                        findNodesPassed(oldNodeBelowHaploParent,(Node)upToWhere.get(0),nodesPassed2,true);
//                    nodesPassed2.trimToSize();
//                    // correct QS start times at these nodes
//                    logHastingsRatio += moveQSInWayAfterMovingHaplotype((double)upToWhere.get(1),tmin,tnewQSstart,toldQSstart,nodesPassed2,parentHaploArray);
                    // at the same time change internal node haploName - both remove and add if necessary
                    // change the internal nodes' haploName, where necessary
//                    if (oldNodeBelowHaploMoved.getHaploAboveName()==node.getNr())
                        oldNodeBelowHaploMoved.setHaploAboveName(-1);
                    nodeBelowHaploMoved.setHaploAboveName(haplo);

                    if (oldNodeBelowHaploParent.getHaploAboveName()==haplotypesParentHaplo)
                        oldNodeBelowHaploParent.setHaploAboveName(-1);
                    nodeBelowHaploParent.setHaploAboveName(haplotypesParentHaplo);

                    recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode)qsTree.getRoot());

//                    // recalculate parentHaplo array, re-assign continuingHaploName
//                    if ((double) upToWhere.get(1)==origin.getValue()){
//                        recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode)qsTree.getRoot());
//                    }
//                    else if (((QuasiSpeciesNode) upToWhere.get(0)).getHaploAboveName() != -1)
//                        recalculateParentHaploAndCorrectContinuingHaploName(parentHaploArray[(int) upToWhere.get(2)], (QuasiSpeciesNode) upToWhere.get(0));
//                    else {
//                        QuasiSpeciesNode nodetemp = findNodeBelowThisHaplo((QuasiSpeciesNode) upToWhere.get(0),(int) upToWhere.get(2));
//                        recalculateParentHaploAndCorrectContinuingHaploName((int) upToWhere.get(2), nodetemp);
//                    }
//
//                    // mark passed nodes as dirty
//                    if (tnewQSstart>toldQSstart){
//                        for (Node nodePassed : nodesPassed)
//                            nodePassed.makeDirty(QuasiSpeciesTree.IS_FILTHY);
//                    }
                //}

            }
            // if the tnew is not above the lastParentCommonAncestorNode, then just reposition the current haplo start
            else {

                // add probability that tempqstimes[1] (and then also tnewQSstart) < tmaxparent
//                logHastingsRatio += Math.log(tmaxparent-tmin) - Math.log(tmax-tmin);
//                logHastingsRatio += Math.log(tmaxparent/tnew-tmin/tnew) - Math.log(tmax/tnew-tmin/tnew);
//                logHastingsRatio -= (Math.log(tmaxparent-tmin) - Math.log(tmax-tmin));
//                logHastingsRatio -= (Math.log(tmaxparent/told-tmin/told) - Math.log(tmax/told-tmin/told));

                // scale all QS start points for all haplotypes on a way from toldQSstart to tnewQSstart
                // find out what nodes have been passed by changing the haplo height
//                ArrayList<QuasiSpeciesNode> nodesPassed = new ArrayList<>(qsTree.getInternalNodeCount());
                // find the new node above which the haplo arises
                QuasiSpeciesNode nodeBelowHaploMoved;
                if (tnewQSstart > toldQSstart) {
//                    nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(oldNodeBelowHaploMoved,tnewQSstart,nodesPassed);
                    nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(oldNodeBelowHaploMoved,tnewQSstart);
                }
                else {
                    // first find the new node below
                    nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(node,tnewQSstart);
                    // then check out which nodes were passed
//                    QuasiSpeciesNode intnodehelper=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(nodeBelowHaploMoved,toldQSstart,nodesPassed);
                }
//                nodesPassed.trimToSize();

                // only shuffle haplotypes start times, if the current haplo moved such that there is a new node below it
                if (oldNodeBelowHaploMoved.getNr()!=nodeBelowHaploMoved.getNr()){
                    // reposition QS start points of the QS on the way and add contribution to the HR for QS start moves
                    // calculate up to where the haplotypes from the passed nodes can max attach to
                    // also get the node up to which they can go to
//                    ArrayList upToWhere;
//
//                    if (oldNodeBelowHaploMoved.getHeight() > nodeBelowHaploMoved.getHeight())
//                        upToWhere = getMaxPossibleHaploAttachTimeForQSStart(oldNodeBelowHaploMoved, haplo, 0);
//                    else if (nodeBelowHaploMoved.getContinuingHaploName() != -1)
//                        upToWhere = getMaxPossibleHaploAttachTimeForQSStart(nodeBelowHaploMoved, nodeBelowHaploMoved.getContinuingHaploName(), 0);
//                    else
//                        // here the haplo is a dummy thing... we just need to find out up to where the QS start from the passed haplotypes can move
//                        upToWhere = getMaxPossibleHaploAttachTimeForQSStart(nodeBelowHaploMoved, haplo, 0);
//
//                    logHastingsRatio += moveQSInWayAfterMovingHaplotype((double)upToWhere.get(1),tmin,tnewQSstart,toldQSstart,nodesPassed,parentHaploArray);
                    // at the same time change internal node haploName - both remove and add if necessary
                    // change the internal nodes' haploName, where necessary
//                    if (oldNodeBelowHaploMoved.getHaploAboveName()==node.getNr())
                        oldNodeBelowHaploMoved.setHaploAboveName(-1);
                    nodeBelowHaploMoved.setHaploAboveName(haplo);

                    recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode)qsTree.getRoot());

//                    // recalculate parentHaplo array, re-assign continuingHaploName
//                    if ((double) upToWhere.get(1)==origin.getValue()){
//                        recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode)qsTree.getRoot());
//                    }
//                    else if (((QuasiSpeciesNode) upToWhere.get(0)).getHaploAboveName() != -1)
//                        recalculateParentHaploAndCorrectContinuingHaploName(parentHaploArray[(int) upToWhere.get(2)], (QuasiSpeciesNode) upToWhere.get(0));
//                    else {
//                        QuasiSpeciesNode nodetemp = findNodeBelowThisHaplo((QuasiSpeciesNode) upToWhere.get(0),(int) upToWhere.get(2));
//                        recalculateParentHaploAndCorrectContinuingHaploName((int) upToWhere.get(2), nodetemp);
//                    }
//
//                    // mark passed nodes as dirty
//                    if (tnewQSstart>toldQSstart){
//                        for (Node nodePassed : nodesPassed)
//                            nodePassed.makeDirty(QuasiSpeciesTree.IS_FILTHY);
//                    }
                }
            }
        }
// NODE ITSELF
// FOR ALL CASES: reposition the HAPLOTYPE START ITSELF

        // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
        int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
        qsTree.setStartBranchCounts(startBranchCountsArray);

    // RETURN log(HASTINGS RATIO)
    return logHastingsRatio; // proper hastings ratio!!!
    }

}
