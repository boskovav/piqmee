package quasispeciestree.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesNode;
import beast.util.Randomizer;


/**
 *  @author Veronika Boskova created on 06/08/2015 finished on 23/04/2016
 */
@Description("Chooses a haplotype at random and moves"
        + "its first attachment time uniformly in interval "
        + "restricted by the grand-parent haplotype"
        +" and tip time of the same haplotype,"
        +" and scales all the attachment times accordingly."
        +" If there is another haplotype in the way, this will"
        +" be scalde down below the MRCA of the two haplotypes.")
public class QuasiSpeciesHaplotypeStartSwap extends QuasiSpeciesTreeOperator{

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
        int event = Randomizer.nextInt(qsTree.getLeafNodeCount());

        QuasiSpeciesNode node = (QuasiSpeciesNode) qsTree.getNode(event);

        // if we did not assign the node at this stage, throw exception
        if (node == null)
            throw new IllegalStateException("Event selection loop fell through!");

        // get the haplotype from which the current haplotype arises
        int[] parentHaploArray = qsTree.getParentHaplo();
        int haplotypesParentHaplo = parentHaploArray[node.getNr()];
        int haplotypesGrandParentHaplo = -1;
        Double[] tempqstimes=qsTree.getAttachmentTimesList(node).clone();

        // get a random number deciding where the current haplo will be moved
        double u = Randomizer.nextDouble();
        double scalefactor=0;
        double toldQSstart=0;
        double tnewQSstart=0;

        // propose a new start time
        // reposition the event (i.e. haplotype's first attachment time) uniformly at random between tmin and tmax
        double tmin, tmax;
        double told=0;
        double tnew=0;
        tmin = node.getHeight();
        // get origin height
        double originheight = origin.getValue();

        // get a node above which the current haplotype arises
        QuasiSpeciesNode currentHaploNodeBelow = null;
        // get a node above which the parent haplotype arises
        QuasiSpeciesNode parentHaploNodeBelow = null;
        // get a node above which the grandparent haplotype arises
        QuasiSpeciesNode grandparentHaploNodeBelow = null;
        // if the parent haplotypes is null, tmax is the origin of the tree
        // We need to only change first position in AttachmentTimesList for the current node & scale the others
        //   startBranchCountsArray
        //   and possibly the haploName for 2 nodes
        //              +  the parentHaplo array
        //
// NO PARENT
        if (haplotypesParentHaplo == -1){

            // get a node above which the current haplotype arises
            QuasiSpeciesNode inode = node;
            while(currentHaploNodeBelow == null){
//            while (inode.getHaploAboveName() != node.getNr()){
                int currentnodeHaplo = inode.getHaploAboveName();
                if (currentnodeHaplo == node.getNr()){
                    currentHaploNodeBelow = inode;
                }
                inode = (QuasiSpeciesNode) inode.getParent();
            }
//            currentHaploNodeBelow = inode;

            // find out tmax
            tmax = originheight;
            // choose new time to attach
            if (tempqstimes.length >1)
                tnew = u*tmin + (1-u)*tmax;

            //logHastingsRatio = 0.0;

            // reposition attachment times: attach ((time - tmin) * (tnew/told)) + tmin
            // get the haplotype's starting time
            toldQSstart = tempqstimes[0];
            if (tempqstimes.length > 1){
                told = tempqstimes[1];
                // scale all the other positions in the array but the 0 position (haplo start time)
                // reposition attachment times: attach ((time - tmin) * (tnew/told)) + tmin
                scalefactor = (tnew - tmin)/(told - tmin);
                for (int i=1; i<tempqstimes.length; i++) {
                    tempqstimes[i] = ((tempqstimes[i] - tmin) * scalefactor) + tmin;
                }
                // set the haplotype's starting time to the new time
                double x = Randomizer.nextDouble();
                tnewQSstart = x*tempqstimes[1] + (1-x)*tmax;
                tempqstimes[0] = tnewQSstart;
                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                logHastingsRatio -= Math.log(tmax - told);
//                logHastingsRatio += Math.log(tmax - tempqstimes[1]);
            }
            else {
                // set the haplotype's starting time to the new time
                double x = Randomizer.nextDouble();
                tnewQSstart = x*tmin + (1-x)*tmax;
                tempqstimes[0] = tnewQSstart;
                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
                //              logHastingsRatio -= Math.log(tmax - tmin);
                //              logHastingsRatio += Math.log(tmax - tmin);
                // since they are identical, no need to change Hastings ratio
            }
            // rewrite the attachment times array
            qsTree.setAttachmentTimesList(node, tempqstimes);

            if (tempqstimes.length > 1){
                // assign contribution to the Hastings ratio for having different possible scales for told
                logHastingsRatio += Math.log(tmax/told - tmin/told);
                logHastingsRatio -= Math.log(tmax/tnew - tmin/tnew);
            }
        }
// PARENT
        // if there is a parent haplotype, we need to change more stuff (also for the parent haplotype)
        else{
            // check what is the grand parent haplotype
            haplotypesGrandParentHaplo = parentHaploArray[haplotypesParentHaplo];
// NO GRANDPARENT
            // if none, change only the arrays for the parent
            if (haplotypesGrandParentHaplo == -1){

                // get a node above which the parent haplotype arises
                QuasiSpeciesNode inode = node;
                while (parentHaploNodeBelow == null || currentHaploNodeBelow == null){
                    int currentnodeHaplo = inode.getHaploAboveName();
                    if (currentnodeHaplo == node.getNr()){
                        currentHaploNodeBelow = inode;
                    }
                    else if (currentnodeHaplo == haplotypesParentHaplo){
                        parentHaploNodeBelow = inode;
                    }
                    inode = (QuasiSpeciesNode) inode.getParent();
                }

                // find out tmax: grandparent haplotype is null, then tmax=origin
                tmax = originheight;
                // choose new time to attach
                if (tempqstimes.length >1)
                    tnew = u*tmin + (1-u)*tmax;
            }
// GRANDPARENT
            // otherwise change stuff also for the grandparent
            else{
                QuasiSpeciesNode inode = node;
                while (grandparentHaploNodeBelow== null || parentHaploNodeBelow == null || currentHaploNodeBelow == null){
                    int currentnodeHaplo = inode.getHaploAboveName();
                    if (currentnodeHaplo == node.getNr()){
                        currentHaploNodeBelow = inode;
                    }
                    else if (currentnodeHaplo == haplotypesParentHaplo){
                        parentHaploNodeBelow = inode;
                    }
                    else if (currentnodeHaplo == haplotypesGrandParentHaplo){
                        grandparentHaploNodeBelow = inode;
                    }
                    inode = (QuasiSpeciesNode) inode.getParent();
                }

                // keep track of common ancestor of parent, grandparent and current haplotype
                QuasiSpeciesNode lastGrandCommonAncestorNode;

                // find the last common ancestor of the chosen haplotype and its grandparent
                lastGrandCommonAncestorNode = findLastCommonAncestor(
                        node, (QuasiSpeciesNode) qsTree.getNode(haplotypesGrandParentHaplo), grandparentHaploNodeBelow);

                // find out tmax: restricted by the last common ancestor node of grandparent haplotype and chosen haplotype
                // the grandparent haplotype can be from the same subtree or from the other subtree
                // in either case (same or other subtree) the last common ancestor node is where the current haplo has tmax
                tmax = lastGrandCommonAncestorNode.getHeight();
                // choose new time to attach
                tnew = u*tmin + (1-u)*tmax;
            }
// PARENT AND GRANDPARENT

            // keep track of common ancestor of parent and current haplotype
            QuasiSpeciesNode lastParentCommonAncestorNode;

            // find the last common ancestor of the chosen haplotype and its parent
            lastParentCommonAncestorNode = findLastCommonAncestor(
                    node, (QuasiSpeciesNode) qsTree.getNode(haplotypesParentHaplo), parentHaploNodeBelow);

            // propose a new start time for the parent node uniformly at random between tminparent and tmaxparent
            double tminparent, tmaxparent;
            tminparent = qsTree.getNode(haplotypesParentHaplo).getHeight();
            tmaxparent = lastParentCommonAncestorNode.getHeight();

            // reposition attachment times: attach ((time - tmin) * (tnew/told)) + tmin
            // get the haplotype's starting time
            toldQSstart = tempqstimes[0];
            if (tempqstimes.length > 1){
                told = tempqstimes[1];
                // scale all the other positions in the array but the 0 position (haplo start time)
                // reposition attachment times: attach ((time - tmin) * (tnew/told)) + tmin
                scalefactor = (tnew - tmin)/(told - tmin);
                for (int i=1; i<tempqstimes.length; i++) {
                    tempqstimes[i] = ((tempqstimes[i] - tmin) * scalefactor) + tmin;
                }
                // set the haplotype's starting time to the new time
                double x = Randomizer.nextDouble();
                tnewQSstart = x*tempqstimes[1] + (1-x)*tmax;
                tempqstimes[0] = tnewQSstart;
                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                logHastingsRatio -= Math.log(tmaxparent - told);
//                logHastingsRatio += Math.log(tmax - tempqstimes[1]);
            }
            else {
                // set the haplotype's starting time to the new time
                double x = Randomizer.nextDouble();
                tnewQSstart = x*tmin + (1-x)*tmax;
                tempqstimes[0] = tnewQSstart;
                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
                //              logHastingsRatio -= Math.log(tmaxparent - tmin);
                //              logHastingsRatio += Math.log(tmax - tmin);
                // since they are identical, no need to change Hastings ratio
            }
            // rewrite the attachment times array
            qsTree.setAttachmentTimesList(node, tempqstimes);

            double tnewparent=0;
            double toldparent=0;
            // move the parent haplotype at random below the common ancestor node
            // ! only if chosen haplo moves above last common ancestor of parent and current haplo

            if (tnewQSstart > tmaxparent) {
                // choose new time to attach for the moved parent node and scale all attachment times of this haplo
                double v = Randomizer.nextDouble();
                // choose new time to attach
                tnewparent = v*tminparent + (1-v)*tmaxparent;

                // reposition parent haplo attachment times: attach ((time - tminparent) * (tnew/told)) + tminparent
                Double[] tempqstimesparent=qsTree.getAttachmentTimesList(haplotypesParentHaplo).clone();
                // get the haplotype's starting time
                double toldQSstartparent = tempqstimesparent[0];
                double grandscalefactor=0;
                double tnewQSstartparent;
                if (tempqstimesparent.length > 1){
                    toldparent = tempqstimesparent[1];
                    // scale all the other positions in the array but the 0 position (haplo start time)
                    grandscalefactor = (tnewparent - tminparent)/(toldparent - tminparent);
                    for (int i=1; i<tempqstimesparent.length; i++) {
                        tempqstimesparent[i] = ((tempqstimesparent[i] - tminparent) * grandscalefactor) + tminparent;
                    }
                    // set the haplotype's starting time to the new time
                    double y = Randomizer.nextDouble();
                    tnewQSstartparent = y*tempqstimesparent[1] + (1-y)*tmaxparent;
                    tempqstimesparent[0] = tnewQSstartparent;
                    // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                    logHastingsRatio -= Math.log(tmax - toldparent);
//                    logHastingsRatio += Math.log(tmaxparent - tempqstimesparent[1]);
                }
                else {
                    // set the haplotype's starting time to the new time
                    double x = Randomizer.nextDouble();
                    tnewQSstartparent = x*tminparent + (1-x)*tmaxparent;
                    tempqstimesparent[0] = tnewQSstartparent;
                    // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
                    //              logHastingsRatio -= Math.log(tmax - tminparent);
                    //              logHastingsRatio += Math.log(tmaxparent - tminparent);
                    // since they are identical, no need to change Hastings ratio
                }

                // rewrite the attachment times array
                qsTree.setAttachmentTimesList(haplotypesParentHaplo, tempqstimesparent);

                // set the log(Hastings ratio) Pr(backwardmove)/Pr(forwardmove)=(1/sthback)/(1/sthforth)=sthforth/sthback
//                if (tempqstimes.length > 1 && tempqstimesparent.length > 1)
//                    logHastingsRatio += Math.log(tmax - tmin) + Math.log(tmaxparent - tminparent) -
//                                        Math.log(tmax - tminparent) - Math.log(tmaxparent - tmin);
//                else if (tempqstimes.length > 1)
//                    logHastingsRatio += Math.log(tmax - tmin) - Math.log(tmaxparent - tmin);
//                else if (tempqstimesparent.length > 1)
//                    logHastingsRatio += Math.log(tmaxparent - tminparent) - Math.log(tmax - tminparent);
                // Incorporate the probability of scaling all the attachment times
                if (tempqstimes.length > 1){
                    // assign contribution to the Hastings ratio for having different possible scales for told
                    logHastingsRatio += Math.log(tmax/told - tmin/told);
                    logHastingsRatio -= Math.log(tmaxparent/tnew - tmin/tnew);
                }

                if (tempqstimesparent.length > 1){
                    // assign contribution to the Hastings ratio for having different possible scales for toldparent
                    logHastingsRatio += Math.log(tmaxparent/toldparent - tminparent/toldparent);
                    logHastingsRatio -= Math.log(tmax/tnewparent - tminparent/tnewparent);
                    // assign contribution of each scaled attachment time
                    logHastingsRatio -= 2 * (Math.log (grandscalefactor * (toldparent - tminparent) + tminparent) - Math.log(toldparent));
                    logHastingsRatio += (tempqstimesparent.length-1) * Math.log(grandscalefactor);
                }

                // change internal node haploName - both remove and add if necessary
                Node intnode = parentHaploNodeBelow;
                Node intnodeparent;
                if (tnewQSstartparent > toldQSstartparent){
//                    while (intnode.getHeight() < tmax && intnode.getHeight() < tnewparent){
                    intnodeparent = intnode.getParent();
                    while (intnodeparent != null && intnodeparent.getHeight() < tnewQSstartparent){
                        intnode = intnodeparent;
                        intnodeparent = intnodeparent.getParent();
                    }
                }
                else {
                    intnode = qsTree.getNode(haplotypesParentHaplo);
                    intnodeparent = intnode.getParent();
                    while (intnodeparent.getHeight() < tnewQSstartparent){
                        intnode = intnodeparent;
                        intnodeparent = intnodeparent.getParent();
                    }
                }
                if (intnode != parentHaploNodeBelow){
                    // change the internal nodes' haploName, where necessary
                    parentHaploNodeBelow.setHaploAboveName(-1);

                    // change the internal nodes' haploName, where necessary
                    ((QuasiSpeciesNode) intnode).setHaploAboveName(haplotypesParentHaplo);

                    // recalculate parentHaplo array, and re-set continuingHaploName -- > done if node itself moved
                }
                // in any case (if changed or not the parentHaplo array) recalculate countPossibleStartBranches -- > done if node itself moved
                // add probability that tnewQSstart > tmaxparent
                if (tempqstimesparent.length > 1)
                logHastingsRatio += Math.log(tmax-tmaxparent) - Math.log(tmax-tmin)
                                  - (Math.log(tmax-tmaxparent) - Math.log(tmax-tminparent));
            }
            // if the tnew is not above the lastParentCommonAncestorNode, then just reposition the current haplo start
            else {
                if (tempqstimes.length > 1)
                    logHastingsRatio += Math.log(tmaxparent-tmin) - Math.log(tmax-tmin)
                                      - (Math.log(tmaxparent-tminparent) - Math.log(tmax-tminparent));
            //    //logHastingsRatio = 0.0;
            }
        }
// NODE ITSELF
// FOR ALL CASES: reposition the HAPLOTYPE START ITSELF

        // Incorporate the probability of scaling all the attachment times
        if (tempqstimes.length > 1){
            // assign contribution of each scaled attachment time
            logHastingsRatio -= 2 * (Math.log (scalefactor * (told - tmin) + tmin) - Math.log(told));
            logHastingsRatio += (tempqstimes.length-1) * Math.log(scalefactor);
        }

        // change internal node haploName - both remove and add if necessary
        Node intnode = currentHaploNodeBelow;
        Node intnodeparent;
        if (tnewQSstart > toldQSstart) {
            intnodeparent = intnode.getParent();
            while (intnodeparent!=null && intnodeparent.getHeight() < tnewQSstart){
                intnode = intnodeparent;
                intnodeparent = intnodeparent.getParent();
            }
        }
        else {
            intnode = node;
            intnodeparent = intnode.getParent();
            while (intnodeparent != null && intnodeparent.getHeight() < tnewQSstart){
                intnode = intnodeparent;
                intnodeparent = intnodeparent.getParent();
            }
        }
        if (intnode != currentHaploNodeBelow){
            // change the internal nodes' haploName, where necessary
            currentHaploNodeBelow.setHaploAboveName(-1);

            ((QuasiSpeciesNode) intnode).setHaploAboveName(node.getNr());

            // recalculate parentHaplo array, re-assign continuingHaploName
            if (grandparentHaploNodeBelow!=null){
                recalculateParentHaploAndCorrectContinuingHaploName(parentHaploArray[haplotypesGrandParentHaplo], grandparentHaploNodeBelow);
            }
            else {
                recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) qsTree.getRoot());
            }
         }
        // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
        int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
        qsTree.setStartBranchCounts(startBranchCountsArray);

    // RETURN log(HASTINGS RATIO)
    return logHastingsRatio; // proper hastings ratio!!!

   }

}
