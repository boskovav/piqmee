package quasispeciestree.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;
import beast.util.Randomizer;

import java.util.ArrayList;
import java.util.Arrays;

/**
 *  @author Veronika Boskova created on 06/08/2015 finished on 31.01.2016
 */
@Description("Chooses a haplotype at random and moves"
        + "its start time uniformly in interval "
        + "restricted by the grand-parent haplotype"
        +" and tip time of the same haplotype.")
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

        // get a random number deciding where the current haplo will be moved
        double u = Randomizer.nextDouble();
        // weight the attachment height by the total haplo count wrt to all others
//        double weight = qsTree.getHaplotypeCounts(node)/qsTree.getTotalAttachmentCounts();
//        u = u*weight;

        // propose a new start time
        // reposition the event (i.e. haplotype starting time) uniformly at random between tmin and tmax
        double tmin, tmax, tnew;
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
            for (Node inode : qsTree.getNodesAsArray()){
                if (((QuasiSpeciesNode) inode).getHaploAboveName() == node.getNr()) {
                    currentHaploNodeBelow = (QuasiSpeciesNode) inode;
                    break;
                }
            }

            // find out tmax
            tmax = originheight;
            // choose new time to attach
            tnew = u*tmin + (1-u)*tmax;

            //logHastingsRatio = 0.0;
        }
// PARENT
        // if there is a parent haplotype, we need to change more stuff (also for the parent haplotype)
        else{
            // check what is the grand parent haplotype
            haplotypesGrandParentHaplo = parentHaploArray[haplotypesParentHaplo];
// NO GRANDPARENT
            // if none, change only the arrays for the parent
            if (haplotypesGrandParentHaplo == -1){

                // find out tmax: grandparent haplotype is null, then tmax=origin
                tmax = originheight;
                // choose new time to attach
                tnew = u*tmin + (1-u)*tmax;

                // get a node above which the parent haplotype arises
                for (Node inode : qsTree.getNodesAsArray()){
                    int currentnodeHaplo = ((QuasiSpeciesNode) inode).getHaploAboveName();
                    if (currentnodeHaplo == node.getNr()) {
                        currentHaploNodeBelow = (QuasiSpeciesNode) inode;
                    }
                    if (currentnodeHaplo == haplotypesParentHaplo){
                        parentHaploNodeBelow = (QuasiSpeciesNode) inode;
                    }
                    if (parentHaploNodeBelow != null && currentHaploNodeBelow != null){
                        break;
                    }
                }
            }
// GRANDPARENT
            // otherwise change stuff also for the grandparent
            else{
                for (Node inode : qsTree.getNodesAsArray()){
                    int currentnodeHaplo = ((QuasiSpeciesNode) inode).getHaploAboveName();
                    if (currentnodeHaplo == node.getNr()) {
                        currentHaploNodeBelow = (QuasiSpeciesNode) inode;
                    }
                    else if  (currentnodeHaplo == haplotypesParentHaplo){
                        parentHaploNodeBelow = (QuasiSpeciesNode) inode;
                    }
                    else if (currentnodeHaplo == haplotypesGrandParentHaplo){
                        grandparentHaploNodeBelow = (QuasiSpeciesNode) inode;
                    }
                    if (grandparentHaploNodeBelow!= null && parentHaploNodeBelow != null && currentHaploNodeBelow != null){
                        break;
                    }
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

            // move the parent haplotype at random below the common ancestor node
            // ! only if chosen haplo moves above last common ancestor of parent and current haplo
            if (tnew > tmaxparent) {
                // choose new time to attach for the moved parent node and scale all attachment times of this haplo
                double v = Randomizer.nextDouble();
                // choose new time to attach
                double tnewparent = v*tminparent + (1-v)*tmaxparent;

                // set the log(Hastings ratio) Pr(backwardmove)/Pr(forwardmove)=(1/sthback)/(1/sthforth)=sthforth/sthback
                logHastingsRatio =  Math.log(
                        ((tmax - tmin)*(tmaxparent - tminparent)) /
                        ((tmax - tminparent)*(tmaxparent - tmin)) );

                // reposition parent haplo attachment times: attach ((time - tminparent) * (tnew/told)) + tminparent
                Double[] tempqstimes=qsTree.getAttachmentTimesList(haplotypesParentHaplo).clone();
                // get the haplotype's starting time
                double toldparent = tempqstimes[0];
                // scale all the other positions in the array but the 0 position (haplo start time)
                double scalefactor = (tnewparent - tminparent)/(toldparent - tminparent);
                for (int i=1; i<tempqstimes.length; i++) {
                    tempqstimes[i] = ((tempqstimes[i] - tminparent) * scalefactor) + tminparent;
                }
                // set the haplotype's starting time to the new time
                tempqstimes[0] = tnewparent;
                // rewrite the attachment times array
                qsTree.setAttachmentTimesList(haplotypesParentHaplo, tempqstimes);


                // change internal node haploName - both remove and add if necessary
                Node intnode;
                Node intnodeparent;
//                if (tnewparent > tmaxparent){
////                    while (intnode.getHeight() < tmax && intnode.getHeight() < tnewparent){
//                    intnodeparent = intnode.getParent();
//                    while (intnodeparent != null && intnodeparent.getHeight() < tnewparent){
//                        intnode = intnodeparent;
//                        intnodeparent = intnodeparent.getParent();
//                    }
//                }
//                else {
                    intnode = qsTree.getNode(haplotypesParentHaplo);
                    intnodeparent = intnode.getParent();
                    while (intnodeparent.getHeight() < tnewparent){
                        intnode = intnodeparent;
                        intnodeparent = intnodeparent.getParent();
                    }
//                }
                if (intnode != parentHaploNodeBelow){
                    // change the internal nodes' haploName, where necessary
                    parentHaploNodeBelow.setHaploAboveName(-1);

                    // change the internal nodes' haploName, where necessary
                    ((QuasiSpeciesNode) intnode).setHaploAboveName(haplotypesParentHaplo);

                    // recalculate parentHaplo array, and re-set continuingHaploName -- > done if node itself moved
                }
               // in any case (if changed or not the parentHaplo array) recalculate countPossibleStartBranches -- > done if node itself moved
            }
            // if the tnew is not above the lastParentCommonAncestorNode, then just reposition the current haplo start
            //else {
            //    //logHastingsRatio = 0.0;
            //}
        }
// NODE ITSELF
// FOR ALL CASES: reposition the HAPLOTYPE START ITSELF
        // reposition attachment times: attach ((time - tmin) * (tnew/told)) + tmin
        Double[] tempqstimes=qsTree.getAttachmentTimesList(node).clone();
        // get the haplotype's starting time
        double told = tempqstimes[0];
        // scale all the other positions in the array but the 0 position (haplo start time)
        double scalefactor = (tnew - tmin)/(told - tmin);
        for (int i=1; i<tempqstimes.length; i++) {
            tempqstimes[i] = ((tempqstimes[i] - tmin) * scalefactor) + tmin;
        }
        // set the haplotype's starting time to the new time
        tempqstimes[0] = tnew;
        // rewrite the attachment times array
        qsTree.setAttachmentTimesList(node, tempqstimes);


        // change internal node haploName - both remove and add if necessary
        Node intnode = currentHaploNodeBelow;
        Node intnodeparent;
        if (tnew > told) {
//                while (intnode.getHeight() < tmax && intnode.getHeight() < tnew){
            intnodeparent = intnode.getParent();
            while (intnodeparent!=null && intnodeparent.getHeight() < tnew){
                intnode = intnodeparent;
                intnodeparent = intnodeparent.getParent();
            }
        }
        else {
            intnode = node;
            intnodeparent = intnode.getParent();
            while (intnodeparent != null && intnodeparent.getHeight() < tnew){
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
