package quasispeciestree.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;
import beast.util.Randomizer;

import java.util.ArrayList;

/**
 *  @author Veronika Boskova created on 06/08/2015
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
        final QuasiSpeciesTree qsTree = quasiSpeciesTreeInput.get();

        // mark the tree as dirty (startEditing)
        qsTree.startEditing(null);

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
        ArrayList<QuasiSpeciesNode> parentHaploArray = qsTree.getParentHaplo();
        QuasiSpeciesNode haplotypesParentHaplo = parentHaploArray.get(node.getNr());
//        QuasiSpeciesNode haplotypesParentHaplo = qsTree.getParentHaplo(node.getNr());
        QuasiSpeciesNode haplotypesGrandParentHaplo = null;

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
        QuasiSpeciesNode nodeBelow = null;

        // if the parent haplotypes is null, tmax is the origin of the tree
        // We need to only change first position in AttachmentTimesList for the current node & scale the others
        //   startBranchCountsArray
        //   and possibly the haploName for 2 nodes
        //                the parentHaplo/aboveNodeHaplo arrays
        //
// NO PARENT
        if (haplotypesParentHaplo == null){

            // get a node above which the current haplotype arises
            for (Node inode : qsTree.getNodesAsArray()){
                if (((QuasiSpeciesNode) inode).getHaploName() == node.getID()) {
                    nodeBelow = (QuasiSpeciesNode) inode;
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
            // get a node above which the parent haplotype arises
            QuasiSpeciesNode parentnodeBelow = null;

            // check what is the grand parent haplotype
            haplotypesGrandParentHaplo = parentHaploArray.get(haplotypesParentHaplo.getNr());
//            haplotypesGrandParentHaplo = qsTree.getParentHaplo(haplotypesParentHaplo.getNr());
// NO GRANDPARENT
            // if none, change only the arrays for the parent
            if (haplotypesGrandParentHaplo == null){

                // find out tmax: grandparent haplotype is null, then tmax=origin
                tmax = originheight;
                // choose new time to attach
                tnew = u*tmin + (1-u)*tmax;

                // get a node above which the parent haplotype arises
                for (Node inode : qsTree.getNodesAsArray()){
                    String currentnodeHaplo = ((QuasiSpeciesNode) inode).getHaploName();
                    if (currentnodeHaplo == node.getID()) {
                        nodeBelow = (QuasiSpeciesNode) inode;
                    }
                    if (currentnodeHaplo == haplotypesParentHaplo.getID()){
                        parentnodeBelow = (QuasiSpeciesNode) inode;
                    }
                    if (parentnodeBelow != null && nodeBelow != null){
                        break;
                    }
                }
            }
// GRANDPARENT
            // otherwise change stuff also for the grandparent
            else{
                // get a node above which the grandparent haplotype arises
                QuasiSpeciesNode grandnodeBelow = null;
                for (Node inode : qsTree.getNodesAsArray()){
                    String currentnodeHaplo = ((QuasiSpeciesNode) inode).getHaploName();
                    if (currentnodeHaplo == node.getID()) {
                        nodeBelow = (QuasiSpeciesNode) inode;
                    }
                    else if  (currentnodeHaplo == haplotypesParentHaplo.getID()){
                        parentnodeBelow = (QuasiSpeciesNode) inode;
                    }
                    else if (currentnodeHaplo == haplotypesGrandParentHaplo.getID()){
                        grandnodeBelow = (QuasiSpeciesNode) inode;
                    }
                    if (grandnodeBelow!= null && parentnodeBelow != null && nodeBelow != null){
                        break;
                    }
                }

                // keep track of common ancestor of parent, grandparent and current haplotype
                QuasiSpeciesNode lastGrandCommonAncestorNode;

                // find the last common ancestor of the chosen haplotype and its grandparent
                lastGrandCommonAncestorNode = findLastCommonAncestor(node, haplotypesGrandParentHaplo, grandnodeBelow);

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
            lastParentCommonAncestorNode = findLastCommonAncestor(node, haplotypesParentHaplo, parentnodeBelow);

            // propose a new start time for the parent node uniformly at random between tminparent and tmaxparent
            double tminparent, tmaxparent;
            tminparent = haplotypesParentHaplo.getHeight();
            tmaxparent = lastParentCommonAncestorNode.getHeight();

            // move the parent haplotype at random below the common ancestor node
            // ! only if chosen haplo moves above last common ancestor of parent and current haplo
            if (tnew > tmaxparent) {
                // choose new time to attach for the moved parent node and scale all attachment times of this haplo
                double v = Randomizer.nextDouble();
                // choose new time to attach
                double tnewparent = v*tminparent + (1-v)*tmaxparent;

                // set the log(Hastings ratio)
                logHastingsRatio =  Math.log(((tmax - tminparent)*(tmaxparent - tmin)) /
                        ((tmax - tmin)*(tmaxparent - tminparent)));

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
                    intnode = haplotypesParentHaplo;
                    intnodeparent = intnode.getParent();
                    while (intnodeparent.getHeight() < tnewparent){
                        intnode = intnodeparent;
                        intnodeparent = intnodeparent.getParent();
                    }
//                }
                if (intnode != parentnodeBelow){
                    // change the internal nodes' haploName, where necessary
                    parentnodeBelow.setHaploName(null);

                    // change the internal nodes' haploName, where necessary
                    ((QuasiSpeciesNode) intnode).setHaploName(haplotypesParentHaplo.getID());

                    // recalculate parentHaplo array, aboveNodeHaplo array -- > done if node itself moved
                }
               // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches -- > done if node itself moved
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
        Node intnode = nodeBelow;
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
        if (intnode != nodeBelow){
            // change the internal nodes' haploName, where necessary
            nodeBelow.setHaploName(null);

            ((QuasiSpeciesNode) intnode).setHaploName(node.getID());

            // recalculate parentHaplo array, aboveNodeHaplo array
            ArrayList<QuasiSpeciesNode> newParentHaplo = new ArrayList<>(qsTree.getLeafNodeCount());
            for (int i=0; i<qsTree.getLeafNodeCount();i++){
                newParentHaplo.add(i,null);
            }
            ArrayList<QuasiSpeciesNode> newAboveNodeHaplo = new ArrayList<>(qsTree.getInternalNodeCount());
            for (int i=0; i<qsTree.getInternalNodeCount();i++){
                newAboveNodeHaplo.add(i,null);
            }

            qsTree.findParentHaploAndAboveNodeHaplo(null, (QuasiSpeciesNode) qsTree.getRoot(), newParentHaplo, newAboveNodeHaplo);

            // clear and re-set the original arrays
//            parentHaploArray.clear();
//            parentHaploArray.addAll(newParentHaplo);
            qsTree.getParentHaplo().clear();
            qsTree.getParentHaplo().addAll(newParentHaplo);
//            ArrayList<QuasiSpeciesNode> aboveNodeHaploArray = qsTree.getAboveNodeHaplo();
//            aboveNodeHaploArray.clear();
//            aboveNodeHaploArray.addAll(newAboveNodeHaplo);
            qsTree.getAboveNodeHaplo().clear();
            qsTree.getAboveNodeHaplo().addAll(newAboveNodeHaplo);
        }
        // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
        int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
        qsTree.setStartBranchCounts(startBranchCountsArray);

    // RETURN log(HASTINGS RATIO)
    return logHastingsRatio; // proper hastings ratio!!!

   }

}
