package piqmee.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import piqmee.tree.QuasiSpeciesNode;
import piqmee.tree.QuasiSpeciesTree;

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
            if (((QuasiSpeciesNode)qsTree.getNode(i)).getParentHaplo() != -1 && ((QuasiSpeciesNode)qsTree.getNode(i)).getAttachmentTimesList().length > 1) {
                haplowithparents += 1;
                break;
            }
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

        // get a node above which the current/parent/grand parent haplotype arises
        QuasiSpeciesNode oldNodeBelowCurrentHaplo = findNodeBelowThisHaplo(node,haplo);
        QuasiSpeciesNode oldNodeBelowParentHaplo = findNodeBelowThisHaplo(node,parentHaplo.getNr());

        // reposition the event (i.e. haplotype's first attachment time) uniformly at random between tmin and tmax
        // what is the parent haplo start time?
        ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTime(oldNodeBelowCurrentHaplo, haplo, 0);
        double maxTimeParent = (double) haploStartMaxNewArray.get(1);

        // get the maximum possible attachment time for the haplo
        double maxTime;
        // NO GRANDPARENT
        // grandparent haplotype is null, then tmax=origin
        if (grandParentHaplo == null){
            maxTime = origin.getValue();
        }
        // GRANDPARENT
        // otherwise tmax is the common ancestor of grandParentHaplo and parentHaplo
        else{
            // get a node above which the grandparent haplotype arises
            maxTime = (double) getMaxPossibleHaploAttachTime(
                    (QuasiSpeciesNode)haploStartMaxNewArray.get(0),parentHaplo.getNr(),0).get(1);
        }

        // find the mrca to determine the min boundary for the haplo
        QuasiSpeciesNode mrca = findLastCommonAncestor(node, parentHaplo, (QuasiSpeciesNode) qsTree.getRoot());

        // scale haplo
        double logHastingsRatioContribution = scaleThisHaplo(node, maxTime, mrca.getHeight(), node.getAttachmentTimesList()[1],
                                                                   maxTimeParent, node.getHeight(), 0);
        if (logHastingsRatioContribution == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;
        else logHastingsRatio += logHastingsRatioContribution;

        // check that the newly scaled attachment times do not go over the boundaries
        if (node.getAttachmentTimesList()[0] < (double) haploStartMaxNewArray.get(1))
            return Double.NEGATIVE_INFINITY;


        // the scale haplo of does not know what the min is on the way back... is it the mrca or the node height?



        // scale also parent haplo
        logHastingsRatioContribution = scaleThisHaplo(parentHaplo, maxTimeParent, parentHaplo.getHeight(), 0,
                                                                   maxTime, mrca.getHeight(), -1);
        if (logHastingsRatioContribution == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;
        else logHastingsRatio += logHastingsRatioContribution;


        // find out what nodes have been passed by changing the haplo height
        QuasiSpeciesNode nodeBelowHaploMoved = (QuasiSpeciesNode) findNodeBelowAfterRepositioningHaploStart(oldNodeBelowCurrentHaplo,node.getAttachmentTimesList()[0]);

        // find the new node above which the haplo arises
        QuasiSpeciesNode nodeBelowHaploParent = (QuasiSpeciesNode) findNodeBelowAfterRepositioningHaploStart(parentHaplo,parentHaplo.getAttachmentTimesList()[0]);


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
        //qsTree.countAndSetPossibleStartBranches();

        // Ensure BEAST knows to recalculate affected likelihood:
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        parentHaplo.makeDirty(QuasiSpeciesTree.IS_FILTHY);

    // RETURN log(HASTINGS RATIO)
    return logHastingsRatio; // proper hastings ratio!!!
    }

}
