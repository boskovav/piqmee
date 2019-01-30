package piqmee.operators;

import beast.core.Description;
import beast.util.Randomizer;
import piqmee.tree.QuasiSpeciesNode;
import piqmee.tree.QuasiSpeciesTree;

import java.util.ArrayList;


/**
 *  @author Veronika Boskova created on 09/06/2016 finished on 27/03/2017
 */
@Description("Chooses a haplotype at random and moves"
        + "its first attachment time uniformly in interval "
        + "restricted by the parent haplotype (start)"
        +" and tip time of the moved haplotype,"
        +" and scales all the attachment times accordingly.")
public class QuasiSpeciesHaplotypeScale extends QuasiSpeciesTreeOperator{

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        if (qsTree.getTotalAttachmentCounts()==0){
            System.out.println("In QuasiSpeciesHaploScale operator --- "
                    +"there are no QS duplicates. The QuasiSpeciesHaplotypeScale "
                    +"operator cannot be used. Please remove it from your xml file.");
        }
    }

    /**
     * Change the start time and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {

        if (qsTree.getTotalAttachmentCounts()==0){
            return 0.0;
        }

        // keep track of the Hastings ratio
        double logHastingsRatio = 0.0;
        // Randomly select event on tree:
        // weighted by the number of events (i.e. count of each haplotype??)
        // for now do random uniform from just the haplotype count (disregarding the counts)
        QuasiSpeciesNode node = null;
        do {
            node = (QuasiSpeciesNode) qsTree.getNode(Randomizer.nextInt(qsTree.getLeafNodeCount()));

        } while (node.getAttachmentTimesList().length < 2);

        // if we did not assign the node at this stage, throw exception
        if (node == null)
            throw new IllegalStateException("Event selection loop fell through!");

        int haplo = node.getNr();

        // get the old first attachment time
        double toldQSstart = node.getAttachmentTimesList()[0];

        // get a node above which the current haplotype arises
        QuasiSpeciesNode oldNodeBelowHaploMoved = findNodeBelowThisHaplo(node,haplo);

        // reposition the event (i.e. haplotype's first attachment time) uniformly at random between tmin and tmax
        // what is the parent first attachment time?
        ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTime(oldNodeBelowHaploMoved, haplo, 0);
        int parentHaplo = (int) haploStartMaxNewArray.get(2);

        // scale the haplotype
        double maxTime = (double)haploStartMaxNewArray.get(1);
        double logHastingsRatioContribution = scaleThisHaplo(node, maxTime, node.getHeight(), 0,
                                                                    maxTime, node.getHeight(),0);
        if (logHastingsRatioContribution == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;
        else logHastingsRatio += logHastingsRatioContribution;

        // get the new first attachment time
        double tnewQSstart = node.getAttachmentTimesList()[0];

        // scale all QS start points for all haplotypes on a way from toldQSstart to tnewQSstart
        // find the new node above which the haplo arises
        QuasiSpeciesNode nodeBelowHaploMoved;
        if (tnewQSstart > toldQSstart) {
            nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(oldNodeBelowHaploMoved,tnewQSstart);
        }
        else {
            nodeBelowHaploMoved=(QuasiSpeciesNode)findNodeBelowAfterRepositioningHaploStart(node,tnewQSstart);
        }

        if (oldNodeBelowHaploMoved.getNr()!=nodeBelowHaploMoved.getNr()){
            // change internal node haploName - both remove and add if necessary
            // change the internal nodes' haploName, where necessary
            oldNodeBelowHaploMoved.setHaploAboveName(-1);
            nodeBelowHaploMoved.setHaploAboveName(node.getNr());
            //recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode)qsTree.getRoot());

            if (oldNodeBelowHaploMoved.getHeight() > nodeBelowHaploMoved.getHeight()){
                if (parentHaplo == -1)
                    recalculateParentHaploAndCorrectContinuingHaploName(parentHaplo, oldNodeBelowHaploMoved);
                else {
                    int grandParentHaplo = ((QuasiSpeciesNode) qsTree.getNode(parentHaplo)).getParentHaplo();
                    QuasiSpeciesNode oldNode = findNodeBelowThisHaplo(oldNodeBelowHaploMoved, parentHaplo);
                    recalculateParentHaploAndCorrectContinuingHaploName(grandParentHaplo, oldNode);
                }
            }
            else
                recalculateParentHaploAndCorrectContinuingHaploName(parentHaplo, nodeBelowHaploMoved);
        }

        // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
        qsTree.countAndSetPossibleStartBranches();

        // Ensure BEAST knows to recalculate affected likelihood:
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);

    // RETURN log(HASTINGS RATIO)
    return logHastingsRatio; // proper hastings ratio!!!
    }

}
