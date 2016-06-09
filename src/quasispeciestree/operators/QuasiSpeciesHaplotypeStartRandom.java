package quasispeciestree.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;

import java.util.ArrayList;

/**
 *  @author Veronika Boskova created on 07/04/2016
 */

@Description("For a given haplotype, randomly and uniformly selects start point "
        + "restricted by the parent haplotype (or origin if parent haplo -1)"
        + "and the first haplotype attachment time.")
public class QuasiSpeciesHaplotypeStartRandom extends QuasiSpeciesTreeOperator{

    // THIS MOVE WILL ALWAYS BE ACCEPTED IF ONLY TREE PRIOR IS INCLUDED

    /**
     * Change the attachment time and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {
        // Randomly select haplotype:
        // weighted by the number of events (i.e. count of each haplotype)
        int event = Randomizer.nextInt(qsTree.getTotalAttachmentCounts()+qsTree.getLeafNodeCount());


        QuasiSpeciesNode node = null;

        // find the haplotype and the sequence position corresponding to the event number
        // backward search...
        for (Node thisNode : qsTree.getExternalNodes()) {
            double tempHaploCount = qsTree.getHaplotypeCounts((QuasiSpeciesNode) thisNode)+1;
            if (event<tempHaploCount) {
                node = (QuasiSpeciesNode)thisNode;
                break;
            }
            event -= tempHaploCount;
        }


        // if we did not assign the node at this stage, throw exception
        if (node == null)
            throw new IllegalStateException("Event selection loop fell through!");

        // reposition the event (i.e. haplotype sequence changeIdx attachment time)
        double tmin, tmax;
        Double[] tempqstimes=qsTree.getAttachmentTimesList(node).clone();

        // choose new max index between 0 - #sequences of this haplotype
        if (qsTree.getHaplotypeCounts(node)>0)
            tmin = tempqstimes[1];
        else
            tmin = node.getHeight();
        ArrayList haploStartMaxNewArray = getMaxPossibleHaploAttachTimeForQSStart(node, node.getNr(), 0);
        tmax = (double) haploStartMaxNewArray.get(1);
//        tmax = getMaxPossibleHaploAttachTime(node.getNr(), node);
        double u = Randomizer.nextDouble();
        double tnewQSstart = u*tmin + (1-u)*tmax;

        double toldQSstart = tempqstimes[0];

        tempqstimes[0]=tnewQSstart;

        qsTree.setAttachmentTimesList(node, tempqstimes);

        QuasiSpeciesNode nodeBelowHaploMoved = node;
        while (nodeBelowHaploMoved != qsTree.getRoot() && tnewQSstart > nodeBelowHaploMoved.getParent().getHeight()){
            nodeBelowHaploMoved = (QuasiSpeciesNode) nodeBelowHaploMoved.getParent();
        }
        QuasiSpeciesNode oldNodeBelowHaploMoved = node;
        while (oldNodeBelowHaploMoved.getHaploAboveName()!=node.getNr()){
            oldNodeBelowHaploMoved = (QuasiSpeciesNode) oldNodeBelowHaploMoved.getParent();
        }
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

        // recalculate countPossibleStartBranches
        int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
        qsTree.setStartBranchCounts(startBranchCountsArray);

        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);

//        if (tnew > told)
            return 0.0;
//        else
//            return Math.log(0.05);

    }


}
