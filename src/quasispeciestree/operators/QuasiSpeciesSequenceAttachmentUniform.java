package quasispeciestree.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;
import beast.util.Randomizer;

/**
 *  @author Veronika Boskova created on 27/07/2015
 */
@Description("Within given haplotype, randomly selects one sequence "
        + "attachment time and moves it uniformly in interval "
        + "restricted by the closest previous and next attachment times "
        + "of another sequence from the same haplotype.")
public class QuasiSpeciesSequenceAttachmentUniform extends QuasiSpeciesTreeOperator{

//    public Input<Boolean> includeRootInput = new Input<>("includeRoot",
//            "Allow modification of root node.", false);

//    public Input<Double> rootScaleFactorInput = new Input<>("rootScaleFactor",
//            "Root scale factor.", 0.9);

    /**
     * Change the attachment time and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {
//        final QuasiSpeciesTree qsTree = quasiSpeciesTreeInput.get(this);
        // Randomly select event on tree:
        // weighted by the number of events (i.e. count of each haplotype)
        int event = Randomizer.nextInt(qsTree.getTotalAttachmentCounts());


        QuasiSpeciesNode node = null;

        // index for the attachment time chosen to change
        int changeIdx = -1;

        // find the haplotype and the sequence position corresponding to the event number
        // backward search...
        for (Node thisNode : qsTree.getExternalNodes()) {
            double tempHaploCount = qsTree.getHaplotypeCounts((QuasiSpeciesNode) thisNode);
            if (event<tempHaploCount) {
                node = (QuasiSpeciesNode)thisNode;
                // change index is event +1 as our arrays are 1-#haplotype repetition instances
                // and position 0 in the array is the haplotype starting point
                // TODO haplotype starting time changes by another operator...
                changeIdx = event+1;
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

        // as we choose index changeIdx between 1 - #sequences of this haplotype then changeIdx-1=0 at minimum
        // TODO change to only tmax = tempqstimes[changeIdx-1]; when testing is done
//        if (changeIdx-1==0){
//            tmax = 16.83;
//        } else {
            tmax = tempqstimes[changeIdx-1];
//        }
//        if (changeIdx+1>qsTree.getHaplotypeCounts((QuasiSpeciesNode) node))
        if (changeIdx+1>(tempqstimes.length-1))
            tmin = node.getHeight();
//        else if (changeIdx-1==0)
//            tmin=16.83;
        else
            tmin = tempqstimes[changeIdx+1];

        double u = Randomizer.nextDouble();
        double tnew = u*tmin + (1-u)*tmax; // = u*(tmin-tmax)+tmax =u*(tmin-(tmin+tmax-tmin))+tmin+tmax-tmin
                                            // invert u and 1-u and have the same as u(tmax-tmin)+tmin
                                            // (1-u)*tmin + u*tmax-u*(tmax-tmin)-tmin=0 ???
                                            // indeed tmin-u*tmin+u*tmin-u*tmax+u*tmin-tmin=0

        tempqstimes[changeIdx]=tnew;
        qsTree.setAttachmentTimesList(node, tempqstimes);

        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);

        return 0.0;

    }

}
