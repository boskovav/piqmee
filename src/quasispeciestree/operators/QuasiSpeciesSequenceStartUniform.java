package quasispeciestree.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;
import beast.util.Randomizer;

/**
 *  @author Veronika Boskova created on 06/08/2015
 */
@Description("Choose a halotype at random and moves"
        + "its start time uniformly in interval "
        + "restricted by the closest previous haplotype"
        +" and next attachment times "
        + "of another sequence from the same haplotype.")
public class QuasiSpeciesSequenceStartUniform extends QuasiSpeciesTreeOperator{

    /**
     * Change the start time and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {
        final QuasiSpeciesTree qsTree = quasiSpeciesTreeInput.get(this);
        // Randomly select event on tree:
        // weighted by the number of events (i.e. count of each haplotype)
// ...  change from here

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

        // as we choose index changeIdx between 1 - #haplo then changeIdx-1=0 at minimum

        tmax = tempqstimes[changeIdx-1];

        if (changeIdx+1>(tempqstimes.length-1))
            tmin = node.getHeight();
        else
            tmin = tempqstimes[changeIdx+1];

        double u = Randomizer.nextDouble();
        double tnew = u*tmin + (1-u)*tmax; // = u*(tmin-tmax)+tmax =u*(tmin-(tmin+tmax-tmin))+tmin+tmax-tmin
        // invert u and 1-u and have the same as u(tmax-tmin)+tmin
        // (1-u)*tmin + u*tmax-u*(tmax-tmin)-tmin=0 ???
        // indeed tmin-u*tmin+u*tmin-u*tmax+u*tmin-tmin=0

        tempqstimes[changeIdx]=tnew;
        qsTree.setAttachmentTimesList(node, tempqstimes);

        return 0.0;

    }






}
