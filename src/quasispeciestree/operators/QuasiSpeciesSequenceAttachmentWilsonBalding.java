package quasispeciestree.operators;

import beast.core.Description;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;
import beast.util.Randomizer;

/**
 *  @author Veronika Boskova created on 07/08/2015
 */

@Description("Within given haplotype, randomly selects one sequence "
        + "attachment time, selects a new interval, restricted by the "
        + "clone start and clone tip times, where this attachment "
        + "time will be added to and attaches it uniformly in that interval.")
public class QuasiSpeciesSequenceAttachmentWilsonBalding extends QuasiSpeciesTreeOperator{

    /**
     * Change the attachment time and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {
        final QuasiSpeciesTree qsTree = quasiSpeciesTreeInput.get(this);
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
        int tminIdx, tmaxIdx;
        Double[] tempqstimes=qsTree.getAttachmentTimesList(node).clone();

        // choose new max index between 0 - #sequences of this haplotype
        tmaxIdx = Randomizer.nextInt(tempqstimes.length);
        tmax = tempqstimes[tmaxIdx];
        tminIdx = tmaxIdx + 1;
        // here we do not allow the event to be repositioned to the interval delimited by itself
        if (tmaxIdx == changeIdx-1 || tmaxIdx == changeIdx){
            return Double.NEGATIVE_INFINITY;
        }

        if (tminIdx == tempqstimes.length)
            tmin = node.getHeight();
        else
            tmin = tempqstimes[tminIdx];

        double newRange = tmax - tmin;
        double oldRange;
        if (changeIdx+1 == tempqstimes.length)
            oldRange = tempqstimes[changeIdx-1] - node.getHeight();
        else
            oldRange = tempqstimes[changeIdx-1] - tempqstimes[changeIdx+1];

        double u = Randomizer.nextDouble();
        double tnew = u*tmin + (1-u)*tmax; // = u*(tmin-tmax)+tmax =u*(tmin-(tmin+tmax-tmin))+tmin+tmax-tmin
        // invert u and 1-u and have the same as u(tmax-tmin)+tmin
        // (1-u)*tmin + u*tmax-u*(tmax-tmin)-tmin=0 ???
        // indeed tmin-u*tmin+u*tmax-u*tmax+u*tmin-tmin=0

        tempqstimes[changeIdx]=tnew;
        if (changeIdx > tmaxIdx){
            for (int i=changeIdx; i>tmaxIdx+1; i--){
                double swaptime = tempqstimes[i];
                tempqstimes[i] = tempqstimes[i-1];
                tempqstimes[i-1] = swaptime;
            }
        } else {
            for (int i=changeIdx; i<tmaxIdx; i++){
                double swaptime = tempqstimes[i];
                tempqstimes[i] = tempqstimes[i+1];
                tempqstimes[i+1] = swaptime;
            }
        }
        qsTree.setAttachmentTimesList(node, tempqstimes);

        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);

        return Math.log(newRange/oldRange);

    }

}
