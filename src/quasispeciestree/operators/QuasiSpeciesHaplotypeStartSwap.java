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
        final QuasiSpeciesTree qsTree = quasiSpeciesTreeInput.get(this);
        // Randomly select event on tree:
        // weighted by the number of events (i.e. count of each haplotype??)
        // for now do random uniform from just the haplotype count (disregarding the counts)
        int event = Randomizer.nextInt(qsTree.getLeafNodeCount());

        QuasiSpeciesNode node = (QuasiSpeciesNode) qsTree.getNode(event);

        // if we did not assign the node at this stage, throw exception
        if (node == null)
            throw new IllegalStateException("Event selection loop fell through!");

        // get the haplotype from which the current haplotype arises
        ArrayList<QuasiSpeciesNode> parentHaplo = qsTree.getParentHaplo();
        QuasiSpeciesNode haplotypesParentHaplo = parentHaplo.get(node.getNr());
        QuasiSpeciesNode haplotypesGrandParentHaplo = null;
        if (haplotypesParentHaplo != null) {
            haplotypesGrandParentHaplo = parentHaplo.get(haplotypesParentHaplo.getNr());
        }

        // propose a new start time
        // reposition the event (i.e. haplotype starting time)
        double tmin, tmax;
        tmin = node.getHeight();
        // find what the max time can be ... restricted by the grandparent haplotype start or internal node
        if (haplotypesParentHaplo == null){
            tmax = qsTree.originInput.get().getValue();
        }
        else {
            if (haplotypesGrandParentHaplo == null){
                tmax = qsTree.originInput.get().getValue();
            }
            else {
                // find the internal node just below the grandparent haplotype
                QuasiSpeciesNode intnode = null;
                for (Node nodeI : qsTree.getInternalNodes()){
                    if (((QuasiSpeciesNode) nodeI).getHaploName() == haplotypesGrandParentHaplo.getID()){
                        intnode = (QuasiSpeciesNode) nodeI;
                        break;
                    }
                }
                // the grandparent haplotype can be from the same subtree or from the other subtree
                // find the last common ancestor of the chosen haplotype and its grandparent
                QuasiSpeciesNode lastCommonAncestorNode = findLastCommonAncestor(node, haplotypesGrandParentHaplo, (QuasiSpeciesNode)intnode.getParent());
                // in either case (same or other subtree) the last common ancestor node is where the current haplo has tmax
                tmax = lastCommonAncestorNode.getHeight();
            }
        }
        // choose new time to attach
        double u = Randomizer.nextDouble();

        // weight the attachment height by the total haplo count wrt to all others
//        double weight = qsTree.getHaplotypeCounts(node)/qsTree.getTotalAttachmentCounts();
//        u = u*weight;
        double tnew = u*tmin + (1-u)*tmax;











        // restructure the attachment time array -- push all the haplotypes in the way below the common internal node (scaling)
        Double[] tempqstimes=qsTree.getAttachmentTimesList(node).clone();


//        qsTree.setAttachmentTimesList(QuasiSpeciesNode node, Double[] tempqstimes)



        /// check consistency of the attachment time array
        // check the consistency of the haplotype start times -- by findParentHaplo?
        // recalculate the countAttachment times array!

        // TODO need to check whether we have to change the start-from-haplotype-array
        // check start time... check the array








//
//        tempqstimes[changeIdx]=tnew;
//        qsTree.setAttachmentTimesList(node, tempqstimes);

        return 0.0; // proper hastings ratio!!!

    }






}
