package quasispeciestree.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import quasispeciestree.tree.QuasiSpeciesNode;

import java.util.ArrayList;

/**
 *  @author Veronika Boskova created on 23/04/17
 */
@Description("Subtree/branch exchange operations for quasispecies tree")
public class QuasiSpeciesSubtreeExchange extends QuasiSpeciesTreeOperator{

    final public Input<Boolean> isNarrowInput = new Input<>("isNarrow",
            "if true (default) a narrow exchange is performed, otherwise a wide exchange", true);

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        if (qsTree.getLeafNodeCount() < 4)
            throw new IllegalArgumentException("QuasiSpeciesSubtreeExchange operator cannot be " +
                    "used since there are only 4 tips in the tree! Remove this " +
                    "operator from your XML file.");
    }

    @Override
    public double proposal() {

        double logHastingsRatio = 0.0;

        // Select source and destination nodes:

        Node srcNode, srcNodeParent, destNode, destNodeParent;

        // Narrow exchange selection:

        if (isNarrowInput.get()) {
            do {
                srcNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
            } while (srcNode.isRoot() || srcNode.getParent().isRoot());
            srcNodeParent = srcNode.getParent();
            destNodeParent = srcNodeParent.getParent();
            destNode = getOtherChild(destNodeParent, srcNodeParent);
        }

        // Wide exchange selection:

        else {
            do {
                srcNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
            } while (srcNode.isRoot());
            srcNodeParent = srcNode.getParent();

            do {
                destNode = qsTree.getNode(Randomizer.nextInt(qsTree.getNodeCount()));
            } while(destNode == srcNode || destNode.isRoot() || destNode.getParent() == srcNodeParent
                    || srcNodeParent == destNode || destNode.getParent() == srcNode);
            destNodeParent = destNode.getParent();
        }


        // Reject if substitution would result in negative branch lengths:
        if (destNode.getHeight() >= srcNodeParent.getHeight() || srcNode.getHeight() >= destNodeParent.getHeight())
            return Double.NEGATIVE_INFINITY;

        // Record probability of forward moves:
        int srcHaplo = ((QuasiSpeciesNode) srcNode).getContinuingHaploName();
        int destHaplo = ((QuasiSpeciesNode) destNode).getContinuingHaploName();
        // if there is no haplotype passing the current chosen node,
        //  get all the possible haplotypes that can be moved up
        ArrayList<Integer> possibleSrcHaplo = new ArrayList<>();
        if (srcHaplo == -1){
            checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, possibleSrcHaplo);
            // choose randomly from the array of haplotypes one to move up
            srcHaplo = possibleSrcHaplo.get(Randomizer.nextInt(possibleSrcHaplo.size()));
        }
        else
            possibleSrcHaplo.add(srcHaplo);

        ArrayList<Integer> possibleDestHaplo = new ArrayList<>();
        if (destHaplo == -1){
            checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) destNode, possibleDestHaplo);
            // choose randomly from the array of haplotypes one to move up
            destHaplo = possibleDestHaplo.get(Randomizer.nextInt(possibleDestHaplo.size()));
        }
        else
            possibleDestHaplo.add(destHaplo);

        // Find current boundaries for haplotype start times
        QuasiSpeciesNode srcnodehaplo = (QuasiSpeciesNode) qsTree.getNode(srcHaplo);
        QuasiSpeciesNode destnodehaplo = (QuasiSpeciesNode) qsTree.getNode(destHaplo);

        double srcHaploStartMin = srcnodehaplo.getHeight();
        double srcHaploStartMax = getMaxPossibleHaploAttachTime((QuasiSpeciesNode) srcNode, srcHaplo);

        double destHaploStartMin = destnodehaplo.getHeight();
        double destHaploStartMax = getMaxPossibleHaploAttachTime((QuasiSpeciesNode) destNode, destHaplo);

        if (destHaploStartMax == srcHaploStartMax){
            QuasiSpeciesNode mrca = findLastCommonAncestor(srcnodehaplo, destnodehaplo, (QuasiSpeciesNode) qsTree.getRoot());
            destHaploStartMax = mrca.getHeight();
            srcHaploStartMax = mrca.getHeight();
        }

        // set all the nodes where src or dest haplo pass to -1 and set above haplo for a node where they arise to -1
        QuasiSpeciesNode srcNodeSister = (QuasiSpeciesNode) getOtherChild(srcNodeParent, srcNode);
        QuasiSpeciesNode correctTreeFromThisNodeSrc1 = getQuasiSpeciesNodeBelowHaploDetached((QuasiSpeciesNode) srcNode,
                srcHaplo, (QuasiSpeciesNode) srcNodeParent, srcNodeSister, srcNodeParent.isRoot());
        QuasiSpeciesNode destNodeSister = (QuasiSpeciesNode) getOtherChild(destNodeParent, destNode);
        QuasiSpeciesNode correctTreeFromThisNodeDest1 = getQuasiSpeciesNodeBelowHaploDetached((QuasiSpeciesNode) destNode,
                destHaplo, (QuasiSpeciesNode) destNodeParent, destNodeSister, destNodeParent.isRoot());

        // Make changes to tree topology:
        replace(srcNodeParent, srcNode, destNode);
        replace(destNodeParent, destNode, srcNode);

        // scale src haplo
        double logHastingsRatioContribution = scaleThisHaplo(srcnodehaplo, destHaploStartMax, srcHaploStartMax, srcHaploStartMin, 0);
        if (logHastingsRatioContribution == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;
        else logHastingsRatio += logHastingsRatioContribution;

        // also scale dest haplo
        logHastingsRatioContribution = scaleThisHaplo(destnodehaplo, srcHaploStartMax, destHaploStartMax, destHaploStartMin, 0);
        if (logHastingsRatioContribution == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;
        else logHastingsRatio += logHastingsRatioContribution;

        // set above haplo for a node where src and dest haplo arise to srcHaplo and destHaplo respectively
        QuasiSpeciesNode correctTreeFromThisNodeSrc2 = getQuasiSpeciesNodeBelowHaploAttached(destNodeParent.getHeight(),
                srcHaplo, srcnodehaplo.getAttachmentTimesList()[0], (QuasiSpeciesNode) destNodeParent, destNodeParent.isRoot());
        QuasiSpeciesNode correctTreeFromThisNodeDest2 = getQuasiSpeciesNodeBelowHaploAttached(srcNodeParent.getHeight(),
                destHaplo, destnodehaplo.getAttachmentTimesList()[0], (QuasiSpeciesNode) srcNodeParent, srcNodeParent.isRoot());

        // Recalculate continuingHaplo and HaploAbove arrays
        recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) qsTree.getRoot());

        // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
        qsTree.countAndSetPossibleStartBranches();

        // Incorporate probability of choosing current haplotype to move
        logHastingsRatio += Math.log(possibleSrcHaplo.size());
        logHastingsRatio += Math.log(possibleDestHaplo.size());

        // Incorporate probability of choosing current haplotype to move back
        if (((QuasiSpeciesNode) srcNode).getContinuingHaploName() != srcHaplo){
            // check how many haplotypes can be pulled up
            ArrayList<Integer> backPossibleSrcHaplo = new ArrayList<>();
            checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, backPossibleSrcHaplo);
            // account for this in the hastings ratio
            logHastingsRatio -= Math.log(backPossibleSrcHaplo.size());
        }

        if (((QuasiSpeciesNode) destNode).getContinuingHaploName() != destHaplo){
            // check how many haplotypes can be pulled up
            ArrayList<Integer> backPossibleDestHaplo = new ArrayList<>();
            checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) destNode, backPossibleDestHaplo);
            // account for this in the hastings ratio
            logHastingsRatio -= Math.log(backPossibleDestHaplo.size());
        }

        // RETURN log(HASTINGS RATIO)
        return logHastingsRatio;
    }
}
