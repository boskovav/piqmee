package piqmee.operators;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import piqmee.tree.QuasiSpeciesNode;

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
        if (qsTree.getLeafNodeCount() < 3)
            throw new IllegalArgumentException("QuasiSpeciesSubtreeExchange operator cannot be " +
                    "used since there are only 2 tips in the tree! Remove this " +
                    "operator from your XML file.");

        if (qsTree.getLeafNodeCount() < 4 && isNarrowInput.get() == false)
            throw new IllegalArgumentException("QuasiSpeciesSubtreeExchange wide version operator cannot be " +
                    "used since there are only 3 tips in the tree! Remove this " +
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
        double srcHaploStartMinBack = srcnodehaplo.getHeight();
        double srcHaploStartMax = getMaxPossibleHaploAttachTime((QuasiSpeciesNode) srcNode, srcHaplo);

        double destHaploStartMin = destnodehaplo.getHeight();
        double destHaploStartMinBack = destnodehaplo.getHeight();
        double destHaploStartMax = getMaxPossibleHaploAttachTime((QuasiSpeciesNode) destNode, destHaplo);

        // check if the two haplotypes from the exchanged subtrees do not coincide
        QuasiSpeciesNode mrca = findLastCommonAncestor(srcnodehaplo, destnodehaplo, (QuasiSpeciesNode) qsTree.getRoot());
        double tbottomsrcforward = 0;
        double tbottomsrcback = 0;
        double tbottomdestforward = 0;
        double tbottomdestback = 0;
        boolean doscale = true;

        // if same max, this means they can only be scaled up the their MRCA
        if (destHaploStartMax == srcHaploStartMax && destHaploStartMax >= mrca.getHeight()){
            destHaploStartMax = mrca.getHeight();
            srcHaploStartMax = mrca.getHeight();
        }
        // else it is possible that one haplo is the parent of the second
        else if (mrca.getContinuingHaploName() != -1){
            if (mrca.getContinuingHaploName() == srcHaplo && destHaploStartMax == mrca.getHeight()){
                if (destnodehaplo.getParentHaplo() != srcHaplo)
                    throw new IllegalStateException("QuasiSpeciesSubtreeExchange: We found the src haplo is a parent " +
                            "of the dest haplo, but this is not annotated as desthaplonode.getParentHaplo().");
                if (destnodehaplo.getAttachmentTimesList().length > 1){
                    destHaploStartMin = mrca.getHeight();
                    tbottomdestforward = destnodehaplo.getAttachmentTimesList()[1];
                    srcHaploStartMinBack = mrca.getHeight();
                    // tbottomsrcback means that the minimum scaling will be srcHaploStartMinBack/attachtime[1] where
                    // attachtime[1] will be taken from the newly scaled attachment time array
                    tbottomsrcback = -1;
                }
                // else we have a src haplotype above mrca and dest haplo below but without duplicate
                // so keep everything as it is - do not scale haplo, only change the tree
                else
                    doscale = false;
            }
            else if (mrca.getContinuingHaploName() == destHaplo && srcHaploStartMax == mrca.getHeight()){
                if (srcnodehaplo.getParentHaplo() != destHaplo)
                    throw new IllegalStateException("QuasiSpeciesSubtreeExchange: We found the dest haplo is a parent " +
                            "of the src haplo, but this is not annotated as srchaplonode.getParentHaplo().");
                if (srcnodehaplo.getAttachmentTimesList().length > 1){
                    srcHaploStartMin = mrca.getHeight();
                    tbottomsrcforward = srcnodehaplo.getAttachmentTimesList()[1];
                    destHaploStartMinBack = mrca.getHeight();
                    // tbottomdestback means that the minimum scaling will be destHaploStartMinBack/attachtime[1] where
                    // attachtime[1] will be taken from the newly scaled attachment time array
                    tbottomdestback = -1;
                }
                // else we have a dest haplotype above mrca and src haplo below but without duplicate
                // so keep everything as it is - do not scale haplo, only change the tree
                else
                    doscale = false;
            }
        }

        if (doscale) {
            // set all the nodes where src or dest haplo pass to -1 and set above haplo for a node where they arise to -1
            if (srcnodehaplo.getAttachmentTimesList().length > 1) {
                QuasiSpeciesNode srcNodeSister = (QuasiSpeciesNode) getOtherChild(srcNodeParent, srcNode);
                QuasiSpeciesNode correctTreeFromThisNodeSrc1 = getQuasiSpeciesNodeBelowHaploDetached((QuasiSpeciesNode) srcNode,
                    srcHaplo, (QuasiSpeciesNode) srcNodeParent, srcNodeSister, srcNodeParent.isRoot());
            }
            if (destnodehaplo.getAttachmentTimesList().length > 1) {
                QuasiSpeciesNode destNodeSister = (QuasiSpeciesNode) getOtherChild(destNodeParent, destNode);
                QuasiSpeciesNode correctTreeFromThisNodeDest1 = getQuasiSpeciesNodeBelowHaploDetached((QuasiSpeciesNode) destNode,
                    destHaplo, (QuasiSpeciesNode) destNodeParent, destNodeSister, destNodeParent.isRoot());
            }
        }

        // Make changes to tree topology:
        replace(srcNodeParent, srcNode, destNode);
        replace(destNodeParent, destNode, srcNode);

        if (doscale) {
            if (srcnodehaplo.getAttachmentTimesList().length > 1) {
                // scale src haplo
                double logHastingsRatioContribution = scaleThisHaplo(srcnodehaplo, destHaploStartMax, srcHaploStartMin, tbottomsrcforward, srcHaploStartMax, srcHaploStartMinBack, tbottomsrcback);
                if (logHastingsRatioContribution == Double.NEGATIVE_INFINITY)
                    return Double.NEGATIVE_INFINITY;
                else
                    logHastingsRatio += logHastingsRatioContribution;

                // set above haplo for a node where src and dest haplo arise to srcHaplo and destHaplo respectively
                QuasiSpeciesNode correctTreeFromThisNodeSrc2 = getQuasiSpeciesNodeBelowHaploAttached(destNodeParent.getHeight(),
                        srcHaplo, srcnodehaplo.getAttachmentTimesList()[0], (QuasiSpeciesNode) destNodeParent, destNodeParent.isRoot());
            }
            if (destnodehaplo.getAttachmentTimesList().length > 1) {
                // also scale dest haplo
                double logHastingsRatioContribution = scaleThisHaplo(destnodehaplo, srcHaploStartMax, destHaploStartMin, tbottomdestforward, destHaploStartMax, destHaploStartMinBack, tbottomdestback);
                if (logHastingsRatioContribution == Double.NEGATIVE_INFINITY)
                    return Double.NEGATIVE_INFINITY;
                else
                    logHastingsRatio += logHastingsRatioContribution;

                // set above haplo for a node where src and dest haplo arise to srcHaplo and destHaplo respectively
                QuasiSpeciesNode correctTreeFromThisNodeDest2 = getQuasiSpeciesNodeBelowHaploAttached(srcNodeParent.getHeight(),
                        destHaplo, destnodehaplo.getAttachmentTimesList()[0], (QuasiSpeciesNode) srcNodeParent, srcNodeParent.isRoot());
            }
        }

        // Recalculate continuingHaplo and HaploAbove arrays
        recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) qsTree.getRoot());

        // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
        qsTree.countAndSetPossibleStartBranches(srcnodehaplo);
        qsTree.countAndSetPossibleStartBranches(destnodehaplo);

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
