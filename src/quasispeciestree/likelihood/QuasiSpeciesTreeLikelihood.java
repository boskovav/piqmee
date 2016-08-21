/*
* File TreeLikelihood.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/


package quasispeciestree.likelihood;

import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.*;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import quasispeciestree.sitemodel.QuasiSpeciesSiteModel;
import quasispeciestree.substitutionmodel.QuasiSpeciesSubstitutionModel;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;

import java.util.*;


@Description("Calculates the probability of sequence data on a beast.quasispeciestree.tree given a site and substitution model using " +
        "a variant of the 'peeling algorithm'. For details, see" +
        "Felsenstein, Joseph (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. J Mol Evol 17 (6): 368-376.")
public class QuasiSpeciesTreeLikelihood extends QuasiSpeciesGenericTreeLikelihood {

    final public Input<Boolean> m_useAmbiguities = new Input<>("useAmbiguities", "flag to indicate that sites containing ambiguous states should be handled instead of ignored (the default)", false);
    final public Input<Boolean> m_useTipLikelihoods = new Input<>("useTipLikelihoods", "flag to indicate that partial likelihoods are provided at the tips", false);


    public static enum Scaling {none, always, _default};
    final public Input<Scaling> scaling = new Input<>("scaling", "type of scaling to use, one of " + Arrays.toString(Scaling.values()) + ". If not specified, the -beagle_scaling flag is used.", Scaling._default, Scaling.values());


    /**
     * calculation engine *
     */
    protected QuasiSpeciesLikelihoodCore likelihoodCore;
    BeagleTreeLikelihood beagle;

    /**
     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    QuasiSpeciesSubstitutionModel substitutionModel;
    protected QuasiSpeciesSiteModel.Base m_siteModel;
    protected BranchRateModel.Base branchRateModel;

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;

    /**
     * memory allocation for likelihoods for each of the patterns *
     */
    protected double[] patternLogLikelihoods;
    /**
     * memory allocation for the root partials *
     */
    protected double[] m_fRootPartials;
    /**
     * memory allocation for the origin partials *
     */
    protected double[] originPartials;
    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    double[] probabilities;

    int matrixSize;

    /**
     * flag to indicate ascertainment correction should be applied *
     */
    boolean useAscertainedSitePatterns = false;

    /**
     * dealing with proportion of site being invariant *
     */
    double proportionInvariant = 0;
    List<Integer> constantPattern = null;

    int nodeCount;
    int leafNodeCount;

    @Override
    public void initAndValidate(){
        // sanity check: alignment should have same #taxa as tree
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }
//        beagle = null;
//        beagle = new BeagleTreeLikelihood();
//        try {
//            beagle.initByName(
//                    "data", dataInput.get(), "tree", treeInput.get(), "siteModel", siteModelInput.get(),
//                    "branchRateModel", branchRateModelInput.get(), "useAmbiguities", m_useAmbiguities.get(),
//                    "useTipLikelihoods", m_useTipLikelihoods.get(),"scaling", scaling.get().toString());
//            if (beagle.beagle != null) {
//                //a Beagle instance was found, so we use it
//                return;
//            }
//        } catch (Exception e) {
//            // ignore
//        }
        // No Beagle instance was found, so we use the good old java likelihood core
        beagle = null;

        nodeCount = treeInput.get().getNodeCount();
        leafNodeCount = treeInput.get().getLeafNodeCount();
        if (!(siteModelInput.get() instanceof QuasiSpeciesSiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (QuasiSpeciesSiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = m_siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        // the entries corresponding to node number, store the branch length above the node
        // the entries corresponding to the node number + node count store the branch lengths above the QS origin, wherever this may be
        m_branchLengths = new double[nodeCount+leafNodeCount];
        storedBranchLengths = new double[nodeCount+leafNodeCount];

        int stateCount = dataInput.get().getMaxStateCount();
        int patterns = dataInput.get().getPatternCount();
//        if (stateCount == 4) {
//            likelihoodCore = new BeerLikelihoodCore4();
//        } else {
            likelihoodCore = new QuasiSpeciesBeerLikelihoodCore(stateCount);
//        }

        String className = getClass().getSimpleName();

        Alignment alignment = dataInput.get();

        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));
        // print startup messages via Log.print*

        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        if (proportionInvariant > 0) {
            calcConstantPatternIndices(patterns, stateCount);
        }

        initCore();

        patternLogLikelihoods = new double[patterns];
        m_fRootPartials = new double[patterns * stateCount];
        originPartials = new double[patterns * stateCount];
        matrixSize = (stateCount + 1) * (stateCount + 1);
        probabilities = new double[(stateCount + 1) * (stateCount + 1)];
        Arrays.fill(probabilities, 1.0);

        if (dataInput.get().isAscertained) {
            useAscertainedSitePatterns = true;
        }
    }

    /**
     * Determine indices of m_fRootProbabilities that need to be updates
     * // due to sites being invariant. If none of the sites are invariant,
     * // the 'site invariant' category does not contribute anything to the
     * // root probability. If the site IS invariant for a certain character,
     * // taking ambiguities in account, there is a contribution of 1 from
     * // the 'site invariant' category.
     */
    void calcConstantPatternIndices(final int patterns, final int stateCount) {
        constantPattern = new ArrayList<>();
        for (int i = 0; i < patterns; i++) {
            final int[] pattern = dataInput.get().getPattern(i);
            final boolean[] isInvariant = new boolean[stateCount];
            Arrays.fill(isInvariant, true);
            for (final int state : pattern) {
                final boolean[] isStateSet = dataInput.get().getStateSet(state);
                if (m_useAmbiguities.get() || !dataInput.get().getDataType().isAmbiguousState(state)) {
                    for (int k = 0; k < stateCount; k++) {
                        isInvariant[k] &= isStateSet[k];
                    }
                }
            }
            for (int k = 0; k < stateCount; k++) {
                if (isInvariant[k]) {
                    constantPattern.add(i * stateCount + k);
                }
            }
        }
    }

    protected void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        likelihoodCore.initialize(
                nodeCount+leafNodeCount,
                dataInput.get().getPatternCount(),
                m_siteModel.getCategoryCount(),
                true, m_useAmbiguities.get()
        );

        final int extNodeCount = nodeCount / 2 + 1;
//        final int intNodeCount = nodeCount / 2;
        // the intNodeCount includes the true internal nodes and the QS start "nodes"
        final int intNodeCount = nodeCount / 2 + extNodeCount;

        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
            setPartials(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        } else {
            setStates(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        }
        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i);
        }
    }

    /**
     * This method samples the sequences based on the tree and site model.
     */
    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Can't sample a fixed alignment!");
    }

    /**
     * set leaf states in likelihood core *
     */
    protected void setStates(Node node, int patternCount) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int i;
            int[] states = new int[patternCount];
            int taxonIndex = getTaxonIndex(node.getID(), data);
            for (i = 0; i < patternCount; i++) {
                int code = data.getPattern(taxonIndex, i);
                int[] statesForCode = data.getDataType().getStatesForCode(code);
                if (statesForCode.length==1)
                    states[i] = statesForCode[0];
                else
                    states[i] = code; // Causes ambiguous states to be ignored.
            }
            likelihoodCore.setNodeStates(node.getNr(), states);

        } else {
            setStates(node.getLeft(), patternCount);
            setStates(node.getRight(), patternCount);
        }
    }

    /**
     *
     * @param taxon the taxon name as a string
     * @param data the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     *         or -1 if the taxon is not in the alignment.
     */
    private int getTaxonIndex(String taxon, Alignment data) {
        int taxonIndex = data.getTaxonIndex(taxon);
        if (taxonIndex == -1) {
            if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
            }
            if (taxonIndex == -1) {
                throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
            }
        }
        return taxonIndex;
    }

    /**
     * set leaf partials in likelihood core *
     */
    protected void setPartials(Node node, int patternCount) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int states = data.getDataType().getStateCount();
            double[] partials = new double[patternCount * states];
            int k = 0;
            int taxonIndex = getTaxonIndex(node.getID(), data);
            for (int patternIndex_ = 0; patternIndex_ < patternCount; patternIndex_++) {
                double[] tipLikelihoods = data.getTipLikelihoods(taxonIndex,patternIndex_);
                if (tipLikelihoods != null) {
                    for (int state = 0; state < states; state++) {
                        partials[k++] = tipLikelihoods[state];
                    }
                }
                else {
                    int stateCount = data.getPattern(taxonIndex, patternIndex_);
                    boolean[] stateSet = data.getStateSet(stateCount);
                    for (int state = 0; state < states; state++) {
                        partials[k++] = (stateSet[state] ? 1.0 : 0.0);
                    }
                }
            }
            likelihoodCore.setNodePartials(node.getNr(), partials);

        } else {
            setPartials(node.getLeft(), patternCount);
            setPartials(node.getRight(), patternCount);
        }
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    double m_fScale = 1.01;
    int m_nScale = 0;
    int X = 100;

    @Override
    public double calculateLogP() {
//        if (beagle != null) {
//            logP = beagle.calculateLogP();
//            return logP;
//        }
        final TreeInterface tree = treeInput.get();

        try {
            if (traverse((QuasiSpeciesNode) tree.getRoot()) != Tree.IS_CLEAN)
                calcLogP();
        }
        catch (ArithmeticException e) {
            return Double.NEGATIVE_INFINITY;
        }
        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
//            System.err.println("Switch off scaling");
//            m_likelihoodCore.setUseScaling(1.0);
//            m_likelihoodCore.unstore();
//            m_nHasDirt = Tree.IS_FILTHY;
//            X *= 2;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
            m_nScale = 0;
            m_fScale *= 1.01;
            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
            likelihoodCore.setUseScaling(m_fScale);
            likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
            traverse((QuasiSpeciesNode) tree.getRoot());
            calcLogP();
            return logP;
        }
        return logP;
    }

    void calcLogP() {
        logP = 0.0;
        if (useAscertainedSitePatterns) {
            final double ascertainmentCorrection = dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * dataInput.get().getPatternWeight(i);
            }
        } else {
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
            }
        }
    }

//    /* Assumes there IS a branch rate model as opposed to traverse() */
//    int traverse(final Node node) {
//
//        int update = (node.isDirty() | hasDirt);
//
//        final int nodeIndex = node.getNr();
//
//        final double branchRate = branchRateModel.getRateForBranch(node);
//        final double branchTime = node.getLength() * branchRate;
//
//        // First update the transition probability matrix(ices) for this branch
//        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
//        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
//            m_branchLengths[nodeIndex] = branchTime;
//            final Node parent = node.getParent();
//            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
//            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
//                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
//                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
//                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
//                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
//            }
//            update |= Tree.IS_DIRTY;
//        }
//
//        // If the node is internal, update the partial likelihoods.
//        if (!node.isLeaf()) {
//
//            // Traverse down the two child nodes
//            final Node child1 = node.getLeft(); //Two children
//            final int update1 = traverse(child1);
//
//            final Node child2 = node.getRight();
//            final int update2 = traverse(child2);
//
//            // If either child node was updated then update this node too
//            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {
//
//                final int childNum1 = child1.getNr();
//                final int childNum2 = child2.getNr();
//
//                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
//                update |= (update1 | update2);
//                if (update >= Tree.IS_FILTHY) {
//                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
//                }
//
//                if (m_siteModel.integrateAcrossCategories()) {
//                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
//                } else {
//                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
//                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
//                }
//
//                if (node.isRoot()) {
//                    // No parent this is the root of the beast.tree -
//                    // calculate the pattern likelihoods
//                    final double[] frequencies = //m_pFreqs.get().
//                            substitutionModel.getFrequencies();
//
//                    final double[] proportions = m_siteModel.getCategoryProportions(node);
//                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);
//
//                    if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
//                        proportionInvariant = m_siteModel.getProportionInvariant();
//                        // some portion of sites is invariant, so adjust root partials for this
//                        for (final int i : constantPattern) {
//                            m_fRootPartials[i] += proportionInvariant;
//                        }
//                    }
//
//                    likelihoodCore.calculateLogLikelihoods(m_fRootPartials, frequencies, patternLogLikelihoods);
//                }
//
//            }
//        }
//        return update;
//    } // traverseWithBRM


    /** CalculationNode methods **/

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
//        if (beagle != null) {
//            return beagle.requiresRecalculation();
//        }
        hasDirt = Tree.IS_CLEAN;

        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            //m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        return treeInput.get().somethingIsDirty();
    }

    @Override
    public void store() {
//        if (beagle != null) {
//            beagle.store();
//            super.store();
//            return;
//        }
        if (likelihoodCore != null) {
            likelihoodCore.store();
        }
        super.store();
        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
    }

    @Override
    public void restore() {
//        if (beagle != null) {
//            beagle.restore();
//            super.restore();
//            return;
//        }
        if (likelihoodCore != null) {
            likelihoodCore.restore();
        }
        super.restore();
        double[] tmp = m_branchLengths;
        m_branchLengths = storedBranchLengths;
        storedBranchLengths = tmp;
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    @Override
    public List<String> getArguments() {
        return Collections.singletonList(dataInput.get().getID());
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    @Override
    public List<String> getConditions() {
        return m_siteModel.getConditions();
    }



    /**
     * QS OWN FUNCTIONS
     */

    /* Assumes there IS a branch rate model as opposed to traverse() */
    int traverse(final QuasiSpeciesNode node){

        QuasiSpeciesTree Tree = (QuasiSpeciesTree) treeInput.get();
        Double originHeight = Tree.originInput.get().getValue();

        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);


        // hacks around unaccessible variables from TreeLikelihood
        QuasiSpeciesSubstitutionModel substitutionModel = m_siteModel.substModelInput.get();
        int stateCount = dataInput.get().getMaxStateCount();
        double[] probabilities = new double[(stateCount + 1) * (stateCount + 1)];

        // get the branch length, if the node is a tip, the total branch length above is the sum of the
        // branch lengths from the origin/attachment time to tip
        final double totalBranchTime;
        if (node.isLeaf())
            totalBranchTime = Tree.getTotalBranchLenghts(node);
        else if (node.isRoot())
            totalBranchTime = originHeight - node.getHeight();
//        else if (node.getContinuingHaploName() != -1)
//            totalBranchTime = 0;
        else
            totalBranchTime = node.getLength();

        final double branchTime =  totalBranchTime * branchRate;

        // also check if the haplotype starts just above the node
        //  if YES, then have to split the branch into part that evolves normally and a part that does not evolve
        //  Note that we create a new variable since the QS can also start directly above the tip and we want the
        //    tip node to hold the probability of no change in the QS sequence for the whole duration of the QS
        double partBranchTime = 0;
        if (node.getHaploAboveName() != -1)
            partBranchTime = ( node.getLength() - (Tree.getAttachmentTimesList(node.getHaploAboveName())[0]-node.getHeight()) ) * branchRate;

        // First update the transition probability matrix(ices) for this branch
        // Update the transition probability for the branches that do not evolve
        // if the node is at tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
        if (node.isLeaf() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])){
            m_branchLengths[nodeIndex] = branchTime;
//            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, branchTime, jointBranchRate, probabilities);
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }
        // Update the transition probability for the partial branches (internal node to QS start)
        if (node.getHaploAboveName() != -1 && (update != Tree.IS_CLEAN || partBranchTime != m_branchLengths[nodeCount + node.getHaploAboveName()])) {
            m_branchLengths[nodeCount + node.getHaploAboveName()] = partBranchTime;
            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(nodeCount + node.getHaploAboveName());
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                // TODO in case we use something else than JC69, check whether the node is used here for something!!! since we do not have full
                // TODO node parent length of the branch
                if (parent != null)
                    substitutionModel.getTransitionProbabilities(null, parent.getHeight(), Tree.getAttachmentTimesList(node.getHaploAboveName())[0], jointBranchRate, probabilities);
                else
                    substitutionModel.getTransitionProbabilities(null, originHeight, Tree.getAttachmentTimesList(node.getHaploAboveName())[0], jointBranchRate, probabilities);
                likelihoodCore.setNodeMatrix(nodeCount + node.getHaploAboveName(), i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }
        //Update the transition probability matrix(ices) for all other branches
        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
        else if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
            m_branchLengths[nodeIndex] = branchTime;
            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }


        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse((QuasiSpeciesNode) child1);

            final Node child2 = node.getRight();
            final int update2 = traverse((QuasiSpeciesNode) child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int child1QS = ((QuasiSpeciesNode) child1).getContinuingHaploName();
                final int child1parentQS = ((QuasiSpeciesNode) child1.getParent()).getContinuingHaploName();
                final int childNum2 = child2.getNr();
                final int child2QS = ((QuasiSpeciesNode) child2).getContinuingHaploName();
                final int child2parentQS = ((QuasiSpeciesNode) child2.getParent()).getContinuingHaploName();
                if (child1parentQS != child2parentQS){
                    System.out.println("In QuasiSpeciesTreeLikelihood - QS of parent of child 1 ne to QS of parent of child 2");
                    System.exit(0);
                }

                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
//                if (node.isLeaf())
//                    likelihoodCore.setNodePartialsForUpdate(nodeCount + nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
//                    if (node.isLeaf())
//                        likelihoodCore.setNodeStatesForUpdate(nodeCount + nodeIndex);
                }

                if (m_siteModel.integrateAcrossCategories()) {
                    likelihoodCore.calculateQSPartials(childNum1, childNum2, nodeIndex, child1QS, child2QS, child1parentQS, nodeCount);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // This is the root of the beast.tree - there can be a root-origin branch

//                    // if NO root-origin branch do the old likelihood calculation
//                    if (Tree.originInput.get().getValue()==null){
//                        // integrate over all possible site categories the sites can be in
//                        final double[] proportions = m_siteModel.getCategoryProportions(node);
//                        likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);
//
//                        if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
//                            double proportionInvariant = m_siteModel.getProportionInvariant();
//                            // some portion of sites is invariant, so adjust root partials for this
//                            for (final int i : constantPattern) {
//                                m_fRootPartials[i] += proportionInvariant;
//                            }
//                        }
//
//                        // calculate the pattern likelihoods
//                        // integrate over all possible starting state
//                        final double[] frequencies = //m_pFreqs.get().
//                                substitutionModel.getFrequencies();
//                        likelihoodCore.calculateLogLikelihoods(m_fRootPartials, frequencies, patternLogLikelihoods);
//                    }
//                    // else continue the likelihood calculation up to the origin
//                    else {
                        // include the root-orig branch!!!

                        likelihoodCore.calculateOriginRootPartials(nodeIndex, child1parentQS, nodeCount, originPartials);

                        // integrate over all possible site categories the sites can be in
                        final double[] proportions = m_siteModel.getCategoryProportions(node);
                        likelihoodCore.integratePartials(node.getNr(), proportions, originPartials);

                        if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
                            double proportionInvariant = m_siteModel.getProportionInvariant();
                            // some portion of sites is invariant, so adjust root partials for this
                            for (final int i : constantPattern) {
                                originPartials[i] += proportionInvariant;
                            }
                        }

                        // calculate the pattern likelihoods
                        // integrate over all possible starting state
                        final double[] frequencies = //m_pFreqs.get().
                                substitutionModel.getFrequencies();
                        likelihoodCore.calculateLogLikelihoods(originPartials, frequencies, patternLogLikelihoods);
//                    }
                }

            }
        }
        // if the tree has only one child
        else {
            if (Tree.getLeafNodeCount()==1){

                final int child1parentQS = ((QuasiSpeciesNode) node).getContinuingHaploName();

                likelihoodCore.calculateOriginTipPartials(nodeIndex, child1parentQS, nodeCount, originPartials);

                // integrate over all possible site categories the sites can be in
                final double[] proportions = m_siteModel.getCategoryProportions(node);
                likelihoodCore.integratePartials(node.getNr(), proportions, originPartials);

                if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
                    double proportionInvariant = m_siteModel.getProportionInvariant();
                    // some portion of sites is invariant, so adjust root partials for this
                    for (final int i : constantPattern) {
                        originPartials[i] += proportionInvariant;
                    }
                }

                // calculate the pattern likelihoods
                // integrate over all possible starting state
                final double[] frequencies = //m_pFreqs.get().
                        substitutionModel.getFrequencies();
                likelihoodCore.calculateLogLikelihoods(originPartials, frequencies, patternLogLikelihoods);
            }
        }
        return update;
    } // traverseWithBRM

} // class TreeLikelihood
