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

    public Input<Boolean> m_useAmbiguities = new Input<Boolean>("useAmbiguities", "flag to indicate that sites containing ambiguous states should be handled instead of ignored (the default)", false);
    public Input<Boolean> m_useTipLikelihoods = new Input<Boolean>("useTipLikelihoods", "flag to indicate that partial likelihoods are provided at the tips", false);


    public static enum Scaling {none, always, _default};
    public Input<Scaling> scaling = new Input<QuasiSpeciesTreeLikelihood.Scaling>(
            "scaling", "type of scaling to use, one of " + Arrays.toString(Scaling.values()) +
            ". If not specified, the -beagle_scaling flag is used.", Scaling._default, Scaling.values());


    /**
     * calculation engine *
     */
//    protected LikelihoodCore likelihoodCore;
    protected QuasiSpeciesBeerLikelihoodCore likelihoodCore;
    BeagleTreeLikelihood beagle;

    /**
     * Plugin associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    QuasiSpeciesSubstitutionModel.Base substitutionModel;
//    SubstitutionModel.Base substitutionModel;
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

//    /**
//     * Lengths of the branches in the tree associated with each of the nodes
//     * in the tree through their node  numbers. By comparing whether the
//     * current branch length differs from stored branch lengths, it is tested
//     * whether a node is dirty and needs to be recomputed (there may be other
//     * reasons as well...).
//     * These lengths take branch rate models in account.
//     */
//    protected double[] m_branchLengths;
//    protected double[] storedBranchLengths;
    protected double[] m_qsTotalBranchLenghts;
    protected double[] storedQsTotalBranchLenghts;

    /**
     * memory allocation for likelihoods for each of the patterns *
     */
    protected double[] patternLogLikelihoods;
    /**
     * memory allocation for the root partials *
     */
    protected double[] m_fRootPartials;
    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
//    double[] probabilities;
    double[] qsprobabilities;

    int arraySize;

    /**
     * flag to indicate ascertainment correction should be applied *
     */
    boolean useAscertainedSitePatterns = false;

//    /**
//     * dealing with proportion of site being invariant *
//     */
//    double proportionInvariant = 0;
//    List<Integer> constantPattern = null;

    @Override
    public void initAndValidate() throws Exception {
        // sanity check: alignment should have same #taxa as tree
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new Exception("The number of nodes in the tree does not match the number of sequences");
        }
//        beagle = null;
//        beagle = new BeagleTreeLikelihood();
//        try {
//	        beagle.initByName(
//                    "data", dataInput.get(), "tree", treeInput.get(), "siteModel", siteModelInput.get(),
//                    "branchRateModel", branchRateModelInput.get(), "useAmbiguities", m_useAmbiguities.get(),
//                    "useTipLikelihoods", m_useTipLikelihoods.get(),"scaling", scaling.get().toString());
//	        if (beagle.beagle != null) {
//	            //a Beagle instance was found, so we use it
//	            return;
//	        }
//        } catch (Exception e) {
//			// ignore
//		}
        // No Beagle instance was found, so we use the good old java likelihood core
        beagle = null;

        int nodeCount = treeInput.get().getNodeCount();
        int leafNodeCount = treeInput.get().getLeafNodeCount();
        if (!(siteModelInput.get() instanceof QuasiSpeciesSiteModel.Base)) {
        	throw new Exception ("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (QuasiSpeciesSiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
//        substitutionModel = m_siteModel.substModelInput.get();
        substitutionModel = m_siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
//        m_branchLengths = new double[nodeCount];
//        storedBranchLengths = new double[nodeCount];

        m_qsTotalBranchLenghts = new double[leafNodeCount];
        storedQsTotalBranchLenghts = new double[leafNodeCount];

        int nStateCount = dataInput.get().getMaxStateCount();
        int nPatterns = dataInput.get().getPatternCount();
//        if (nStateCount == 4) {
//            likelihoodCore = new BeerLikelihoodCore4();
//        } else {
//            likelihoodCore = new BeerLikelihoodCore(nStateCount);
            likelihoodCore = new QuasiSpeciesBeerLikelihoodCore(nStateCount);
//        }

        String className = getClass().getSimpleName();

        Alignment alignment = dataInput.get();

        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));
        // print startup messages via Log.print*

//        proportionInvariant = m_siteModel.getProportionInvariant();
//        m_siteModel.setPropInvariantIsCategory(false);
//        if (proportionInvariant > 0) {
//            calcConstantPatternIndices(nPatterns, nStateCount);
//        }

        initCore();

        patternLogLikelihoods = new double[nPatterns];
        m_fRootPartials = new double[nPatterns * nStateCount];
//        matrixSize = (nStateCount + 1) * (nStateCount + 1);
        arraySize = (nStateCount + 1);
//        probabilities = new double[(nStateCount + 1) * (nStateCount + 1)];
        qsprobabilities = new double[(nStateCount + 1) * (nStateCount + 1)];
//        Arrays.fill(probabilities, 1.0);
        Arrays.fill(qsprobabilities,1.0);

        if (dataInput.get().isAscertained) {
            useAscertainedSitePatterns = true;
        }
    }


//    /**
//     * Determine indices of m_fRootProbabilities that need to be updates
//     * // due to sites being invariant. If none of the sites are invariant,
//     * // the 'site invariant' category does not contribute anything to the
//     * // root probability. If the site IS invariant for a certain character,
//     * // taking ambiguities in account, there is a contribution of 1 from
//     * // the 'site invariant' category.
//     */
//    void calcConstantPatternIndices(final int nPatterns, final int nStateCount) {
//        constantPattern = new ArrayList<Integer>();
//        for (int i = 0; i < nPatterns; i++) {
//            final int[] pattern = dataInput.get().getPattern(i);
//            final boolean[] bIsInvariant = new boolean[nStateCount];
//            Arrays.fill(bIsInvariant, true);
//            for (final int state : pattern) {
//                final boolean[] bStateSet = dataInput.get().getStateSet(state);
//                if (m_useAmbiguities.get() || !dataInput.get().getDataType().isAmbiguousState(state)) {
//                    for (int k = 0; k < nStateCount; k++) {
//                        bIsInvariant[k] &= bStateSet[k];
//                    }
//                }
//            }
//            for (int k = 0; k < nStateCount; k++) {
//                if (bIsInvariant[k]) {
//                    constantPattern.add(i * nStateCount + k);
//                }
//            }
//        }
//    }

    void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        likelihoodCore.initialize(
                nodeCount,
                dataInput.get().getPatternCount(),
                m_siteModel.getCategoryCount(),
                true, m_useAmbiguities.get()
        );

        final int extNodeCount = nodeCount / 2 + 1;
        final int intNodeCount = nodeCount / 2;

        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
            setPartials(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        } else {
            setStates(treeInput.get().getRoot(), dataInput.get().getPatternCount());
        }
        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i);
        }
        //TODO just for now
        likelihoodCore.createNodePartials(0);
    }

    /**
     * This method samples the sequences based on the tree and site model.
     */
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Can't sample a fixed alignment!");
    }

    /**
     * set leaf states in likelihood core *
     */
    void setStates(Node node, int patternCount) {
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
    void setPartials(Node node, int patternCount) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int nStates = data.getDataType().getStateCount();
            double[] partials = new double[patternCount * nStates];
            int k = 0;
            int iTaxon = getTaxonIndex(node.getID(), data);
            for (int iPattern = 0; iPattern < patternCount; iPattern++) {                
                double[] tipLikelihoods = data.getTipLikelihoods(iTaxon,iPattern);
                if (tipLikelihoods != null) {
                	for (int iState = 0; iState < nStates; iState++) {
                		partials[k++] = tipLikelihoods[iState];
                	}
                }
                else {
                	int nState = data.getPattern(iTaxon, iPattern);
	                boolean[] stateSet = data.getStateSet(nState);
	                for (int iState = 0; iState < nStates; iState++) {
	                	 partials[k++] = (stateSet[iState] ? 1.0 : 0.0);                
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
    public double calculateLogP() throws Exception {
//        if (beagle != null) {
//            logP = beagle.calculateLogP();
//            return logP;
//        }
        final TreeInterface tree = treeInput.get();

        if (traverse((QuasiSpeciesNode) tree.getRoot()) != Tree.IS_CLEAN)
            calcLogP();

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
            System.err.println("Turning on scaling to prevent numeric instability " + m_fScale);
            likelihoodCore.setUseScaling(m_fScale);
            likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
            traverse((QuasiSpeciesNode) tree.getRoot());
            calcLogP();
            return logP;
        }
        return logP;
    }

    void calcLogP() throws Exception {
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
//    int traverse(final Node node) throws Exception {
//
//        int update = (node.isDirty() | hasDirt);
//
//        final int iNode = node.getNr();
//
//        final double branchRate = branchRateModel.getRateForBranch(node);
//        final double branchTime = node.getLength() * branchRate;
//
//        // First update the transition probability matrix(ices) for this branch
//        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[iNode])) {
//        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[iNode])) {
//            m_branchLengths[iNode] = branchTime;
//            final Node parent = node.getParent();
//            likelihoodCore.setNodeMatrixForUpdate(iNode);
//            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
//                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
//                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
//                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
//                likelihoodCore.setNodeMatrix(iNode, i, probabilities);
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
//                likelihoodCore.setNodePartialsForUpdate(iNode);
//                update |= (update1 | update2);
//                if (update >= Tree.IS_FILTHY) {
//                    likelihoodCore.setNodeStatesForUpdate(iNode);
//                }
//
//                if (m_siteModel.integrateAcrossCategories()) {
//                    likelihoodCore.calculatePartials(childNum1, childNum2, iNode);
//                } else {
//                    throw new Exception("Error TreeLikelihood 201: Site categories not supported");
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


    /* Assumes there IS a branch rate model as opposed to traverse() */
    int traverse(final QuasiSpeciesNode node) throws Exception {

        QuasiSpeciesTree tree = (QuasiSpeciesTree) treeInput.get();

        int update = (node.isDirty() | hasDirt);

        final int iNode = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);

        final double totalBranchTime = tree.getTotalAttachmentTimes(node);

        final double branchTime =  totalBranchTime * branchRate;

        // First update the transition probability matrix(ices) for this branch
        if (update != Tree.IS_CLEAN || branchTime != m_qsTotalBranchLenghts[iNode]) {
            m_qsTotalBranchLenghts[iNode] = branchTime;
            likelihoodCore.setNodeMatrixForUpdate(iNode);
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getQSTransitionProbabilities(node, totalBranchTime, jointBranchRate, qsprobabilities);
                likelihoodCore.setNodeMatrix(iNode, i, qsprobabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // update the partial likelihoods.
        likelihoodCore.setNodePartialsForUpdate(iNode);
        if (update >= Tree.IS_FILTHY) {
            likelihoodCore.setNodeStatesForUpdate(iNode);
        }

//        if (m_siteModel.integrateAcrossCategories()) {
            likelihoodCore.calculatePartials(iNode);
//        } else {
//            throw new Exception("Error TreeLikelihood 201: Site categories not supported");
//        }

        // this node acts as the root of the beast.tree -
        // calculate the pattern likelihoods
//        final double[] frequencies = //m_pFreqs.get().
//        substitutionModel.getFrequencies();

        final double[] proportions = m_siteModel.getCategoryProportions(node);
        likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

        likelihoodCore.calculateQSLogLikelihoods(m_fRootPartials, patternLogLikelihoods);

        return update;
    } // traverseWithBRM

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
        if (beagle != null) {
            beagle.store();
            super.store();
            return;
        }
        if (likelihoodCore != null) {
            likelihoodCore.store();
        }
        super.store();
//        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
        System.arraycopy(m_qsTotalBranchLenghts, 0, storedQsTotalBranchLenghts, 0, m_qsTotalBranchLenghts.length);
    }

    @Override
    public void restore() {
        if (beagle != null) {
            beagle.restore();
            super.restore();
            return;
        }
        if (likelihoodCore != null) {
            likelihoodCore.restore();
        }
        super.restore();
//        double[] tmp = m_branchLengths;
//        m_branchLengths = storedBranchLengths;
//        storedBranchLengths = tmp;
        double[] tmp = m_qsTotalBranchLenghts;
        m_qsTotalBranchLenghts = storedQsTotalBranchLenghts;
        storedQsTotalBranchLenghts = tmp;
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    public List<String> getArguments() {
        return Collections.singletonList(dataInput.get().getID());
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    public List<String> getConditions() {
        return m_siteModel.getConditions();
    }

} // class TreeLikelihood
