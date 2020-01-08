package piqmee.likelihood;

import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.util.Log;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.branchratemodel.UCRelaxedClockModel;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.evolution.likelihood.LikelihoodCore;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import piqmee.tree.QuasiSpeciesNode;
import piqmee.tree.QuasiSpeciesTree;

import java.util.*;

/**
 *  @author Veronika Boskova created on summer 2016 finished 21/04/2017
 */
@Description("Calculates the probability of sequence data on a beast.piqmee.tree " +
        "given a site and substitution model using a variant of the 'peeling algorithm'. " +
        "For details, see Felsenstein, Joseph (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. J Mol Evol 17 (6): 368-376.")
public class QuasiSpeciesTreeLikelihood extends GenericTreeLikelihood {

    final public Input<Boolean> useAmbiguities = new Input<>("useAmbiguities", "flag to indicate that sites containing ambiguous states should be handled instead of ignored (the default)", false);
    final public Input<Boolean> useTipLikelihoods = new Input<>("useTipLikelihoods", "flag to indicate that partial likelihoods are provided at the tips", false);
    public static enum Scaling {none, always, _default};
    final public Input<Scaling> scaling = new Input<>("scaling", "type of scaling to use, one of " + Arrays.toString(Scaling.values()) + ". If not specified, the -beagle_scaling flag is used.", Scaling._default, Scaling.values());

//    public Input<RealParameter> origin =
//            new Input<RealParameter>("origin", "The time from origin to last sample (must be larger than tree height)", (RealParameter) null, Input.Validate.REQUIRED);

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    protected LikelihoodCore likelihoodCore;
//    BeagleTreeLikelihood beagle;
    protected SubstitutionModel substitutionModel;
    protected SiteModel.Base siteModel;
    protected BranchRateModel.Base branchRateModel;
    protected double[] branchLengths;
    protected double[] storedBranchLengths;
    protected double[] patternLogLikelihoods;
    protected double[] rootPartials;
    protected double[] originPartials;
    int nodeCount;
    int leafNodeCount;
    int matrixSize;
    int nStates;
    Alignment alignment;

    /**
     * Memory for transition probabilities.
     */
    double[] probabilities;
    /**
     * Memory for substitution rates QS.
     */
    double[] rates;
    double[] storedRates;
    double[] tmpevectimesevals;
    /**
     * flag to indicate ascertainment correction should be applied *
     */
    boolean useAscertainedSitePatterns = false;
    /**
     * dealing with proportion of site being invariant *
     */
    double proportionInvariant = 0;
    List<Integer> constantPattern = null;

    Node toyNode = new Node();


    @Override
    public void initAndValidate(){
        // sanity check: alignment should have same #taxa as tree -- not the case for QS trees
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
//            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
            // subset the alignment to match the taxa in the tree
            alignment = subset(dataInput,treeInput.get().getExternalNodes());
        }
        else
            alignment = dataInput.get();
//        beagle = null;
//        beagle = new BeagleTreeLikelihood();
//        try {
//            beagle.initByName(
//                    "data", dataInput.get(), "tree", treeInput.get(), "siteModel", siteModelInput.get(),
//                    "branchRateModel", branchRateModelInput.get(), "useAmbiguities", useAmbiguities.get(),
//                    "useTipLikelihoods", useTipLikelihoods.get(),"scaling", scaling.get().toString());
//            if (beagle.beagle != null) {
//                //a Beagle instance was found, so we use it
//                return;
//            }
//        } catch (Exception e) {
//            // ignore
//        }
        // No Beagle instance was found, so we use the good old java likelihood core
//        beagle = null;

        nodeCount = treeInput.get().getNodeCount();
        leafNodeCount = treeInput.get().getLeafNodeCount();

        nStates = alignment.getMaxStateCount();

        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        siteModel = (SiteModel.Base) siteModelInput.get();
        siteModel.setDataType(alignment.getDataType());
        substitutionModel = siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
//            if (!(branchRateModel instanceof StrictClockModel))
//                throw new IllegalArgumentException("PIQMEE currently only" +
//                        "supports strict clock model.");
        } else {
            branchRateModel = new StrictClockModel();
        }
        // the entries corresponding to node number, store the branch length above the node
        // the entries corresponding to the node number + node count store the branch lengths above the QS origin, wherever this may be
        branchLengths = new double[nodeCount+leafNodeCount];
        storedBranchLengths = new double[nodeCount+leafNodeCount];

        int patterns = alignment.getPatternCount();
//        if (nStates == 4) {
//            likelihoodCore = new BeerLikelihoodCore4();
//        } else {
            likelihoodCore = new QuasiSpeciesBeerLikelihoodCore(nStates);
//        }

        String className = getClass().getSimpleName();

        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));
        // print startup messages via Log.print*

        proportionInvariant = siteModel.getProportionInvariant();
        siteModel.setPropInvariantIsCategory(false);
        if (proportionInvariant > 0) {
            calcConstantPatternIndices(patterns, nStates);
        }

        initCore();

        patternLogLikelihoods = new double[patterns];
        rootPartials = new double[patterns * nStates * siteModel.getCategoryCount()];
        originPartials = new double[patterns * nStates];
        matrixSize = (nStates + 1) * (nStates + 1);
        probabilities = new double[(nStates + 1) * (nStates + 1)];
        Arrays.fill(probabilities, 1.0);

        rates = new double[nStates];
        storedRates = new double [nStates];
        tmpevectimesevals = new double[nStates * nStates];
        getNoChangeRates(rates);

        if (alignment.isAscertained) {
            useAscertainedSitePatterns = true;
        }
    }

    /**
     * Get the transition rates for the QS branches,
     * i.e. the rate of no change, i.e. the diagonal entries in the rate matix
     *
     * @param rates the matrix where the new rates should be stored
     *
     */
    public void getNoChangeRates(double[] rates) {
        EigenDecomposition eigenDecomposition = substitutionModel.getEigenDecomposition(null);
        double[] evec = eigenDecomposition.getEigenVectors();
        double[] eval = eigenDecomposition.getEigenValues();
        double[] ievc = eigenDecomposition.getInverseEigenVectors();
        for (int i = 0; i < nStates; i++) {
            for (int j = 0; j < nStates; j++) {
                tmpevectimesevals[i * nStates + j] = evec[i * nStates + j] * eval[j];
            }
        }
        for (int i = 0; i < nStates; i++) {
            rates[i] = 0;
            for (int j = 0; j < nStates; j++) {
                rates[i] += tmpevectimesevals[i * nStates + j] * ievc[j * nStates + i];
            }
            //            System.out.println(rates[i]);
        }
    }

    /**
     * Determine indices of rootProbabilities that need to be updated
     * // due to sites being invariant. If none of the sites are invariant,
     * // the 'site invariant' category does not contribute anything to the
     * // root probability. If the site IS invariant for a certain character,
     * // taking ambiguities in account, there is a contribution of 1 from
     * // the 'site invariant' category.
     */
    void calcConstantPatternIndices(final int patterns, final int nStates) {
        constantPattern = new ArrayList<>();
        for (int i = 0; i < patterns; i++) {
            final int[] pattern = alignment.getPattern(i);
            final boolean[] isInvariant = new boolean[nStates];
            Arrays.fill(isInvariant, true);
            for (final int state : pattern) {
                final boolean[] isStateSet = alignment.getStateSet(state);
                if (useAmbiguities.get() || !alignment.getDataType().isAmbiguousState(state)) {
                    for (int k = 0; k < nStates; k++) {
                        isInvariant[k] &= isStateSet[k];
                    }
                }
            }
            for (int k = 0; k < nStates; k++) {
                if (isInvariant[k]) {
                    constantPattern.add(i * nStates + k);
                }
            }
        }
    }

    protected void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        likelihoodCore.initialize(
                nodeCount+leafNodeCount,
                alignment.getPatternCount(),
                siteModel.getCategoryCount(),
                true, useAmbiguities.get()
        );

        final int extNodeCount = nodeCount / 2 + 1;
        // the intNodeCount includes the true internal nodes and the QS start "nodes"
        final int intNodeCount = nodeCount / 2 + extNodeCount;

        if (useAmbiguities.get() || useTipLikelihoods.get()) {
            setPartials(treeInput.get().getRoot(), alignment.getPatternCount());
        } else {
            setStates(treeInput.get().getRoot(), alignment.getPatternCount());
        }
        hasDirt = QuasiSpeciesTree.IS_FILTHY;
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
            Alignment data = alignment;
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
            Alignment data = alignment;
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

        if (siteModel.isDirtyCalculation())
            getNoChangeRates(rates);

        try {
            if (traverse((QuasiSpeciesNode) tree.getRoot()) != QuasiSpeciesTree.IS_CLEAN)
                calcLogP();
        }
        catch (ArithmeticException e) {
            return Double.NEGATIVE_INFINITY;
        }
        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
            System.err.println("Switch off scaling");
            likelihoodCore.setUseScaling(1.0);
            likelihoodCore.unstore();
            hasDirt = QuasiSpeciesTree.IS_FILTHY;
            X *= 2;
            traverse((QuasiSpeciesNode) tree.getRoot());
            calcLogP();
            return logP;
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
            m_nScale = 0;
            m_fScale *= 1.01;
            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
            likelihoodCore.setUseScaling(m_fScale);
            likelihoodCore.unstore();
            hasDirt = QuasiSpeciesTree.IS_FILTHY;
            traverse((QuasiSpeciesNode) tree.getRoot());
            calcLogP();
            return logP;
        }
        return logP;
    }

    void calcLogP() {
        logP = 0.0;
        if (useAscertainedSitePatterns) {
            final double ascertainmentCorrection = alignment.getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < alignment.getPatternCount(); i++) {
                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * alignment.getPatternWeight(i);
            }
        } else {
            for (int i = 0; i < alignment.getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * alignment.getPatternWeight(i);
            }
        }
    }


    /** CalculationNode methods **/

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
//        if (beagle != null) {
//            return beagle.requiresRecalculation();
//        }
        hasDirt = QuasiSpeciesTree.IS_CLEAN;

        if (alignment.isDirtyCalculation()) {
            hasDirt = QuasiSpeciesTree.IS_FILTHY;
            return true;
        }
        if (siteModel.isDirtyCalculation()) {
            hasDirt = QuasiSpeciesTree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            //m_nHasDirt = QuasiSpeciesTree.IS_DIRTY;
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
        System.arraycopy(branchLengths, 0, storedBranchLengths, 0, branchLengths.length);
        System.arraycopy(rates, 0, storedRates, 0, rates.length);
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
        double[] tmp = branchLengths;
        branchLengths = storedBranchLengths;
        storedBranchLengths = tmp;
        tmp = rates;
        rates = storedRates;
        storedRates = tmp;
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    @Override
    public List<String> getArguments() {
        return Collections.singletonList(alignment.getID());
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    @Override
    public List<String> getConditions() {
        return siteModel.getConditions();
    }



    /**
     * QS OWN FUNCTIONS
     */

    /* Assumes there IS a branch rate model as opposed to traverse() */
    int traverse(final QuasiSpeciesNode node){

        QuasiSpeciesTree Tree = (QuasiSpeciesTree) treeInput.get();
//        Double originHeight = origin.get().getValue();

        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);


        // get the branch length, if the node is a tip, the total branch length above is the sum of the
        // branch lengths from the origin/attachment time to tip
        final double totalBranchTime;
        if (node.isLeaf())
            totalBranchTime = node.getTotalBranchLengths();
        else if (node.isRoot())
//            totalBranchTime = originHeight - node.getHeight();
            totalBranchTime = 0.0;
        else
            totalBranchTime = node.getLength();

        final double branchTime =  totalBranchTime * branchRate;

        // also check if the haplotype starts just above the node
        //  if YES, then have to split the branch into part that evolves normally and a part that does not evolve
        //  Note that we create a new variable since the QS can also start directly above the tip and we want the
        //    tip node to hold the probability of no change in the QS sequence for the whole duration of the QS
        //  For the relaxed clock model, this also means we need to create a new "branch rate category" for when haplo
        //    starts just above the tip by splitting the rate for the branch before and after the the haplo start
        double partBranchTime = 0;
        double partBranchRate = 0;
        if (node.getHaploAboveName() != -1) {
            if (node.isLeaf()){
                toyNode.setNr(nodeCount+node.getNr());
                partBranchRate = branchRateModel.getRateForBranch(toyNode);
                partBranchTime = (node.getLength() - (((QuasiSpeciesNode)
                        Tree.getNode(node.getHaploAboveName())).getAttachmentTimesList()[0] - node.getHeight())) * partBranchRate;
            }
            else {
                partBranchRate = branchRate;
                partBranchTime = (node.getLength() - (((QuasiSpeciesNode)
                        Tree.getNode(node.getHaploAboveName())).getAttachmentTimesList()[0] - node.getHeight())) * partBranchRate;
            }
        }
        // First update the transition probability matrix(ices) for this branch
        // Update the transition probability for the branches that do not evolve
        // if the node is at tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
        if (node.isLeaf() && (update != Tree.IS_CLEAN  || branchTime != branchLengths[nodeIndex])){
            branchLengths[nodeIndex] = branchTime;
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRate;
                // fill the transition probability matrix with move probabilities
                Arrays.fill(probabilities, 0);
                for (int j = 0; j < nStates; j++) {
                    probabilities[j * (nStates + 1)] = Math.exp(totalBranchTime * jointBranchRate * rates[j]);
                }
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }
        // Update the transition probability for the partial branches (internal node to QS start)
        if (node.getHaploAboveName() != -1 && (update != Tree.IS_CLEAN || partBranchTime != branchLengths[nodeCount + node.getHaploAboveName()])) {
            branchLengths[nodeCount + node.getHaploAboveName()] = partBranchTime;
            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(nodeCount + node.getHaploAboveName());
            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = siteModel.getRateForCategory(i, node) * partBranchRate;
                if (parent != null)
                    substitutionModel.getTransitionProbabilities(null, parent.getHeight(), ((QuasiSpeciesNode) Tree.getNode(node.getHaploAboveName())).getAttachmentTimesList()[0], jointBranchRate, probabilities);
                else
//                    substitutionModel.getTransitionProbabilities(null, originHeight, ((QuasiSpeciesNode) Tree.getNode(node.getHaploAboveName())).getAttachmentTimesList()[0], jointBranchRate, probabilities);
                    substitutionModel.getTransitionProbabilities(null, ((QuasiSpeciesNode) Tree.getNode(node.getHaploAboveName())).getAttachmentTimesList()[0], ((QuasiSpeciesNode) Tree.getNode(node.getHaploAboveName())).getAttachmentTimesList()[0], jointBranchRate, probabilities);
                likelihoodCore.setNodeMatrix(nodeCount + node.getHaploAboveName(), i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }
        //Update the transition probability matrix(ices) for all other branches
        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
        if (!node.isRoot() && ! node.isLeaf() && (update != Tree.IS_CLEAN || branchTime != branchLengths[nodeIndex])) {
            branchLengths[nodeIndex] = branchTime;
            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }
        //Update the transition probability matrix(ices) for root-origin branch
        else if (node.isRoot() && (update != Tree.IS_CLEAN || branchTime != branchLengths[nodeIndex])) {
            branchLengths[nodeIndex] = branchTime;
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRate;
//                substitutionModel.getTransitionProbabilities(node, originHeight, node.getHeight(), jointBranchRate, probabilities);
                substitutionModel.getTransitionProbabilities(node, node.getHeight(), node.getHeight(), jointBranchRate, probabilities);
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
                if (child1parentQS != child2parentQS)
                    throw new IllegalStateException("In QuasiSpeciesTreeLikelihood - QS of parent of child 1 ne to QS of parent of child 2");

                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY)
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);

                if (siteModel.integrateAcrossCategories()) {
                    //if (nStates == 4)
                    //    ((QuasiSpeciesBeerLikelihoodCore4) likelihoodCore).calculateQSPartials(childNum1, childNum2, nodeIndex, child1QS, child2QS, child1parentQS, nodeCount);
                    //else
                        ((QuasiSpeciesBeerLikelihoodCore) likelihoodCore).calculateQSPartials(childNum1, childNum2, nodeIndex, child1QS, child2QS, child1parentQS, nodeCount);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // This is the root of the beast.tree - there can be a root-origin branch

//                    // if NO root-origin branch do the old likelihood calculation
//                    if (Tree.originInput.get().getValue()==null){
//                        // integrate over all possible site categories the sites can be in
//                        final double[] proportions = siteModel.getCategoryProportions(node);
//                        likelihoodCore.integratePartials(node.getNr(), proportions, rootPartials);
//
//                        if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
//                            double proportionInvariant = siteModel.getProportionInvariant();
//                            // some portion of sites is invariant, so adjust root partials for this
//                            for (final int i : constantPattern) {
//                                rootPartials[i] += proportionInvariant;
//                            }
//                        }
//
//                        // calculate the pattern likelihoods
//                        // integrate over all possible starting state
//                        final double[] frequencies = //m_pFreqs.get().
//                                substitutionModel.getFrequencies();
//                        likelihoodCore.calculateLogLikelihoods(rootPartials, frequencies, patternLogLikelihoods);
//                    }
//                    // else continue the likelihood calculation up to the origin
//                    else {
                        // include the root-orig branch!!!

                        //if (nStates == 4)
                        //    ((QuasiSpeciesBeerLikelihoodCore4)likelihoodCore).calculateOriginRootPartials(nodeIndex, child1parentQS, nodeCount, rootPartials);
                        //else
                            ((QuasiSpeciesBeerLikelihoodCore)likelihoodCore).calculateOriginRootPartials(nodeIndex, child1parentQS, nodeCount, rootPartials);

                        // integrate over all possible site categories the sites can be in
                        final double[] proportions = siteModel.getCategoryProportions(node);
                        //if (nStates == 4)
                        //    ((QuasiSpeciesBeerLikelihoodCore4)likelihoodCore).integratePartials(rootPartials, proportions, originPartials);
                        //else
                            ((QuasiSpeciesBeerLikelihoodCore)likelihoodCore).integratePartials(rootPartials, proportions, originPartials);

                        if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
                            proportionInvariant = siteModel.getProportionInvariant();
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

                //if (nStates == 4)
                //    ((QuasiSpeciesBeerLikelihoodCore4)likelihoodCore).calculateOriginTipPartials(nodeIndex, child1parentQS, nodeCount, rootPartials);
                //else
                    ((QuasiSpeciesBeerLikelihoodCore)likelihoodCore).calculateOriginTipPartials(nodeIndex, child1parentQS, nodeCount, rootPartials);

                // integrate over all possible site categories the sites can be in
                final double[] proportions = siteModel.getCategoryProportions(node);
                //if (nStates == 4)
                //    ((QuasiSpeciesBeerLikelihoodCore4)likelihoodCore).integratePartials(rootPartials, proportions, originPartials);
                //else
                    ((QuasiSpeciesBeerLikelihoodCore)likelihoodCore).integratePartials(rootPartials, proportions, originPartials);

                if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
                    proportionInvariant = siteModel.getProportionInvariant();
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
    }

    /**
     * Function to subset full alignment to only sequences carried at tips in the
     * unique sequence tree
     *
     * @param data  full alignment
     * @param leafs leaf nodes of the unique-sequence tree
     * @return subsetted alignment
     */
    public Alignment subset(Input<Alignment> data, List<Node> leafs){
        Alignment fullData = data.get();
        int tipCount = leafs.size();
        ArrayList sequences = new ArrayList(tipCount);
        for (int i = 0; i < tipCount; i++){
            sequences.add(leafs.get(i).getNr(),fullData.sequenceInput.get().get(fullData.getTaxonIndex(leafs.get(i).getID())));
        }

        Alignment subsetData = new Alignment(sequences, fullData.dataTypeInput.get());
        return subsetData;
    }
}