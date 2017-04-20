
package quasispeciestree.likelihood;

import beast.evolution.likelihood.LikelihoodCore;
/**
 * The likelihood core is the class that performs the calculations
 * in the peeling algorithm (see Felsenstein, Joseph (1981).
 * Evolutionary trees from DNA sequences: a maximum likelihood approach.
 * J Mol Evol 17 (6): 368-376.). It does this by calculating the partial
 * results for all sites, node by node. The order in which the nodes
 * are visited is controlled by the TreeLikelihood. T
 * <p/>
 * In order to reuse computations of previous likelihood calculations,
 * a current state, and a stored state are maintained. Again, the TreeLikelihood
 * controls when to update from current to stored and vice versa. So, in
 * LikelihoodCore implementations, duplicates need to be kept for all partials.
 * Also, a set of indices to indicate which of the data is stored state and which
 * is current state is commonly the most efficient way to sort out which is which.
 */

abstract public class QuasiSpeciesLikelihoodCore extends LikelihoodCore{

    /**
     * integrate partials over categories at the origin (if any). *
     */
    abstract public void integratePartials(double[] inPartials, double[] proportions, double[] outPartials);

    /**
     * Calculate partials for node node3, with children node1 and node2Index.
     * NB Depending on whether the child nodes contain states or partials, the
     * calculation differs-*
     */
    abstract public void calculateQSPartials(int nodeIndex1, int nodeIndex2, int nodeIndex3, int child1QS, int child2QS, int child1parentQS, int nodeCount);

    /**
     * Calculates partial likelihoods at a node.
     *
     * @param rootNodeIndex     the root node
     * @param rootQS            QS passing through the root node
     * @param nodeCount         QS passing through child 2 - if none this is -1
     * @param originPartials   probability vector at origin (of length nrOfStates * nrOfPatterns)
     */
    abstract public void calculateOriginRootPartials(int rootNodeIndex, int rootQS, int nodeCount, double[] originPartials);


    /**
     * Calculates partial likelihoods at origin when coming from the root.
     *
     * @param nodeIndex1       the tip (child 1) node
     * @param child1QS         QS passing through the tip (child 1)
     * @param nodeCount        total count of the true nodes in the tree
     * @param originPartials   probability vector at origin (of length nrOfStates * nrOfPatterns)
     */
    abstract public void calculateOriginTipPartials(int nodeIndex1, int child1QS, int nodeCount, double[] originPartials);

}
