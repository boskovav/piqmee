package piqmee.likelihood;

import beast.core.Description;
import beast.evolution.likelihood.BeerLikelihoodCore;

/**
 *  @author Veronika Boskova created on 03/03/17
 */
@Description("quasispecies likelihood core")
public class QuasiSpeciesBeerLikelihoodCore extends BeerLikelihoodCore {

    public QuasiSpeciesBeerLikelihoodCore(int nrOfStates) {
        super (nrOfStates);
    } // c'tor

    /**
     * QS OWN FUNCTION
     */

    /**
     * Integrates partial likelihoods at origin.
     *
     * @param inPartials        partials from the origin
     * @param proportions
     * @param outPartials       partials to be outputted
     */
    public void integratePartials(double[] inPartials, double[] proportions, double[] outPartials) {
        calculateIntegratePartials(inPartials, proportions, outPartials);
    }

    /**
     * Calculates partial likelihoods at origin when coming from the root.
     *
     * @param nodeIndex1       the tip (child 1) node
     * @param child1QS         QS passing through the tip (child 1)
     * @param nodeCount        total count of the true nodes in the tree
     * @param originPartials   probability vector at origin (of length nrOfStates * nrOfPatterns)
     */
    public void calculateOriginTipPartials(int nodeIndex1, int child1QS, int nodeCount, double[] originPartials) {
        // tip has a state
        if (states[nodeIndex1] != null){
            calculateOriginTipPruning(
                    states[nodeIndex1],
                    matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                    matrices[currentMatrixIndex[nodeCount+nodeIndex1]][nodeCount+nodeIndex1],
                    originPartials,child1QS);
        }
        // tip has partials
        else {
            calculateOriginTipPruning(
                    partials[currentPartialsIndex[nodeIndex1]][nodeIndex1],
                    matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],
                    matrices[currentMatrixIndex[nodeCount+nodeIndex1]][nodeCount+nodeIndex1],
                    originPartials, child1QS);
        }
    }

    /**
     * Calculates partial likelihood at origin if tree has only one tip with states.
     *
     * @param stateIndex            alignment at the tip (child 1)
     * @param matricesQS1           transition probability matrix from parent to the tip (child 1) - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through the tip (child 1)
     * @param originPartials        probability vector at origin (of length nrOfStates * nrOfPatterns)
     * @param child1QS              QS passing through parent tip (child 1)
     */
    protected void calculateOriginTipPruning(int[] stateIndex, double[] matricesQS1, double[] matrices1aboveQSstart,
                                             double[] originPartials, int child1QS){

        double tmp, sum;

        // v keeps track of the pattern we are about to calculate
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {

            // w keeps track of the state the internal node evolves from
            int w = l * matrixSize;

            for (int k = 0; k < nrOfPatterns; k++) {
                // note down the state at the tip
                int state = stateIndex[k];

                // the child has a state
                if (state < nrOfStates) {
                    // note down the transition probabilities (of no change) on the sum of the QS branch lengths
                    // P(QS start -> QS tip)
                    tmp = matricesQS1[w + nrOfStates * state + state];
                    // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
                    for (int i = 0; i < nrOfStates; i++) {
                        originPartials[v + i] = tmp * matrices1aboveQSstart[w + nrOfStates * i + state];
                    }

                } else {
                    // the alignment at node has a gap or unknown state so treat it as unknown
                    // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
                    for (int i = 0; i < nrOfStates; i++) {
                        sum = 0.0;
                        for (int j = 0; j < nrOfStates; j++){
                            tmp = matricesQS1[w + nrOfStates * j + j];
                            sum += tmp * matrices1aboveQSstart[w + nrOfStates * i + j];
                        }
                        originPartials[v + i] = sum;
                    }
                }

                v += nrOfStates;

            }
        }
    }

    /**
     * Calculates partial likelihood at origin if tree has only one tip with partials.
     *
     * @param partials1             partials at the tip (child 1)
     * @param matricesQS1           transition probability matrix from parent to the tip (child 1) - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through the tip (child 1)
     * @param originPartials        probability vector at origin (of length nrOfStates * nrOfPatterns)
     * @param child1QS              QS passing through parent tip (child 1)
     */
    protected void calculateOriginTipPruning(double[] partials1, double[] matricesQS1, double[] matrices1aboveQSstart,
                                             double[] originPartials, int child1QS){

        double tmp, sum;

        // v keeps track of the pattern we are about to calculate
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {

            // w keeps track of the state the internal node evolves from
            int w = l * matrixSize;

            for (int k = 0; k < nrOfPatterns; k++) {
                // note down the state at the tip

                for (int i = 0; i < nrOfStates; i++) {
                    // since state at QS tip is unknown, take into account all the possibilities
                    // for each state at the origin calculate prob. of going to the state of the tip * P(QS start -> QS tip)
                    sum = 0.0;
                    for (int j = 0; j < nrOfStates; j++) {
                        // note down the partial at the child 1
                        tmp = partials1[v + j] * matricesQS1[w + nrOfStates * j + j];
                        // calculate the partials at the origin
                        sum += tmp * matrices1aboveQSstart[w + nrOfStates * i + j];
                    }
                    originPartials[v + i] = sum;
                }

                v += nrOfStates;

            }
        }
    }


    /**
     * Calculates partial likelihoods at origin when coming from the root.
     *
     * @param rootNodeIndex     the root node
     * @param rootQS            QS passing through the root node
     * @param nodeCount         total count of the true nodes in the tree
     * @param originPartials    probability vector at origin (of length nrOfStates * nrOfPatterns)
     */
    public void calculateOriginRootPartials(int rootNodeIndex, int rootQS, int nodeCount, double[] originPartials) {
        if (rootQS == -1) {
            calculateOriginRootPruning(null, partials[currentPartialsIndex[rootNodeIndex]][rootNodeIndex],
                    matrices[currentMatrixIndex[rootNodeIndex]][rootNodeIndex],
                    null,
                    originPartials, rootQS);
        } else {
            calculateOriginRootPruning(states[rootQS], partials[currentPartialsIndex[rootNodeIndex]][rootNodeIndex],
                    matrices[currentMatrixIndex[rootNodeIndex]][rootNodeIndex],
                    matrices[currentMatrixIndex[nodeCount + rootQS]][nodeCount + rootQS],
                    originPartials, rootQS);
        }
    }

    /**
     * Calculates partial likelihood at origin.
     *
     * @param stateIndexRoot    alignment at node corresponding to the QS passing through the root
     * @param partialsRoot      probability vector at root node (of length nrOfStates * nrOfPatterns)
     * @param matricesRoot      transition probability matrix from origin to root
     * @param originPartials    probability vector at origin (of length nrOfStates * nrOfPatterns)
     * @param rootQS            QS passing through parent node
     */
    protected void calculateOriginRootPruning(int[] stateIndexRoot, double[] partialsRoot, double[] matricesRoot, double[] matricesRootaboveQSstart,
                                              double[] originPartials, int rootQS) {

        // v keeps track of the pattern we are about to calculate
        int v = 0;

        // there is no QS arising above the root
        if (rootQS==-1 || stateIndexRoot == null){

            for (int l = 0; l < nrOfMatrices; l++) {

                // w keeps track of the state the internal node evolves from
                int w = l * matrixSize;

                for (int k = 0; k < nrOfPatterns; k++) {

                    calculateOriginRootPruningHelperNoQSorUnknownState(partialsRoot, matricesRoot, originPartials, w, v);

                    v += nrOfStates;

                }
            }
        }
        // there is QS passing through the root node
        else {

            for (int l = 0; l < nrOfMatrices; l++) {

                // w keeps track of the state the internal node evolves from
                int w = l * matrixSize;

                for (int k = 0; k < nrOfPatterns; k++) {
                    // note down the states at the tips belonging to QS passing through root
                    int rootState = stateIndexRoot[k];

                    // root's QS has a state
                    if (rootState < nrOfStates){
                        calculateOriginRootPruningHelperWithQS(partialsRoot, matricesRootaboveQSstart, originPartials, w, v, rootState);
                    // root's QS has a gap or unknown state
                    } else {
                        calculateOriginRootPruningHelperNoQSorUnknownState(partialsRoot, matricesRootaboveQSstart, originPartials, w, v);
                    }

                    v += nrOfStates;

                }
            }
        }
    }


    /**
     * Helper function to calculateOriginRootPruning for case of no QS passing or unknown state/gap at QS passing the root.
     *
     * @param partialsRoot                      probability vector at root node (of length nrOfStates * nrOfPatterns)
     * @param matricesRootOrRootaboveQSstart    transition probability matrix from origin to root
     * @param originPartials                    probability vector at origin (of length nrOfStates * nrOfPatterns)
     */
    protected void calculateOriginRootPruningHelperNoQSorUnknownState(double[] partialsRoot, double[] matricesRootOrRootaboveQSstart,
                                                               double[] originPartials, int w, int v) {

        double tmp, sum;

        for (int i = 0; i < nrOfStates; i++) {
            // since state at root is unknown, take into account all the possibilities
            // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
            sum = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the root
                tmp = partialsRoot[v + j];
                sum += tmp * matricesRootOrRootaboveQSstart[w + nrOfStates * i + j];
            }
            originPartials[v + i] = sum;
        }
    }

    /**
     * Helper function to calculateOriginRootPruning for case a QS passing.
     *
     * @param rootState                     states at the tip belonging to QS passing through root
     * @param partialsRoot                  probability vector at root node (of length nrOfStates * nrOfPatterns)
     * @param matricesRootaboveQSstart      transition probability matrix from origin to root
     * @param originPartials                probability vector at origin (of length nrOfStates * nrOfPatterns)
     */
    protected void calculateOriginRootPruningHelperWithQS(double[] partialsRoot, double[] matricesRootaboveQSstart,
                                                          double[] originPartials, int w, int v, int rootState) {

        double tmp;

        // note down the partial at the child 1
        tmp = partialsRoot[v + rootState];
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            originPartials[v + i] = tmp * matricesRootaboveQSstart[w + nrOfStates * i + rootState];
        }
    }

    /**
     * Calculates partial likelihoods at a node when both children have states.
     *
     * @param stateIndex1           alignment at child 1
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param stateIndex2           alignment at child 2
     * @param matricesQS2           transition probability matrix from parent to child 2 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices2aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 2
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param child1QS              QS passing through child 1
     * @param child2QS              QS passing through child 2
     * @param parentQS              QS passing through parent node
     */
    protected void calculateStatesStatesPruning(int[] stateIndex1, double[] matricesQS1, double[] matrices1aboveQSstart,
                                                int[] stateIndex2, double[] matricesQS2, double[] matrices2aboveQSstart,
                                                double[] partials3, int child1QS, int child2QS, int parentQS) {

        // v keeps track of the pattern we are about to calculate
        int v = 0;

        // both children have the QS start on the branch leading to the tip
        if (parentQS==-1){

            for (int l = 0; l < nrOfMatrices; l++) {

                // w keeps track of the state the internal node evolves from
                int w = l * matrixSize;

                for (int k = 0; k < nrOfPatterns; k++) {
                    // note down the states at the tips
                    int state1 = stateIndex1[k];
                    int state2 = stateIndex2[k];

                    // both children have a state
                    if (state1 < nrOfStates && state2 < nrOfStates) {
                        calculateStatesStatesPruningHelperBothQSbelow(matricesQS1, matrices1aboveQSstart,
                                                                      matricesQS2, matrices2aboveQSstart,
                                                                      partials3, w, v, state1, state2);
                    // child 2 has a gap or unknown state so treat it as unknown
                    } else if (state1 < nrOfStates) {
                        calculateStatesStatesPruningHelperBothQSbelowOneUnknownState(matricesQS1, matrices1aboveQSstart,
                                                                                     matricesQS2, matrices2aboveQSstart,
                                                                                     partials3, w, v, state1);
                    // child 1 has a gap or unknown state so treat it as unknown
                    } else if (state2 < nrOfStates) {
                        calculateStatesStatesPruningHelperBothQSbelowOneUnknownState(matricesQS2, matrices2aboveQSstart,
                                                                                     matricesQS1, matrices1aboveQSstart,
                                                                                     partials3, w, v, state2);
                    // both children have a gap or unknown state
                    } else {
                        calculateStatesStatesPruningHelperBothQSbelowBothUnknownState(matricesQS1, matrices1aboveQSstart,
                                                                                      matricesQS2, matrices2aboveQSstart,
                                                                                      partials3, w, v);
                    }

                    v += nrOfStates;

                }
            }
        }
        // there is QS passing through the parent node
        else {
            // child 1 has QS start on the branch leading to the tip
            if (child1QS!=parentQS){

                for (int l = 0; l < nrOfMatrices; l++) {

                    // w keeps track of the state the internal node evolves from
                    int w = l * matrixSize;

                    for (int k = 0; k < nrOfPatterns; k++) {
                        // note down the states at the tips
                        int state1 = stateIndex1[k];
                        int state2 = stateIndex2[k];

                        // both children have a state
                        if (state1 < nrOfStates && state2 < nrOfStates) {
                            calculateStatesStatesPruningHelperOneQSbelow(matricesQS1, matrices1aboveQSstart,
                                                                         matricesQS2,
                                                                         partials3, w, v, state1, state2);
                        // child 2 has a gap or unknown state so treat it as unknown
                        } else if (state1 < nrOfStates) {
                            calculateStatesStatesPruningHelperOneQSbelowAndUnknownState(matricesQS1, matrices1aboveQSstart,
                                                                                        matricesQS2,
                                                                                        partials3, w, v, state1);
                        // child 1 has a gap or unknown state so treat it as unknown
                        } else if (state2 < nrOfStates) {
                            calculateStatesStatesPruningHelperOneQSbelowWithUnknownState(matricesQS1, matrices1aboveQSstart,
                                                                                         matricesQS2,
                                                                                         partials3, w, v, state2);
                        // both children have a gap or unknown state
                        } else {
                            calculateStatesStatesPruningHelperOneQSbelowBothUnknownState(matricesQS1, matrices1aboveQSstart,
                                                                                         matricesQS2,
                                                                                         partials3, w, v);
                        }

                        v += nrOfStates;

                    }
                }
            }
            // child 2 has QS start on the branch leading to the tip
            else if (child2QS!=parentQS){

                for (int l = 0; l < nrOfMatrices; l++) {

                    // w keeps track of the state the internal node evolves from
                    int w = l * matrixSize;

                    for (int k = 0; k < nrOfPatterns; k++) {
                        // note down the states at the tips
                        int state1 = stateIndex1[k];
                        int state2 = stateIndex2[k];

                        // both children have a state
                        if (state1 < nrOfStates && state2 < nrOfStates) {
                            calculateStatesStatesPruningHelperOneQSbelow(matricesQS2, matrices2aboveQSstart,
                                                                         matricesQS1,
                                                                         partials3, w, v, state2, state1);
                        // child 2 has a gap or unknown state so treat it as unknown
                        } else if (state1 < nrOfStates) {
                            calculateStatesStatesPruningHelperOneQSbelowWithUnknownState(matricesQS2, matrices2aboveQSstart,
                                                                                         matricesQS1,
                                                                                         partials3, w, v, state1);
                        // child 1 has a gap or unknown state so treat it as unknown
                        } else if (state2 < nrOfStates) {
                            calculateStatesStatesPruningHelperOneQSbelowAndUnknownState(matricesQS2, matrices2aboveQSstart,
                                                                                        matricesQS1,
                                                                                        partials3, w, v, state2);
                        // both children have a gap or unknown state
                        } else {
                            calculateStatesStatesPruningHelperOneQSbelowBothUnknownState(matricesQS2, matrices2aboveQSstart,
                                                                                         matricesQS1,
                                                                                         partials3, w, v);
                        }

                        v += nrOfStates;

                    }
                }
            }
        }
    }

    /**
     * Helper function to calculateStatesStatesPruning for case of no QS passing through the parent node.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param matricesQS2           transition probability matrix from parent to child 2 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices2aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 2
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                state at child 1
     * @param state2                state at child 2
     */
    protected void calculateStatesStatesPruningHelperBothQSbelow(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] matricesQS2, double[] matrices2aboveQSstart,
            double[] partials3, int w, int v, int state1, int state2) {

        double tmp1, tmp2;

        // note down the transition probabilities (of no change) on the sum of the QS branch lengths
        // P(QS start -> QS tip)
        tmp1 = matricesQS1[w + nrOfStates * state1 + state1];
        tmp2 = matricesQS2[w + nrOfStates * state2 + state2];
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++) {
            partials3[v + i] = tmp1 * matrices1aboveQSstart[w + nrOfStates * i + state1]
                             * tmp2 * matrices2aboveQSstart[w + nrOfStates * i + state2];
        }
    }

    /**
     * Helper function to calculateStatesStatesPruning for case of no QS passing through the parent node and one node with a gap.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param matricesQS2           transition probability matrix from parent to child 2 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices2aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 2
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                state at child 1
     */
    protected void calculateStatesStatesPruningHelperBothQSbelowOneUnknownState(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] matricesQS2, double[] matrices2aboveQSstart,
            double[] partials3, int w, int v, int state1) {

        double tmp1, tmp2, sum2;

        // note down the transition probabilities (of no change) on the sum of the QS branch lengths
        // P(QS start -> QS tip)
        tmp1 = matricesQS1[w + nrOfStates * state1 + state1];
        // since state at tip 2 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            sum2 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                tmp2 = matricesQS2[w + nrOfStates * j + j];
                sum2 += tmp2 * matrices2aboveQSstart[w + nrOfStates * i + j];
            }
            partials3[v + i] = tmp1 * matrices1aboveQSstart[w + nrOfStates * i + state1]
                             * sum2;
        }
    }

    /**
     * Helper function to calculateStatesStatesPruning for case of no QS passing through the parent node and both nodes with a gap.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param matricesQS2           transition probability matrix from parent to child 2 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices2aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 2
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     */
    protected void calculateStatesStatesPruningHelperBothQSbelowBothUnknownState(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] matricesQS2, double[] matrices2aboveQSstart,
            double[] partials3, int w, int v) {

        double tmp1, tmp2, sum1, sum2;

        // note down the transition probabilities (of no change) on the sum of the QS branch lengths
        // P(QS start -> QS tip)
        // since states at tips are unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            sum1 = sum2 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                tmp1 = matricesQS1[w + nrOfStates * j + j];
                tmp2 = matricesQS2[w + nrOfStates * j + j];
                sum1 += tmp1 * matrices1aboveQSstart[w + nrOfStates * i + j];
                sum2 += tmp2 * matrices2aboveQSstart[w + nrOfStates * i + j];
            }
            partials3[v + i] = sum1 * sum2;
        }
    }

    /**
     * Helper function to calculateStatesStatesPruning for case of a QS passing through the parent node.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param matricesQS2           transition probability matrix from parent to child 2 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                state at child 1
     * @param state2                state at child 2
     */
    protected void calculateStatesStatesPruningHelperOneQSbelow(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] matricesQS2,
            double[] partials3, int w, int v, int state1, int state2) {

        double tmp1, tmp2;

        // note down the transition probabilities (of no change) on the sum of the QS branch lengths
        // P(QS start -> QS tip)
        tmp1 = matricesQS1[w + nrOfStates * state1 + state1];
        tmp2 = matricesQS2[w + nrOfStates * state2 + state2];
        // for the state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        // 0 elsewhere (if not the state at the parent)
        for (int i = 0; i < nrOfStates; i++) {
            partials3[v + i] = 0;
        }
        partials3[v + state2] = tmp1 * matrices1aboveQSstart[w + nrOfStates * state2 + state1]
                              * tmp2;
    }

    /**
     * Helper function to calculateStatesStatesPruning for case of a QS passing through the parent node and the corresponding node has a gap.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param matricesQS2           transition probability matrix from parent to child 2 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                state at child 1
     */
    protected void calculateStatesStatesPruningHelperOneQSbelowAndUnknownState(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] matricesQS2,
            double[] partials3, int w, int v, int state1) {

        double tmp1, tmp2;

        // note down the transition probabilities (of no change) on the sum of the QS branch lengths
        // P(QS start -> QS tip)
        tmp1 = matricesQS1[w + nrOfStates * state1 + state1];
        // since state at tip 2 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++) {
            tmp2 = matricesQS2[w + nrOfStates * i + i];
            partials3[v + i] = tmp1 * matrices1aboveQSstart[w + nrOfStates * i + state1]
                             * tmp2;
        }
    }

    /**
     * Helper function to calculateStatesStatesPruning for case of a QS passing through the parent node and the other node has a gap.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param matricesQS2           transition probability matrix from parent to child 2 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state2                state at child 2
     */
    protected void calculateStatesStatesPruningHelperOneQSbelowWithUnknownState(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] matricesQS2,
            double[] partials3, int w, int v, int state2) {

        double tmp1, tmp2, sum1;

        // note down the transition probabilities (of no change) on the sum of the QS branch lengths
        // P(QS start -> QS tip)
        tmp2 = matricesQS2[w + nrOfStates * state2 + state2];
        // since state at tip 1 is unknown, take into account all the possibilities
        // for the state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        sum1 = 0.0;
        for (int i = 0; i < nrOfStates; i++) {
            tmp1 = matricesQS1[w + nrOfStates * i + i];
            sum1 += tmp1 * matrices1aboveQSstart[w + nrOfStates * state2 + i];
            partials3[v + i] = 0;
        }
        partials3[v + state2] = sum1 * tmp2;
    }

    /**
     * Helper function to calculateStatesStatesPruning for case of a QS passing through the parent node both nodes with a gap.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param matricesQS2           transition probability matrix from parent to child 2 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     */
    protected void calculateStatesStatesPruningHelperOneQSbelowBothUnknownState(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] matricesQS2,
            double[] partials3, int w, int v) {

        double tmp1, tmp2, sum1;

        // note down the transition probabilities (of no change) on the sum of the QS branch lengths
        // P(QS start -> QS tip)
        // since states at tips are unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            tmp2 = matricesQS2[w + nrOfStates * i + i];
            sum1 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                tmp1 = matricesQS1[w + nrOfStates * j + j];
                sum1 += tmp1 * matrices1aboveQSstart[w + nrOfStates * i + j];
            }
            partials3[v + i] = sum1 * tmp2;
        }
    }

    /**
     * Calculates partial likelihoods at a node when one child has states and one has partials.
     *
     * @param stateIndex1           alignment at child 1
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param stateIndex2           alignment at node corresponding to the QS passing through child 2
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2             transition probability matrix from parent to child 2
     * @param matrices2aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 2
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param child1QS              QS passing through child 1
     * @param child2QS              QS passing through child 2
     * @param parentQS              QS passing through parent node
     */
    protected void calculateStatesPartialsPruning(int[] stateIndex1, double[] matricesQS1, double[] matrices1aboveQSstart,
                                                  int[] stateIndex2, double[] partials2, double[] matrices2, double[] matrices2aboveQSstart,
                                                  double[] partials3, int child1QS, int child2QS, int parentQS) {

        // v keeps track of the pattern we are about to calculate
        int v = 0;

        // both children have the QS start on the branch leading from the parent to the child or below
        if (parentQS==-1){
            // both children have the QS start on the branch leading from the parent to the child
            if (child2QS!=-1){

                for (int l = 0; l < nrOfMatrices; l++) {

                    // w keeps track of the state the internal node evolves from
                    int w = l * matrixSize;

                    for (int k = 0; k < nrOfPatterns; k++) {
                        // note down the states at the tips (tip being true tip or tip belonging to QS passing through child2)
                        int state1 = stateIndex1[k];
                        int state2 = stateIndex2[k];

                        // child 1 has a state and child2's QS has a state
                        if (state1 < nrOfStates && state2 < nrOfStates) {
                            calculateStatesPartialsPruningHelperBothQSbelow(matricesQS1, matrices1aboveQSstart,
                                                                            partials2, matrices2aboveQSstart,
                                                                            partials3, w, v, state1, state2);
                        // child 1 has a state but child2's QS has a gap or unknown sequence
                        } else if (state1 < nrOfStates){
                            calculateStatesPartialsPruningHelperBothQSbelowQSPartialsUnknownState(matricesQS1, matrices1aboveQSstart,
                                                                                                  partials2, matrices2aboveQSstart,
                                                                                                  partials3, w, v, state1);
                        // child 1 has a gap or unknown state but child2's QS has a state
                        } else if (state2 < nrOfStates){
                            calculateStatesPartialsPruningHelperBothQSbelowTipUnknownState(matricesQS1, matrices1aboveQSstart,
                                                                                           partials2, matrices2aboveQSstart,
                                                                                           partials3, w, v, state2);
                        // both children have a gap or unknown state
                        } else {
                            calculateStatesPartialsPruningHelperBothQSbelowBothUnknownState(matricesQS1, matrices1aboveQSstart,
                                                                                            partials2, matrices2aboveQSstart,
                                                                                            partials3, w, v);
                        }

                        v += nrOfStates;

                    }
                }
            }
            // child 2 (partials) has the QS start below the branch leading from the parent to the child (i.e. no QS on current branch)
            else {

                for (int l = 0; l < nrOfMatrices; l++) {

                    // w keeps track of the state the internal node evolves from
                    int w = l * matrixSize;

                    for (int k = 0; k < nrOfPatterns; k++) {
                        // note down the states at the tips (tip being true tip or tip belonging to QS passing through child2)
                        int state1 = stateIndex1[k];

                        // child 1 has a state
                        if (state1 < nrOfStates){
                            calculateStatesPartialsPruningHelperBothQSbelowQSPartialsUnknownState(matricesQS1, matrices1aboveQSstart,
                                                                                                  partials2, matrices2,
                                                                                                  partials3, w, v, state1);
                        // child1 has a gap or unknown state
                        } else {
                            calculateStatesPartialsPruningHelperBothQSbelowBothUnknownState(matricesQS1, matrices1aboveQSstart,
                                                                                            partials2, matrices2,
                                                                                            partials3, w, v);
                        }

                        v += nrOfStates;

                    }
                }
            }
        }
        // there is QS passing through the parent node
        else {
            // child 1 has QS start on the branch leading to the tip
            if (child1QS!=parentQS){

                for (int l = 0; l < nrOfMatrices; l++) {

                    // w keeps track of the state the internal node evolves from
                    int w = l * matrixSize;

                    for (int k = 0; k < nrOfPatterns; k++) {
                        // note down the states at the tips (tip being true tip or tip belonging to QS passing through child2)
                        int state1 = stateIndex1[k];
                        int state2 = stateIndex2[k];

                        // child 1 has a state and child2's QS has a state
                        if (state1 < nrOfStates && state2 < nrOfStates) {
                            calculateStatesPartialsPruningHelperTipQSbelow(matricesQS1, matrices1aboveQSstart,
                                                                           partials2,
                                                                           partials3, w, v, state1, state2);
                        // child 1 has a state but child2's QS has a gap or unknown sequence
                        } else if (state1 < nrOfStates){
                            calculateStatesPartialsPruningHelperTipQSbelowQSPartialsUnknownState(matricesQS1, matrices1aboveQSstart,
                                                                                                 partials2,
                                                                                                 partials3, w, v, state1);
                        // child 1 has a gap or unknown state but child2's QS has a state
                        } else if (state2 < nrOfStates){
                            calculateStatesPartialsPruningHelperTipQSbelowTipUnknownState(matricesQS1, matrices1aboveQSstart,
                                                                                          partials2,
                                                                                          partials3, w, v, state2);
                        // both children have a gap or unknown state
                        } else {
                            calculateStatesPartialsPruningHelperTipQSbelowBothUnknownState(matricesQS1, matrices1aboveQSstart,
                                                                                           partials2,
                                                                                           partials3, w, v);
                        }

                        v += nrOfStates;

                    }
                }
            }
            // child 2 has QS start on the branch leading from the parent to the child or below
            else if (child2QS!=parentQS){
                // child2 has the QS start on the branch leading from the parent to the child
                if (child2QS!=-1){

                    for (int l = 0; l < nrOfMatrices; l++) {

                        // w keeps track of the state the internal node evolves from
                        int w = l * matrixSize;

                        for (int k = 0; k < nrOfPatterns; k++) {
                            // note down the states at the tips (tip being true tip or tip belonging to QS passing through child2)
                            int state1 = stateIndex1[k];
                            int state2 = stateIndex2[k];

                            // child 1 has a state and child2's QS has a state
                            if (state1 < nrOfStates && state2 < nrOfStates) {
                                calculateStatesPartialsPruningHelperQSPartialsQSbelow(matricesQS1,
                                                                                      partials2, matrices2aboveQSstart,
                                                                                      partials3, w, v, state1, state2);
                            // child 1 has a state but child2's QS has a gap or unknown sequence
                            } else if (state1 < nrOfStates){
                                calculateStatesPartialsPruningHelperQSPartialsQSbelowTipUnknownState(matricesQS1,
                                                                                                     partials2, matrices2aboveQSstart,
                                                                                                     partials3, w, v, state1);
                            // child 1 has a gap or unknown state but child2's QS has a state
                            } else if (state2 < nrOfStates){
                                calculateStatesPartialsPruningHelperQSPartialsQSbelowQSPartialsUnknownState(matricesQS1,
                                                                                                            partials2, matrices2aboveQSstart,
                                                                                                            partials3, w, v, state2);
                            // both children have a gap or unknown state
                            } else {
                                calculateStatesPartialsPruningHelperQSPartialsQSbelowBothUnknownState(matricesQS1,
                                                                                                      partials2, matrices2aboveQSstart,
                                                                                                      partials3, w, v);
                            }

                            v += nrOfStates;

                        }
                    }
                }
                // child 2 (partials) has the QS start below the branch leading from the parent to the child (i.e. no QS on current branch)
                else {

                    for (int l = 0; l < nrOfMatrices; l++) {

                        // w keeps track of the state the internal node evolves from
                        int w = l * matrixSize;

                        for (int k = 0; k < nrOfPatterns; k++) {
                            // note down the states at the tips (tip being true tip or tip belonging to QS passing through child2)
                            int state1 = stateIndex1[k];

                            // child 1 has a state
                            if (state1 < nrOfStates){
                                calculateStatesPartialsPruningHelperQSPartialsQSbelowTipUnknownState(matricesQS1,
                                                                                                     partials2, matrices2,
                                                                                                     partials3, w, v, state1);
                            // child1 has a gap or unknown state
                            } else {
                                calculateStatesPartialsPruningHelperQSPartialsQSbelowBothUnknownState(matricesQS1,
                                                                                                      partials2, matrices2,
                                                                                                      partials3, w, v);
                            }

                            v += nrOfStates;

                        }
                    }
                }
            }
            else{
                throw new IllegalStateException("What case did we miss in QuasiSpeciesBeerLikelihoodCore function calculateStatesPartialsPruning");
            }
        }
    }

    /**
     * Helper function to calculateStatesPartialsPruning for case of no QS passing through the parent node.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 2
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                state at child 1
     * @param state2                state at child 2
     */
    protected void calculateStatesPartialsPruningHelperBothQSbelow(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] partials2, double[] matrices2aboveQSstart,
            double[] partials3, int w, int v, int state1, int state2) {

        double tmp1, tmp2;

        // note down the transition probability (of no change) on the sum of the QS branch lengths
        // P(QS start -> QS tip)
        tmp1 = matricesQS1[w + nrOfStates * state1 + state1];
        // note down the partial at the child 2
        tmp2 = partials2[v + state2];
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++) {
            partials3[v + i] = tmp1 * matrices1aboveQSstart[w + nrOfStates * i + state1]
                             * tmp2 * matrices2aboveQSstart[w + nrOfStates * i + state2];
        }
    }

    /**
     * Helper function to calculateStatesPartialssPruning for case of no QS passing through the parent node and the (QS from partials) node has a gap.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2aboveQSstartOrMatrices2 transition probability matrix from node above QS start to QS start for QS passing through child 2 / transition probability matrix between two nodes
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                state at child 1
     */
    protected void calculateStatesPartialsPruningHelperBothQSbelowQSPartialsUnknownState(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] partials2, double[] matrices2aboveQSstartOrMatrices2,
            double[] partials3, int w, int v, int state1) {

        double tmp1, tmp2, sum2;

        // note down the transition probability (of no change) on the sum of the QS branch lengths
        // P(QS start -> QS tip)
        tmp1 = matricesQS1[w + nrOfStates * state1 + state1];
        // since state at child 2 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            sum2 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 2
                tmp2 = partials2[v + j];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[w + nrOfStates * i + j];
            }
            partials3[v + i] = tmp1 * matrices1aboveQSstart[w + nrOfStates * i + state1]
                             * sum2;
        }
    }

    /**
     * Helper function to calculateStatesPartialsPruning for case of no QS passing through the parent node and the (tip) node has a gap.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 2
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state2                state at child 2
     */
    protected void calculateStatesPartialsPruningHelperBothQSbelowTipUnknownState(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] partials2, double[] matrices2aboveQSstart,
            double[] partials3, int w, int v, int state2) {

        double tmp1, tmp2, sum1;

        // note down the partial at the child 2
        tmp2 = partials2[v + state2];
        // since state at tip 1 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            sum1 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the transition probability (of no change) on the sum of the QS branch lengths
                // P(QS start -> QS tip)
                tmp1 = matricesQS1[w + nrOfStates * j + j];
                sum1 += tmp1 * matrices1aboveQSstart[w + nrOfStates * i + j];
            }
            partials3[v + i] = sum1
                             * tmp2 * matrices2aboveQSstart[w + nrOfStates * i + state2];
        }
    }

    /**
     * Helper function to calculateStatesPartialsPruning for case of no QS passing through the parent node and both nodes have a gap.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2aboveQSstartOrMatrices2 transition probability matrix from node above QS start to QS start for QS passing through child 2 / transition probability matrix between two nodes
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     */
    protected void calculateStatesPartialsPruningHelperBothQSbelowBothUnknownState(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] partials2, double[] matrices2aboveQSstartOrMatrices2,
            double[] partials3, int w, int v) {

        double tmp1, tmp2, sum1, sum2;

        // since states at tip/QS tip are unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            sum1 = sum2 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the transition probability (of no change) on the sum of the QS branch lengths
                // P(QS start -> QS tip)
                tmp1 = matricesQS1[w + nrOfStates * j + j];
                // note down the partial at the child 2
                tmp2 = partials2[v + j];
                sum1 += tmp1 * matrices1aboveQSstart[w + nrOfStates * i + j];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[w + nrOfStates * i + j];
            }
            partials3[v + i] = sum1 * sum2;
        }
    }

    /**
     * Helper function to calculateStatesPartialsPruning for case of a QS passing through the parent node.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                state at child 1
     * @param state2                state at child 2
     */
    protected void calculateStatesPartialsPruningHelperTipQSbelow(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] partials2,
            double[] partials3, int w, int v, int state1, int state2) {

        double tmp1, tmp2;

        // note down the transition probability (of no change) on the sum of the QS branch lengths
        // P(QS start -> QS tip)
        tmp1 = matricesQS1[w + nrOfStates * state1 + state1];
        // note down the partial at the child 2
        tmp2 = partials2[v + state2];
        // for the state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        // 0 elsewhere (if not the state at the parent)
        for (int i = 0; i < nrOfStates; i++) {
            partials3[v + i] = 0;
        }
        partials3[v + state2] = tmp1 * matrices1aboveQSstart[w + nrOfStates * state2 + state1]
                              * tmp2;
    }

    /**
     * Helper function to calculateStatesPartialsPruning for case of a QS passing through the parent node and the (QS from partials) node has a gap.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                state at child 1
     */
    protected void calculateStatesPartialsPruningHelperTipQSbelowQSPartialsUnknownState(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] partials2,
            double[] partials3, int w, int v, int state1) {

        double tmp1, tmp2;

        // note down the transition probability (of no change) on the sum of the QS branch lengths
        // P(QS start -> QS tip)
        tmp1 = matricesQS1[w + nrOfStates * state1 + state1];
        // since state at child 2 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            // note down the partial at the child 2
            tmp2 = partials2[v + i];
            partials3[v + i] = tmp1 * matrices1aboveQSstart[w + nrOfStates * i + state1]
                             * tmp2;
        }
    }

    /**
     * Helper function to calculateStatesPartialsPruning for case of a QS passing through the parent node and the (tip) node has a gap.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state2                state at child 2
     */
    protected void calculateStatesPartialsPruningHelperTipQSbelowTipUnknownState(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] partials2,
            double[] partials3, int w, int v, int state2) {

        double tmp1, tmp2, sum1;

        // note down the partial at the child 2
        tmp2 = partials2[v + state2];
        // since state at tip 1 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        sum1 = 0.0;
        for (int i = 0; i < nrOfStates; i++){
            partials3[v + i] = 0;
            // note down the transition probability (of no change) on the sum of the QS branch lengths
            // P(QS start -> QS tip)
            tmp1 = matricesQS1[w + nrOfStates * i + i];
            sum1 += tmp1 * matrices1aboveQSstart[w + nrOfStates * state2 + i];
        }
        partials3[v + state2] = sum1 * tmp2;
    }

    /**
     * Helper function to calculateStatesPartialsPruning for case of a QS passing through the parent node and both nodes have a gap.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     */
    protected void calculateStatesPartialsPruningHelperTipQSbelowBothUnknownState(
            double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] partials2,
            double[] partials3, int w, int v) {

        double tmp1, tmp2, sum1;

        // since states at tip/QS tip are unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            // note down the partial at the child 2
            tmp2 = partials2[v + i];
            sum1 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the transition probability (of no change) on the sum of the QS branch lengths
                // P(QS start -> QS tip)
                tmp1 = matricesQS1[w + nrOfStates * j + j];
                sum1 += tmp1 * matrices1aboveQSstart[w + nrOfStates * i + j];
            }
            partials3[v + i] = sum1 * tmp2;
        }
    }

    /**
     * Helper function to calculateStatesPartialsPruning for case of a QS passing through the parent node.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 2
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                state at child 1
     * @param state2                state at child 2
     */
    protected void calculateStatesPartialsPruningHelperQSPartialsQSbelow(
            double[] matricesQS1,
            double[] partials2, double[] matrices2aboveQSstart,
            double[] partials3, int w, int v, int state1, int state2) {

        double tmp1, tmp2;

        // note down the transition probability (of no change) on the sum of the QS branch lengths
        // P(QS start -> QS tip)
        tmp1 = matricesQS1[w + nrOfStates * state1 + state1];
        // note down the partial at the child 2
        tmp2 = partials2[v + state2];
        // for the state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        // 0 elsewhere (if not the state at the parent)
        for (int i = 0; i < nrOfStates; i++) {
            partials3[v + i] = 0;
        }
        partials3[v + state1] = tmp1
                              * tmp2 * matrices2aboveQSstart[w + nrOfStates * state1 + state2];

    }

    /**
     * Helper function to calculateStatesPartialsPruning for case of a QS passing through the parent node and the (QS from partials) node has a gap.
     *
     * @param matricesQS1                       transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param partials2                         probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2aboveQSstartOrMatrices2  transition probability matrix from node above QS start to QS start for QS passing through child 2 / transition probability matrix between two nodes
     * @param partials3                         probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                            state at child 1
     */
    protected void calculateStatesPartialsPruningHelperQSPartialsQSbelowTipUnknownState(
            double[] matricesQS1,
            double[] partials2, double[] matrices2aboveQSstartOrMatrices2,
            double[] partials3, int w, int v, int state1) {

        double tmp1, tmp2, sum2;

        // note down the transition probability (of no change) on the sum of the QS branch lengths
        // P(QS start -> QS tip)
        tmp1 = matricesQS1[w + nrOfStates * state1 + state1];
        // since state at child 2 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        sum2 = 0.0;
        for (int i = 0; i < nrOfStates; i++){
            partials3[v + i] = 0;
            // note down the partial at the child 2
            tmp2 = partials2[v + i];
            sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[w + nrOfStates * state1 + i];
        }
        partials3[v + state1] = tmp1 * sum2;
    }

    /**
     * Helper function to calculateStatesPartialsPruning for case of a QS passing through the parent node and the (tip) node has a gap.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 2
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state2                state at child 2
     */
    protected void calculateStatesPartialsPruningHelperQSPartialsQSbelowQSPartialsUnknownState(
            double[] matricesQS1,
            double[] partials2, double[] matrices2aboveQSstart,
            double[] partials3, int w, int v, int state2) {

        double tmp1, tmp2;

        // note down the partial at the child 2
        tmp2 = partials2[v + state2];
        // since state at tip 1 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            // note down the transition probability (of no change) on the sum of the QS branch lengths
            // P(QS start -> QS tip)
            tmp1 = matricesQS1[w + nrOfStates * i + i];
            partials3[v + i] = tmp1
                             * tmp2 * matrices2aboveQSstart[w + nrOfStates * i + state2];
        }
    }

    /**
     * Helper function to calculateStatesPartialsPruning for case of a QS passing through the parent node and both nodes have a gap.
     *
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2aboveQSstartOrMatrices2 transition probability matrix from node above QS start to QS start for QS passing through child 2 / transition probability matrix between two nodes
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     */
    protected void calculateStatesPartialsPruningHelperQSPartialsQSbelowBothUnknownState(
            double[] matricesQS1,
            double[] partials2, double[] matrices2aboveQSstartOrMatrices2,
            double[] partials3, int w, int v) {

        double tmp1, tmp2, sum2;

        // since states at tip/QS tip are unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            // note down the transition probability (of no change) on the sum of the QS branch lengths
            // P(QS start -> QS tip)
            tmp1 = matricesQS1[w + nrOfStates * i + i];
            sum2 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 2
                tmp2 = partials2[v + j];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[w + nrOfStates * i + j];
            }
            partials3[v + i] = tmp1 * sum2;
        }
    }







    /**
     * Helper function to calculatePartialsPartialsPruning for case of no QS passing through the parent node.
     *
     * @param partials1             probability vector at child 1 (of length nrOfStates * nrOfPatterns)
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 2
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                state at child 1
     * @param state2                state at child 2
     */
    protected void calculatePartialsPartialsPruningHelperBothQSbelow(
            double[] partials1, double[] matrices1aboveQSstart,
            double[] partials2, double[] matrices2aboveQSstart,
            double[] partials3, int w, int v, int state1, int state2) {

        double tmp1, tmp2;

        // note down the partial at the child1/child 2
        tmp1 = partials1[v + state1];
        tmp2 = partials2[v + state2];
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++) {
            partials3[v + i] = tmp1 * matrices1aboveQSstart[w + nrOfStates * i + state1]
                             * tmp2 * matrices2aboveQSstart[w + nrOfStates * i + state2];
        }
    }

    /**
     * Helper function to calculatePartialsPartialsPruning for case of no QS passing through the parent node and one node has a gap.
     *
     * @param partials1                         probability vector at child 1 (of length nrOfStates * nrOfPatterns)
     * @param matrices1aboveQSstart             transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2                         probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2aboveQSstartOrMatrices2  transition probability matrix from node above QS start to QS start for QS passing through child 2 / transition probability matrix between two nodes
     * @param partials3                         probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                            state at child 1
     */
    protected void calculatePartialsPartialsPruningHelperBothQSbelowOnePartialsUnknownState(
            double[] partials1, double[] matrices1aboveQSstart,
            double[] partials2, double[] matrices2aboveQSstartOrMatrices2,
            double[] partials3, int w, int v, int state1) {

        double tmp1, tmp2, sum2;

        // note down the partial at the child 1
        tmp1 = partials1[v + state1];
        // since state at child 2 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            sum2 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 2
                tmp2 = partials2[v + j];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[w + nrOfStates * i + j];
            }
            partials3[v + i] = tmp1 * matrices1aboveQSstart[w + nrOfStates * i + state1]
                             * sum2;
        }
    }

    /**
     * Helper function to calculatePartialsPartialsPruning for case of no QS passing through the parent node and both nodes have a gap.
     *
     * @param partials1                         probability vector at child 1 (of length nrOfStates * nrOfPatterns)
     * @param matrices1aboveQSstartOrMatrices1  transition probability matrix from node above QS start to QS start for QS passing through child 1 / transition probability matrix between two nodes
     * @param partials2                         probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2aboveQSstartOrMatrices2  transition probability matrix from node above QS start to QS start for QS passing through child 2 / transition probability matrix between two nodes
     * @param partials3                         probability vector at parent node (of length nrOfStates * nrOfPatterns)
     */
    protected void calculatePartialsPartialsPruningHelperBothQSbelowBothUnknownState(
            double[] partials1, double[] matrices1aboveQSstartOrMatrices1,
            double[] partials2, double[] matrices2aboveQSstartOrMatrices2,
            double[] partials3, int w, int v) {

        double tmp1, tmp2, sum1, sum2;

        // since states at QS tips are unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            sum1 = sum2 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 1
                tmp1 = partials1[v + j];
                // note down the partial at the child 2
                tmp2 = partials2[v + j];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[w + nrOfStates * i + j];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[w + nrOfStates * i + j];
            }
            partials3[v + i] = sum1 * sum2;
        }
    }

    /**
     * Helper function to calculatePartialsPartialsPruning for case of no QS passing through the parent node and both nodes have a gap.
     *
     * @param partials1             probability vector at child 1 (of length nrOfStates * nrOfPatterns)
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matricesQS2           transition probability matrix from parent to child 2 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices2aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 2
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     */
    protected void calculatePartialsPartialsPruningHelperBothQSbelowBothQSTipsWithPartials(
            double[] partials1, double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] partials2, double[] matricesQS2, double[] matrices2aboveQSstart,
            double[] partials3, int w, int v){

        double tmp1, tmp2, sum1, sum2;

        // since states at QS tips are partials, take into account all viable possibilities
        // for each state at the parent node calculate prob. of going to the state (partial) of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            sum1 = sum2 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 1
                tmp1 = partials1[v + j] * matricesQS1[w + nrOfStates * j + j];
                // note down the partial at the child 2
                tmp2 = partials2[v + j] * matricesQS2[w + nrOfStates * j + j];
                sum1 += tmp1 * matrices1aboveQSstart[w + nrOfStates * i + j];
                sum2 += tmp2 * matrices2aboveQSstart[w + nrOfStates * i + j];
            }
            partials3[v + i] = sum1 * sum2;
        }
    }

    /**
     * Helper function to calculatePartialsPartialsPruning for case of no QS passing through the parent node and both nodes have a gap.
     *
     * @param partials1             probability vector at child 1 (of length nrOfStates * nrOfPatterns)
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 2
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     */
    protected void calculatePartialsPartialsPruningHelperBothQSbelowOneQSTipWithPartials(
            double[] partials1, double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] partials2, double[] matrices2aboveQSstart,
            double[] partials3, int w, int v){

        double tmp1, tmp2, sum1, sum2;

        // since states at QS tips are partials, take into account all viable possibilities
        // for each state at the parent node calculate prob. of going to the state (partial) of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            sum1 = sum2 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 1
                tmp1 = partials1[v + j] * matricesQS1[w + nrOfStates * j + j];
                // note down the partial at the child 2
                tmp2 = partials2[v + j];
                sum1 += tmp1 * matrices1aboveQSstart[w + nrOfStates * i + j];
                sum2 += tmp2 * matrices2aboveQSstart[w + nrOfStates * i + j];
            }
            partials3[v + i] = sum1 * sum2;
        }
    }

    /**
     * Helper function to calculatePartialsPartialsPruning for case of a QS passing through the parent node.
     *
     * @param partials1             probability vector at child 1 (of length nrOfStates * nrOfPatterns)
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                state at child 1
     * @param state2                state at child 2
     */
    protected void calculatePartialsPartialsPruningHelperOneQSbelow(
            double[] partials1, double[] matrices1aboveQSstart,
            double[] partials2,
            double[] partials3, int w, int v, int state1, int state2) {

        double tmp1, tmp2;

        // note down the partial at the child 1/child 2
        tmp1 = partials1[v + state1];
        tmp2 = partials2[v + state2];
        // for the state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        // 0 elsewhere (if not the state at the parent)
        for (int i = 0; i < nrOfStates; i++) {
            partials3[v + i] = 0;
        }
        partials3[v + state2] = tmp1 * matrices1aboveQSstart[w + nrOfStates * state2 + state1]
                              * tmp2;
    }

    /**
     * Helper function to calculatePartialsPartialsPruning for case of a QS passing through the parent node and the corresponding node has a gap.
     *
     * @param partials1                         probability vector at child 1 (of length nrOfStates * nrOfPatterns)
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2                         probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state1                state at child 1
     */
    protected void calculatePartialsPartialsPruningHelperOneQSbelowAndPartialsUnknownState(
            double[] partials1, double[] matrices1aboveQSstart,
            double[] partials2,
            double[] partials3, int w, int v, int state1) {

        double tmp1, tmp2;

        // note down the partial at the child 1
        tmp1 = partials1[v + state1];
        // since state at child 2 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            // note down the partial at the child 2
            tmp2 = partials2[v + i];
            partials3[v + i] = tmp1 * matrices1aboveQSstart[w + nrOfStates * i + state1]
                             * tmp2;
        }
    }

    /**
     * Helper function to calculatePartialsPartialsPruning for case of a QS passing through the parent node and the other node has a gap.
     *
     * @param partials1                         probability vector at child 1 (of length nrOfStates * nrOfPatterns)
     * @param matrices1aboveQSstartOrMatrices1  transition probability matrix from node above QS start to QS start for QS passing through child 1 / transition probability matrix between two nodes
     * @param partials2                         probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param partials3                         probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param state2                            state at child 2
     */
    protected void calculatePartialsPartialsPruningHelperOneQSbelowWithPartialsUnknownState(
            double[] partials1, double[] matrices1aboveQSstartOrMatrices1,
            double[] partials2,
            double[] partials3, int w, int v, int state2) {

        double tmp1, tmp2, sum1;

        // note down the partial at the child 2
        tmp2 = partials2[v + state2];
        // since state at child 1 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        sum1 = 0.0;
        for (int i = 0; i < nrOfStates; i++){
            partials3[v + i] = 0;
            // note down the partial at the child 1
            tmp1 = partials1[v + i];
            sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[w + nrOfStates * state2 + i];
        }
        partials3[v + state2] = sum1 * tmp2;
    }

    /**
     * Helper function to calculatePartialsPartialsPruning for case of a QS passing through the parent node both nodes with a gap.
     *
     * @param partials1                         probability vector at child 1 (of length nrOfStates * nrOfPatterns)
     * @param matrices1aboveQSstartOrMatrices1  transition probability matrix from node above QS start to QS start for QS passing through child 1 / transition probability matrix between two nodes
     * @param partials2                         probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param partials3                         probability vector at parent node (of length nrOfStates * nrOfPatterns)
     */
    protected void calculatePartialsPartialsPruningHelperOneQSbelowBothUnknownState(
            double[] partials1, double[] matrices1aboveQSstartOrMatrices1,
            double[] partials2,
            double[] partials3, int w, int v) {

        double tmp1, tmp2, sum1;

        // since states at QS tips are unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            // note down the partial at the child 2
            tmp2 = partials2[v + i];
            sum1 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 1
                tmp1 = partials1[v + j];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[w + nrOfStates * i + j];
            }
            partials3[v + i] = sum1 * tmp2;
        }
    }

    /**
     * Helper function to calculatePartialsPartialsPruning for case of a QS passing through the parent node and both nodes have a gap.
     *
     * @param partials1             probability vector at child 1 (of length nrOfStates * nrOfPatterns)
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matricesQS2           transition probability matrix from parent to child 2 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     */
    protected void calculatePartialsPartialsPruningHelperOneQSbelowBothQSTipsWithPartials(
            double[] partials1, double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] partials2, double[] matricesQS2,
            double[] partials3, int w, int v){


        double tmp1, tmp2, sum1;

        // since states at QS tips are partials, take into account all viable possibilities
        // for each state at the parent node calculate prob. of going to the state (partial) of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            // note down the partial at the child 2
            tmp2 = partials2[v + i] * matricesQS2[w + nrOfStates * i + i];
            sum1 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 1
                tmp1 = partials1[v + j] * matricesQS1[w + nrOfStates * j + j];
                sum1 += tmp1 * matrices1aboveQSstart[w + nrOfStates * i + j];
            }
            partials3[v + i] = sum1 * tmp2;
        }
    }

    /**
     * Helper function to calculatePartialsPartialsPruning for case of a QS passing through the parent node and both nodes have a gap.
     *
     * @param partials1             probability vector at child 1 (of length nrOfStates * nrOfPatterns)
     * @param matricesQS1           transition probability matrix from parent to child 1 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     */
    protected void calculatePartialsPartialsPruningHelperOneQSbelowOneQSTipWithPartials(
            double[] partials1, double[] matricesQS1, double[] matrices1aboveQSstart,
            double[] partials2,
            double[] partials3, int w, int v){

        double tmp1, tmp2, sum1;

        // since states at QS tips are partials, take into account all viable possibilities
        // for each state at the parent node calculate prob. of going to the state (partial) of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            // note down the partial at the child 2
            tmp2 = partials2[v + i];
            sum1 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 1
                tmp1 = partials1[v + j] * matricesQS1[w + nrOfStates * j + j];
                sum1 += tmp1 * matrices1aboveQSstart[w + nrOfStates * i + j];
            }
            partials3[v + i] = sum1 * tmp2;
        }
    }

    /**
     * Helper function to calculatePartialsPartialsPruning for case of a QS passing through the parent node and both nodes have a gap.
     *
     * @param partials1                         probability vector at child 1 (of length nrOfStates * nrOfPatterns)
     * @param matrices1aboveQSstartOrMatrices1  transition probability matrix from node above QS start to QS start for QS passing through child 1 / transition probability matrix between two nodes
     * @param partials2                         probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matricesQS2                       transition probability matrix from parent to child 2 - if the node is a tip, it holds the probability that the sequence does not change from the tip to the start of the haplo
     * @param partials3                         probability vector at parent node (of length nrOfStates * nrOfPatterns)
     */
    protected void calculatePartialsPartialsPruningHelperOneQSbelowOtherQSTipWithPartials(
            double[] partials1, double[] matrices1aboveQSstartOrMatrices1,
            double[] partials2, double[] matricesQS2,
            double[] partials3, int w, int v){


        double tmp1, tmp2, sum1;

        // since states at QS tips are partials, take into account all viable possibilities
        // for each state at the parent node calculate prob. of going to the state (partial) of the tip * P(QS start -> QS tip)
        for (int i = 0; i < nrOfStates; i++){
            // note down the partial at the child 2
            tmp2 = partials2[v + i] * matricesQS2[w + nrOfStates * i + i];
            sum1 = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 1
                tmp1 = partials1[v + j];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[w + nrOfStates * i + j];
            }
            partials3[v + i] = sum1 * tmp2;
        }
    }




    // TODO continue here  --- include the case where the QS tips have partials and not states





    /**
     * Calculates partial likelihoods at a node when both children have partials.
     *
     * @param stateIndex1           alignment at child 1
     * @param partials1             probability vector at child 1 (of length nrOfStates * nrOfPatterns)
     * @param matrices1             transition probability matrix from parent to child 1
     * @param matrices1aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 1
     * @param stateIndex2           alignment at node corresponding to the QS passing through child 2
     * @param partials2             probability vector at child 2 (of length nrOfStates * nrOfPatterns)
     * @param matrices2             transition probability matrix from parent to child 2
     * @param matrices2aboveQSstart transition probability matrix from node above QS start to QS start for QS passing through child 2
     * @param partials3             probability vector at parent node (of length nrOfStates * nrOfPatterns)
     * @param child1QS              QS passing through child 1
     * @param child2QS              QS passing through child 2
     * @param parentQS              QS passing through parent node
     */
    protected void calculatePartialsPartialsPruning(int[] stateIndex1, double[] partials1, double[] matrices1, double[] matrices1aboveQSstart,
                                                    int[] stateIndex2, double[] partials2, double[] matrices2, double[] matrices2aboveQSstart,
                                                    double[] partials3, int child1QS, int child2QS, int parentQS, int nodeIndex1, int nodeIndex2) {

        double sum1, sum2, tmp1, tmp2;

        // v keeps track of the pattern we are about to calculate
        int v = 0;

        // both children have the QS start on the branch leading from the parent to the child or below
        if (parentQS==-1){
            // both children have the QS start on the branch leading from the parent to the child
            if (child1QS!=-1 && child2QS!=-1) {

                // the QS passing the nodes below have a state (as opposed to uncertain state = partials)
                if (stateIndex1 != null && stateIndex2 != null){

                    for (int l = 0; l < nrOfMatrices; l++) {

                        // w keeps track of the state the internal node evolves from
                        int w = l * matrixSize;

                        for (int k = 0; k < nrOfPatterns; k++) {
                            // note down the states at the tips belonging to QS passing through child1/child2
                            int state1 = stateIndex1[k];
                            int state2 = stateIndex2[k];

                            // child1's QS and child2's QS have a state
                            if (state1 < nrOfStates && state2 < nrOfStates) {
                                calculatePartialsPartialsPruningHelperBothQSbelow(partials1, matrices1aboveQSstart,
                                                                                  partials2, matrices2aboveQSstart,
                                                                                  partials3, w, v, state1, state2);
                            // child1's QS has a state but child2's QS has a gap or unknown sequence
                            } else if (state1 < nrOfStates){
                                calculatePartialsPartialsPruningHelperBothQSbelowOnePartialsUnknownState(partials1, matrices1aboveQSstart,
                                                                                                         partials2, matrices2aboveQSstart,
                                                                                                         partials3, w, v, state1);
                            // child1's QS has a gap or unknown state but child2's QS has a state
                            } else if (state2 < nrOfStates){
                                calculatePartialsPartialsPruningHelperBothQSbelowOnePartialsUnknownState(partials2, matrices2aboveQSstart,
                                                                                                         partials1, matrices1aboveQSstart,
                                                                                                         partials3, w, v, state2);
                            // both children have a gap or unknown state
                            } else {
                                calculatePartialsPartialsPruningHelperBothQSbelowBothUnknownState(partials1, matrices1aboveQSstart,
                                                                                                  partials2, matrices2aboveQSstart,
                                                                                                  partials3, w, v);
                            }

                            v += nrOfStates;

                        }
                    }
                }
                // the QS do not have a state so they have partials
                else {

                    for (int l = 0; l < nrOfMatrices; l++) {

                        // w keeps track of the state the internal node evolves from
                        int w = l * matrixSize;

                        for (int k = 0; k < nrOfPatterns; k++) {

                            // check if both of the child nodes is a tip
                            if (nodeIndex1 == child1QS && nodeIndex2 == child2QS) {
                                // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                calculatePartialsPartialsPruningHelperBothQSbelowBothQSTipsWithPartials(partials1, matrices1,
                                                                                                        matrices1aboveQSstart,
                                                                                                        partials2, matrices2,
                                                                                                        matrices2aboveQSstart,
                                                                                                        partials3, w, v);
                            }
                            // or just child node 1 is a tip
                            else if (nodeIndex1 == child1QS) {
                                // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                calculatePartialsPartialsPruningHelperBothQSbelowOneQSTipWithPartials(partials1, matrices1,
                                                                                                      matrices1aboveQSstart,
                                                                                                      partials2,
                                                                                                      matrices2aboveQSstart,
                                                                                                      partials3, w, v);
                            }
                            // or just child node 2 is a tip
                            else if (nodeIndex2 == child2QS) {
                                // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                calculatePartialsPartialsPruningHelperBothQSbelowOneQSTipWithPartials(partials2, matrices2,
                                                                                                      matrices2aboveQSstart,
                                                                                                      partials1,
                                                                                                      matrices1aboveQSstart,
                                                                                                      partials3, w, v);
                            }
                            // or no child node is a tip
                            else {
                                // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                calculatePartialsPartialsPruningHelperBothQSbelowBothUnknownState(partials1, matrices1aboveQSstart,
                                                                                                  partials2, matrices2aboveQSstart,
                                                                                                  partials3, w, v);
                            }

                            v += nrOfStates;

                        }
                    }
                }
            }
            // child 2 (partials) has the QS start below the branch leading from the parent to the child (i.e. no QS on current branch)
            else if (child1QS!=-1) {

                // the QS passing the child 1 below has a state (as opposed to uncertain state = partials)
                if (stateIndex1 != null) {

                    for (int l = 0; l < nrOfMatrices; l++) {

                        // w keeps track of the state the internal node evolves from
                        int w = l * matrixSize;

                        for (int k = 0; k < nrOfPatterns; k++) {
                            // note down the states at the tips belonging to QS passing through child1/child2
                            int state1 = stateIndex1[k];

                            // child1's QS has a state
                            if (state1 < nrOfStates){
                                calculatePartialsPartialsPruningHelperBothQSbelowOnePartialsUnknownState(partials1, matrices1aboveQSstart,
                                                                                                         partials2, matrices2,
                                                                                                         partials3, w, v, state1);
                            // child1's QS has a gap or unknown state
                            } else {
                                calculatePartialsPartialsPruningHelperBothQSbelowBothUnknownState(partials1, matrices1aboveQSstart,
                                                                                                  partials2, matrices2,
                                                                                                  partials3, w, v);
                            }

                            v += nrOfStates;

                        }
                    }
                }
                // the QS do not have a state so they have partials
                else {

                    for (int l = 0; l < nrOfMatrices; l++) {

                        // w keeps track of the state the internal node evolves from
                        int w = l * matrixSize;

                        for (int k = 0; k < nrOfPatterns; k++) {

                            // check if the child node 1 is a tip
                            if (nodeIndex1 == child1QS) {
                                // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                calculatePartialsPartialsPruningHelperBothQSbelowOneQSTipWithPartials(partials1, matrices1,
                                                                                                      matrices1aboveQSstart,
                                                                                                      partials2,
                                                                                                      matrices2aboveQSstart,
                                                                                                      partials3, w, v);
                            }
                            // or child 1 node is not a tip
                            else {
                                // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                calculatePartialsPartialsPruningHelperBothQSbelowBothUnknownState(partials1, matrices1aboveQSstart,
                                                                                                  partials2, matrices2aboveQSstart,
                                                                                                  partials3, w, v);
                            }

                            v += nrOfStates;

                        }
                    }
                }
            }
            // child 1 (partials) has the QS start below the branch leading from the parent to the child (i.e. no QS on current branch)
            else if (child2QS!=-1) {

                // the QS passing the child 2 below has a state (as opposed to uncertain state = partials)
                if (stateIndex2 != null) {

                    for (int l = 0; l < nrOfMatrices; l++) {

                        // w keeps track of the state the internal node evolves from
                        int w = l * matrixSize;

                        for (int k = 0; k < nrOfPatterns; k++) {
                            // note down the states at the tips belonging to QS passing through child1/child2
                            int state2 = stateIndex2[k];

                            // child2's QS has a state
                            if (state2 < nrOfStates) {
                                calculatePartialsPartialsPruningHelperBothQSbelowOnePartialsUnknownState(partials2, matrices2aboveQSstart,
                                                                                                         partials1, matrices1,
                                                                                                         partials3, w, v, state2);
                            // child2's QS has a gap or unknown state
                            } else {
                                calculatePartialsPartialsPruningHelperBothQSbelowBothUnknownState(partials2, matrices2aboveQSstart,
                                                                                                  partials1, matrices1,
                                                                                                  partials3, w, v);
                            }

                            v += nrOfStates;

                        }
                    }
                }
                // the QS do not have a state so they have partials
                else {

                    for (int l = 0; l < nrOfMatrices; l++) {

                        // w keeps track of the state the internal node evolves from
                        int w = l * matrixSize;

                        for (int k = 0; k < nrOfPatterns; k++) {

                            // check if the child node 2 is a tip
                            if (nodeIndex2 == child2QS) {
                                // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                calculatePartialsPartialsPruningHelperBothQSbelowOneQSTipWithPartials(partials2, matrices2,
                                                                                                      matrices2aboveQSstart,
                                                                                                      partials1,
                                                                                                      matrices1aboveQSstart,
                                                                                                      partials3, w, v);
                            }
                            // or child 2 node is not a tip
                            else {
                                // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                calculatePartialsPartialsPruningHelperBothQSbelowBothUnknownState(partials1, matrices1aboveQSstart,
                                                                                                  partials2, matrices2aboveQSstart,
                                                                                                  partials3, w, v);
                            }

                            v += nrOfStates;

                        }
                    }
                }
            }
            // both children the QS start below the branch leading from the parent to the child (i.e. no QS on current branch)
            else {
                for (int l = 0; l < nrOfMatrices; l++) {

                    // w keeps track of the state the internal node evolves from
                    int w = l * matrixSize;

                    for (int k = 0; k < nrOfPatterns; k++) {
                        // note down the states at the tips belonging to QS passing through child1/child2
                        calculatePartialsPartialsPruningHelperBothQSbelowBothUnknownState(partials1, matrices1,
                                                                                          partials2, matrices2,
                                                                                          partials3, w, v);

                        v += nrOfStates;

                    }
                }
            }
        }
        // there is QS passing through the parent node
        else {
            // child 1 has QS start on the branch leading from the parent to the child or below
            if (child1QS!=parentQS){

                // child1 has the QS start on the branch leading from the parent to the child
                if (child1QS!=-1) {

                    // the QS passing the nodes below have a state (as opposed to uncertain state = partials)
                    if (stateIndex1 != null && stateIndex2 != null){

                        for (int l = 0; l < nrOfMatrices; l++) {

                            // w keeps track of the state the internal node evolves from
                            int w = l * matrixSize;

                            for (int k = 0; k < nrOfPatterns; k++) {
                                // note down the states at the tips belonging to QS passing through child1/child2
                                int state1 = stateIndex1[k];
                                int state2 = stateIndex2[k];

                                // child1's QS and child2's QS have a state
                                if (state1 < nrOfStates && state2 < nrOfStates) {
                                    calculatePartialsPartialsPruningHelperOneQSbelow(partials1, matrices1aboveQSstart,
                                                                                     partials2,
                                                                                     partials3, w, v, state1, state2);
                                // child1's QS has a state but child2's QS has a gap or unknown sequence
                                } else if (state1 < nrOfStates){
                                    calculatePartialsPartialsPruningHelperOneQSbelowAndPartialsUnknownState(partials1, matrices1aboveQSstart,
                                                                                                            partials2,
                                                                                                            partials3, w, v, state1);

                                // child1's QS has a gap or unknown state but child2's QS has a state
                                } else if (state2 < nrOfStates){
                                    calculatePartialsPartialsPruningHelperOneQSbelowWithPartialsUnknownState(partials1, matrices1aboveQSstart,
                                                                                                             partials2,
                                                                                                             partials3, w, v, state2);
                                // both children have a gap or unknown state
                                } else {
                                    calculatePartialsPartialsPruningHelperOneQSbelowBothUnknownState(partials1, matrices1aboveQSstart,
                                                                                                     partials2,
                                                                                                     partials3, w, v);
                                }

                                v += nrOfStates;

                            }
                        }

                    }
                    // the QS do not have a state so they have partials
                    else {

                        for (int l = 0; l < nrOfMatrices; l++) {

                            // w keeps track of the state the internal node evolves from
                            int w = l * matrixSize;

                            for (int k = 0; k < nrOfPatterns; k++) {

                                // check if both of the child nodes is a tip
                                if (nodeIndex1 == child1QS && nodeIndex2 == child2QS) {
                                    // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                    calculatePartialsPartialsPruningHelperOneQSbelowBothQSTipsWithPartials(partials1, matrices1,
                                                                                                           matrices1aboveQSstart,
                                                                                                           partials2, matrices2,
                                                                                                           partials3, w, v);
                                }
                                // or just child node 1 is a tip
                                else if (nodeIndex1 == child1QS) {
                                    // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                    calculatePartialsPartialsPruningHelperOneQSbelowOneQSTipWithPartials(partials1, matrices1,
                                                                                                         matrices1aboveQSstart,
                                                                                                         partials2,
                                                                                                         partials3, w, v);
                                }
                                // or just child node 2 is a tip
                                else if (nodeIndex2 == child2QS) {
                                    // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                    calculatePartialsPartialsPruningHelperOneQSbelowOtherQSTipWithPartials(partials1, matrices1aboveQSstart,
                                                                                                           partials2, matrices2,
                                                                                                           partials3, w, v);
                                }
                                // or no child node is a tip
                                else {
                                    // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                    calculatePartialsPartialsPruningHelperOneQSbelowBothUnknownState(partials1, matrices1aboveQSstart,
                                                                                                     partials2,
                                                                                                     partials3, w, v);
                                }

                                v += nrOfStates;

                            }
                        }
                    }
                }
                // child 1 (partials) has the QS start below the branch leading from the parent to the child (i.e. no QS on current branch)
                else {

                    // the QS passing the child 2 below has a state (as opposed to uncertain state = partials)
                    if (stateIndex2 != null) {

                        for (int l = 0; l < nrOfMatrices; l++) {

                            // w keeps track of the state the internal node evolves from
                            int w = l * matrixSize;

                            for (int k = 0; k < nrOfPatterns; k++) {
                                // note down the states at the tips belonging to QS passing through child1/child2
                                int state2 = stateIndex2[k];

                                // child 2 has a state
                                if (state2 < nrOfStates) {
                                    calculatePartialsPartialsPruningHelperOneQSbelowWithPartialsUnknownState(partials1, matrices1, partials2, partials3, w, v, state2);
                                // child 2 has a gap or unknown state
                                } else {
                                    calculatePartialsPartialsPruningHelperOneQSbelowBothUnknownState(partials1, matrices1, partials2, partials3, w, v);
                                }

                                v += nrOfStates;

                            }
                        }
                    }
                    // the QS do not have a state so they have partials
                    else {

                        for (int l = 0; l < nrOfMatrices; l++) {

                            // w keeps track of the state the internal node evolves from
                            int w = l * matrixSize;

                            for (int k = 0; k < nrOfPatterns; k++) {

                                // check if the child node 2 is a tip
                                if (nodeIndex2 == child2QS) {
                                    // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                    calculatePartialsPartialsPruningHelperOneQSbelowOtherQSTipWithPartials(partials1, matrices1,
                                                                                                           partials2, matrices2,
                                                                                                           partials3, w, v);
                                }
                                // or child 2 node is not a tip
                                else {
                                    // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                    calculatePartialsPartialsPruningHelperOneQSbelowBothUnknownState(partials1, matrices1,
                                                                                                     partials2,
                                                                                                     partials3, w, v);
                                }

                                v += nrOfStates;

                            }
                        }
                    }
                }
            // child 2 has QS start on the branch leading from the parent to the child or below
            } else if (child2QS!=parentQS){
                // child2 has the QS start on the branch leading from the parent to the child
                if (child2QS!=-1){

                    // the QS passing the nodes below have a state (as opposed to uncertain state = partials)
                    if (stateIndex1 != null && stateIndex2 != null) {

                        for (int l = 0; l < nrOfMatrices; l++) {

                            // w keeps track of the state the internal node evolves from
                            int w = l * matrixSize;

                            for (int k = 0; k < nrOfPatterns; k++) {
                                // note down the states at the tips belonging to QS passing through child1/child2
                                int state1 = stateIndex1[k];
                                int state2 = stateIndex2[k];

                                // child1's QS and child2's QS have a state
                                if (state1 < nrOfStates && state2 < nrOfStates) {
                                    calculatePartialsPartialsPruningHelperOneQSbelow(partials2, matrices2aboveQSstart,
                                                                                     partials1,
                                                                                     partials3, w, v, state2, state1);
                                // child1's QS has a state but child2's QS has a gap or unknown sequence
                                } else if (state1 < nrOfStates) {
                                    calculatePartialsPartialsPruningHelperOneQSbelowWithPartialsUnknownState(partials2, matrices2aboveQSstart,
                                                                                                             partials1,
                                                                                                             partials3, w, v, state1);
                                // child1's QS has a gap or unknown state but child2's QS has a state
                                } else if (state2 < nrOfStates) {
                                    calculatePartialsPartialsPruningHelperOneQSbelowAndPartialsUnknownState(partials2, matrices2aboveQSstart,
                                                                                                            partials1,
                                                                                                            partials3, w, v, state2);
                                // both children have a gap or unknown state
                                } else {
                                    calculatePartialsPartialsPruningHelperOneQSbelowBothUnknownState(partials2, matrices2aboveQSstart,
                                                                                                     partials1,
                                                                                                     partials3, w, v);
                                }

                                v += nrOfStates;

                            }
                        }
                    }
                    // the QS do not have a state so they have partials
                    else {

                        for (int l = 0; l < nrOfMatrices; l++) {

                            // w keeps track of the state the internal node evolves from
                            int w = l * matrixSize;

                            for (int k = 0; k < nrOfPatterns; k++) {

                                // check if both of the child nodes is a tip
                                if (nodeIndex1 == child1QS && nodeIndex2 == child2QS) {
                                    // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                    calculatePartialsPartialsPruningHelperOneQSbelowBothQSTipsWithPartials(partials2, matrices2,
                                                                                                           matrices2aboveQSstart,
                                                                                                           partials1, matrices1,
                                                                                                           partials3, w, v);
                                }
                                // or just child node 2 is a tip
                                else if (nodeIndex2 == child1QS) {
                                    // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                    calculatePartialsPartialsPruningHelperOneQSbelowOneQSTipWithPartials(partials2, matrices2,
                                                                                                         matrices2aboveQSstart,
                                                                                                         partials1,
                                                                                                         partials3, w, v);
                                }
                                // or just child node 1 is a tip
                                else if (nodeIndex1 == child1QS) {
                                    // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                    calculatePartialsPartialsPruningHelperOneQSbelowOtherQSTipWithPartials(partials2, matrices2aboveQSstart,
                                                                                                           partials1, matrices1,
                                                                                                           partials3, w, v);
                                }
                                // or no child node is a tip
                                else {
                                    // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                    calculatePartialsPartialsPruningHelperOneQSbelowBothUnknownState(partials2, matrices2aboveQSstart,
                                                                                                     partials1,
                                                                                                     partials3, w, v);
                                }

                                v += nrOfStates;

                            }
                        }
                    }
                }
                // child 2 (partials) has the QS start below the branch leading2 from the parent to the child (i.e. no QS on current branch)
                else {

                    // the QS passing the child 2 below has a state (as opposed to uncertain state = partials)
                    if (stateIndex2 != null) {

                        for (int l = 0; l < nrOfMatrices; l++) {

                            // w keeps track of the state the internal node evolves from
                            int w = l * matrixSize;

                            for (int k = 0; k < nrOfPatterns; k++) {
                                // note down the states at the tips belonging to QS passing through child1/child2
                                int state1 = stateIndex1[k];

                                // child 1 has a state
                                if (state1 < nrOfStates){
                                    calculatePartialsPartialsPruningHelperOneQSbelowWithPartialsUnknownState(partials2, matrices2,
                                                                                                             partials1,
                                                                                                             partials3, w, v, state1);
                                // child 1 has a gap or unknown state
                                } else {
                                    calculatePartialsPartialsPruningHelperOneQSbelowBothUnknownState(partials2, matrices2,
                                                                                                     partials1,
                                                                                                     partials3, w, v);
                                }

                                v += nrOfStates;

                            }
                        }
                    }
                    // the QS do not have a state so they have partials
                    else {

                        for (int l = 0; l < nrOfMatrices; l++) {

                            // w keeps track of the state the internal node evolves from
                            int w = l * matrixSize;

                            for (int k = 0; k < nrOfPatterns; k++) {

                                // check if the child node 1 is a tip
                                if (nodeIndex1 == child1QS) {
                                    // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                    calculatePartialsPartialsPruningHelperOneQSbelowOtherQSTipWithPartials(partials2, matrices2,
                                                                                                           partials1, matrices1,
                                                                                                           partials3, w, v);
                                }
                                // or child 1 node is not a tip
                                else {
                                    // note that the tips belonging to QS passing through child1/child2 have partials (no states) or gap (does not matter)
                                    calculatePartialsPartialsPruningHelperOneQSbelowBothUnknownState(partials2, matrices2,
                                                                                                     partials1,
                                                                                                     partials3, w, v);
                                }

                                v += nrOfStates;

                            }
                        }
                    }
                }
            }
            else{
                throw new IllegalStateException("What case did we miss in QuasiSpeciesBeerLikelihoodCore function calculateStatesPartialsPruning");
            }
        }
    }

    /**
     * Calculates partial likelihoods at a node.
     *
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     * @param child1QS   QS passing through child 1 - if none this is -1
     * @param child2QS   QS passing through child 2 - if none this is -1
     * @param parentQS   QS passing through the parent node - if none this is -1
     * @param nodeCount  total count of the true nodes in the tree
     */
    public void calculateQSPartials(int nodeIndex1, int nodeIndex2, int nodeIndex3, int child1QS, int child2QS, int parentQS, int nodeCount) {
        if (states[nodeIndex1] != null) {
            if (states[nodeIndex2] != null) {
                calculateStatesStatesPruning(
                        states[nodeIndex1],matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],matrices[currentMatrixIndex[nodeCount+nodeIndex1]][nodeCount+nodeIndex1],
                        states[nodeIndex2],matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],matrices[currentMatrixIndex[nodeCount+nodeIndex2]][nodeCount+nodeIndex2],
                        partials[currentPartialsIndex[nodeIndex3]][nodeIndex3],child1QS,child2QS,parentQS);
            } else {
                if (child2QS == -1){
                    calculateStatesPartialsPruning(
                        states[nodeIndex1],matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],matrices[currentMatrixIndex[nodeCount+nodeIndex1]][nodeCount+nodeIndex1],
                        null,partials[currentPartialsIndex[nodeIndex2]][nodeIndex2],matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],null,
                        partials[currentPartialsIndex[nodeIndex3]][nodeIndex3],child1QS,child2QS,parentQS);
                }
                else {
                    calculateStatesPartialsPruning(
                        states[nodeIndex1],matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],matrices[currentMatrixIndex[nodeCount+nodeIndex1]][nodeCount+nodeIndex1],
                        states[child2QS],partials[currentPartialsIndex[nodeIndex2]][nodeIndex2],matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],matrices[currentMatrixIndex[nodeCount+child2QS]][nodeCount+child2QS],
                        partials[currentPartialsIndex[nodeIndex3]][nodeIndex3],child1QS,child2QS,parentQS);
                }
            }
        } else {
            if (states[nodeIndex2] != null) {
                if (child1QS == -1){
                    calculateStatesPartialsPruning(
                        states[nodeIndex2],matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],matrices[currentMatrixIndex[nodeCount+nodeIndex2]][nodeCount+nodeIndex2],
                        null,partials[currentPartialsIndex[nodeIndex1]][nodeIndex1],matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],null,
                        partials[currentPartialsIndex[nodeIndex3]][nodeIndex3],child2QS,child1QS,parentQS);
                }
                else{
                    calculateStatesPartialsPruning(
                        states[nodeIndex2],matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],matrices[currentMatrixIndex[nodeCount+nodeIndex2]][nodeCount+nodeIndex2],
                        states[child1QS],partials[currentPartialsIndex[nodeIndex1]][nodeIndex1],matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],matrices[currentMatrixIndex[nodeCount+child1QS]][nodeCount+child1QS],
                        partials[currentPartialsIndex[nodeIndex3]][nodeIndex3],child2QS,child1QS,parentQS);
                }
            } else {
                if (child1QS == -1 && child2QS == -1){
                    calculatePartialsPartialsPruning(
                        null,partials[currentPartialsIndex[nodeIndex1]][nodeIndex1],matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],null,
                        null,partials[currentPartialsIndex[nodeIndex2]][nodeIndex2],matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],null,
                        partials[currentPartialsIndex[nodeIndex3]][nodeIndex3],child1QS,child2QS,parentQS,nodeIndex1,nodeIndex2);
                }
                else if (child1QS == -1){
                    calculatePartialsPartialsPruning(
                        null,partials[currentPartialsIndex[nodeIndex1]][nodeIndex1],matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],null,
                        states[child2QS],partials[currentPartialsIndex[nodeIndex2]][nodeIndex2],matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],matrices[currentMatrixIndex[nodeCount+child2QS]][nodeCount+child2QS],
                        partials[currentPartialsIndex[nodeIndex3]][nodeIndex3],child1QS,child2QS,parentQS,nodeIndex1,nodeIndex2);
                }
                else if (child2QS == -1){
                    calculatePartialsPartialsPruning(
                        states[child1QS],partials[currentPartialsIndex[nodeIndex1]][nodeIndex1],matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],matrices[currentMatrixIndex[nodeCount+child1QS]][nodeCount+child1QS],
                        null,partials[currentPartialsIndex[nodeIndex2]][nodeIndex2],matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],null,
                        partials[currentPartialsIndex[nodeIndex3]][nodeIndex3],child1QS,child2QS,parentQS,nodeIndex1,nodeIndex2);
                }
                else{
                    calculatePartialsPartialsPruning(
                        states[child1QS],partials[currentPartialsIndex[nodeIndex1]][nodeIndex1],matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1],matrices[currentMatrixIndex[nodeCount+child1QS]][nodeCount+child1QS],
                        states[child2QS],partials[currentPartialsIndex[nodeIndex2]][nodeIndex2],matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2],matrices[currentMatrixIndex[nodeCount+child2QS]][nodeCount+child2QS],
                        partials[currentPartialsIndex[nodeIndex3]][nodeIndex3],child1QS,child2QS,parentQS,nodeIndex1,nodeIndex2);
                }
            }
        }

        if (useScaling) {
            scalePartials(nodeIndex3);
        }
    }

} // class QuasiSpeciesBeerLikelihoodCore
