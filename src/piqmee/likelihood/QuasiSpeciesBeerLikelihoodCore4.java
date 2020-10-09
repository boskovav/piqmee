package piqmee.likelihood;

public class QuasiSpeciesBeerLikelihoodCore4 extends QuasiSpeciesBeerLikelihoodCore {

	public QuasiSpeciesBeerLikelihoodCore4(int nrOfStates) {
		super(nrOfStates);
		if (nrOfStates !=4) {
			throw new IllegalArgumentException("This only works with 4 states");
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

        for (int i = 0; i < 4; i++) {
            // since state at root is unknown, take into account all the possibilities
            // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
            sum = 0.0;
//            for (int j = 0; j < nrOfStates; j++) {
            	final int s = w + 4*i;
                // note down the partial at the root
                tmp = partialsRoot[v];
                sum += tmp * matricesRootOrRootaboveQSstart[s];
                tmp = partialsRoot[v + 1];
                sum += tmp * matricesRootOrRootaboveQSstart[s + 1];
                tmp = partialsRoot[v + 2];
                sum += tmp * matricesRootOrRootaboveQSstart[s + 2];
                tmp = partialsRoot[v + 3];
                sum += tmp * matricesRootOrRootaboveQSstart[s + 3];
                
//            }
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
        final int s = w + rootState;
//        for (int i = 0; i < nrOfStates; i++){
            originPartials[v    ] = tmp * matricesRootaboveQSstart[s];
            originPartials[v + 1] = tmp * matricesRootaboveQSstart[s + 4];
            originPartials[v + 2] = tmp * matricesRootaboveQSstart[s + 4 * 2];
            originPartials[v + 3] = tmp * matricesRootaboveQSstart[s + 4 * 3];
//        }
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
        tmp1 = matricesQS1[w + 4 * state1 + state1];
        tmp2 = matricesQS2[w + 4 * state2 + state2];
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        final int s1 = w  + state1;
        final int s2 = w  + state2;
        //for (int i = 0; i < nrOfStates; i++) {
            partials3[v]     = tmp1 * matrices1aboveQSstart[s1]
                             * tmp2 * matrices2aboveQSstart[s2];
            partials3[v + 1] = tmp1 * matrices1aboveQSstart[s1 + 4]
                    * tmp2 * matrices2aboveQSstart[s2 + 4];
            partials3[v + 2] = tmp1 * matrices1aboveQSstart[s1 + 4 * 2]
                    * tmp2 * matrices2aboveQSstart[s2 + 4 * 2];
            partials3[v + 3] = tmp1 * matrices1aboveQSstart[s1 + 4 * 3]
                    * tmp2 * matrices2aboveQSstart[s2 + 4 * 3];
            
        //}
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
        tmp1 = matricesQS1[w + 4 * state1 + state1];
        // since state at tip 2 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < 4; i++){
            sum2 = 0.0;
            final int s = w + 4 * i;
            //for (int j = 0; j < nrOfStates; j++) {
                tmp2 = matricesQS2[w];
                sum2 += tmp2 * matrices2aboveQSstart[s];

                tmp2 = matricesQS2[w + 4 + 1];
                sum2 += tmp2 * matrices2aboveQSstart[s + 1];
                
                tmp2 = matricesQS2[w + 4 * 2 + 2];
                sum2 += tmp2 * matrices2aboveQSstart[s + 2];
                
                tmp2 = matricesQS2[w + 4 * 3 + 3];
                sum2 += tmp2 * matrices2aboveQSstart[s + 3];
                
            //}
            partials3[v + i] = tmp1 * matrices1aboveQSstart[w + 4 * i + state1]
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
        for (int i = 0; i < 4; i++){
            sum1 = sum2 = 0.0;
            final int s = w + 4 * i;
            //for (int j = 0; j < nrOfStates; j++) {
                tmp1 = matricesQS1[w];
                tmp2 = matricesQS2[w];
                sum1 += tmp1 * matrices1aboveQSstart[s];
                sum2 += tmp2 * matrices2aboveQSstart[s];
                
                tmp1 = matricesQS1[w + 4 * 1 + 1];
                tmp2 = matricesQS2[w + 4 * 1 + 1];
                sum1 += tmp1 * matrices1aboveQSstart[s + 1];
                sum2 += tmp2 * matrices2aboveQSstart[s + 1];

                tmp1 = matricesQS1[w + 4 * 2 + 2];
                tmp2 = matricesQS2[w + 4 * 2 + 2];
                sum1 += tmp1 * matrices1aboveQSstart[s + 2];
                sum2 += tmp2 * matrices2aboveQSstart[s + 2];

                tmp1 = matricesQS1[w + 4 * 3 + 3];
                tmp2 = matricesQS2[w + 4 * 3 + 3];
                sum1 += tmp1 * matrices1aboveQSstart[s + 3];
                sum2 += tmp2 * matrices2aboveQSstart[s + 3];
                //}
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
        tmp1 = matricesQS1[w + 4 * state1 + state1];
        tmp2 = matricesQS2[w + 4 * state2 + state2];
        // for the state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        // 0 elsewhere (if not the state at the parent)
//        for (int i = 0; i < nrOfStates; i++) {
            partials3[v] = 0;
            partials3[v + 1] = 0;
            partials3[v + 2] = 0;
            partials3[v + 3] = 0;
//        }
        partials3[v + state2] = tmp1 * matrices1aboveQSstart[w + 4 * state2 + state1]
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
        tmp1 = matricesQS1[w + 4 * state1 + state1];
        // since state at tip 2 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        final int s = w + state1;
        //for (int i = 0; i < nrOfStates; i++) {
            tmp2 = matricesQS2[w];
            partials3[v    ] = tmp1 * matrices1aboveQSstart[s]
                             * tmp2;
            
            tmp2 = matricesQS2[w + 4 * 1 + 1];
            partials3[v + 1] = tmp1 * matrices1aboveQSstart[s + 4 * 1]
                             * tmp2;
            
            tmp2 = matricesQS2[w + 4 * 2 + 2];
            partials3[v + 2] = tmp1 * matrices1aboveQSstart[s + 4 * 2]
                             * tmp2;
            
            tmp2 = matricesQS2[w + 4 * 3 + 3];
            partials3[v + 3] = tmp1 * matrices1aboveQSstart[s + 4 * 3]
                             * tmp2;
            
            
        //}
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
        tmp2 = matricesQS2[w + 4 * state2 + state2];
        // since state at tip 1 is unknown, take into account all the possibilities
        // for the state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        sum1 = 0.0;
        final int s = w + 4 * state2;
//        for (int i = 0; i < nrOfStates; i++) {
            tmp1 = matricesQS1[w];
            sum1 += tmp1 * matrices1aboveQSstart[s];
            partials3[v] = 0;

            tmp1 = matricesQS1[w + 4 * 1 + 1];
            sum1 += tmp1 * matrices1aboveQSstart[s + 1];
            partials3[v + 1] = 0;

            tmp1 = matricesQS1[w + 4 * 2 + 2];
            sum1 += tmp1 * matrices1aboveQSstart[s + 2];
            partials3[v + 2] = 0;

            tmp1 = matricesQS1[w + 4 * 3 + 3];
            sum1 += tmp1 * matrices1aboveQSstart[s + 3];
            partials3[v + 3] = 0;

//        }
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
        for (int i = 0; i < 4; i++){
            final int s = w + 4 * i;
            tmp2 = matricesQS2[s + i];
            sum1 = 0.0;
//            for (int j = 0; j < nrOfStates; j++) {
                tmp1 = matricesQS1[w];
                sum1 += tmp1 * matrices1aboveQSstart[s];
                
                tmp1 = matricesQS1[w + 4 * 1 + 1];
                sum1 += tmp1 * matrices1aboveQSstart[s + 1];
                
                tmp1 = matricesQS1[w + 4 * 2 + 2];
                sum1 += tmp1 * matrices1aboveQSstart[s + 2];
                
                tmp1 = matricesQS1[w + 4 * 3 + 3];
                sum1 += tmp1 * matrices1aboveQSstart[s + 3];
                
                
//            }
            partials3[v + i] = sum1 * tmp2;
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
        tmp1 = matricesQS1[w + 4 * state1 + state1];
        // note down the partial at the child 2
        tmp2 = partials2[v + state2];
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        //for (int i = 0; i < nrOfStates; i++) {
        final int s1 = w + state1;
        final int s2 = w + state2;
        
            partials3[v    ] = tmp1 * matrices1aboveQSstart[s1]
                             * tmp2 * matrices2aboveQSstart[s2];

            partials3[v + 1] = tmp1 * matrices1aboveQSstart[s1 + 4]
                    * tmp2 * matrices2aboveQSstart[s2 + 4];
            
            partials3[v + 2] = tmp1 * matrices1aboveQSstart[s1 + 4 * 2]
                    * tmp2 * matrices2aboveQSstart[s2 + 4 * 2];
            
            partials3[v + 3] = tmp1 * matrices1aboveQSstart[s1 + 4 * 3]
                    * tmp2 * matrices2aboveQSstart[s2 + 4 * 3];
        //}
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
        tmp1 = matricesQS1[w + 4 * state1 + state1];
        // since state at child 2 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        for (int i = 0; i < 4; i++){
            sum2 = 0.0;
            final int s = w + 4 * i;
//            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 2
                tmp2 = partials2[v];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s];
                tmp2 = partials2[v + 1];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 1];
                tmp2 = partials2[v + 2];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 2];
                tmp2 = partials2[v + 3];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 3];
                
//            }
            partials3[v + i] = tmp1 * matrices1aboveQSstart[s + state1]
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
        for (int i = 0; i < 4; i++){
            sum1 = 0.0;
            final int s = w + 4 * i;
//            for (int j = 0; j < nrOfStates; j++) {
                // note down the transition probability (of no change) on the sum of the QS branch lengths
                // P(QS start -> QS tip)
                tmp1 = matricesQS1[w];
                sum1 += tmp1 * matrices1aboveQSstart[s];
                tmp1 = matricesQS1[w + 4 + 1];
                sum1 += tmp1 * matrices1aboveQSstart[s + 1];
                tmp1 = matricesQS1[w + 4 * 2 + 2];
                sum1 += tmp1 * matrices1aboveQSstart[s + 2];
                tmp1 = matricesQS1[w + 4 * 3 + 3];
                sum1 += tmp1 * matrices1aboveQSstart[s + 3];
//            }
            partials3[v + i] = sum1
                             * tmp2 * matrices2aboveQSstart[s + state2];
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
        for (int i = 0; i < 4; i++){
            sum1 = sum2 = 0.0;
            //for (int j = 0; j < nrOfStates; j++) {
                // note down the transition probability (of no change) on the sum of the QS branch lengths
                // P(QS start -> QS tip)
            	final int s = w + 4 * i;
            	
                tmp1 = matricesQS1[w];
                // note down the partial at the child 2
                tmp2 = partials2[v];
                sum1 += tmp1 * matrices1aboveQSstart[s];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s];
 
                tmp1 = matricesQS1[w + 4 + 1];
                tmp2 = partials2[v + 1];
                sum1 += tmp1 * matrices1aboveQSstart[s + 1];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 1];

                tmp1 = matricesQS1[w + 4 * 2 + 2];
                tmp2 = partials2[v + 2];
                sum1 += tmp1 * matrices1aboveQSstart[s + 2];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 2];

                tmp1 = matricesQS1[w + 4 * 3 + 3];
                tmp2 = partials2[v + 3];
                sum1 += tmp1 * matrices1aboveQSstart[s + 3];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 3];
            //}
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
        tmp1 = matricesQS1[w + 4 * state1 + state1];
        // note down the partial at the child 2
        tmp2 = partials2[v + state2];
        // for the state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        // 0 elsewhere (if not the state at the parent)
//        for (int i = 0; i < nrOfStates; i++) {
            partials3[v] = 0;
            partials3[v + 1] = 0;
            partials3[v + 2] = 0;
            partials3[v + 3] = 0;
//        }
        partials3[v + state2] = tmp1 * matrices1aboveQSstart[w + 4 * state2 + state1]
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
        tmp1 = matricesQS1[w + 4 * state1 + state1];
        // since state at child 2 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        //for (int i = 0; i < nrOfStates; i++){
            // note down the partial at the child 2
        final int s = w + state1;
            tmp2 = partials2[v];
            partials3[v    ] = tmp1 * matrices1aboveQSstart[s]
                             * tmp2;
            tmp2 = partials2[v + 1];
            partials3[v + 1] = tmp1 * matrices1aboveQSstart[s + 4 * 1]
                             * tmp2;
            tmp2 = partials2[v + 2];
            partials3[v + 2] = tmp1 * matrices1aboveQSstart[s + 4 * 2]
                             * tmp2;
            tmp2 = partials2[v + 3];
            partials3[v + 3] = tmp1 * matrices1aboveQSstart[s + 4 * 3]
                             * tmp2;
        //}
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
//        for (int i = 0; i < nrOfStates; i++){
            partials3[v] = 0;
            // note down the transition probability (of no change) on the sum of the QS branch lengths
            // P(QS start -> QS tip)
            final int s = w + 4 * state2;
            tmp1 = matricesQS1[w];
            sum1 += tmp1 * matrices1aboveQSstart[s];
            tmp1 = matricesQS1[w + 4 + 1];
            sum1 += tmp1 * matrices1aboveQSstart[s + 1];
            tmp1 = matricesQS1[w + 4 * 2 + 2];
            sum1 += tmp1 * matrices1aboveQSstart[s + 2];
            tmp1 = matricesQS1[w + 4 * 3 + 3];
            sum1 += tmp1 * matrices1aboveQSstart[s + 3];
//        }
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
        for (int i = 0; i < 4; i++){
            // note down the partial at the child 2
            tmp2 = partials2[v + i];
            sum1 = 0.0;
//            for (int j = 0; j < nrOfStates; j++) {
                // note down the transition probability (of no change) on the sum of the QS branch lengths
                // P(QS start -> QS tip)
            final int s = w + 4 * i;
                tmp1 = matricesQS1[w];
                sum1 += tmp1 * matrices1aboveQSstart[s];
                tmp1 = matricesQS1[w + 4 + 1];
                sum1 += tmp1 * matrices1aboveQSstart[s + 1];
                tmp1 = matricesQS1[w + 4 * 2 + 2];
                sum1 += tmp1 * matrices1aboveQSstart[s + 2];
                tmp1 = matricesQS1[w + 4 * 3 + 3];
                sum1 += tmp1 * matrices1aboveQSstart[s + 3];
//            }
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
        tmp1 = matricesQS1[w + 4 * state1 + state1];
        // note down the partial at the child 2
        tmp2 = partials2[v + state2];
        // for the state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        // 0 elsewhere (if not the state at the parent)
//        for (int i = 0; i < nrOfStates; i++) {
            partials3[v] = 0;
            partials3[v + 1] = 0;
            partials3[v + 2] = 0;
            partials3[v + 3] = 0;
//        }
        partials3[v + state1] = tmp1
                              * tmp2 * matrices2aboveQSstart[w + 4 * state1 + state2];

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
        tmp1 = matricesQS1[w + 4 * state1 + state1];
        // since state at child 2 is unknown, take into account all the possibilities
        // for each state at the parent node calculate prob. of going to the state of the tip * P(QS start -> QS tip)
        sum2 = 0.0;
        final int s = w + 4 * state1;
        //for (int i = 0; i < nrOfStates; i++){
            partials3[v] = 0;
            // note down the partial at the child 2
            tmp2 = partials2[v];
            sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s];
            partials3[v + 1] = 0;
            // note down the partial at the child 2
            tmp2 = partials2[v + 1];
            sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 1];
            partials3[v + 2] = 0;
            // note down the partial at the child 2
            tmp2 = partials2[v + 2];
            sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 2];
            partials3[v + 3] = 0;
            // note down the partial at the child 2
            tmp2 = partials2[v + 3];
            sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 3];
        
        //}
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
        //for (int i = 0; i < nrOfStates; i++){
            // note down the transition probability (of no change) on the sum of the QS branch lengths
            // P(QS start -> QS tip)
        final int s = w + state2;
            tmp1 = matricesQS1[w];
            partials3[v] = tmp1
                             * tmp2 * matrices2aboveQSstart[s];
            tmp1 = matricesQS1[w + 4 + 1];
            partials3[v + 1] = tmp1
                             * tmp2 * matrices2aboveQSstart[s + 4];
            tmp1 = matricesQS1[w + 4*2 + 2];
            partials3[v + 2] = tmp1
                             * tmp2 * matrices2aboveQSstart[s + 4 * 2];
            tmp1 = matricesQS1[w + 4*3 + 3];
            partials3[v + 3] = tmp1
                             * tmp2 * matrices2aboveQSstart[s + 4 * 3];
        //}
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
        for (int i = 0; i < 4; i++){
            // note down the transition probability (of no change) on the sum of the QS branch lengths
            // P(QS start -> QS tip)
        	final int s = w + 4 * i;
            tmp1 = matricesQS1[s + i];
            sum2 = 0.0;
 //           for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 2
                tmp2 = partials2[v];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s];
                tmp2 = partials2[v + 1];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 1];
                tmp2 = partials2[v + 2];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 2];
                tmp2 = partials2[v + 3];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 3];
 //           }
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
//        for (int i = 0; i < nrOfStates; i++) {
  
            partials3[v    ] = tmp1 * matrices1aboveQSstart[w  + state1]
                    * tmp2 * matrices2aboveQSstart[w + state2];
            partials3[v + 1] = tmp1 * matrices1aboveQSstart[w + 4 * 1 + state1]
                    * tmp2 * matrices2aboveQSstart[w + 4 * 1 + state2];
            partials3[v + 2] = tmp1 * matrices1aboveQSstart[w + 4 * 2 + state1]
                    * tmp2 * matrices2aboveQSstart[w + 4 * 2 + state2];
            partials3[v + 3] = tmp1 * matrices1aboveQSstart[w + 4 * 3 + state1]
                    * tmp2 * matrices2aboveQSstart[w + 4 * 3 + state2];
//        }
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
        for (int i = 0; i < 4; i++){
            sum2 = 0.0;
            final int s = w + 4 * i;

//            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 2
                tmp2 = partials2[v];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s];
                tmp2 = partials2[v + 1];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 1];
                tmp2 = partials2[v + 2];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 2];
                tmp2 = partials2[v + 3];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 3];
//            }
            partials3[v + i] = tmp1 * matrices1aboveQSstart[s + state1]
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
        for (int i = 0; i < 4; i++){
            sum1 = sum2 = 0.0;
//            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 1
            final int s = w + 4 * i;

            	tmp1 = partials1[v];
                // note down the partial at the child 2
                tmp2 = partials2[v];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s];

            	tmp1 = partials1[v + 1];
                // note down the partial at the child 2
                tmp2 = partials2[v + 1];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s + 1];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 1];

            	tmp1 = partials1[v + 2];
                // note down the partial at the child 2
                tmp2 = partials2[v + 2];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s + 2];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 2];

            	tmp1 = partials1[v + 3];
                // note down the partial at the child 2
                tmp2 = partials2[v + 3];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s + 3];
                sum2 += tmp2 * matrices2aboveQSstartOrMatrices2[s + 3];

                //            }
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
        for (int i = 0; i < 4; i++){
            sum1 = sum2 = 0.0;
            //for (int j = 0; j < nrOfStates; j++) {
            final int s = w + 4 * i;

            
                // note down the partial at the child 1
                tmp1 = partials1[v    ] * matricesQS1[w];
                // note down the partial at the child 2
                tmp2 = partials2[v    ] * matricesQS2[w];
                sum1 += tmp1 * matrices1aboveQSstart[s];
                sum2 += tmp2 * matrices2aboveQSstart[s];

                // note down the partial at the child 1
                tmp1 = partials1[v + 1] * matricesQS1[w + 4 * 1 + 1];
                // note down the partial at the child 2
                tmp2 = partials2[v + 1] * matricesQS2[w + 4 * 1 + 1];
                sum1 += tmp1 * matrices1aboveQSstart[s + 1];
                sum2 += tmp2 * matrices2aboveQSstart[s + 1];

                //                // note down the partial at the child 1
                tmp1 = partials1[v + 2] * matricesQS1[w + 4 * 2 + 2];
                // note down the partial at the child 2
                tmp2 = partials2[v + 2] * matricesQS2[w + 4 * 2 + 2];
                sum1 += tmp1 * matrices1aboveQSstart[s + 2];
                sum2 += tmp2 * matrices2aboveQSstart[s + 2];

                // note down the partial at the child 1
                tmp1 = partials1[v + 3] * matricesQS1[w + 4 * 3 + 3];
                // note down the partial at the child 2
                tmp2 = partials2[v + 3] * matricesQS2[w + 4 * 3 + 3];
                sum1 += tmp1 * matrices1aboveQSstart[s + 3];
                sum2 += tmp2 * matrices2aboveQSstart[s + 3];
        	
        
         //}
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
        for (int i = 0; i < 4; i++){
            sum1 = sum2 = 0.0;
            //for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 1
            final int s = w + 4 * i;
                tmp1 = partials1[v    ] * matricesQS1[w];
                // note down the partial at the child 2
                tmp2 = partials2[v    ];
                sum1 += tmp1 * matrices1aboveQSstart[s];
                sum2 += tmp2 * matrices2aboveQSstart[s];
                
                tmp1 = partials1[v + 1] * matricesQS1[w + 4 + 1];
                // note down the partial at the child 2
                tmp2 = partials2[v + 1];
                sum1 += tmp1 * matrices1aboveQSstart[s + 1];
                sum2 += tmp2 * matrices2aboveQSstart[s + 1];
                
                tmp1 = partials1[v + 2] * matricesQS1[w + 4 * 2 + 2];
                // note down the partial at the child 2
                tmp2 = partials2[v + 2];
                sum1 += tmp1 * matrices1aboveQSstart[s + 2];
                sum2 += tmp2 * matrices2aboveQSstart[s + 2];
                
                tmp1 = partials1[v + 3] * matricesQS1[w + 4 * 3 + 3];
                // note down the partial at the child 2
                tmp2 = partials2[v + 3];
                sum1 += tmp1 * matrices1aboveQSstart[s + 3];
                sum2 += tmp2 * matrices2aboveQSstart[s + 3];
            // }
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
//        for (int i = 0; i < nrOfStates; i++) {
            partials3[v] = 0;
            partials3[v+1] = 0;
            partials3[v+2] = 0;
            partials3[v+3] = 0;
//        }
        partials3[v + state2] = tmp1 * matrices1aboveQSstart[w + 4 * state2 + state1]
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
//        for (int i = 0; i < nrOfStates; i++){
            // note down the partial at the child 2
            w += state1;
            tmp2 = partials2[v];
            partials3[v    ] = tmp1 * matrices1aboveQSstart[w]
                    * tmp2;
            tmp2 = partials2[v + 1];
            w += 4;
            partials3[v + 1] = tmp1 * matrices1aboveQSstart[w]
                    * tmp2;
            tmp2 = partials2[v + 2];
            w += 4;
            partials3[v + 2] = tmp1 * matrices1aboveQSstart[w]
                    * tmp2;
            tmp2 = partials2[v + 3];
            w += 4;
            partials3[v + 3] = tmp1 * matrices1aboveQSstart[w]
                    * tmp2;
//        }
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
        final int s = w + 4 * state2;
//        for (int i = 0; i < nrOfStates; i++){
            partials3[v ] = 0;
            // note down the partial at the child 1
            tmp1 = partials1[v];
            sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s];
            partials3[v + 1] = 0;
            // note down the partial at the child 1
            tmp1 = partials1[v + 1];
            sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s + 1];
            partials3[v + 2] = 0;
            // note down the partial at the child 1
            tmp1 = partials1[v + 2];
            sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s + 2];
            partials3[v + 3] = 0;
            // note down the partial at the child 1
            tmp1 = partials1[v + 3];
            sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s + 3];
//        }
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
        for (int i = 0; i < 4; i++){
            // note down the partial at the child 2
            tmp2 = partials2[v + i];
            sum1 = 0.0;
            final int s = w + 4 * i;
//            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 1
                tmp1 = partials1[v];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s];
                tmp1 = partials1[v + 1];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s + 1];
                tmp1 = partials1[v + 2];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s + 2];
                tmp1 = partials1[v + 3];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s + 3];
//            }
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
        for (int i = 0; i < 4; i++){
            // note down the partial at the child 2
            tmp2 = partials2[v + i] * matricesQS2[w + 4 * i + i];
            sum1 = 0.0;
            final int s = 4 * i;
//            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 1
                tmp1 = partials1[v    ] * matricesQS1[w];
                sum1 += tmp1 * matrices1aboveQSstart[w + s];
                tmp1 = partials1[v + 1] * matricesQS1[w + 4 + 1];
                sum1 += tmp1 * matrices1aboveQSstart[w + s + 1];
                tmp1 = partials1[v + 2] * matricesQS1[w + 4 * 2 + 2];
                sum1 += tmp1 * matrices1aboveQSstart[w + s + 2];
                tmp1 = partials1[v + 3] * matricesQS1[w + 4 * 3 + 3];
                sum1 += tmp1 * matrices1aboveQSstart[w + s + 3];
//            }
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
        for (int i = 0; i < 4; i++){
            // note down the partial at the child 2
            tmp2 = partials2[v + i];
            sum1 = 0.0;
            final int s = 4 * i;
//            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 1
                tmp1 = partials1[v    ] * matricesQS1[w];
                sum1 += tmp1 * matrices1aboveQSstart[w + s];
                tmp1 = partials1[v + 1] * matricesQS1[w + 4 + 1];
                sum1 += tmp1 * matrices1aboveQSstart[w + s + 1];
                tmp1 = partials1[v + 2] * matricesQS1[w + 4 * 2 + 2];
                sum1 += tmp1 * matrices1aboveQSstart[w + s + 2];
                tmp1 = partials1[v + 3] * matricesQS1[w + 4 * 3 + 3];
                sum1 += tmp1 * matrices1aboveQSstart[w + s + 3];
//            }
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
        for (int i = 0; i < 4; i++){
            // note down the partial at the child 2
            tmp2 = partials2[v + i] * matricesQS2[w + 4 * i + i];
            final int s = w + 4 * i;
            sum1 = 0.0;
//            for (int j = 0; j < nrOfStates; j++) {
                // note down the partial at the child 1
                tmp1 = partials1[v];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s];
                tmp1 = partials1[v + 1];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s + 1];
                tmp1 = partials1[v + 2];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s + 2];
                tmp1 = partials1[v + 3];
                sum1 += tmp1 * matrices1aboveQSstartOrMatrices1[s + 3];
//            }
            partials3[v + i] = sum1 * tmp2;
        }
    }

	
	
	
}
