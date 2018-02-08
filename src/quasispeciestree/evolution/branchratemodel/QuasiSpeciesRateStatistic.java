package quasispeciestree.evolution.branchratemodel;

import beast.core.Description;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.RateStatistic;
import beast.evolution.tree.Node;
import beast.math.statistic.DiscreteStatistics;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;

/**
 *  @author Veronika Boskova created on 08/02/2018
 */
@Description("Makes sure the summary statistic is calculated properly also for QS trees ")
public class QuasiSpeciesRateStatistic extends RateStatistic {

    protected QuasiSpeciesTree tree = null;
    protected BranchRateModel branchRateModel = null;
    protected boolean internal = true;
    protected boolean external = true;

    final static int MEAN = 0;
    final static int VARIANCE = 1;
    final static int COEFFICIENT_OF_VARIATION = 2;

    // Empty constructor as required:
    public QuasiSpeciesRateStatistic() { };

    /**
     * calculate the three statistics from scratch *
     */
    @Override
    public double[] calcValues() {
        // literal copy from RateStatistic besides lines TODO & TODO //
        int length = 0;
        int offset = 0;
        int offsetForBranchLengthsAndRates = 0;

        final int nrOfLeafs = tree.getLeafNodeCount();
        final int nrOfQSDupl = tree.getTotalAttachmentCounts();
        final int nrOfLeafsAndQSDupl = nrOfQSDupl + nrOfLeafs;

        if (external) {
            length += nrOfLeafsAndQSDupl;
        }
        if (internal) {
            length += tree.getInternalNodeCount() - 1 + nrOfQSDupl;
        }

        // these are the rates for each branch as in the full tree when reconstructed from the QS tree
        final double[] rates = new double[length];
        // these are the rates for each branch length, where branch length for the tip is the total
        //  branch length spanned by the haplotype
        final double[] ratesForBranchLengths = new double[length];
        final double[] internalQSBranchLengths = new double[nrOfLeafs];
        final double[] internalQSBranchRates = new double[nrOfLeafs];
        // need those only for mean
        final double[] branchLengths = new double[length];

        final Node[] nodes = tree.getNodesAsArray();

        /** handle leaf nodes **/
        if (external) {
            int k = 0;
            for (int i = 0; i < nrOfLeafs; i++) {
                final QuasiSpeciesNode child = (QuasiSpeciesNode) nodes[i];
                double rate = branchRateModel.getRateForBranch(child);
                final Node parent = child.getParent();

                int haploAboveChild = child.getHaploAboveName();
                double[] attachTimes = child.getAttachmentTimesList();
                int attachTimeArrayLength = attachTimes.length;
                // count here only the QS subtrees without the partial branch above it!!!
                // it can be that there are no duplicates, so the branch is from the child to the parent
                if (haploAboveChild != -1 && attachTimeArrayLength == 1) {
                    branchLengths[i] = parent.getHeight() - child.getHeight();
                    rates[i] = rate;
                }
                // and now for each "external branch" from QS duplicates
                else {
                    double[] tipTimes = child.getTipTimesList();
                    branchLengths[i] = child.getTotalBranchLengths() - attachTimes[0] - tipTimes[0];

                    // the first external branch can be till the real internal node or till the next attachment time
                    if (parent.getHeight() > attachTimes[attachTimeArrayLength - 1]) {
                        branchLengths[i] += attachTimes[attachTimeArrayLength - 1] - tipTimes[0];
                        internalQSBranchLengths[i] = attachTimes[0] - attachTimes[attachTimeArrayLength - 1];
                    }
                    else {
                        branchLengths[i] += parent.getHeight() - tipTimes[0];
                        internalQSBranchLengths[i] = attachTimes[0] - parent.getHeight();
                    }
                    internalQSBranchRates[i] = rate;
                    for (int q = 0; q < attachTimeArrayLength; q++) {
                        rates[k + q] = rate;
                    }
                    k += attachTimeArrayLength;
                }
                ratesForBranchLengths[i] = rate;
            }
            offset = nrOfLeafsAndQSDupl;
            offsetForBranchLengthsAndRates = nrOfLeafs;
        }

        /** handle internal nodes **/
        if (internal) {
            final int n = tree.getNodeCount();
            int k = offset;
            int l = offsetForBranchLengthsAndRates;
            // contribution from the real internal nodes only
            for (int i = nrOfLeafs; i < n; i++) {
                final QuasiSpeciesNode child = (QuasiSpeciesNode) nodes[i];
                //orig-root branch has rate 1 --- defined in the relaxed clock models
                if (!child.isRoot()) {
                    double parentHeight = child.getParent().getHeight();
                    int haploAboveChild = child.getHaploAboveName();
                    int haploContChild = child.getContinuingHaploName();
                    // Case 1: branch can be free - i.e. no haplo passing, branch has its rate
                    if (haploContChild == -1) {
                        branchLengths[k] = parentHeight - child.getHeight();
                    }
                    // Case 2: it can be that a haplo starts above this node, in that case
                    //          only the partial branch has the branch's rate
                    else if (haploAboveChild != -1) {
                        branchLengths[k] = parentHeight - ((QuasiSpeciesNode) tree.getNode(haploAboveChild)).getAttachmentTimesList()[0];
                    }
                    rates[k] = branchRateModel.getRateForBranch(child);
                    k++;
                    // Case 3: it can be that a haplo passes this branch, in that case the rate is the same as for the
                    //          corresponding tip branch rate -- treated separately below
                }
            }
            // contribution from the "internal nodes" of the QS subtree only
            for (int i = 0; i < nrOfLeafs; i++){
                final QuasiSpeciesNode child = (QuasiSpeciesNode) nodes[i];
                Node parent = child.getParent();

                int haploAboveChild = child.getHaploAboveName();
                double[] attachTimes = child.getAttachmentTimesList();
                int attachTimeArrayLength = attachTimes.length;

                double rate = branchRateModel.getRateForBranch(child);
                if (attachTimeArrayLength > 1) {
                    // take care of internal branches that are partial branches above the tip
                    if (haploAboveChild != -1) {
                        branchLengths[k] = parent.getHeight() - attachTimes[0];
                        rates[k] = rate;
                        k++;
                    }
                    // else taken care of above

                    // for each internal branch of the QS subtree
                    // first and any next attach time can be above the parent/grandparent/... already
                    // so check which is the first node below the current attach time
                    int position = attachTimeArrayLength - 1;
                    for (int q = position; q > 1; q--) {
                        if (attachTimes[position] > parent.getHeight()) {
                            while (attachTimes[position] > parent.getHeight()) {
                                // is it also above the grandparent?
                                if (!parent.isRoot()) {
                                    if (attachTimes[position] > parent.getParent().getHeight()) {
                                        branchLengths[k] = parent.getParent().getHeight() - parent.getHeight();
                                        rates[k] = rate;
                                        k++;
                                        parent = parent.getParent();
                                    }
                                } else {
                                    for (int o = position; o > 1; o--) {
                                        branchLengths[k] = attachTimes[o - 1] - attachTimes[o];
                                        rates[k] = rate;
                                        k++;
                                    }
                                }
                            }
                        }
                        else {
                            parent = parent.getParent();
                        }
                    }
                    // now check for the rest of attach times
                    else {
                        for (int q = position; q > 1; q--) {
                            // watch our for the origin
                            if (parent.getHeight() > attachTimes[q]) {
                                branchLengths[k] = attachTimes[q] - parent.getHeight();
                                rates[k] = rate;
                                k++;
                            } else if ()

                                rates[k] = rate;


                            branchLengths[k] = attachTimes[q] - height;
                            rates[k] = rate;
                            k++;
                        }
                    else{
                            branchLengths[k] = parent.getHeight() - child.getHeight();
                            rates[k] = branchRateModel.getRateForBranch(child);
                            k++;
                        }


                    else{
                            double height = child.getHeight();
                            branchLengths[k] = attachTimes[attachTimeArrayLength - 1] - height;
                            rates[k] = rate;
                            k++;
                            for (int q = 1; q < attachTimeArrayLength; q++) {
                                branchLengths[k] = attachTimes[q] - height;
                                rates[k] = rate;
                                k++;
                            }
                        }
                    }
                }
            }
        }

        // from here as in original RateStatistic class //
        final double[] values = new double[3];
        double totalWeightedRate = 0.0;
        double totalTreeLength = 0.0;

        for (int i = 0; i < rates.length; i++) {

            totalWeightedRate += rates[i] * branchLengths[i];
            totalTreeLength += branchLengths[i];
        }

        // from here as in original RateStatistic class //
        values[MEAN] = totalWeightedRate / totalTreeLength;
        // Q2R why not?
        //  final double mean = DiscreteStatistics.mean(rates);
        //        values[VARIANCE] = DiscreteStatistics.variance(rates, mean);
        //        values[COEFFICIENT_OF_VARIATION] = Math.sqrt(D values[VARIANCE]) / mean;
        values[VARIANCE] = DiscreteStatistics.variance(rates);
        final double mean = DiscreteStatistics.mean(rates);
        values[COEFFICIENT_OF_VARIATION] = Math.sqrt(DiscreteStatistics.variance(rates, mean)) / mean;
        return values;
    }

}
