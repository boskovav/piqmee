package piqmee.evolution.branchratemodel;

import beast.base.core.Description;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.RateStatistic;
import beast.base.evolution.tree.Node;
import beast.base.util.DiscreteStatistics;
import piqmee.tree.QuasiSpeciesNode;
import piqmee.tree.QuasiSpeciesTree;
import java.util.Arrays;

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

    Node toyNode = new Node();

    // Empty constructor as required:
//    public QuasiSpeciesRateStatistic() { };

    @Override
    public void initAndValidate() {
        tree = (QuasiSpeciesTree) treeInput.get();
        branchRateModel = branchRateModelInput.get();
        if (branchRateModel == null) {
            branchRateModel = likelihoodInput.get().branchRateModelInput.get();
        }
        this.internal = internalInput.get();
        this.external = externalInput.get();
    }



    /**
     * calculate the three statistics from scratch *
     */
    @Override
    public double[] calcValues() {
        int length = 0;
        int lengthBL = 0;
        int offset = 0;
        int offsetForBranchLengthsAndRates = 0;

        final int nrOfLeafs = tree.getLeafNodeCount();
        final int nrOfQSDupl = tree.getTotalAttachmentCounts();
        final int nrOfLeafsAndQSDupl = nrOfQSDupl + nrOfLeafs;

        if (external) {
            length += nrOfLeafsAndQSDupl;
            lengthBL += nrOfLeafs;
        }
        if (internal) {
            length += tree.getInternalNodeCount() - 1 + nrOfQSDupl;
            // nrOfLeafs is 2 * because
            //      1) for the branches above haplo, when haplo starts above a tip and
            //      2) for internal QS branches
            lengthBL += tree.getInternalNodeCount() - 1 + 2 * nrOfLeafs;
        }

        // these are the rates for each branch as in the full tree when reconstructed from the QS tree
        final double[] rates = new double[length];
        // these are the rates for each branch as in the full tree when reconstructed from the QS tree in case one haplo
        //  is above the root --- and so we have one less rate (no partial branch rate above that haplo)
        final double[] ratestmp = new double[length-1];
        // for internal QS branches and lengths...need to split from external QS branch lengths
//        final double[] internalQSBranchLengths = new double[nrOfLeafs];
//        final double[] internalQSBranchRates = new double[nrOfLeafs];
        // these are the rates for each branch length, where branch length for the tip is the total
        //  branch length spanned by the haplotype (minus the internal QS branch lengths -- in internalQSBranchLengths)
        final double[] ratesForBranchLengths = new double[lengthBL];
        Arrays.fill(ratesForBranchLengths, Double.NaN );
        // need those only for mean
        final double[] branchLengths = new double[lengthBL];

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
                // 1) it can be that there are no duplicates, so the branch is from the child to the parent
                //      and has a rate of tree.getNodeCount()+child.getNr()
                if (haploAboveChild != -1 && attachTimeArrayLength == 1) {
                    branchLengths[i] = parent.getHeight() - child.getHeight();
                    toyNode.setNr(tree.getNodeCount()+child.getNr());
                    rate = branchRateModel.getRateForBranch(toyNode);
                    rates[k] = rate;
                    k++;
                }
                // 2) and now for QS with duplicates
                else {
                    branchLengths[i] = child.getTotalBranchLengths() - (attachTimes[0]-attachTimes[attachTimeArrayLength-1]);

                    // for each branch of the haplotype add contribution to rates array (also for 0 - first = "true" branch)
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
            // contribution from the partial QS branches
            for (int i = 0; i < nrOfLeafs; i++){
                final QuasiSpeciesNode child = (QuasiSpeciesNode) nodes[i];
                double rate = branchRateModel.getRateForBranch(child);
                double[] attachTimes = child.getAttachmentTimesList();
                int attachTimeArrayLength = attachTimes.length;
                // check if the partial branch is above the tip
                //      (else above another internal node --- treated in internal node loop)
                if (attachTimeArrayLength > 1){
                    // this is for the QS duplicate branches and respective internal QS branches
                    for (int q = 0; q < attachTimeArrayLength - 2; q++) {
                        rates[k + q] = rate;
                    }
                    k += attachTimeArrayLength - 2;

                    if (attachTimeArrayLength > 2) {
                        branchLengths[l] = attachTimes[0] - attachTimes[attachTimeArrayLength - 1];
                        ratesForBranchLengths[l] = rate;
                        l++;
                    }
                    // it always has toyNode branchRate
                    // but check what length to count for the branch rate
                    //  is haplo above tip count it here, if haplo above internal node count it below
                    if (child.getHaploAboveName() != -1) {
                        toyNode.setNr(tree.getNodeCount()+child.getNr());
                        rate = branchRateModel.getRateForBranch(toyNode);
                        rates[k] = rate;
                        k++;

                        branchLengths[l] = child.getParent().getHeight() - attachTimes[0];
                        ratesForBranchLengths[l] = rate;
                        l++;
                    }
                }
            }

            // contribution from the real internal nodes only
            for (int i = nrOfLeafs; i < n; i++) {
                final QuasiSpeciesNode child = (QuasiSpeciesNode) nodes[i];
                //orig-root branch has rate 1 --- defined in the relaxed clock models --- so do not include here
                if (!child.isRoot()) {
                    double parentHeight = child.getParent().getHeight();
                    int haploAboveChild = child.getHaploAboveName();
                    int haploContChild = child.getContinuingHaploName();
                    double rate = branchRateModel.getRateForBranch(child);
                    // Case 1: branch can be free - i.e. no haplo passing, branch has its rate
                    if (haploContChild == -1) {
                        rates[k] = rate;

                        branchLengths[l] = parentHeight - child.getHeight();
                        ratesForBranchLengths[l] = rate;
                        l++;
                    }
                    // Case 2: it can be that a haplo starts above this node, in that case
                    //          count only the rate for the partial branch which has its own toyNode branchRate
                    else if (haploAboveChild != -1) {
                        toyNode.setNr(tree.getNodeCount()+haploAboveChild);
                        rate = branchRateModel.getRateForBranch(toyNode);
                        rates[k] = rate;
                        k++;
                        // and one of the haplo's branches was split in two
                        rates[k] = branchRateModel.getRateForBranch(tree.getNode(haploAboveChild));

                        branchLengths[l] = parentHeight - ((QuasiSpeciesNode) tree.getNode(haploAboveChild)).getAttachmentTimesList()[0];
                        ratesForBranchLengths[l] = rate;
                        l++;
                    }
                    // Case 3: it can be that a haplo passes this branch, in that case
                    //         one of the haplo's branches was split in two and the rate
                    //         is the same as for the corresponding tip branch rate
                    else {
                        rates[k] = branchRateModel.getRateForBranch(tree.getNode(haploContChild));
                    }
                    k++;
                }
            }
        }

        final double[] values = new double[3];
        double totalWeightedRate = 0.0;
        double totalTreeLength = 0.0;

        for (int i = 0; i < ratesForBranchLengths.length ; i++) {
            if (! Double.isNaN(ratesForBranchLengths[i])) {
                totalWeightedRate += ratesForBranchLengths[i] * branchLengths[i];
                totalTreeLength += branchLengths[i];
            }
        }

        // from here as in original RateStatistic class //
        values[MEAN] = totalWeightedRate / totalTreeLength;
        // Q2R why not?
        //  final double mean = DiscreteStatistics.mean(rates);
        //        values[VARIANCE] = DiscreteStatistics.variance(rates, mean);
        //        values[COEFFICIENT_OF_VARIATION] = Math.sqrt(D values[VARIANCE]) / mean;

        // if a haplo is above a root, we do not have a partial branch rate for the branch above that haplo
        if (((QuasiSpeciesNode)tree.getRoot()).getHaploAboveName() != -1){
            System.arraycopy(rates,0,ratestmp, 0, length -1);
            values[VARIANCE] = DiscreteStatistics.variance(ratestmp);
            final double mean = DiscreteStatistics.mean(ratestmp);
            values[COEFFICIENT_OF_VARIATION] = Math.sqrt(DiscreteStatistics.variance(ratestmp, mean)) / mean;
        } else {
            values[VARIANCE] = DiscreteStatistics.variance(rates);
            final double mean = DiscreteStatistics.mean(rates);
            values[COEFFICIENT_OF_VARIATION] = Math.sqrt(DiscreteStatistics.variance(rates, mean)) / mean;
        }
        return values;
    }

}
