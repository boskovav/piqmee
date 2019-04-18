package piqmee.evolution.branchratemodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.UCRelaxedClockModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;
import org.apache.commons.math.MathException;
import piqmee.tree.QuasiSpeciesNode;

import java.util.Arrays;

/**
 *  @author Veronika Boskova created on 28/03/2019 based on A.Drummond's UCRelaxedClockModel class
 */
@Description("Defines an uncorrelated relaxed molecular clock for PIQMEE model.")
public class QuasiSpeciesUCRelaxedClockModel extends BranchRateModel.Base {

    final public Input<ParametricDistribution> rateDistInput = new Input<>("distr", "the distribution governing the rates among branches. Must have mean of 1. The clock.rate parameter can be used to change the mean rate.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> categoryInput = new Input<>("rateCategories", "the rate categories associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);

    final public Input<Integer> numberOfDiscreteRates = new Input<>("numberOfDiscreteRates", "the number of discrete rate categories to approximate the rate distribution by. A value <= 0 will cause the number of categories to be set equal to the number of branches in the tree. (default = -1)", -1);

    final public Input<RealParameter> quantileInput = new Input<>("rateQuantiles", "the rate quantiles associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.XOR, categoryInput);

    final public Input<Tree> treeInput = new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
    final public Input<Boolean> normalizeInput = new Input<>("normalize", "Whether to normalize the average rate (default false).", false);

    RealParameter meanRate;

    int LATTICE_SIZE_FOR_DISCRETIZED_RATES = 100;

    // true if quantiles are used, false if discrete rate categories are used.
    boolean usingQuantiles;

    private int branchCount;

    private Node toyNode = new Node();

    @Override
    public void initAndValidate() {

        tree = treeInput.get();
        // for purposes of relaxed clock in piqmee, we have to have number of branches equal to number of nodes+tips
        // (1 for root) but this should always be 1
        branchCount = tree.getNodeCount() + tree.getLeafNodeCount();

        categories = categoryInput.get();
        usingQuantiles = (categories == null);

        if (!usingQuantiles) {
            LATTICE_SIZE_FOR_DISCRETIZED_RATES = numberOfDiscreteRates.get();
            if (LATTICE_SIZE_FOR_DISCRETIZED_RATES <= 0) LATTICE_SIZE_FOR_DISCRETIZED_RATES = branchCount;
            Log.info.println("  UCRelaxedClockModel: using " + LATTICE_SIZE_FOR_DISCRETIZED_RATES + " rate " +
                    "categories to approximate rate distribution across branches.");
        } else {
            if (numberOfDiscreteRates.get() != -1) {
                throw new RuntimeException("Can't specify both numberOfDiscreteRates and rateQuantiles inputs.");
            }
            Log.info.println("  UCRelaxedClockModel: using quantiles for rate distribution across branches.");
        }

        if (usingQuantiles) {
            quantiles = quantileInput.get();
            quantiles.setDimension(branchCount);
            Double[] initialQuantiles = new Double[branchCount];
            for (int i = 0; i < branchCount; i++) {
                initialQuantiles[i] = Randomizer.nextDouble();
            }
            RealParameter other = new RealParameter(initialQuantiles);
            quantiles.assignFromWithoutID(other);
            quantiles.setLower(0.0);
            quantiles.setUpper(1.0);
        } else {
            categories.setDimension(branchCount);
            Integer[] initialCategories = new Integer[branchCount];
            for (int i = 0; i < branchCount; i++) {
                initialCategories[i] = Randomizer.nextInt(LATTICE_SIZE_FOR_DISCRETIZED_RATES);
            }
            // set initial values of rate categories
            IntegerParameter other = new IntegerParameter(initialCategories);
            categories.assignFromWithoutID(other);
            categories.setLower(0);
            categories.setUpper(LATTICE_SIZE_FOR_DISCRETIZED_RATES - 1);
        }

        distribution = rateDistInput.get();

        if (!usingQuantiles) {
            // rates are initially zero and are computed by getRawRate(int i) as needed
            rates = new double[LATTICE_SIZE_FOR_DISCRETIZED_RATES];
            storedRates = new double[LATTICE_SIZE_FOR_DISCRETIZED_RATES];

            //System.arraycopy(rates, 0, storedRates, 0, rates.length);
        }
        normalize = normalizeInput.get();

        meanRate = meanRateInput.get();
        if (meanRate == null) {
            meanRate = new RealParameter("1.0");
        }

        try {
            double mean = rateDistInput.get().getMean();
            if (Math.abs(mean - 1.0) > 1e-6) {
                Log.warning.println("WARNING: mean of distribution for relaxed clock model is not 1.0.");
            }
        } catch (RuntimeException e) {
            // ignore
        }
    }

    @Override
    public double getRateForBranch(Node node) {
        if (node.isRoot() && node.getNr() == (tree.getNodeCount()-1)) {
            // root has no rate
            return 1;
        }

        if (recompute) {
            // this must be synchronized to avoid being called simultaneously by
            // two different likelihood threads
            synchronized (this) {
                prepare();
                recompute = false;
            }
        }

        if (renormalize) {
            if (normalize) {
                synchronized (this) {
                    computeFactor();
                }
            }
            renormalize = false;
        }

        return getRawRate(node) * scaleFactor * meanRate.getValue();
    }

    /**
     * Computes a scale factor for normalization. Only called if normalize=true.
     */
    private void computeFactor() {

        //scale mean rate to 1.0 or separate parameter

        double treeRate = 0.0;
        double treeTime = 0.0;

        if (!usingQuantiles) {
            for (int i = 0; i < (tree.getNodeCount() + tree.getLeafNodeCount()); i++) {
                Node node;
                if (i >= tree.getNodeCount()) {
                    toyNode.setNr(i);
                    node = tree.getNode(i-tree.getNodeCount());
                    treeRate += getRawRateForCategory(toyNode) *
                            (node.getLength() - ((QuasiSpeciesNode) node).getAttachmentTimesList()[0]);
                    treeTime += node.getLength() - ((QuasiSpeciesNode) node).getAttachmentTimesList()[0];
                } else {
                    node = tree.getNode(i);
                    if (!node.isRoot()) {
                        if (node.isLeaf()) {
                            treeRate += getRawRateForCategory(node) * ((QuasiSpeciesNode) node).getTotalBranchLengths();
                            treeTime += ((QuasiSpeciesNode) node).getTotalBranchLengths();
                        } else {
                            treeRate += getRawRateForCategory(node) * node.getLength();
                            treeTime += node.getLength();
                        }
                    }
                }
            }
        } else {
            for (int i = 0; i < (tree.getNodeCount() + tree.getLeafNodeCount()); i++) {
                Node node;
                if (i >= tree.getNodeCount()) {
                    toyNode.setNr(i);
                    node = tree.getNode(i-tree.getNodeCount());
                    treeRate += getRawRateForQuantile(toyNode) *
                            (node.getLength() - ((QuasiSpeciesNode) node).getAttachmentTimesList()[0]);
                    treeTime += node.getLength() - ((QuasiSpeciesNode) node).getAttachmentTimesList()[0];
                } else {
                    node = tree.getNode(i);
                    if (!node.isRoot()) {
                        if (node.isLeaf()) {
                            treeRate += getRawRateForQuantile(node) * ((QuasiSpeciesNode) node).getTotalBranchLengths();
                            treeTime += ((QuasiSpeciesNode) node).getTotalBranchLengths();
                        } else {
                            treeRate += getRawRateForQuantile(node) * node.getLength();
                            treeTime += node.getLength();
                        }
                    }
                }
            }
        }

        scaleFactor = 1.0 / (treeRate / treeTime);
    }

    private double getRawRate(Node node) {
        if (usingQuantiles) {
            return getRawRateForQuantile(node);
        }
        return getRawRateForCategory(node);
    }

    /**
     * @param node the node to get the rate of
     * @return the rate of the branch
     */
    private double getRawRateForCategory(Node node) {

        int nodeNumber = node.getNr();
        // VB: what the hell is this??? when does this happen that NodeNr > branchCount?
        // if (nodeNumber == branchCount) {
        if (nodeNumber == (tree.getNodeCount() - 1)) {
            // root node has nr less than #categories, so use that nr
            nodeNumber = node.getTree().getRoot().getNr();
        }

        int category = categories.getValue(nodeNumber);

        if (rates[category] == 0.0) {
            try {
                rates[category] = distribution.inverseCumulativeProbability((category + 0.5) / rates.length);
            } catch (MathException e) {
                throw new RuntimeException("Failed to compute inverse cumulative probability!");
            }
        }
        return rates[category];
    }

    private double getRawRateForQuantile(Node node) {

        int nodeNumber = node.getNr();
        // VB: what the hell is this??? when does this happen that NodeNr > branchCount?
        // if (nodeNumber == branchCount) {
        if (nodeNumber == (tree.getNodeCount() - 1)) {
            // root node has nr less than #categories, so use that nr
            nodeNumber = node.getTree().getRoot().getNr();
        }

        try {
            return distribution.inverseCumulativeProbability(quantiles.getValue(nodeNumber));
        } catch (MathException e) {
            throw new RuntimeException("Failed to compute inverse cumulative probability!");
        }
    }

    private void prepare() {

        categories = categoryInput.get();

        usingQuantiles = (categories == null);

        distribution = rateDistInput.get();

        tree = treeInput.get();

        if (!usingQuantiles) {
            // rates array initialized to correct length in initAndValidate
            // here we just reset rates to zero and they are computed by getRawRate(int i) as needed
            Arrays.fill(rates, 0.0);
        }
    }


    @Override
    protected boolean requiresRecalculation() {
        recompute = false;
        renormalize = true;

        if (rateDistInput.get().isDirtyCalculation()) {
            recompute = true;
            return true;
        }
        // NOT processed as trait on the tree, so DO mark as dirty
        if (categoryInput.get() != null && categoryInput.get().somethingIsDirty()) {
            //recompute = true;
            return true;
        }

        if (quantileInput.get() != null && quantileInput.get().somethingIsDirty()) {
            return true;
        }

        if (meanRate.somethingIsDirty()) {
            return true;
        }

        return recompute;
    }

    @Override
    public void store() {
        if (!usingQuantiles) System.arraycopy(rates, 0, storedRates, 0, rates.length);

        storedScaleFactor = scaleFactor;
        super.store();
    }

    @Override
    public void restore() {
        if (!usingQuantiles) {
            double[] tmp = rates;
            rates = storedRates;
            storedRates = tmp;
        }
        scaleFactor = storedScaleFactor;
        super.restore();
    }

    ParametricDistribution distribution;
    IntegerParameter categories;
    RealParameter quantiles;
    Tree tree;

    private boolean normalize = false;
    private boolean recompute = true;
    private boolean renormalize = true;

    private double[] rates;
    private double[] storedRates;
    private double scaleFactor = 1.0;
    private double storedScaleFactor = 1.0;

}
