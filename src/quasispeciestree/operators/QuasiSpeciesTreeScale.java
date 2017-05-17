package quasispeciestree.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import quasispeciestree.tree.QuasiSpeciesNode;

import java.util.ArrayList;
import java.util.List;

/**
 *  @author Veronika Boskova created on 21/04/2016 finished on (24/04/2016) 13/01/2017
 */
@Description("Scale operator for quasispecies trees. Also allows additional "
        + "scalar parameters to be rescaled (either forward or inversely) "
        + "at the same time.")
public class QuasiSpeciesTreeScale extends QuasiSpeciesTreeOperator{

    public Input<List<RealParameter>> parametersInput =
        new Input<>("parameter",
            "Scale this scalar parameter by the same amount as tree.",
            new ArrayList<RealParameter>());

    public Input<List<BooleanParameter>> indicatorsInput =
        new Input<>("indicator",
            "If provided, used to specify a subset of parameter elements to scale.",
            new ArrayList<BooleanParameter>());

    public Input<List<RealParameter>> parametersInverseInput =
        new Input<>("parameterInverse",
            "Scale this scalar parameter inversely.",
            new ArrayList<RealParameter>());

    public Input<List<BooleanParameter>> indicatorsInverseInput =
        new Input<>("indicatorInverse",
            "If provided, used to specify a subset of parameter elements to scale "
            + "inversely.",
            new ArrayList<BooleanParameter>());

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
           "Scaling is restricted to the range [1/scaleFactor, scaleFactor]");

    final public Input<Boolean> rootOnlyInput = new Input<>("rootOnly",
            "scale root of a tree only, ignored if tree is not specified (default false)", false);

    boolean indicatorsUsed, indicatorsInverseUsed;

    @Override
    public void initAndValidate() {

        super.initAndValidate();

        if (indicatorsInput.get().size() > 0) {
            if (indicatorsInput.get().size() != parametersInput.get().size())
                throw new IllegalArgumentException("If an indicator element "
                        + "exists, the number of such elements must equal "
                        + "the number of parameter elements.");

            for (int pidx = 0; pidx < parametersInput.get().size(); pidx++) {
                if (parametersInput.get().get(pidx).getDimension() !=
                        indicatorsInput.get().get(pidx).getDimension()) {
                    throw new IllegalArgumentException("The number of boolean "
                            + "values in indicator element "
                            + String.valueOf(pidx+1)
                            + " doesn't match the dimension of the "
                            + "corresponding parameter element.");
                }
            }
            indicatorsUsed = true;
        } else
            indicatorsUsed = false;

        if (indicatorsInverseInput.get().size() > 0) {
            if (indicatorsInverseInput.get().size() != parametersInverseInput.get().size())
                throw new IllegalArgumentException("If an indicatorInverse element "
                        + "exists, the number of such elements must equal "
                        + "the number of parameterInverse elements.");

            for (int pidx = 0; pidx < parametersInverseInput.get().size(); pidx++) {
                if (parametersInverseInput.get().get(pidx).getDimension() !=
                        indicatorsInverseInput.get().get(pidx).getDimension()) {
                    throw new IllegalArgumentException("The number of boolean "
                            + "values in indicatorInverse element "
                            + String.valueOf(pidx+1)
                            + " doesn't match the dimension of the "
                            + "corresponding parameterInverse element.");
                }
            }
            indicatorsInverseUsed = true;
        } else
            indicatorsInverseUsed = false;
    }

    @Override
    public double proposal() {

        try {

            // Choose scale factor:
            double u = Randomizer.nextDouble();
            double f = u * scaleFactorInput.get() + (1.0 - u) / scaleFactorInput.get();

            // Keep track of Hastings ratio:
            double logf = Math.log(f);
            double logHastingsRatio = 0.0;

            final Node root = qsTree.getRoot();

            // scaling only the root
            if (rootOnlyInput.get()) {
                double oldHeight = root.getHeight();
                final double newHeight = oldHeight * f;
                logHastingsRatio -= logf;

                if (newHeight < Math.max(root.getLeft().getHeight(), root.getRight().getHeight())) {
                    return Double.NEGATIVE_INFINITY;
                }

                // There are in total 5 possibilities when moving root down and 3 possibilities when moving it up
                //  down:   1) no haplo passing the root before, nor after the move
                //          2) no haplo passing the root before, but one passing it after
                //          3) no haplo passing the root before, but 2 passing it after -- this has to be rejected
                //          4) one haplo passing the root before, the same single haplo passing it after
                //          5) one haplo passing the root before, but 2 passing it after -- this has to be rejected
                //          cases 3 and 5 need to be rejected
                //  up:     1) no haplo passing the root before, nor after the move
                //          2) one haplo passing the root before, none after
                //          3) one haplo passing the root before, the same single haplo passing it after
                //          in theory, none of the up moves need extra modifications BUT
                //
                //          NO FURTHER MOVES:       case "down 1" is a reciprocal of "up 1"
                //                                  case "down 4" is a reciprocal of "up 3"
                //          NEED ANNOTATION CHANGES:case "down 2" is a reciprocal of "up 2"
                //
                //
                //

                int haplo = ((QuasiSpeciesNode) root).getHaploAboveName();
                QuasiSpeciesNode left = (QuasiSpeciesNode) root.getLeft();
                int haploleft = ((QuasiSpeciesNode) root.getLeft()).getHaploAboveName();
                double[] lefttemqstimes = null;
                QuasiSpeciesNode right = (QuasiSpeciesNode) root.getRight();
                int haploright = ((QuasiSpeciesNode) root.getRight()).getHaploAboveName();
                double[] righttemqstimes = null;

                // check if the new height is above at least one haplo (or both haplo from left and right)
                boolean leftabove = false;
                boolean rightabove = false;
                if (haplo != -1 && haplo == left.getContinuingHaploName()) {
                    lefttemqstimes = ((QuasiSpeciesNode) qsTree.getNode(haplo)).getAttachmentTimesList();
                    if (newHeight < lefttemqstimes[0]) {
                        leftabove = true;
                    }
                } else if (haploleft != -1) {
                    lefttemqstimes = ((QuasiSpeciesNode) qsTree.getNode(haploleft)).getAttachmentTimesList();
                    if (newHeight < lefttemqstimes[0]) {
                        leftabove = true;
                    }
                }
                if (haplo != -1 && haplo == right.getContinuingHaploName()) {
                    righttemqstimes = ((QuasiSpeciesNode) qsTree.getNode(haplo)).getAttachmentTimesList();
                    if (newHeight < righttemqstimes[0]) {
                        rightabove = true;
                    }
                } else if (haploright != -1) {
                    righttemqstimes = ((QuasiSpeciesNode) qsTree.getNode(haploright)).getAttachmentTimesList();
                    if (newHeight < righttemqstimes[0]) {
                        rightabove = true;
                    }
                }
                // if both haplo are above the new height, reject the move
                if (leftabove && rightabove) {
                    return Double.NEGATIVE_INFINITY;
                }
                root.setHeight(newHeight);

                // check if annotations need to be changed
                // the haplo at the root can come from the left child
                if (haplo != -1 && haplo == left.getContinuingHaploName()) {
                    if (haploleft != -1) {
                        throw new IllegalStateException("Seems like we found a haplotype at the root " + "that should come from the left subtree but yet there is another haplotype " + "arising at the left child of the root. This should not happen!");
                    } else
                        haploleft = haplo;

                    // change annotations if the left haplo is no longer above the root
                    if (newHeight > oldHeight && haploleft != -1 && lefttemqstimes[0] < newHeight) {
                        left.setHaploAboveName(haploleft);
                        ((QuasiSpeciesNode) root).setHaploAboveName(-1);
                        // correct continuing haplo and parent haplo
                        recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) root);
                    }
                }
                // the haplo at the root can come from the right child
                else if (haplo != -1 && haplo == right.getContinuingHaploName()) {
                    if (haploright != -1) {
                        throw new IllegalStateException("Seems like we found a haplotype at the root " + "that should come from the right subtree but yet there is another haplotype " + "arising at the right child of the root. This should not happen!");
                    } else
                        haploright = haplo;

                    // change annotations if the right haplo is no longer above the root
                    if (newHeight > oldHeight && righttemqstimes[0] < newHeight) {
                        right.setHaploAboveName(haploright);
                        ((QuasiSpeciesNode) root).setHaploAboveName(-1);
                        // correct continuing haplo and parent haplo
                        recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) root);
                    }
                } else if (haplo == -1 && newHeight < oldHeight && (leftabove || rightabove)) {
                    if (leftabove) {
                        left.setHaploAboveName(-1);
                        ((QuasiSpeciesNode) root).setHaploAboveName(haploleft);
                        // correct continuing haplo and parent haplo
                        recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) root);
                    } else {
                        right.setHaploAboveName(-1);
                        ((QuasiSpeciesNode) root).setHaploAboveName(haploright);
                        // correct continuing haplo and parent haplo
                        recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) root);
                    }
                }
                // recalculate countPossibleStartBranches
                qsTree.countAndSetPossibleStartBranches();

                // scaling the entire tree
            } else {
                // scale the quasispecies.tree
                final int totalNodes = qsTree.scale(f);
                logHastingsRatio += (totalNodes - 2) * logf;
                //return Math.log(f) * (totalNodes - 2);
            }

            // Reject invalid tree/root scaling: - done in BDSKY
            //if (f > 1.0 && root.getHeight() > origin.getValue())
            //        return Double.NEGATIVE_INFINITY;

            // Scale parameters:
            for (int pidx = 0; pidx < parametersInput.get().size(); pidx++) {
                RealParameter param = parametersInput.get().get(pidx);
                for (int i = 0; i < param.getDimension(); i++) {
                    if (!indicatorsUsed || indicatorsInput.get().get(pidx).getValue(i)) {
                        double oldValue = param.getValue(i);
                        double newValue = oldValue * f;
                        if (newValue < param.getLower() || newValue > param.getUpper())
                            return Double.NEGATIVE_INFINITY;

                        param.setValue(i, newValue);
                        logHastingsRatio += logf;
                    }
                }
            }

            // Scale parameters inversely:
            for (int pidx = 0; pidx < parametersInverseInput.get().size(); pidx++) {
                RealParameter param = parametersInverseInput.get().get(pidx);
                for (int i = 0; i < param.getDimension(); i++) {
                    if (!indicatorsInverseUsed || indicatorsInverseInput.get().get(pidx).getValue(i)) {
                        double oldValue = param.getValue(i);
                        double newValue = oldValue / f;
                        if (newValue < param.getLower() || newValue > param.getUpper())
                            return Double.NEGATIVE_INFINITY;

                        param.setValue(i, newValue);
                        logHastingsRatio -= logf;
                    }
                }
            }

            // Return Hastings ratio:
            return logHastingsRatio;
        } catch (Exception e) {
            // whatever went wrong, we want to abort this operation...
            return Double.NEGATIVE_INFINITY;
        }
    }
}
