package quasispeciestree.distributions;

import java.util.List;
import java.util.Random;

import beast.core.*;
//import beast.evolution.tree.coalescent.TreeIntervals;
import quasispeciestree.tree.QuasiSpeciesTree;


/**
 * @author Veronika Boskova created on 26/06/2015
 */

@Description("A quasi-species phylogenetic tree base class.")
@Citation("Boskova, V., Stadler, T., Quasispecies algorithm")
public abstract class QuasiSpeciesBaseTreeDistribution extends Distribution {

    final public Input<QuasiSpeciesTree> qsTreeInput = new Input<>(
            "quasiSpeciesTree", "Quasi-Species tree over which to calculate a prior or likelihood",
            Input.Validate.REQUIRED);

//    final public Input<TreeIntervals> treeIntervalsInput = new Input<>("treeIntervals", "Intervals for a phylogenetic beast tree", Input.Validate.XOR, qsTreeInput);

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }

//    @Override
//    protected boolean requiresRecalculation() {
//        final TreeIntervals ti = treeIntervalsInput.get();
//        if (ti != null) {
//            //boolean d = ti.isDirtyCalculation();
//            //assert d;
//            assert ti.isDirtyCalculation();
//            return true;
//        }
//        return treeInput.get().somethingIsDirty();
//    }

    /** Indicate that the tree distribution can deal with dated tips in the tree
     * Some tree distributions like the Yule prior cannot handle this.
     * @return true by default
     */
    public boolean canHandleTipDates() {
        return true;
    }
}