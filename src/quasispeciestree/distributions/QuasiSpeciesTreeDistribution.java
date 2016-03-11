package quasispeciestree.distributions;

import beast.core.*;
import quasispeciestree.tree.QuasiSpeciesTree;

import java.util.List;
import java.util.Random;


/**
 * @author Veronika Boskova created on 26/06/2015
 */

@Description("A quasi-species phylogenetic tree.")
@Citation("Boskova, V., Stadler, T., Quasispecies algorithm")

//public class QuasiSpeciesTreeDistribution extends SpeciesTreeDistribution {
public abstract class QuasiSpeciesTreeDistribution extends Distribution {

    public Input<QuasiSpeciesTree> qsTreeInput = new Input<>(
            "quasiSpeciesTree", "Quasi-Species tree over which to calculate a prior or likelihood",
            Input.Validate.REQUIRED);

    protected QuasiSpeciesTree qsTree;

    @Override
    public void initAndValidate(){
        qsTree = qsTreeInput.get();
    }

    @Override
    public double calculateLogP() {
        final QuasiSpeciesTree tree = qsTree;
        logP = calculateTreeLogLikelihood(tree);
        return logP;
    }

    public abstract double calculateTreeLogLikelihood(QuasiSpeciesTree qsTree);



    // Interface requirements:


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
}