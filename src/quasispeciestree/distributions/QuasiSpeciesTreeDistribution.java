package quasispeciestree.distributions;

import java.util.List;
import java.util.Random;

import beast.core.*;
import quasispeciestree.tree.QuasiSpeciesTree;
import quasispeciestree.distributions.QuasiSpeciesBaseTreeDistribution;

/**
 * Created by Veronika Boskova on 28/02/2017.
 */

@Description("A quasi-species phylogenetic tree likelihood class.")
@Citation("Boskova, V., Stadler, T., Quasispecies algorithm")
abstract public class QuasiSpeciesTreeDistribution extends QuasiSpeciesBaseTreeDistribution {

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     *
     * @return the log likelihood
     */
    @Override
    public double calculateLogP() {
        final QuasiSpeciesTree tree = qsTreeInput.get();
        logP = calculateTreeLogLikelihood(tree);
        return logP;
    } // calculateLogP


    /**
     * Generic likelihood calculation
     *
     * @param qsTree
     * @return log-likelihood of density
     */
    public abstract double calculateTreeLogLikelihood(QuasiSpeciesTree qsTree);

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************


    /*****************************************/
    /** Distribution implementation follows **/
    /**
     * *************************************
     */
    @Override
    public List<String> getArguments() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public List<String> getConditions() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("This should eventually sample a tree conditional on provided speciation model.");
    }
}