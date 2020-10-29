package piqmee.tree;

import beast.core.*;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.math.distributions.MRCAPrior;
import beast.evolution.tree.RandomTree;

import java.util.*;

/**
 * @author Veronika Boskova created on 25/08/2020.
 */

@Description("Class to initialize a QuasiSpeciesTree from the provided alignment" +
             " using coalescent to determine provided sequence topology." +
             " If duplicate sequence counts are provided in haplotypeCountsInput," +
             " their branching times will be evenly spread between the MRCA" +
             " of identical sequence sub-tree and respective tip times.")
public class QuasiSpeciesRandomTree extends QuasiSpeciesTree implements StateNodeInitialiser{

    final public Input<PopulationFunction> populationFunctionInput = new Input<>("populationModel",
            "population function for generating coalescent???", Input.Validate.REQUIRED);
    final public Input<Double> rootHeightInput = new Input<>("rootHeight",
            "If specified the tree will be scaled to match the root height, if constraints allow this");
    public Input<Boolean> collapseIdenticalSequencesInput = new Input<>("collapseIdenticalSequences",
            "Should nodes that have identical sequences be collapsed to one haplotype? " +
                    "Default true.", true);


    public QuasiSpeciesRandomTree() {
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // make sure to use date and haploCount traits
        if (m_initial.get() != null)
            processTraits(m_initial.get().m_traitList.get());
        else
            processTraits(m_traitList.get());

        // initialize the tree
        // get the input alignment
        Alignment data = dataInput.get();
        if (data instanceof FilteredAlignment) {
            data = ((FilteredAlignment) data).alignmentInput.get();
        }
        if (data == null)
            throw new RuntimeException("The data input needs to be specified");

        // create a toy tree for calculation of distances
        RandomTree toyRandomTree = new RandomTree();
        toyRandomTree.setDateTrait(timeTraitSet);
        toyRandomTree.initByName(
                "taxa", data,
                "populationModel", populationFunctionInput.get());

        // get monophyletic constraints necessary to cluster identical sequences
        // specify monophyletic clusters from distance matrix
        List<MRCAPrior> monophyleticGroups = new ArrayList();
        // Get the distances for the sequences:
        double[][] distanceMatrix = getDistanceMatrix(data, toyRandomTree, collapseSequencesWithMissingDataInput.get());
        // specify monophyletic constraints
        for (int i = 0; i < distanceMatrix.length; i++){
            List<Taxon> identical = new ArrayList<>();
            Taxon currentTaxon = new Taxon(data.getTaxaNames().get(i));
            identical.add(currentTaxon);
            for (int j = 0; j < i; j++){
                if (distanceMatrix[j][i] == 0)
                    continue;
            }
            for (int k = i+1; k < distanceMatrix.length; k++ ) {
                if (distanceMatrix[k][i] == 0) {
                    Taxon newTaxon = new Taxon(data.getTaxaNames().get(k));
                    identical.add(newTaxon);
                }
            }
            if (identical.size() == 1)
                continue;

            MRCAPrior group = new MRCAPrior();
            group.initByName(
                    "tree", toyRandomTree,
                    "taxonset", new TaxonSet(identical),
                    "monophyletic", "true");
            monophyleticGroups.add(group);
        }

        // initialize random tree with constraints specified such that all identical sequences form monophyletic clusters
        RandomTree inputTree = new RandomTree();
        inputTree.setDateTrait(timeTraitSet);
        inputTree.initByName(
                "taxa", data,
                "populationModel", populationFunctionInput.get(),
                "constraint", monophyleticGroups,
                "rootHeight", rootHeightInput.get());

        // initialize the quasispecies tree - and collapse identical sequences, if necessary
        if (haplotypeCountsSet != null && !haplotypeCountIsAll1(haplotypeCountsSet))
            initFromUniqueHaploTree(inputTree, data,
                    collapseIdenticalSequencesInput.get(), collapseSequencesWithMissingDataInput.get(),
                    haplotypeCountsSet);
        else
            initFromFullTree(inputTree, data,
                    collapseIdenticalSequencesInput.get(), collapseSequencesWithMissingDataInput.get());

        initStateNodes();
    }

    @Override
    public void initStateNodes() {
        if (m_initial.get() != null) {
            m_initial.get().assignFromWithoutID(this);
        }
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        if (m_initial.get() != null) {
            stateNodes.add(m_initial.get());
        }
    }

}
