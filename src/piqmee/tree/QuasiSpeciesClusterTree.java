package piqmee.tree;

import beast.core.*;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.TraitSet;
import beast.util.ClusterTree;

import java.util.List;

/**
 * @author Veronika Boskova created on 09/03/2017.
 */

@Description("Class to initialize a QuasiSpeciesTree from the alignment only")
public class QuasiSpeciesClusterTree extends QuasiSpeciesTree implements StateNodeInitialiser{

    final public Input<Alignment> dataInput = new Input<>("data",
            "Alignment data used for calculating distances for clustering",
            Input.Validate.REQUIRED);
    public Input<Boolean> collapseIdenticalSequencesInput = new Input<>("collapseIdenticalSequences",
            "Should nodes that have identical sequences be collapsed to one haplotype? " +
                    "Default true.", true);

    public QuasiSpeciesClusterTree() {

        // When specifying the input tree with full newick tree, we do not allow for duplicate counts input on the top.
        haplotypeCountsInput.setRule(Input.Validate.OPTIONAL);

    }

    public void initAndValidate() {
        super.initAndValidate();

        ClusterTree inputTree = new ClusterTree();
        inputTree.setDateTrait(timeTraitSet);
        inputTree.initByName(
                "clusterType", "upgma",
                "taxa", dataInput.get());

        if (dataInput.get() == null)
            throw new RuntimeException("The data input needs to be specified");

        if (haplotypeCountsInput.get() != null || !haplotypeCountIsAll1(haplotypeCountsInput.get()))
            initFromUniqueHaploTree(inputTree, dataInput.get(),collapseIdenticalSequencesInput.get(),haplotypeCountsInput.get());
        else
            initFromFullTree(inputTree, dataInput.get(),collapseIdenticalSequencesInput.get());
    }

    @Override
    public void initStateNodes(){
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(this);
    }

}
