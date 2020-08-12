package piqmee.tree;

import beast.core.*;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.TraitSet;
import beast.util.ClusterTree;
import beast.util.ClusterTree.*;
import java.util.*;

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
    final public Input<Type> clusterTypeInput = new Input<>("clusterType", "type of clustering algorithm used for generating initial beast.tree. " +
            "Should be one of " + Arrays.toString(Type.values()) + " (default " + Type.average + ")", Type.average, Type.values());

    public QuasiSpeciesClusterTree() {
    }

    public void initAndValidate() {
        super.initAndValidate();

        if (dataInput.get() == null)
            throw new RuntimeException("The data input needs to be specified");

        ClusterTree inputTree = new ClusterTree();
        inputTree.setDateTrait(timeTraitSet);
        inputTree.initByName(
                "clusterType", clusterTypeInput.get(),
                "taxa", dataInput.get());

        if (haplotypeCountsSet != null && !haplotypeCountIsAll1(haplotypeCountsSet))
            initFromUniqueHaploTree(inputTree, dataInput.get(),collapseIdenticalSequencesInput.get(),haplotypeCountsSet);
        else
            initFromFullTree(inputTree, dataInput.get(),collapseIdenticalSequencesInput.get());

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
        stateNodes.add(this);
    }

}
