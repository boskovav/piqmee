package piqmee.tree;

import beast.core.*;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.util.ClusterTree;
import beast.util.ClusterTree.*;

import java.util.*;

/**
 * @author Veronika Boskova created on 09/03/2017.
 */

@Description("Class to initialize a QuasiSpeciesTree from the alignment only")
public class QuasiSpeciesClusterTree extends QuasiSpeciesTree implements StateNodeInitialiser{

    final public Input<Type> clusterTypeInput = new Input<>("clusterType", "type of clustering algorithm used for generating initial beast.tree. " +
            "Should be one of " + Arrays.toString(Type.values()) + " (default " + Type.average + ")", Type.average, Type.values());
    public Input<Boolean> collapseIdenticalSequencesInput = new Input<>("collapseIdenticalSequences",
            "Should nodes that have identical sequences be collapsed to one haplotype? " +
                    "Default true.", true);
    public Input<Boolean> collapseSequencesWithMissingDataInput = new Input<>("collapseSequencesIfIdenticalUpToMissingParts",
            "Flag to indicate if sequences that have missing data (stretches of N's) should" +
                    "be collapsed with a sequence that is identical to it up the missing data. Default false.",
            false);

    public QuasiSpeciesClusterTree() {
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

        if (collapseSequencesWithMissingDataInput.get() && !collapseIdenticalSequencesInput.get())
            throw new RuntimeException("It seems that you want to collapse sequences that have parts of missing data to" +
                    "their closest resembling sequence. If this is the case, you need to set both " +
                    "collapseIdenticalSequences and collapseSequencesIfIdenticalUpToMissingParts to true.");

        ClusterTree inputTree = new ClusterTree();
        inputTree.setDateTrait(timeTraitSet);
        inputTree.initByName(
                "clusterType", clusterTypeInput.get(),
                "taxa", data);

        // initialize the quasispecies tree - and collapse identical sequences, if necessary
        if (haplotypeCountsSet != null && !haplotypeCountIsAll1(haplotypeCountsSet))
            initFromUniqueHaploTree(inputTree, data,
                    collapseIdenticalSequencesInput.get(), collapseSequencesWithMissingDataInput.get(),
                    haplotypeCountsSet);
        else
            initFromFullTree(inputTree, data,
                    collapseIdenticalSequencesInput.get(),collapseSequencesWithMissingDataInput.get());

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
