package piqmee.tree;

import beast.core.*;
import beast.core.Input.Validate;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.TraitSet;
import beast.util.TreeParser;

import java.util.List;

/**
 * @author Veronika Boskova created on 16/06/2015.
 */

@Description("Class to initialize a QuasiSpeciesTree from newick tree format with quasispecies count trait set")
public class QuasiSpeciesTreeFromNewick extends QuasiSpeciesTree implements StateNodeInitialiser {

    public Input<String> newickStringInput = new Input<>("newick",
            "Tree in Newick format.", Validate.REQUIRED);
    public Input<Boolean> adjustTipHeightsInput = new Input<>("adjustTipHeights",
            "Adjust tip heights in tree? Default true.", true);
    public Input<Boolean> collapseIdenticalSequencesInput = new Input<>("collapseIdenticalSequences",
            "Should nodes that have identical sequences be collapsed to one haplotype? " +
                    "Default true.", true);

    public QuasiSpeciesTreeFromNewick() {
    }

    @Override
    public void initAndValidate(){
        super.initAndValidate();

        // make sure to use date and haploCount traits
        if (m_initial.get() != null)
            processTraits(m_initial.get().m_traitList.get());
        else
            processTraits(m_traitList.get());

        // Ensure tree is compatible with traits.
        if (hasDateTrait())
            adjustTreeNodeHeights(root);

        // initialize the tree
        TreeParser inputTree = new TreeParser();
        if (this.getDateTrait()!=null) {
            TraitSet times = this.getDateTrait();
            inputTree.initByName("IsLabelledNewick", true, "adjustTipHeights", adjustTipHeightsInput.get(), "newick", newickStringInput.get(), "trait", times);
        }
        else {
            inputTree.initByName("IsLabelledNewick", true, "adjustTipHeights", adjustTipHeightsInput.get(), "newick", newickStringInput.get());

        }

        if (dataInput.get() == null)
            throw new RuntimeException("The data input needs to be specified");

        // When specifying the input tree with newick tree, duplicate counts input NEEDS to be specified
        if (haplotypeCountsSet == null || haplotypeCountIsAll1(haplotypeCountsSet)){
            throw new RuntimeException("The haplotypeCounts input was not specified. " +
                    "This is not the proper class to initiate such tree. Use QuasiSpeciesTreeFromFullNewick.");
        }

        initFromUniqueHaploTree(inputTree, dataInput.get(),collapseIdenticalSequencesInput.get(),haplotypeCountsSet);

        initStateNodes();
    }

    @Override
    public void initStateNodes(){
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
