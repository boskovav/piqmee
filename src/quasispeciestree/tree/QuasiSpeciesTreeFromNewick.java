package quasispeciestree.tree;

import beast.core.*;
import beast.core.Input.Validate;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

import java.util.List;

/**
 * @author Veronika Boskova created on 16/06/2015.
 */

@Description("Class to initialize a QuasiSpeciesTree from newick tree format with quasispecies count trait set")
public class QuasiSpeciesTreeFromNewick extends QuasiSpeciesTree implements StateNodeInitialiser{

    public Input<String> newickStringInput = new Input<>("newick",
            "Tree in Newick format.", Validate.REQUIRED);
    public Input<Boolean> adjustTipHeightsInput = new Input<>("adjustTipHeights",
            "Adjust tip heights in tree? Default true.", true);
//    public Input<TraitSet> quasiSpeciesCountsInput =
//            new Input<TraitSet>("quasiSpeciesCounts","QS counts for each haplotype (excluding the one representative of each haplotype in the tree input)", Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() throws Exception {
// gives NULL pointer warning here... what can I do --> nothing this is normal behavior!!!

        TreeParser inputTree = new TreeParser();
        inputTree.initByName(
                "IsLabelledNewick", true,
                "adjustTipHeights", adjustTipHeightsInput.get(),
                "newick", newickStringInput.get());
        Tree regularTree = inputTree;
        initFromRegularTree(regularTree);

// this was used when QuasiSpeciesNode class was not implemented yet
//        root = inputTree.root.copy();
//        nodeCount = inputTree.nodeCount;
//        internalNodeCount = inputTree.internalNodeCount;
//        leafNodeCount = inputTree.leafNodeCount;

        super.initAndValidate();
    }

    @Override
    public void initStateNodes() throws Exception {

    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(this);
    }

}
