package quasispeciestree.tree;

import beast.core.*;
import beast.core.Input.Validate;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

import java.util.List;

/**
 * @author Veronika Boskova created on 02/03/2017.
 */

@Description("Class to initialize a QuasiSpeciesTree from the full tree in newick tree format")
public class QuasiSpeciesTreeFromFullNewick extends QuasiSpeciesTree implements StateNodeInitialiser{

    final public Input<Alignment> dataInput = new Input<>("data",
            "Alignment data used for calculating distances for clustering",
            Input.Validate.REQUIRED);
    public Input<String> newickStringInput = new Input<>("newick",
            "Tree in Newick format.", Input.Validate.REQUIRED);
    public Input<Boolean> adjustTipHeightsInput = new Input<>("adjustTipHeights",
            "Adjust tip heights in tree? Default true.", true);

    public QuasiSpeciesTreeFromFullNewick() {

        // When specifying the input tree with full newick tree, we do not allow for duplicate counts input on the top.
        haplotypeCountsInput.setRule(Input.Validate.FORBIDDEN);

    }

    public void initAndValidate() {
        super.initAndValidate();

        TreeParser inputTree = new TreeParser();
        inputTree.initByName(
                "IsLabelledNewick", true,
                "adjustTipHeights", adjustTipHeightsInput.get(),
                "newick", newickStringInput.get());

        if (dataInput.get() == null)
            throw new RuntimeException("The data input needs to be specified");

        initFromFullTree(inputTree,dataInput.get());

    }

    @Override
    public void initStateNodes(){
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.add(this);
    }
}
