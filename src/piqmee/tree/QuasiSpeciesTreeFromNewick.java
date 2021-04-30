package piqmee.tree;

import beast.core.*;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.util.TreeParser;

import java.util.*;

/**
 * @author Veronika Boskova created on 16/06/2015.
 */

@Description("Class to initialize a QuasiSpeciesTree from newick tree format")
public class QuasiSpeciesTreeFromNewick extends QuasiSpeciesTree implements StateNodeInitialiser {

    public Input<String> newickStringInput = new Input<>("newick",
            "Tree in Newick format.", Input.Validate.REQUIRED);
    public Input<Boolean> adjustTipHeightsInput = new Input<>("adjustTipHeights",
            "Adjust tip heights in tree? Default true.", true);
    public Input<Boolean> collapseIdenticalSequencesInput = new Input<>("collapseIdenticalSequences",
            "Should nodes that have identical sequences be collapsed to one haplotype? " +
                    "Default true.", true);
    public Input<Boolean> collapseSequencesWithMissingDataInput = new Input<>("collapseSequencesIfIdenticalUpToMissingParts",
            "Flag to indicate if sequences that have missing data (stretches of N's) should" +
                    "be collapsed with a sequence that is identical to it up the missing data. Default false.",
            false);

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

        // read in the user input tree
        TreeParser inputTree = new TreeParser();
        inputTree.setDateTrait(timeTraitSet);
        inputTree.initByName(
                "IsLabelledNewick", true,
                "adjustTipHeights", adjustTipHeightsInput.get(),
                "newick", newickStringInput.get());

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
