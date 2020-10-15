package piqmee.tree;

import beast.app.beauti.BeautiDoc;
import beast.core.*;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.distance.Distance;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import piqmee.distance.DifferenceCount;
import beast.evolution.likelihood.GenericTreeLikelihood;

import java.io.PrintStream;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

/**
 * @author Veronika Boskova created on 29/04/2015
 */

public class QuasiSpeciesTree extends Tree {

    final public Input<Alignment> dataInput = new Input<>("data",
            "Alignment data used for calculating distances for clustering",
            Input.Validate.REQUIRED);
    public Input<TraitSet> haplotypeCountsInput =
            new Input<TraitSet>("haplotypeCounts","Count of sequences for each haplotype (including the one representative of each haplotype in the tree input)");//,
           // Input.Validate.REQUIRED);

    protected TraitSet haplotypeCountsSet;
    protected Map<String,Integer> haplotypeCounts;
    protected String qsLabel = "qscounts";

    public QuasiSpeciesTree() { }

    public QuasiSpeciesTree(Node rootNode) {

        if (!(rootNode instanceof QuasiSpeciesNode))
            throw new IllegalArgumentException("Attempted to instantiate "
                    + "quasi-species tree with regular root node.");

        setRoot((QuasiSpeciesNode) rootNode);
        initArrays();
    }

    // init and validate from scratch in order to implement the quasi-species node -- holding the haplotype starting above
    @Override
    public void initAndValidate() {

        if (m_initial.get() != null && !(this instanceof StateNodeInitialiser)) {

            if (!(m_initial.get() instanceof QuasiSpeciesTree)) {
                throw new IllegalArgumentException("Attempted to initialise "
                        + "quasi-species tree with regular tree object.");
            }

            QuasiSpeciesTree other = (QuasiSpeciesTree) m_initial.get();
            root = ((QuasiSpeciesNode)other.root).copy();
            nodeCount = other.nodeCount;
            internalNodeCount = other.internalNodeCount;
            leafNodeCount = other.leafNodeCount;
        }

        if (nodeCount < 0) {
            if (m_taxonset.get() != null) {
                // make a caterpillar
                List<String> sTaxa = m_taxonset.get().asStringList();
                Node left = new QuasiSpeciesNode();
                left.setNr(0);
                left.setHeight(0);
                left.setID(sTaxa.get(0));
                for (int i = 1; i < sTaxa.size(); i++) {
                    Node right = new QuasiSpeciesNode();
                    right.setNr(i);
                    right.setHeight(0);
                    right.setID(sTaxa.get(i));
                    Node parent = new QuasiSpeciesNode();
                    parent.setNr(sTaxa.size() + i - 1);
                    parent.setHeight(i);
                    left.setParent(parent);
                    parent.setLeft(left);
                    right.setParent(parent);
                    parent.setRight(right);
                    left = parent;
                }
                root = left;
                leafNodeCount = sTaxa.size();
                nodeCount = leafNodeCount * 2 - 1;
                internalNodeCount = leafNodeCount - 1;

            } else {
                // make dummy tree with a single root node
                root = new QuasiSpeciesNode();
                root.setNr(0);
                ((QuasiSpeciesNode) root).setqsTree(this);
                nodeCount = 1;
                internalNodeCount = 0;
                leafNodeCount = 1;
            }
        }

        if (nodeCount >= 0) {
            initArrays();
        }

        haplotypeCounts = new HashMap<>();

        processTraits(m_traitList.get());

        // Ensure tree is compatible with traits.
        if (hasDateTrait())
            adjustTreeNodeHeights(root);

        if (dataInput.get() == null)
            throw new RuntimeException("The data input needs to be specified");
    }

    /*
    //
    //
    //          OWN FUNCTIONS
    //
    //
    */

    @Override
    protected void processTraits(List<TraitSet> traitList) {
        super.processTraits(traitList);

        // Record trait set associated with leaf types.
        for (TraitSet traitSet : traitList) {
            if (traitSet.getTraitName().equals(qsLabel)) {
                haplotypeCountsSet = traitSet;
                break;
            }
        }

        // Use explicitly-identified type trait set if available.
        // Seems dumb, but needed for BEAUti as ListInputEditors
        // muck things up...

        if (haplotypeCountsInput.get() != null) {
            haplotypeCountsSet = haplotypeCountsInput.get();
            qsLabel = haplotypeCountsSet.getTraitName();
        } else if (m_initial.get() != null && ((QuasiSpeciesTree)m_initial.get()).haplotypeCountsInput.get() != null) {
            haplotypeCountsInput = ((QuasiSpeciesTree) m_initial.get()).haplotypeCountsInput;
            haplotypeCountsSet = ((QuasiSpeciesTree)m_initial.get()).haplotypeCountsInput.get();
            qsLabel = haplotypeCountsSet.getTraitName();
        } else if (m_initial.get() != null && ((QuasiSpeciesTree)m_initial.get()).haplotypeCountsSet != null) {
            haplotypeCountsSet = ((QuasiSpeciesTree)m_initial.get()).haplotypeCountsSet;
            qsLabel = haplotypeCountsSet.getTraitName();
        }

        if (haplotypeCountsSet == null) {

            if (getTaxonset() != null) {
                TraitSet dummyTraitSet = new TraitSet();

                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < getTaxonset().getTaxonCount(); i++) {
                    if (i > 0)
                        sb.append(",\n");
                    sb.append(getTaxonset().getTaxonId(i)).append("=1");
                }
                try {
                    dummyTraitSet.initByName("traitname", "qscounts", "taxa", getTaxonset(), "value", sb.toString());
                    if (getID() == null)
                        dummyTraitSet.setID("haplotypeCounts.t:dummy");
                    else
                    dummyTraitSet.setID("haplotypeCounts.t:" + BeautiDoc.parsePartition(getID()));
                    setHaplotypeCountsTrait(dummyTraitSet);
                } catch (Exception ex) {
                    System.out.println("Error setting default haplotype count trait.");
                }
            }
        }

        setHaploCounts(haplotypeCountsSet);

    }

    /**
     * Function to initiate the array list of attachment times for each haplotype in quasispecies
     */
    // for those nodes where haplotype arises, change the haploAboveNode to haplotype's (corresponding tip node) number
    //  and for nodes below up to the tip set continuingHaploName to the same haplotype's (tip)
    private void initAttachmentTimes(){
        // for those nodes where haplotype arises, change the haploAboveNode to haplotype's (corresponding tip node) number
        //  and for nodes below up to the tip set continuingHaploName to the same haplotype's (tip)
        for (Node node : this.getExternalNodes()){
            double[] attachmentTimesListOld = ((QuasiSpeciesNode) node).getAttachmentTimesList();
            // check if getNr() always returns the same >>> Node number is guaranteed not to change during an MCMC run.
            //      written in the Node class)
            // assigns to each unique haplotype an array of size haplotypeCount
            double[] tempqstimes = new double[((QuasiSpeciesNode) node).getHaplotypeCountsFromTips()];
            double[] temptiptimes = ((QuasiSpeciesNode) node).getTipTimesList();
            int[] temptiptimescount = ((QuasiSpeciesNode) node).getTipTimesCountList();

            // Assign the first entry to be the same as the first split
            // this one is marking the start of the haplotype
            if (tempqstimes.length > 1){
                // check if the tree has more than one tip
                if (this.getLeafNodeCount() > 1)
                    initAttachmentTimesHelper(tempqstimes,temptiptimes,temptiptimescount,node.getParent().getHeight(),attachmentTimesListOld);
                // or there is just one sequence possibly sampled through time
                else
//                    initAttachmentTimesHelper(tempqstimes,temptiptimes,temptiptimescount,originInput.get().getValue(),attachmentTimesListOld);
                    initAttachmentTimesHelper(tempqstimes,temptiptimes,temptiptimescount,attachmentTimesListOld[0],attachmentTimesListOld);
            }
            else
                tempqstimes[0] = node.getHeight();

            ((QuasiSpeciesNode) node).setHaploAboveName(node.getNr());
            ((QuasiSpeciesNode) node).setContinuingHaploName(node.getNr());
            // attachment time list defined as "node" height, going from present (0) to past (positive height)
            ((QuasiSpeciesNode) node).setAttachmentTimesList(tempqstimes);
            ((QuasiSpeciesNode) node).setFirstEntryAndSortAttachTimeList();
        }
    }

    /**
     * helper to initAttachmentTimes
     */
    private void initAttachmentTimesHelper(double[] tempqstimes, double[] temptiptimes, int[] temptiptimescount,
                                                double maxTime, double[] maxTimeArray) {
        int currentPosition = 0;
        for (int j = 1; j < temptiptimescount[temptiptimes.length - 1]; j++) {
            tempqstimes[currentPosition] = maxTime - (j + 1) * ((maxTime - temptiptimes[temptiptimes.length - 1]) / (1 + temptiptimescount[temptiptimes.length - 1]));
            currentPosition++;
        }
        for (int i = temptiptimes.length - 2; i >= 0; i--) {
            double currentMaxTime = maxTimeArray[maxTimeArray.length-i-2];
            tempqstimes[currentPosition] = currentMaxTime;
            currentPosition++;
            for (int j = 0; j < temptiptimescount[i] - 1; j++) {
                tempqstimes[currentPosition] = currentMaxTime - (j + 1) * ((currentMaxTime - temptiptimes[i]) / (1 + temptiptimescount[i]));
                currentPosition++;
            }
        }
    }

    /**
     * Gets the total count of duplicates for a given haplotype (node)
     *
     */
    public int getHaplotypeCounts(Node node){
        return haplotypeCounts.get(node.getID());
    }

    /**
     * Gets the total count of attachment points for all haplotypes
     *
     */
    public int getTotalAttachmentCounts(){
        int totalCount=0;
        for (Node node : this.getExternalNodes()){
            totalCount += getHaplotypeCounts(node) - 1;
        }
        return totalCount;
    }

    /**
     * Sets the parent haplo name to -1 for all tips
     *
     */
    public void clearParentHaploNames() {
        for (Node node : this.getExternalNodes()){
            ((QuasiSpeciesNode) node).setParentHaplo(-1);
        }
    }

    /**
     * Sets the continuing haplo name to -1 for all nodes
     *
     */
    public void clearContinuingHaploNames() {
        for (Node node : this.getNodesAsArray()){
            ((QuasiSpeciesNode) node).setContinuingHaploName(-1);
        }
    }

    /**
     * Sets the continuing haplo name to -1 for all nodes below given node
     *
     * @param belowThisNode node below which the continuing haplo will be set to -1
     */
    public void clearContinuingHaploNames(Node belowThisNode) {
        ((QuasiSpeciesNode) belowThisNode).setContinuingHaploName(-1);
        for (Node node : belowThisNode.getAllChildNodes()){
            ((QuasiSpeciesNode) node).setContinuingHaploName(-1);
        }
    }

    /**
     * Sets the number of haplotype counts for each haplo
     *
     * @param haploCounts new set of haplotype counts
     */
    private void setHaploCounts(TraitSet haploCounts){
        for (Node node : this.getExternalNodes()) {
//            if (((QuasiSpeciesNode) node).getAttachmentTimesList().length != (int) haploCounts.getValue(node.getID()))
//                throw new RuntimeException("QuasiSpeciesTree class: The new to be set haploCount is not " +
//                                           "equal the the number of attachment times. Why?");
           if (node.getID() != null)
               haplotypeCounts.put(node.getID(), (int) haploCounts.getValue(node.getID()));
        }
    }

    /**
     * Sets the number of haplotype counts for each haplo
     *
     * @param node node for which the haplo count is to be set
     * @param value new haplotype count
     */
    protected void setHaploCounts(Node node, int value){
        haplotypeCounts.put(node.getID(), value);
    }

    /**
     * Sets the number of haplotype counts for each haplo
     *
     * @param haploCounts new set of haplotype counts
     * @param tree a tree whose nodes are the be used for assignment
     */
    protected void setHaploCounts(TraitSet haploCounts, Tree tree){
        for (Node node : tree.getExternalNodes()) {
//            if (((QuasiSpeciesNode) node).getAttachmentTimesList().length != (int) haploCounts.getValue(node.getID()))
//                throw new RuntimeException("QuasiSpeciesTree class: The new to be set haploCount is not " +
//                                           "equal the the number of attachment times. Why?");
            haplotypeCounts.put(node.getID(), (int) haploCounts.getValue(node.getID()));
        }
    }

    /**
     * Clears the number of haplotype counts for each haplo
     *
     */
    protected void clearHaploCounts(Tree tree){
        for (Node node : tree.getExternalNodes()) {
            haplotypeCounts.remove(node.getID());
        }
    }

    /**
     * Method to count and set the possible number of start branches for each internal node
     *  at one pre-order tree pass (prevent energy & resources waste)
     *
     */
    public void countAndSetPossibleStartBranches(){
        // get the counts of possible number of start branches for each haplotype (i.e. internal node attachment branches)
        int nTips = this.getLeafNodeCount();
        for (int i = 0; i < this.getInternalNodeCount(); i++){
//            // once determined the parent haplotype find out from how many branches
//            // of the parent haplotype can it actually branch off
            QuasiSpeciesNode node = (QuasiSpeciesNode) this.getNode(nTips + i);
//            int haplo = node.getContinuingHaploName();
//            if ( haplo != -1 ){
//                QuasiSpeciesNode haploTip = (QuasiSpeciesNode) this.getNode(haplo);
//                node.setStartBranchCounts(haploTip.countPossibleAttachmentBranches(0, this.getNode(nTips + i).getHeight()));
//            } else {
                node.setStartBranchCounts(1);
//            }
        }
    }

    /**
     * Method to determine the parent haplotype for each haplotype
     *  and determine the haplotype above each internal node (if none, set to null)
     * at one pre-order tree pass (prevent energy & resources waste)
     *
     * input and return:
     *         output) Array holding for each haplotype = node (#this.getExternalNodes()) the haplotype
     *         it arises from, held in array at position determined by node.getNr()
     */
    protected void fillParentHaplo(){
        // set all the parent haplotype to -1
        clearParentHaploNames();
        clearContinuingHaploNames();

        // starting at the root of the tree, see what is the order of the haplotype, always the haplotype below,
        // starts from the haplotype above (parent)
        //
        // if there is already haplotype arising at root-origin branch, this will have to arise from NULL (-1) haplotype anyways
        // and by default we set each entry of parentHaplo array to -1
        findParentHaplo(-1, (QuasiSpeciesNode) root);
    }

    /**
     * Helper method used by fillParentHaplo to assign parent to
     * each haplotype within quasi-species. This is a pre-order traversal, meaning that
     * starting at the root of the tree, we track what is the order of the
     * haplotypes arising, the haplotype below always has to start from the haplotype above (parent)
     *
     * @param currentHaploType Node holding the name of haplotype coming from the node directly
     *                      above nextNode (disregarding possible change on the incoming branch)
     * @param nextNode Node whose incoming branch we are testing for haplotype change
     */
    public void findParentHaplo(int currentHaploType, QuasiSpeciesNode nextNode){
        // get haplotype starting at a branch directly above next node --- if no haplotype (-1), check children nodes
        // check whether there is a new haplotype arising
        if (nextNode.getHaploAboveName() !=-1 ){
            // check whether the new haplotype is on the same or different branch than that leading to the parent haplotype
            // TODO can we do kind of binary search for String here instead of for loop???
            int newHaploType = nextNode.getHaploAboveName();
            QuasiSpeciesNode thisNode = (QuasiSpeciesNode) this.getNode(newHaploType);
            // set the haplotype from which the next haplotype arises as a parent
            thisNode.setParentHaplo(currentHaploType);
            // set the continuingHaploName at corresponding internal nodes leading from nextNode to the respective child\
            nextNode.setContinuingHaploName(newHaploType);
            while ( nextNode != thisNode ) {
                thisNode.setContinuingHaploName(newHaploType);
                thisNode = (QuasiSpeciesNode) thisNode.getParent();
            }
            currentHaploType = newHaploType;
        }
        if (nextNode.getChildCount() > 0 ) {
            // apply the same getParentHaplo function to the child nodes, with updated current haplotype
            for (Node childNode : nextNode.getChildren()) {
                findParentHaplo(currentHaploType, (QuasiSpeciesNode) childNode);
            }
        }
    }

    /**
     * Function checking if the haplotype counts array contains other entry than 1
     *
     * @param counts    traitSet with haplotype counts
     * @return          returns true if all the sequences have are present in 1 copy each, i.e. if the haplotype count
     *                      array only contains 1s
     */
    public boolean haplotypeCountIsAll1(TraitSet counts){
        for (int i = 0; i < counts.taxaInput.get().asStringList().size(); i++) {
            if (counts.getValue(i) > 1)
                return false;
        }
        return(true);
    }

    /**
     * @return Haplotype counts trait set if available, null otherwise.
     */
    public TraitSet getHaplotypeCountsTrait() {
        if (!traitsProcessed)
            processTraits(m_traitList.get());

        return haplotypeCountsSet;
    }

    /**
     * Determine whether tree has a haplotype counts trait set associated with it.
     *
     * @return true if so
     */
    public boolean hasHaplotypeCountsTrait() {
        return getHaplotypeCountsTrait() != null;
    }


    /**
     * Specifically set the haplotype counts trait set for this tree. A null value simply
     * removes the existing trait set. -- needed for BEAUti to correctly make XML
     *
     * @param traitSet
     */
    public void setHaplotypeCountsTrait(TraitSet traitSet) {
        if (hasHaplotypeCountsTrait()) {
            m_traitList.get().remove(haplotypeCountsSet);
        }

        if (traitSet != null)
            m_traitList.get().add(traitSet);

        haplotypeCountsSet = traitSet;
    }



    /*
    //
    //
    //          FUNCTIONS COPIED AND ADAPTED FROM MULTI-TYPE TREE CLASS
    //
    //
    */

    /**
     * Initiate node and storedNodes arrays
     *
     */
    @Override
    protected final void initArrays() {
        // initialise tree-as-array representation + its stored variant
        m_nodes = new QuasiSpeciesNode[nodeCount];
        listNodes((QuasiSpeciesNode)root, (QuasiSpeciesNode[]) m_nodes);
        m_storedNodes = new QuasiSpeciesNode[nodeCount];
        Node copy = root.copy();
        listNodes((QuasiSpeciesNode)copy, (QuasiSpeciesNode[])m_storedNodes);
    }

    /**
     * Convert quasi-species tree to array representation.
     *
     * @param node Root of sub-tree to convert.
     * @param nodes Array to populate with tree nodes.
     */
    private void listNodes(QuasiSpeciesNode node, QuasiSpeciesNode[] nodes) {
        nodes[node.getNr()] = node;
        node.setqsTree(this);
        if (!node.isLeaf()) {
            listNodes((QuasiSpeciesNode) node.getLeft(), nodes);
            if (node.getRight() != null)
                listNodes((QuasiSpeciesNode) node.getRight(), nodes);
        }
    }

    /**
     * Deep copy, returns a completely new quasi-species tree.
     *
     * @return a deep copy of this quasi-species tree
     */
    @Override
    public QuasiSpeciesTree copy() {
        QuasiSpeciesTree tree = new QuasiSpeciesTree();
        tree.ID = ID;
        tree.index = index;
        tree.root = ((QuasiSpeciesNode)root).copy();
        tree.nodeCount = nodeCount;
        tree.internalNodeCount = internalNodeCount;
        tree.leafNodeCount = leafNodeCount;
        tree.qsLabel = qsLabel;
        tree.haplotypeCounts = haplotypeCounts;
        return tree;
    }

    /**
     * Copy all values from an existing quasi-species tree.
     *
     * @param other
     */
    @Override
    public void assignFrom(StateNode other) {
        QuasiSpeciesTree qsTree = (QuasiSpeciesTree) other;

        QuasiSpeciesNode[] qsNodes = new QuasiSpeciesNode[qsTree.getNodeCount()];
        for (int i = 0; i < qsTree.getNodeCount(); i++)
            qsNodes[i] = new QuasiSpeciesNode();

        ID = qsTree.ID;
        root = qsNodes[qsTree.root.getNr()];
        ((QuasiSpeciesNode) root).assignFrom(qsNodes, qsTree.root);
        root.setParent(null);

        nodeCount = qsTree.nodeCount;
        internalNodeCount = qsTree.internalNodeCount;
        leafNodeCount = qsTree.leafNodeCount;
        if (qsTree.haplotypeCounts != null)
            haplotypeCounts = qsTree.haplotypeCounts;
        initArrays();
    }

    /**
     * Copy all values aside from IDs from an existing quasi-species tree.
     *  Important for restoration from state file
     *  This function is called after initFromFlatTree so we need to copy EVERYTHING
     *
     * @param other
     */
    @Override
    public void assignFromFragile(StateNode other) {
        QuasiSpeciesTree qsTree = (QuasiSpeciesTree) other;

        haplotypeCounts = qsTree.haplotypeCounts;

        if (m_nodes == null) {
            initArrays();
        }
        root = m_nodes[qsTree.root.getNr()];
        Node[] otherNodes = qsTree.m_nodes;
        int iRoot = root.getNr();
        assignFromFragileHelper(0, iRoot, otherNodes);
        root.setHeight(otherNodes[iRoot].getHeight());
        root.setParent(null);

        QuasiSpeciesNode qsRoot = (QuasiSpeciesNode)root;
        qsRoot.setHaploAboveName(((QuasiSpeciesNode)(otherNodes[iRoot])).getHaploAboveName());
        qsRoot.setContinuingHaploName(((QuasiSpeciesNode)(otherNodes[iRoot])).getContinuingHaploName());
        //qsRoot.setStartBranchCounts(((QuasiSpeciesNode)(otherNodes[iRoot])).getStartBranchCounts());
        qsRoot.setAttachmentTimesList(((QuasiSpeciesNode)(otherNodes[iRoot])).getAttachmentTimesList());
        qsRoot.setTipTimesList(((QuasiSpeciesNode)(otherNodes[iRoot])).getTipTimesList());
        qsRoot.setTipTimesCountList(((QuasiSpeciesNode)(otherNodes[iRoot])).getTipTimesCountList());
        qsRoot.setParentHaplo(((QuasiSpeciesNode)(otherNodes[iRoot])).getParentHaplo());

        if (otherNodes[iRoot].getLeft() != null) {
            root.setLeft(m_nodes[otherNodes[iRoot].getLeft().getNr()]);
        } else {
            root.setLeft(null);
        }
        if (otherNodes[iRoot].getRight() != null) {
            root.setRight(m_nodes[otherNodes[iRoot].getRight().getNr()]);
        } else {
            root.setRight(null);
        }
        // the second assignFromFragile Helper is called because root may not be the node with
        //  the highest number but needs special treatment as opposed to the rest of the nodes
        assignFromFragileHelper(iRoot + 1, nodeCount, otherNodes);

        //make sure to correctly assign the haplotypeCounts array
        haplotypeCounts.clear();

        for (Node node : this.getExternalNodes()){
            setHaploCounts(node, ((QuasiSpeciesNode) node).getHaplotypeCountsFromTips());
        }

        initArrays();
    }

    /**
     * helper to assignFromFragile *
     */
    private void assignFromFragileHelper(int iStart, int iEnd, Node[] otherNodes) {
        for (int i = iStart; i < iEnd; i++) {
            QuasiSpeciesNode sink = (QuasiSpeciesNode) m_nodes[i];
            QuasiSpeciesNode src = (QuasiSpeciesNode) otherNodes[i];
            sink.setHeight(src.getHeight());
            sink.setParent(m_nodes[src.getParent().getNr()]);

            sink.setHaploAboveName(src.getHaploAboveName());
            sink.setContinuingHaploName(src.getContinuingHaploName());
            //sink.setStartBranchCounts(src.getStartBranchCounts());
            sink.setAttachmentTimesList(src.getAttachmentTimesList());
            sink.setTipTimesList(src.getTipTimesList());
            sink.setTipTimesCountList(src.getTipTimesCountList());
            sink.setParentHaplo(src.getParentHaplo());

            if (src.getLeft() != null) {
                sink.setLeft(m_nodes[src.getLeft().getNr()]);
                if (src.getRight() != null) {
                    sink.setRight(m_nodes[src.getRight().getNr()]);
                } else {
                    sink.setRight(null);
                }
            }
        }
    }

    /**
     * Generates a new tree in which the duplicates of each haplotypes are represented
     * as individual tips.
     *
     * This method is useful for logging the quasi-speces trees into state file
     * NOTICE: More useful than getFullTree is the getFlattenedTree
     *         with attachment times and tip times as metadata
     *
     * @return Regular tree.
     */
    public Tree getFullTree() {

        // Create new tree to modify.  Note that copy() doesn't
        // initialise the node array lists, so initArrays() must
        // be called manually.
        QuasiSpeciesTree regularTree = copy();
        regularTree.initArrays();

        int additionalNodes = regularTree.getTotalAttachmentCounts();
        regularTree.nodeCount = regularTree.getNodeCount()+(additionalNodes*2);
        regularTree.leafNodeCount = (regularTree.nodeCount+1)/2;
        regularTree.internalNodeCount = regularTree.nodeCount/2-1;

        Node[] nodestemp = regularTree.m_nodes.clone();
        regularTree.m_nodes = new Node[regularTree.nodeCount];
        for (Node node : nodestemp){
            regularTree.m_nodes[node.getNr()] = node;
        }

        int nextNodeNr = getNodeCount();
        Node nextNode;
        Node nextTip;
        Node maxNode = regularTree.getRoot();

        for (Node node : getExternalNodes()) {
            QuasiSpeciesNode qsNode = (QuasiSpeciesNode) node;
            int nodeNr = node.getNr();
            Node nodeBelow = regularTree.getNode(nodeNr);
            Node nodeAbove = nodeBelow.getParent();
            boolean left = false;
            if (nodeAbove.getLeft() == nodeBelow) left = true;
            double[] tempqstimes = qsNode.getAttachmentTimesList();
            double[] temptiptimes = qsNode.getTipTimesList();
            int[] temptiptimescount = qsNode.getTipTimesCountList();
            // the first element of the TipTimesList is the last tip sampling time
            //  so the count for that entry in TipTimesCountList is # duplicates at that time
            //  plus 1 corresponding to the tip time of the node in the actual tree
            int currentTipTimePosition = temptiptimes.length-1;
            int currentTipArrayPosition = 0;
            for (int i = tempqstimes.length-1; i>0; i--) {

                // Check that the next attachment time is not above the current nodeAbove
                if (nodeAbove != null && tempqstimes[i] > nodeAbove.getHeight()){
                    nodeBelow = nodeAbove;
                    nodeAbove = nodeAbove.getParent();
                    // Check if nodeBelow is left or right child
                    if (nodeAbove != null) {
                        left = (nodeAbove.getLeft() == nodeBelow);
                    }
                }

                // Create and label new node:
                nextNode = new QuasiSpeciesNode();
                nextNode.setNr(nextNodeNr);
                nextNode.setID(String.valueOf(nextNodeNr));
                regularTree.m_nodes[nextNodeNr] = nextNode;
                nextNodeNr += 1;

                // Connect to child and parent:
                nextNode.addChild(nodeBelow);
                if (nodeAbove != null) {
                    if (left)
                        nodeAbove.setLeft(nextNode);
                    else
                        nodeAbove.setRight(nextNode);
                }
                nextNode.setParent(nodeAbove);

                // Ensure height is set:
                nextNode.setHeight(tempqstimes[i]);

                // Create and label new tip:
                nextTip = new QuasiSpeciesNode();
                nextTip.setNr(nextNodeNr);
                nextTip.setID(String.valueOf(nextNodeNr));
                regularTree.m_nodes[nextNodeNr] = nextTip;
                nextNodeNr += 1;

                // Connect to parent:
                nextNode.addChild(nextTip);

                // Ensure height is set:
                nextTip.setHeight(temptiptimes[currentTipTimePosition]);
                currentTipArrayPosition += 1;

                if (currentTipArrayPosition == temptiptimescount[currentTipTimePosition]) {
                    currentTipArrayPosition = 0;
                    currentTipTimePosition -= 1;
                }

                // Adjust variables
                nodeBelow = nextNode;
                if (regularTree.getRoot().getHeight() < nextNode.getHeight())
                    maxNode = nextNode;
            }
        }
        regularTree.setRoot((QuasiSpeciesNode) maxNode);
        return regularTree;
    }

    /**
     * Helper function for getFullTree
     */
    public void setRoot(final QuasiSpeciesNode root){
        super.setRoot(root);
        this.root = root;
    }

    /**
     *
     * Generates a new tree in which the duplicates' attachment times and tip times/counts
     * of each haplotypes are represented as metadata of the respective tips.
     *
     * This method is useful for logging the quasi-speces trees into state file
     *
     * @return Flattened tree.
     */
    protected Tree getFlattenedHaploTree() {
        // Create new tree to modify.  Note that copy() doesn't
        // initialise the node array lists, so initArrays() must
        // be called manually.
        QuasiSpeciesTree flatTree = copy();
        flatTree.initArrays();

        for (Node node : getExternalNodes()) {
            QuasiSpeciesNode qsNode = (QuasiSpeciesNode) node;
            int nodeNr = node.getNr();
            Node thisNode = flatTree.getNode(nodeNr);
            double[] tempqstimes = qsNode.getAttachmentTimesList();
            double[] temptiptimes = qsNode.getTipTimesList();
            int[] temptiptimescount = qsNode.getTipTimesCountList();
            thisNode.setMetaData("AttachTimes", Arrays.toString(tempqstimes));
            thisNode.setMetaData("TipTimes", Arrays.toString(temptiptimes));
            thisNode.setMetaData("TipCounts", Arrays.toString(temptiptimescount));
            thisNode.metaDataString = String.format("AttachTimes={%s},TipTimes={%s},TipCounts={%s}",
                         Arrays.toString(tempqstimes), Arrays.toString(temptiptimes), Arrays.toString(temptiptimescount));
        }
        return flatTree;
    }

    /**
     * Initialise tree topology from Tree object with attachment times and tip times as metadata
     *
     * @param flatHaploTree
     */
    protected void initFromFlatHaploTree(Tree flatHaploTree) {

        // Create quasispecies tree
        QuasiSpeciesNode[] quasiSpeciesNodes = new QuasiSpeciesNode[flatHaploTree.getNodeCount()];

        List<Node> nodesToAssignNext = new ArrayList<>();
        List<Node> nodesToAssign = new ArrayList<>();
        nodesToAssign.add(flatHaploTree.getRoot());
        quasiSpeciesNodes[flatHaploTree.getRoot().getNr()] = new QuasiSpeciesNode();
        quasiSpeciesNodes[flatHaploTree.getRoot().getNr()].setHeight(flatHaploTree.getRoot().getHeight());

        QuasiSpeciesNode newRoot = quasiSpeciesNodes[flatHaploTree.getRoot().getNr()];

        while (!nodesToAssign.isEmpty()) {

            nodesToAssignNext.clear();

            for (Node node : nodesToAssign) {

                QuasiSpeciesNode thisQuasiSpeciesNode = quasiSpeciesNodes[node.getNr()];

                switch (node.getChildCount()) {
                    case 0:
                        // Leaf at base of branch
                        thisQuasiSpeciesNode.setNr(node.getNr());
                        thisQuasiSpeciesNode.setID(String.valueOf(node.getID()));

                        // Set attachmentTimesList for each tip
                        Object typeObject1 = node.getMetaData("AttachTimes");
                        if (typeObject1 instanceof Double[]) {
                            double[] tempqstimes = new double[((Double[]) typeObject1).length];
                            int count1 = 0;
                            for (Double nextEntry : (Double[]) typeObject1) {
                                double attachTime = nextEntry;
                                tempqstimes[count1] = attachTime;
                                count1++;
                            }
                            thisQuasiSpeciesNode.setAttachmentTimesList(tempqstimes);
                        } else if (typeObject1 instanceof double[]) {
                            thisQuasiSpeciesNode.setAttachmentTimesList((double[]) typeObject1);
                        } else
                            throw new IllegalArgumentException("Unrecognised type metadata.");

                        // Set tipTimesList for each tip
                        Object typeObject2 = node.getMetaData("TipTimes");
                        if (typeObject2 instanceof Double[]) {
                            double[] temptiptimes = new double[((Double[]) typeObject2).length];
                            int count2 = 0;
                            for (Double nextEntry : (Double[]) typeObject2) {
                                double tipTime = nextEntry;
                                temptiptimes[count2] = tipTime;
                                count2++;
                            }
                            thisQuasiSpeciesNode.setTipTimesList(temptiptimes);
                        } else if (typeObject2 instanceof double[]) {
                            thisQuasiSpeciesNode.setTipTimesList((double[]) typeObject2);
                        } else
                            throw new IllegalArgumentException("Unrecognised type metadata.");

                        // Set tipTimesCountList for each tip
                        Object typeObject3 = node.getMetaData("TipCounts");
                        if (typeObject3 instanceof Integer[]) {
                            int[] temptiptimescount = new int[((Integer[]) typeObject3).length];
                            int count3 = 0;
                            for (Integer nextEntry : (Integer[]) typeObject3) {
                                int tipTimeCount = nextEntry;
                                temptiptimescount[count3] = tipTimeCount;
                                count3++;
                            }
                            thisQuasiSpeciesNode.setTipTimesCountList(temptiptimescount);
                        } else if (typeObject3 instanceof int[]) {
                            thisQuasiSpeciesNode.setTipTimesCountList((int[]) typeObject3);
                        } else if (typeObject3 instanceof Double[]) {
                            int[] temptiptimescount = new int[((Double[]) typeObject3).length];
                            int count3 = 0;
                            for (Double nextEntry : (Double[]) typeObject3) {
                                int tipTimeCount = nextEntry.intValue();
                                temptiptimescount[count3] = tipTimeCount;
                                count3++;
                            }
                            thisQuasiSpeciesNode.setTipTimesCountList(temptiptimescount);
                        } else if (typeObject3 instanceof double[]) {
                            thisQuasiSpeciesNode.setTipTimesCountList((int[]) typeObject3);
                        } else
                            throw new IllegalArgumentException("Unrecognised type metadata.");

                        //attachmentTimesList and tipTimesList created on a go...sort them now
                        // do not use initAttachmentTimes();
                        thisQuasiSpeciesNode.sortAttachTimeList();
                        thisQuasiSpeciesNode.sortTipTimeAndCountList();

                        // reassign the tip heights from the tiptimes array
                        thisQuasiSpeciesNode.setHeight(thisQuasiSpeciesNode.getTipTimesList()[0]);

                        break;

                    case 2:
                        // Non-leaf at base of branch
                        quasiSpeciesNodes[node.getChild(0).getNr()] = new QuasiSpeciesNode();
                        quasiSpeciesNodes[node.getChild(0).getNr()].setHeight(node.getChild(0).getHeight());

                        quasiSpeciesNodes[node.getChild(1).getNr()] = new QuasiSpeciesNode();
                        quasiSpeciesNodes[node.getChild(1).getNr()].setHeight(node.getChild(1).getHeight());

                        thisQuasiSpeciesNode.addChild(quasiSpeciesNodes[node.getChild(0).getNr()]);
                        thisQuasiSpeciesNode.addChild(quasiSpeciesNodes[node.getChild(1).getNr()]);

                        nodesToAssignNext.add(node.getLeft());
                        nodesToAssignNext.add(node.getRight());

                        break;
                }
            }

            nodesToAssign.clear();
            nodesToAssign.addAll(nodesToAssignNext);
        }

        // Number internal nodes:
        numberInternalNodes(newRoot, newRoot.getLeafNodeCount());

        // Assign tree topology:
        assignFromWithoutID(new QuasiSpeciesTree(newRoot));
        initArrays();

        //traverse a tree and assign nodes above and continuing haplo annotations
        assignContinuingHaploAndHaploAbove();

        fillParentHaplo();

        countAndSetPossibleStartBranches();

        haplotypeCounts.clear();

        for (Node node : this.getExternalNodes()){
            setHaploCounts(node, ((QuasiSpeciesNode) node).getHaplotypeCountsFromTips());
        }
    }

    /**
     * Initialise tree topology from Tree object with duplicate counts (but no attachment times) in trait set array
     *
     * @param uniqueHaploTree
     */
    public void initFromUniqueHaploTree(Tree uniqueHaploTree, Alignment data, boolean collapseIdentical, TraitSet haplotypeCountsTrait){

        // In unique haplo tree, there can still be duplicate sequences, if found at different points in time
        // Get the distances for the sequences:
        int taxaSize = data.getTaxonCount();
        double[][] distanceMatrix = new double[taxaSize][taxaSize];
        double[][] distanceMatrixSum = new double[taxaSize][taxaSize];
        double[][] distanceMatrixTmp;
        // 1) check if there are multiple alignments linked with this tree -- such that unique sequences correctly identified
        Set<BEASTInterface> outputset;
        if (m_initial.get() != null)
            outputset = m_initial.get().getOutputs();
        else
            outputset = this.getOutputs();
        for (BEASTInterface o : outputset) {
            if (o instanceof GenericTreeLikelihood) {
                GenericTreeLikelihood likelihood = (GenericTreeLikelihood) o;
                Alignment odata = likelihood.dataInput.get();
                if (odata instanceof FilteredAlignment) {
                    odata = ((FilteredAlignment) odata).alignmentInput.get();
                }
                if (odata.getTaxaNames() == null){
                    Alignment odatatmp = new Alignment(odata.sequenceInput.get(), odata.dataTypeInput.get());
                    odata = odatatmp;
                }
                // 2) make a distance matrix for each such alignment
                distanceMatrixTmp = getSequenceDistances(odata, uniqueHaploTree);
                // 3) add this distance to distances from other alignments
                for (int i = 0; i < taxaSize - 1; i++) {
                    for (int j = i + 1; j < taxaSize; j++) {
                        distanceMatrixSum[i][j] = distanceMatrix[i][j] + distanceMatrixTmp[i][j];
                        distanceMatrixSum[j][i] = distanceMatrixSum[i][j];
                    }
                }
                // 4)copy sum of distances to distanceMatrix
                System.arraycopy(distanceMatrixSum, 0, distanceMatrix, 0, distanceMatrix.length);
            }
        }
        // 5) it could be, especially in a test case, that the tree is not linked with any output - check for this
        if (outputset.size() == 0)
            distanceMatrix = getSequenceDistances(data, uniqueHaploTree);

        // Build new quasi-species tree:
        ArrayList haplotypesSeen = new ArrayList<>();
        List<QuasiSpeciesNode> qsTips = new ArrayList<>();
        List<QuasiSpeciesNode> qsInternalNodes = new ArrayList<>();

        // set haplo count to those found in input for uniqueHaploTree
        setHaploCounts(haplotypeCountsTrait,uniqueHaploTree);

        ArrayList result = processNextNodeOfFullNewickTree(
                uniqueHaploTree.getRoot(),qsTips,qsInternalNodes,distanceMatrix,haplotypesSeen,collapseIdentical);

        // renumber tips to match the number of tips in the qsTree (so far matching fullTree node numbers)
        // need to match the tip times and attach time and haplo count lists!! -- this should not affect the order
        //  in these arrays, as the values were put in in the order of tips put in qsTips array
        QuasiSpeciesNode[] qsTipsTmp = new QuasiSpeciesNode[qsTips.size()];
        System.arraycopy(qsTips.toArray(),0,qsTipsTmp,0,qsTips.size());

        int[] labelsold = new int[qsTipsTmp.length];
        int[] labelsnew = new int[qsTipsTmp.length];
        int[] index = new int[qsTipsTmp.length];
        List<String> taxanames = m_taxonset.get().asStringList();
        for (int i = 0; i<qsTipsTmp.length; i++){
            for (int j = 0; j<taxanames.size(); j++) {
                if (qsTipsTmp[i].getID().equals(taxanames.get(j))) {
                    labelsold[i] = j;
                    labelsnew[i] = j;
                }
            }
        }
        Arrays.sort(labelsnew);
        for (int i = 0; i<labelsnew.length; i++){
            for (int j = 0; j<labelsold.length; j++){
                if (labelsnew[i]==labelsold[j]){
                    index[i]=j;
                    break;
                }
            }
        }
        for (int i = 0; i<qsTipsTmp.length; i++){
            qsTips.set(i,qsTipsTmp[index[i]]);
            qsTips.get(i).setNr(i);
        }
        QuasiSpeciesNode newRoot = (QuasiSpeciesNode) result.get(0);
        // Number internal nodes:
        numberInternalNodes(newRoot, newRoot.getLeafNodeCount());

        // Assign tree topology:
        assignFromWithoutID(new QuasiSpeciesTree(newRoot));

        // tipTimesList created on a go...
        for (Node node : getExternalNodes()){
            // sort and reverse the rest of the array to start with the largest value
            ((QuasiSpeciesNode) node).sortAttachTimeList();
            // sort tip times
            ((QuasiSpeciesNode) node).sortTipTimeAndCountList();
        }

        // Make sure to properly assign the attachmentTimes
        // Set start of haplotype times as default to belong to the leaf node
        // treeNode.setHaploname(treeNode.getID());
        // done in initAndValidate!!!
        initAttachmentTimes();

        initArrays();

        //traverse a tree and assign nodes above and continuing haplo annotations
        assignContinuingHaploAndHaploAbove();

        fillParentHaplo();

        countAndSetPossibleStartBranches();

        // remove all entries put from the uniqueHaploTree
        clearHaploCounts(uniqueHaploTree);

        for (Node node : this.getExternalNodes()){
            setHaploCounts(node, ((QuasiSpeciesNode) node).getHaplotypeCountsFromTips());
        }
    }

    /**
     * Initialise tree topology from full Tree object - have to check for duplicates in post order (children first) traversal
     *                                                  and collapse duplicates, throw error if tree does not fullfill
     *                                                  recursively monophyletic constraint of our QS model
     *
     * @param fullTree
     */
    public void initFromFullTree(Tree fullTree, Alignment data, boolean collapseIdentical){

        // Get the distances for the sequences:
        int taxaSize = data.getTaxonCount();
        double[][] distanceMatrix = new double[taxaSize][taxaSize];
        double[][] distanceMatrixSum = new double[taxaSize][taxaSize];
        double[][] distanceMatrixTmp;
        // 1) check if there are multiple alignments linked with this tree -- such that unique sequences correctly identified
        Set<BEASTInterface> outputset;
        if (m_initial.get() != null)
            outputset = m_initial.get().getOutputs();
        else
            outputset = this.getOutputs();
        for (BEASTInterface o : outputset) {
            if (o instanceof GenericTreeLikelihood) {
                GenericTreeLikelihood likelihood = (GenericTreeLikelihood) o;
                Alignment odata = likelihood.dataInput.get();
                if (odata instanceof FilteredAlignment) {
                    odata = ((FilteredAlignment) odata).alignmentInput.get();
                }
                if (odata.getTaxaNames() == null){
                    Alignment odatatmp = new Alignment(odata.sequenceInput.get(), odata.dataTypeInput.get());
                    odata = odatatmp;
                }
                // 2) make a distance matrix for each such alignment
                distanceMatrixTmp = getSequenceDistances(odata, fullTree);
                // 3) add this distance to distances from other alignments
                for (int i = 0; i < taxaSize - 1; i++) {
                    for (int j = i + 1; j < taxaSize; j++) {
                        distanceMatrixSum[i][j] = distanceMatrix[i][j] + distanceMatrixTmp[i][j];
                        distanceMatrixSum[j][i] = distanceMatrixSum[i][j];
                    }
                }
                // 4)copy sum of distances to distanceMatrix
                System.arraycopy(distanceMatrixSum, 0, distanceMatrix, 0, distanceMatrix.length);
            }
        }
        // 5) it could be, especially in a test case, that the tree is not linked with any output - check for this
        if (outputset.size() == 0)
            distanceMatrix = getSequenceDistances(data, fullTree);

        // Build new quasi-species tree:
        ArrayList haplotypesSeen = new ArrayList<>();
        List<QuasiSpeciesNode> qsTips = new ArrayList<>();
        List<QuasiSpeciesNode> qsInternalNodes = new ArrayList<>();

        ArrayList result = processNextNodeOfFullNewickTree(fullTree.getRoot(),qsTips,qsInternalNodes,
                                                            distanceMatrix,haplotypesSeen,collapseIdentical);

        // renumber tips to match the number of tips in the qsTree (so far matching fullTree node numbers)
        // need to match the tip times and attach time and haplo count lists!! -- this should not affect the order
        //  in these arrays, as the values were put in in the order of tips put in qsTips array
        QuasiSpeciesNode[] qsTipsTmp = new QuasiSpeciesNode[qsTips.size()];
        System.arraycopy(qsTips.toArray(),0,qsTipsTmp,0,qsTips.size());

        int[] labelsold = new int[qsTipsTmp.length];
        int[] labelsnew = new int[qsTipsTmp.length];
        int[] index = new int[qsTipsTmp.length];
        List<String> taxanames = m_taxonset.get().asStringList();
        for (int i = 0; i<qsTipsTmp.length; i++){
            for (int j = 0; j<taxanames.size(); j++) {
                if (qsTipsTmp[i].getID().equals(taxanames.get(j))) {
                    labelsold[i] = j;
                    labelsnew[i] = j;
                }
            }
        }
        Arrays.sort(labelsnew);
        for (int i = 0; i<labelsnew.length; i++){
            for (int j = 0; j<labelsold.length; j++){
                if (labelsnew[i]==labelsold[j]){
                    index[i]=j;
                    break;
                }
            }
        }
        for (int i = 0; i<qsTipsTmp.length; i++){
            qsTips.set(i,qsTipsTmp[index[i]]);
            qsTips.get(i).setNr(i);
        }
        QuasiSpeciesNode newRoot = (QuasiSpeciesNode) result.get(0);
        // Number internal nodes:
        numberInternalNodes(newRoot, newRoot.getLeafNodeCount());

        // Assign tree topology:
        assignFromWithoutID(new QuasiSpeciesTree(newRoot));

        //attachmentTimesList and tipTimesList created on a go...
        // do not use initAttachmentTimes();
        for (Node node : getExternalNodes()){
            // sort and reverse the rest of the array to start with the largest value
            ((QuasiSpeciesNode) node).setFirstEntryAndSortAttachTimeList();
            // sort tip times
            ((QuasiSpeciesNode) node).sortTipTimeAndCountList();
        }

        initArrays();

        //traverse a tree and assign nodes above and continuing haplo annotations
        assignContinuingHaploAndHaploAbove();

        fillParentHaplo();

        countAndSetPossibleStartBranches();

        // remove all entries put from the fullTree
        clearHaploCounts(fullTree);

        for (Node node : this.getExternalNodes()){
            setHaploCounts(node, ((QuasiSpeciesNode) node).getHaplotypeCountsFromTips());
        }
    }

    /**
     * Helper method used by initFromFullTree/initFromUniqueHaploTree to evaluate
     * which nodes are to be kept as internal nodes and to assign the attachmentTimes array.
     * This is a post-order traversal, meaning the root is given returned
     * as the last true QuasiSpeciesNode.
     *
     * @param node
     * @param qsTips
     * @param qsInternalNodes
     * @param distanceMatrix
     * @param haplotypesSeen list of taxon names that will be tips in the qsTree with unique sequences already seen
     * @return
     */
    private ArrayList processNextNodeOfFullNewickTree(
            Node node, List<QuasiSpeciesNode> qsTips, List<QuasiSpeciesNode> qsInternalNodes,
            double[][] distanceMatrix, ArrayList haplotypesSeen, boolean collapseIdentical){

        QuasiSpeciesNode returnNode = null;
        ArrayList haplotypesAtThisNode = new ArrayList();
        // fakeHaplo the haplotype number (in the qsTree); a node in the regular that is a sequence duplicate of a previously
        //      seen haplotype is called "fakeHaplo" since all the nodes that attach to this branch should attach
        //      to the branch leading to the true tip
        int fakeHaplo = -1;
        // for leaf nodes check if the sequence has been seen at another node already
        // pass on to the parent the info on which haplo is at the tip
        if (node.isLeaf()){
            boolean skip = false;
            // check if the sequence has been seen already
            if (collapseIdentical==true) {
                for (int i = 0; i < haplotypesSeen.size(); i++) {
                    if (distanceMatrix[node.getNr()][(int) haplotypesSeen.get(i)] == 0) {
                        QuasiSpeciesNode seenNode = qsTips.get(i);
                        // check if the time of the tip is less than the uniqueHaploTree tip
                        // if not, rewrite the info on the uniqueHaploTree tip
                        if ( (node.getHeight() - seenNode.getHeight()) < -1e-10) {
                            seenNode.setHeight(node.getHeight());
                        }
                        // to robustly always pick the same node as the unique haplo node even after restart from state file
                        if (!Pattern.compile("^t").matcher(seenNode.getID()).find() && Pattern.compile("^t").matcher(node.getID()).find()) {
                            seenNode.setID(String.valueOf(node.getID()));
                        }
                        // if the tips do not start with t - then keep always the tip with smallest number
                        else if (seenNode.getID().compareTo(node.getID()) > 0) {
                            seenNode.setID(String.valueOf(node.getID()));
                        }
                        // since the sequence has been seen already, assign the tip time to array
                        double[] tipTimesListTmp = seenNode.getTipTimesList();
                        int[] tipTimesCountListTmp = seenNode.getTipTimesCountList();
                        if (haplotypeCounts.size() != 0 && !haplotypeCountIsAll1(haplotypeCountsSet)) {
                            // throw an error if if tip has already been found with same seq and at the same time...
                            // the user wants to use full tree? or did not correctly merge duplicate sequences?
                            for (int j = 0; j < tipTimesListTmp.length; j++) {
                                // check if the tips with the same sequence that have been seen had also the current node's sampling time
                                if (tipTimesListTmp[j] == node.getHeight()) {
                                    throw new IllegalArgumentException(
                                            "There are at least two tips with the same sequence and same sampling time." +
                                            " Please, either input a tree with all sequences as tips, or remove duplicates" +
                                            " and use haplotypeCounts traitset to annotate the duplicate counts.");
                                }
                            }
                            // if with different time, add to tip times and counts
                            // expand the TipTimesList and add a new value
                            addNewTimesAndCountEntry(seenNode, tipTimesListTmp, tipTimesCountListTmp, node.getHeight(), getHaplotypeCounts(node));
                        } else {
                            // since we are assigning from the full tree, we need to check if we have
                            //  already observed the same time for another already processed tip
                            boolean haploSeen = false;
                            for (int j = 0; j < tipTimesListTmp.length; j++) {
                                // if yes, just increase the corresponding timecount array by one
                                if ( Math.abs(tipTimesListTmp[j] - node.getHeight()) < 1e-10) {
                                    tipTimesCountListTmp[j] += 1;
                                    haploSeen = true;
                                    break;
                                }
                            }
                            // if not create a new entry
                            if (!haploSeen) {
                                // expand the TipTimesList and add a new value
                                addNewTimesAndCountEntry(seenNode, tipTimesListTmp, tipTimesCountListTmp, node.getHeight(), 1);
                            }
                        }
                        skip = true;
                        haplotypesAtThisNode.add(haplotypesSeen.get(i));
                        // make this node to be a "fake" node
                        returnNode = null;
                        fakeHaplo = i;
                        break;
                    }
                }
            }
            if (!skip) {
                haplotypesSeen.add(node.getNr());
                haplotypesAtThisNode.add(node.getNr());
                returnNode = new QuasiSpeciesNode();
                qsTips.add(returnNode);
                returnNode.setHeight(node.getHeight());
                returnNode.setID(String.valueOf(node.getID()));
                returnNode.setNr(node.getNr());
                // create a new attachmentTimesList and tipTimesList entry and check how long it needs to be
                int newEntryLength = 0;
                double[] distances = distanceMatrix[node.getNr()].clone();
                Arrays.sort(distances);
                for (int i = 0; i < distances.length; i++){
                    if (distances[i] == 0)
                        newEntryLength += 1;
                    else
                        break;
                }
                returnNode.setAttachmentTimesList(new double[newEntryLength]);
                returnNode.setTipTimesList(new double[1]);
                returnNode.getTipTimesList()[0] = node.getHeight();
                returnNode.setTipTimesCountList(new int[1]);
                if (haplotypeCounts.size() != 0)
                    returnNode.getTipTimesCountList()[0] = getHaplotypeCounts(node);
                else
                    returnNode.getTipTimesCountList()[0] = 1;
            }
        }
        else {
            ArrayList leftOut = processNextNodeOfFullNewickTree(node.getLeft(),qsTips,qsInternalNodes,
                                                                distanceMatrix,haplotypesSeen,collapseIdentical);
            ArrayList rightOut = processNextNodeOfFullNewickTree(node.getRight(),qsTips,qsInternalNodes,
                                                                 distanceMatrix,haplotypesSeen,collapseIdentical);
            QuasiSpeciesNode leftNode = (QuasiSpeciesNode) leftOut.get(0);
            QuasiSpeciesNode rightNode = (QuasiSpeciesNode) rightOut.get(0);
            ArrayList leftHaplo = (ArrayList) leftOut.get(1);
            ArrayList rightHaplo = (ArrayList) rightOut.get(1);
            int leftFakeHaplo = (int) leftOut.get(2);
            int rightFakeHaplo = (int) rightOut.get(2);
            // Case 1:  there is more than one QS identical throw exception = this is not a QS tree
            // Case 2:  both arrays have exactly one same QS this means we found a duplicate and need to record this time
            //            into attachmentTimesList
            // Case 3:  this is a fake internal node where several fake haplo meet -- this option is not implemented yet - throw error
            // Case 4:  this is a real QS tree node where several haplotypes meet -- check if any of haplo has been seen already
            //            if it has, assign the internal node to the existing branch of that haplo
            ArrayList sameHaploLeftRight = new ArrayList();
            for (int i = 0; i < leftHaplo.size(); i++){
                for (int j = 0; j < rightHaplo.size(); j++){
                    if ((int) leftHaplo.get(i) == (int) rightHaplo.get(j)) {
                        sameHaploLeftRight.add(leftHaplo.get(i));
                    }
                }
            }
            //case 1
            if (sameHaploLeftRight.size() > 1 && collapseIdentical){
                throw new IllegalArgumentException("The input tree is not recursively monophyletic and therefore" +
                                                   " cannot be converted to a quasi-species tree. Try to input a" +
                                                   " different tree. Alternatively, input sequences only.");
            }
            //case 2
            else if(sameHaploLeftRight.size() == 1 && collapseIdentical){
                // since the sequence has been seen already, assign the attachment haplo time to array
                int haplo = (int) sameHaploLeftRight.get(0);
                QuasiSpeciesNode tip = null;
                for (int i = 0; i < qsTips.size(); i++){
                    if (qsTips.get(i).getNr() == haplo){
                        tip = qsTips.get(i);
                        break;
                    }
                }
                // set the next entry in the attach time array to the currently processed internal node height
                double[] attachmentTimesListTmp = tip.getAttachmentTimesList();
                for (int i = 0; i < attachmentTimesListTmp.length; i++) {
                    if (attachmentTimesListTmp[i] == 0) {
                        attachmentTimesListTmp[i] = node.getHeight();
                        break;
                    }
                }
                haplotypesAtThisNode = leftHaplo;
                haplotypesAtThisNode.addAll(haplotypesAtThisNode.size(),rightHaplo);
                // return the first index of the recurring element
                int index1 = haplotypesAtThisNode.indexOf(haplo);
                haplotypesAtThisNode.remove(index1);
                //check what is the returnNode : will be null if both nodes are fake
                //  otherwise the higher of the two QS nodes from left and right
                //  the lower node is has already been appended to the branch leading to haplo (see case 3)
                if (leftNode == null && rightNode == null)
                    returnNode = null;
                else if (rightNode == null)
                    returnNode = leftNode;
                else if (leftNode == null)
                    returnNode = rightNode;
                else if (leftNode.getHeight() > rightNode.getHeight())
                    returnNode = leftNode;
                else
                    returnNode = rightNode;
                // if both nodes are fake, set fakeNode to that haplo
                if ((leftFakeHaplo == -1 && qsTips.get(rightFakeHaplo).getNr() == haplo)
                        || (qsTips.get(leftFakeHaplo).getNr() == haplo && rightFakeHaplo == -1))
                    fakeHaplo=-1;
                else if (qsTips.get(leftFakeHaplo).getNr() == haplo && qsTips.get(rightFakeHaplo).getNr() == haplo)
                    fakeHaplo=leftFakeHaplo;
                else if (leftFakeHaplo != -1 && rightFakeHaplo != -1){
                    throw new IllegalArgumentException("There is something seriously wrong with the tree or with this algorithm.");
                }
                else {
                    throw new IllegalArgumentException("What case have I forgotten?");
                }
            }
            //case 3
            else if (leftFakeHaplo != -1 && rightFakeHaplo != -1 && collapseIdentical){
                throw new IllegalArgumentException("The input tree is not recursively monophyletic and therefore" +
                        " cannot be converted to a quasi-species tree. Try to input a" +
                        " different tree. Alternatively, input sequences only.");
            }
            //case 4
            else {
                haplotypesAtThisNode = leftHaplo;
                haplotypesAtThisNode.addAll(haplotypesAtThisNode.size(),rightHaplo);
                // look for fake haplotype (i.e. one of the nodes below this node are not qsTips
                QuasiSpeciesNode nodeToPlace = null;
                // if we have a fake node, connect the left-right haplotypes at a new node corresponding to fake haplo
                if(leftFakeHaplo != -1 || rightFakeHaplo != -1){
                    //find the TRUE tip belonging to the found haplo and the node that should be the parent of the new node
                    QuasiSpeciesNode checkThisNode = null;
                    QuasiSpeciesNode checkThisNodeKid = null;
                    if (leftFakeHaplo != -1){
                        checkThisNode = qsTips.get(leftFakeHaplo);
                        nodeToPlace = rightNode;
                        fakeHaplo = leftFakeHaplo;
                    }
                    else if (rightFakeHaplo != -1){
                        checkThisNode = qsTips.get(rightFakeHaplo);
                        nodeToPlace = leftNode;
                        fakeHaplo = rightFakeHaplo;
                    }

                    // paste all the events on the branch leading to the tip
                    while (checkThisNode != null && checkThisNode.getHeight() < node.getHeight()) {
                        checkThisNodeKid = checkThisNode;
                        checkThisNode = (QuasiSpeciesNode) checkThisNode.getParent();
                    }
                    QuasiSpeciesNode newInternalNode = new QuasiSpeciesNode();
                    qsInternalNodes.add(newInternalNode);
                    newInternalNode.addChild(checkThisNodeKid);
                    newInternalNode.addChild(nodeToPlace);
                    newInternalNode.setHeight(node.getHeight());
                    if (checkThisNode != null) {
                        if (checkThisNode.getLeft().getNr() == checkThisNodeKid.getNr()
                                && checkThisNode.getLeft().getHeight() == checkThisNodeKid.getHeight())
                            checkThisNode.setLeft(newInternalNode);
                        else
                            checkThisNode.setRight(newInternalNode);
                        newInternalNode.setParent(checkThisNode);
                    }
                    else
                        returnNode = newInternalNode;
                    // correct the annotations of haplotypesAtThisNode from new node up
                    //      actually do not correct these, this would mess things up later on, as
                    //      then several haplotypes could meet at one node suddenly
                }
                // otherwise connect the two sides at the current branches
                else{
                    returnNode = new QuasiSpeciesNode();
                    qsInternalNodes.add(returnNode);
                    returnNode.addChild(leftNode);
                    returnNode.addChild(rightNode);
                    returnNode.setHeight(node.getHeight());
                }
            }
        }

        ArrayList output = new ArrayList();
        output.add(returnNode);
        output.add(haplotypesAtThisNode);
        output.add(fakeHaplo);
        return output;
    }

    /**
     * Helper method used by processNextNodeOfFullNewickTree to
     * expand the tipTimes and tipTimesCount arrays by one
     */
    protected void addNewTimesAndCountEntry(QuasiSpeciesNode seenNode, double[] tipTimesListTmp,
                                          int[] tipTimesCountListTmp, double newHeight, int newCount) {
        // expand the TipTimesList and add a new value
        double[] tipTimesTempArray = tipTimesListTmp;
        tipTimesListTmp = new double[tipTimesTempArray.length+1];
        System.arraycopy(tipTimesTempArray,0,tipTimesListTmp,0,tipTimesTempArray.length);
        tipTimesListTmp[tipTimesListTmp.length-1] = newHeight;
        seenNode.setTipTimesList(tipTimesListTmp);

        int[] tipTimesCountTempArray = tipTimesCountListTmp;
        tipTimesCountListTmp = new int[tipTimesCountTempArray.length+1];
        System.arraycopy(tipTimesCountTempArray,0,tipTimesCountListTmp,0,tipTimesCountTempArray.length);
        tipTimesCountListTmp[tipTimesCountListTmp.length-1] = newCount;
        seenNode.setTipTimesCountList(tipTimesCountListTmp);
    }

    /**
     * Helper method used by initFromFullTree/initFromUniqueHaploTree/initFromFlatHaploTree to
     * assign to each node in the tree continuingHaplo and haploAboveName.
     */
    private void assignContinuingHaploAndHaploAbove() {
        //traverse a tree and assign nodes above and continuing haplo annotations
        for (Node thisNode : this.getExternalNodes()){
            int haplo = thisNode.getNr();
            double maxHaploTime = ((QuasiSpeciesNode) thisNode).getAttachmentTimesList()[0];
            while (thisNode.getParent() != null && thisNode.getParent().getHeight() < maxHaploTime){
                ((QuasiSpeciesNode) thisNode).setContinuingHaploName(haplo);
                ((QuasiSpeciesNode) thisNode).setHaploAboveName(-1);
                thisNode = thisNode.getParent();
            }
            ((QuasiSpeciesNode) thisNode).setHaploAboveName(haplo);
            ((QuasiSpeciesNode) thisNode).setContinuingHaploName(haplo);
        }
    }

    /**
     * Helper method used by initFromFullTree/initFromUniqueHaploTree to
     * calculate pairwise sequence distances and return the matrix of them.
     *
     * @param data
     * @param tree
     * @return
     */
    protected double[][] getSequenceDistances(Alignment data, Tree tree) {
    	// collect unique sequences
    	Map<String, List<Integer>> sequenceMap = new HashMap<>();
    	for (Sequence seq : data.sequenceInput.get()) {
    		String sequence = seq.dataInput.get();
    		if (!sequenceMap.containsKey(sequence)) {
    			sequenceMap.put(sequence, new ArrayList<>());
    		}
    		sequenceMap.get(sequence).add(data.getTaxonIndex(seq.getTaxon()));
    	}
    	
    	
        // Get the distances for the sequences:
        Distance distance = new DifferenceCount();
        ((Distance.Base) distance).setPatterns(data);
        // make a tip sequence distance matrix
        int taxaSize = data.getTaxonCount();
        double[][] distanceMatrix = new double[taxaSize][taxaSize];

        String [] uniuqueSequences = sequenceMap.keySet().toArray(new String[]{});
        int n = uniuqueSequences.length;

        Log.info("Found " + n + " unique sequences out of " + taxaSize + " sequences");
        
        if (n > 1000) {
        	Log.warning("\nWARNING: with " + n + " unique sequences you might consider sub-sampling\n");
        }
        
        Log.warning.print("Prepping distance matrix");        
        for (int i = 0; i < n - 1; i++) {
        	List<Integer> taxa1 = sequenceMap.get(uniuqueSequences[i]);
            for (int j = i + 1; j < n; j++) {
            	List<Integer> taxa2 = sequenceMap.get(uniuqueSequences[j]);
             	double dist = distance.pairwiseDistance(taxa1.get(0), taxa2.get(0));
            	for (int k = 0; k < taxa1.size(); k++) {
            		int taxon1 = taxa1.get(k);
            		for (int m = 0; m < taxa2.size(); m++) {
                		int taxon2 = taxa2.get(m);
                		distanceMatrix[taxon1][taxon2] = dist;
                		distanceMatrix[taxon2][taxon1] = dist;
            		}
            	}
            }
            if (i % 100 == 0) {
            	Log.warning.print(".");
            }
        }
        Log.warning.println("Done.");
    	
        return distanceMatrix;
    }

    /**
     * Helper method used by initFromFullTree/initFromUniqueHaploTree to assign sensible node numbers
     * to each internal node.  This is a post-order traversal, meaning the
     * root is given the largest number.
     *
     * @param node
     * @param nextNr
     * @return
     */
    protected int numberInternalNodes(Node node, int nextNr) {
        if (node.isLeaf())
            return nextNr;

        for (Node child : node.getChildren())
            nextNr = numberInternalNodes(child, nextNr);

        node.setNr(nextNr);
        node.setID(String.valueOf(nextNr));

        return nextNr+1;
    }

    /**
     * Return string representation of quasi-species tree.  We use reflection
     * here to determine whether this is being called as part of writing
     * the state file.
     *
     * @return Quasi-species tree string in Newick format.
     */
    @Override
    public String toString() {

        // Behaves differently if writing a state file
        StackTraceElement[] ste = Thread.currentThread().getStackTrace();
        if (ste[2].getMethodName().equals("toXML")) {
            // Use toShortNewick to generate Newick string without taxon labels
            String string = getFlattenedHaploTree().getRoot().toShortNewick(true);

            // Sanitize ampersands if this is destined for a state file.
            return string;//.replaceAll("&", "&amp;");
        } else{
            return getFlattenedHaploTree().getRoot().toSortedNewick(new int[1], true);
        }
    }

    /**
     * Return number of internalnodes and duplicate attachment points
     * for scaling operator to correctly determine the hastings ratio
     *
     * @return number of internal nodes + duplicate attachment times
     */
    @Override
    public int scale(final double scale) {
        int dof = ((QuasiSpeciesNode) root).scale(scale);
        this.countAndSetPossibleStartBranches();
        return getInternalNodeCount() + getTotalAttachmentCounts() - getDirectAncestorNodeCount();
    }

    /////////////////////////////////////////////////
    //           StateNode implementation          //
    /////////////////////////////////////////////////
    /**
     * Store method for storing state of the tree/nodes before the new proposal
     *
     */
    @Override
    protected void store() {
        int iRoot = root.getNr();
        storedRoot = m_storedNodes[iRoot];

        storeNodes(0, iRoot);

        storedRoot.setHeight(m_nodes[iRoot].getHeight());
        if (this.getLeafNodeCount()>1){
            storedRoot.setParent(null);

            if (root.getLeft()!=null)
                storedRoot.setLeft(m_storedNodes[root.getLeft().getNr()]);
            else
                storedRoot.setLeft(null);
            if (root.getRight()!=null)
                storedRoot.setRight(m_storedNodes[root.getRight().getNr()]);
            else
                storedRoot.setRight(null);
        }
        QuasiSpeciesNode qsStoredRoot = (QuasiSpeciesNode)storedRoot;
        qsStoredRoot.setHaploAboveName(((QuasiSpeciesNode)m_nodes[iRoot]).getHaploAboveName());
        qsStoredRoot.setContinuingHaploName(((QuasiSpeciesNode)m_nodes[iRoot]).getContinuingHaploName());

        if (m_nodes[iRoot].isLeaf()) {
            System.arraycopy(((QuasiSpeciesNode) m_nodes[iRoot]).getAttachmentTimesList(),0,
                    qsStoredRoot.getAttachmentTimesList(),0,
                    ((QuasiSpeciesNode) m_nodes[iRoot]).getAttachmentTimesList().length);
            System.arraycopy(((QuasiSpeciesNode) m_nodes[iRoot]).getTipTimesList(),0,
                    qsStoredRoot.getTipTimesList(),0,
                    ((QuasiSpeciesNode) m_nodes[iRoot]).getTipTimesList().length);
            System.arraycopy(((QuasiSpeciesNode) m_nodes[iRoot]).getTipTimesCountList(),0,
                    qsStoredRoot.getTipTimesCountList(),0,
                    ((QuasiSpeciesNode) m_nodes[iRoot]).getTipTimesCountList().length);
            qsStoredRoot.setParentHaplo(((QuasiSpeciesNode) m_nodes[iRoot]).getParentHaplo());
        }
        //else
        //    qsStoredRoot.setStartBranchCounts(((QuasiSpeciesNode)m_nodes[iRoot]).getStartBranchCounts());

        storeNodes(iRoot+1, nodeCount);
    }

    /**
     * helper to store *
     */
    private void storeNodes(int iStart, int iEnd) {
        for (int i = iStart; i < iEnd; i++) {
            QuasiSpeciesNode sink = (QuasiSpeciesNode)m_storedNodes[i];
            QuasiSpeciesNode src = (QuasiSpeciesNode)m_nodes[i];

            sink.setHeight(src.getHeight());
            sink.setParent(m_storedNodes[src.getParent().getNr()]);

            sink.setHaploAboveName(src.getHaploAboveName());
            sink.setContinuingHaploName(src.getContinuingHaploName());

            if (src.isLeaf()){
                System.arraycopy(src.getAttachmentTimesList(),0,sink.getAttachmentTimesList(),0,src.getAttachmentTimesList().length);
                System.arraycopy(src.getTipTimesList(),0,sink.getTipTimesList(),0,src.getTipTimesList().length);
                System.arraycopy(src.getTipTimesCountList(),0,sink.getTipTimesCountList(),0,src.getTipTimesCountList().length);
                sink.setParentHaplo(src.getParentHaplo());
            } //else
              //  sink.setStartBranchCounts(src.getStartBranchCounts());

            if (src.getLeft()!=null) {
                sink.setLeft(m_storedNodes[src.getLeft().getNr()]);
                if (src.getRight()!=null)
                    sink.setRight(m_storedNodes[src.getRight().getNr()]);
                else
                    sink.setRight(null);
            }
        }
    }


    @Override
    public void restore() {

        // necessary for sampled ancestor trees
        nodeCount = m_storedNodes.length;

        final Node[] tmp = m_storedNodes;
        m_storedNodes = m_nodes;
        m_nodes = tmp;
        root = m_nodes[storedRoot.getNr()];

        // necessary for sampled ancestor trees,
        // we have the nodes, no need for expensive recursion
        leafNodeCount = 0;
        for( Node n : m_nodes ) {
            leafNodeCount += n.isLeaf() ? 1 : 0;
        }

        //leafNodeCount = root.getLeafNodeCount();

        hasStartedEditing = false;

        for( Node n : m_nodes ) {
            n.makeDirty(Tree.IS_CLEAN);
        }

        postCache = null;
    }


    /////////////////////////////////////////////////
    // Methods implementing the Loggable interface //
    /////////////////////////////////////////////////
    /**
     * Method for logging the header of the tree
     *
     */
    @Override
    public void init(PrintStream printStream){

        printStream.println("#NEXUS\n");
        printStream.println("Begin taxa;");
        printStream.println("\tDimensions ntax="+getLeafNodeCount()+";");
        printStream.println("\t\tTaxlabels");
        for (int i = 0; i<getLeafNodeCount(); i++)
            printStream.println("\t\t\t"+getNodesAsArray()[i].getID());
        printStream.println("\t\t\t;");
        printStream.println("End;");

        printStream.println("Begin trees;");
        printStream.println("\tTranslate");
        for (int i = 0; i<getLeafNodeCount(); i++) {
            printStream.print("\t\t\t"+(getNodesAsArray()[i].getNr()+1)
                    +" "+getNodesAsArray()[i].getID());
            if (i<getLeafNodeCount()-1)
                printStream.print(",");
            printStream.print("\n");
        }
        printStream.print("\t\t\t;");
    }

    /**
     * Method for logging the actual tree to state file
     *
     */
    @Override
    public void log(int i, PrintStream printStream) {
        printStream.print("tree STATE_"+i+" = ");
        printStream.print(toString());
        printStream.print(";");
    }

    /**
     * Method for finishing logging the tree
     *
     */
    @Override
    public void close(PrintStream printStream) {
        printStream.println("End;");
    }


    /////////////////////////////////////////////////
    // Serialization and deserialization for state //
    /////////////////////////////////////////////////

    /**
     * reconstruct tree from XML fragment in the form of a DOM node *
     *
     * @param node
     */
    @Override
    public void fromXML(org.w3c.dom.Node node) {
        try {
            String sNewick = node.getTextContent().replace("{[","{").replace("]}","}");

            TreeParser parser = new TreeParser();
            parser.initByName(
                    "IsLabelledNewick", false,
                    "offset", 0,
                    "adjustTipHeights", false,
                    "newick", sNewick);
            //parser.m_nThreshold.setValue(1e-10, parser);
            //parser.m_nOffset.setValue(0, parser);

            initFromFlatHaploTree(parser);

            initArrays();
        } catch (Exception ex) {
            Logger.getLogger(QuasiSpeciesTree.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    
    Node [] externalNodeArray = null;
    // fast alternative for getExternalNodes()
    // assumes all leave nodes are numbered 0,...,leafnodecount-1
	public Node [] getExternalNodesArray() {
		if (externalNodeArray == null) {
			externalNodeArray = new Node[getLeafNodeCount()];
		}
		System.arraycopy(m_nodes, 0, externalNodeArray, 0, externalNodeArray.length);
        return externalNodeArray;
    }

}