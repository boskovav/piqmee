package quasispeciestree.tree;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.distance.Distance;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import quasispeciestree.distance.DifferenceCount;

import java.io.PrintStream;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

/**
 * @author Veronika Boskova created on 29/04/2015
 */

public class QuasiSpeciesTree extends Tree {

// TODO the haplotype count can be input as Integer counts only -- implement uncertainty? A:100-120 reads, B: 3200-3700reads, etc
    public Input<TraitSet> haplotypeCountsInput =
            new Input<TraitSet>("haplotypeCounts","Count of sequences for each haplotype (including the one representative of each haplotype in the tree input)",
            Input.Validate.REQUIRED);
    public Input<RealParameter> originInput = new Input<>(
            "origin", "The time from origin to last sample (must be larger than tree height)",
            Input.Validate.REQUIRED);


    private Map<String,Integer> haplotypeCounts;
    private String qsLabel;


    public QuasiSpeciesTree() { }

    protected RealParameter origin;


    public QuasiSpeciesTree(Node rootNode) {

        if (!(rootNode instanceof QuasiSpeciesNode))
            throw new IllegalArgumentException("Attempted to instantiate "
                    + "quasi-species tree with regular root node.");

        setRoot(rootNode);
        initArrays();
    }

    // init and validate from scratch in order to implement the quasi-species node -- holding the haplotype starting above
    public void initAndValidate() {

        if (m_initial.get() != null && !(this instanceof StateNodeInitialiser)) {

            if (!(m_initial.get() instanceof QuasiSpeciesTree)) {
                throw new IllegalArgumentException("Attempted to initialise "
                        + "quasi-species tree with regular tree object.");
            }

            QuasiSpeciesTree other = (QuasiSpeciesTree) m_initial.get();
            root = other.root.copy();
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

        processTraits(m_traitList.get());

        // Ensure tree is compatible with traits.
        if (hasDateTrait())
            adjustTreeNodeHeights(root);


        origin = originInput.get();

        haplotypeCounts = new HashMap<>();

        if (haplotypeCountsInput.get()!=null){
            TraitSet haplotypeCountsTrait= haplotypeCountsInput.get();
            qsLabel = haplotypeCountsTrait.getTraitName();

            setHaploCounts(haplotypeCountsTrait);

//            initAttachmentTimes();
//
//            fillParentHaplo();
////
//            startBranchCounts = countPossibleStartBranches();
        }

    }

    /*
    //
    //
    //          OWN FUNCTIONS
    //
    //
    */

    /**
     * Function to initiate the array list of attachment times for each haplotype in quasispecies
     */
    private void initAttachmentTimes(){
        ArrayList attachmentTimesListTmp = null;
        // for those nodes where haplotype arises, change the haploAboveNode to haplotype's (corresponding tip node) number
        //  and for nodes below up to the tip set continuingHaploName to the same haplotype's (tip)
        for (Node node : this.getExternalNodes()){
            // check if getNr() always returns the same >>> Node number is guaranteed not to change during an MCMC run.
            //      written in the Node class)
            // assigns to each unique haplotype an array of size haplotypeCount
            double[] tempqstimes = new double[((QuasiSpeciesTip) node).getHaplotypeCountsFromTips()];
            double[] temptiptimes = ((QuasiSpeciesTip) node).getTipTimesList();
            int[] temptiptimescount = ((QuasiSpeciesTip) node).getTipTimesCountList();
            // start with the star tree, and define start of haplotype at the multi-furcating node time
            int currentPosition = 0;
            // check if the tree has more than one tip
            if (this.getLeafNodeCount()>1) {
                if (haplotypeCounts != null){
                    for (int j = 1; j < temptiptimescount[0]; j++) {
                        tempqstimes[currentPosition] = node.getParent().getHeight()
                                - (j + 1) * ((node.getParent().getHeight() - temptiptimes[0]) / (1 + temptiptimescount[0]));
                        currentPosition++;
                    }
                    for (int i = 1; i < temptiptimes.length; i++) {
                        tempqstimes[currentPosition]=((double[])attachmentTimesListTmp.get(node.getNr()))[0];
                        currentPosition++;
                        for (int j = 0; j < temptiptimescount[i]; j++) {
                            tempqstimes[currentPosition] = ((double[])attachmentTimesListTmp.get(node.getNr()))[i-1]
                                    - (j + 1) * ((((double[])attachmentTimesListTmp.get(node.getNr()))[i-1] - temptiptimes[i]) / (1 + temptiptimescount[i]));
                            currentPosition++;
                        }
                    }
                }
                else {
                    for (int j = 0; j < temptiptimescount[0]; j++) {
                        tempqstimes[currentPosition] = node.getParent().getHeight()
                                - (j + 1) * ((node.getParent().getHeight() - temptiptimes[0]) / (1 + temptiptimescount[0]));
                        currentPosition++;
                    }
                    for (int i = 1; i < temptiptimes.length; i++) {
                        for (int j = 0; j < temptiptimescount[i]; j++) {
                            tempqstimes[currentPosition] = node.getParent().getHeight()
                                    - (j + 1) * ((node.getParent().getHeight() - temptiptimes[i]) / (1 + temptiptimescount[i]));
                            currentPosition++;
                        }
                    }
                }
            }
            // or there is just one sequence possibly sampled through time
            else {
                if (haplotypeCounts != null){
                    for (int j = 1; j < temptiptimescount[0]; j++) {
                        tempqstimes[currentPosition] = origin.getValue()
                                - (j + 1) * ((origin.getValue() - temptiptimes[0]) / (1 + temptiptimescount[0]));
                        currentPosition++;
                    }
                    for (int i = 1; i < temptiptimes.length; i++) {
                        tempqstimes[currentPosition]=((double[])attachmentTimesListTmp.get(node.getNr()))[0];
                        currentPosition++;
                        for (int j = 0; j < temptiptimescount[i]; j++) {
                            tempqstimes[currentPosition] = ((double[])attachmentTimesListTmp.get(node.getNr()))[i-1]
                                    - (j + 1) * ((((double[])attachmentTimesListTmp.get(node.getNr()))[i-1] - temptiptimes[i]) / (1 + temptiptimescount[i]));
                            currentPosition++;
                        }
                    }
                }
                else {

                    for (int j = 0; j < temptiptimescount[0]; j++) {
                        tempqstimes[currentPosition] = origin.getValue()
                                - (j + 1) * ((origin.getValue() - temptiptimes[0]) / (1 + temptiptimescount[0]));
                        currentPosition++;
                    }
                    for (int i = 1; i < temptiptimes.length; i++) {
                        for (int j = 0; j < temptiptimescount[i]; j++) {
                            tempqstimes[currentPosition] = node.getParent().getHeight()
                                    - (j + 1) * ((node.getParent().getHeight() - temptiptimes[i]) / (1 + temptiptimescount[i]));
                            currentPosition++;
                        }
                    }
                }
            }
            // Assign the first entry to be the same as the first split
            // this one is marking the start of the haplotype
             if (tempqstimes.length > 1){
                 Arrays.sort(tempqstimes);
                 tempqstimes[0] = tempqstimes[tempqstimes.length-1];
             }
            else
                tempqstimes[0]=node.getHeight();

            ((QuasiSpeciesNode) node).setHaploAboveName(node.getNr());
            ((QuasiSpeciesNode) node).setContinuingHaploName(node.getNr());
            // attachment time list defined as "node" height, going from present (0) to past (positive height)
            ((QuasiSpeciesTip) node).setAttachmentTimesList(tempqstimes);
        }
    }

//    private void setTipTimesLists(QuasiSpeciesTip node, double[] uniquetemptiptimes, int[] tiptimescounts){
//        node.setTipTimesList(uniquetemptiptimes);
//        node.setTipTimesCountList(tiptimescounts);
//    }

//    private void setTipTimesLists(int position, double[] uniquetemptiptimes, int[] tiptimescounts){
//        ((QuasiSpeciesTip) this.getNode(position)).setTipTimesList(uniquetemptiptimes);
//        ((QuasiSpeciesTip) this.getNode(position)).setTipTimesCountList(tiptimescounts);
//    }

    public int getHaplotypeCounts(QuasiSpeciesNode node){
        return haplotypeCounts.get(node.getNr());
    }


    /**
     * Gets the total count of attachment points for all haplotypes
     *
     */
    public int getTotalAttachmentCounts(){
        int totalCount=0;
        for (Node node : this.getExternalNodes()){

            //TODO is this working properly.. just transformed the haplotypeCounts to the trait set from the hashmap
            totalCount += haplotypeCounts.get(node.getID()) - 1;
        }
        return totalCount;
    }


    /**
     * Sets the parent haplo name to -1 for all tips
     *
     */
    public void clearParentHaploNames() {
        for (Node node : this.getExternalNodes()){
            ((QuasiSpeciesTip) node).setParentHaplo(-1);
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
            if (((QuasiSpeciesTip) node).getAttachmentTimesList().length + 1 != (int) haploCounts.getValue(node.getID()))
                throw new RuntimeException("QuasiSpeciesTree class: The new to be set haploCount is not " +
                                           "equal the the number of attachmetn times. Why?");
            haplotypeCounts.put(node.getID(), (int) haploCounts.getValue(node.getID()));
        }
    }

    /**
     * Sets the number of haplotype counts for each haplo
     *
     * @param haploCounts new set of haplotype counts
     * @param tree a tree whose nodes are the be used for assignment
     */
    private void setHaploCounts(TraitSet haploCounts, Tree tree){
        for (Node node : tree.getExternalNodes()) {
            if (((QuasiSpeciesTip) node).getAttachmentTimesList().length + 1 != (int) haploCounts.getValue(node.getID()))
                throw new RuntimeException("QuasiSpeciesTree class: The new to be set haploCount is not " +
                        "equal the the number of attachmetn times. Why?");
            haplotypeCounts.put(node.getID(), (int) haploCounts.getValue(node.getID()));
        }
    }

    /**
     * Method to count the possible number of start branches for each internal node
     * at one pre-order tree pass (prevent energy & resources waste)
     *
     * @return Array holding for each internal node the number
     *         of possible start branches, held in array at position determined by node.getNr()-nTips
     */
    public int[] countPossibleStartBranches(){
        // get the counts of possible number of start branches for each haplotype (i.e. internal node attachment branches)
        int nTips = this.getLeafNodeCount();
        // array to store the count of possible start branches for each internal node
        int[] startBranchCountArray = new int[this.getInternalNodeCount()];
        for (int i = 0; i<this.getInternalNodeCount(); i++){
            // once determined the parent haplotype find out from how many branches
            // of the parent haplotype can it actually branch off
            int haplo = ((QuasiSpeciesNode) this.getNode(nTips + i)).getContinuingHaploName();
            if ( haplo != -1 ){
                QuasiSpeciesTip haploTip = (QuasiSpeciesTip) this.getNode(haplo);
                startBranchCountArray[i] = haploTip.countPossibleAttachmentBranches(0, this.getNode(nTips + i).getHeight());
            } else {
                startBranchCountArray[i] = 1;
            }
        }
        return startBranchCountArray;
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
    private void fillParentHaplo(){
        // set all the parent haplotype to -1
        clearParentHaploNames();

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
            QuasiSpeciesNode thisNode = (QuasiSpeciesTip) this.getNode(newHaploType);
            // set the haplotype from which the next haplotype arises as a parent
            ((QuasiSpeciesTip) thisNode).setParentHaplo(currentHaploType);
            // set the continuingHaploName at corresponding internal nodes leading from nextNode to the respective child\
            nextNode.setContinuingHaploName(newHaploType);
            while ( nextNode != thisNode ) {
                thisNode = (QuasiSpeciesNode) thisNode.getParent();
                thisNode.setContinuingHaploName(newHaploType);
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


    /*
    //
    //
    //          FUNCTIONS COPIED AND ADAPTED FROM MULTI-TYPE TREE CLASS
    //
    //
    */

    @Override
    protected final void initArrays() {
        // initialise tree-as-array representation + its stored variant
        m_nodes = new QuasiSpeciesTip[nodeCount];
        listNodes((QuasiSpeciesNode)root, (QuasiSpeciesNode[])m_nodes);
        m_storedNodes = new QuasiSpeciesTip[nodeCount];
        Node copy = root.copy();
        listNodes((QuasiSpeciesNode)copy, (QuasiSpeciesNode[])m_storedNodes);



        // todo DOES IS MATTER HERE THAT WE HAVE ARRAY OF QSnodes and not QStips?????
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
            if (node.getRight()!=null)
                listNodes((QuasiSpeciesNode) node.getRight(), nodes);
        }


        // todo DOES IS MATTER HERE THAT WE HAVE ARRAY OF QSnodes and not QStips?????
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
        tree.root = root.copy();
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

        QuasiSpeciesTip[] qsNodes = new QuasiSpeciesTip[qsTree.getNodeCount()];
        for (int i=0; i<qsTree.getNodeCount(); i++)
            qsNodes[i] = new QuasiSpeciesTip();

        ID = qsTree.ID;
        root = qsNodes[qsTree.root.getNr()];
        root.assignFrom(qsNodes, qsTree.root);
        root.setParent(null);

        nodeCount = qsTree.nodeCount;
        internalNodeCount = qsTree.internalNodeCount;
        leafNodeCount = qsTree.leafNodeCount;
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
        qsRoot.haploAboveName = ((QuasiSpeciesNode)(otherNodes[iRoot])).haploAboveName;
        qsRoot.continuingHaploName = ((QuasiSpeciesNode)(otherNodes[iRoot])).continuingHaploName;

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
        assignFromFragileHelper(iRoot + 1, nodeCount, otherNodes);

        //TODO Why is there twice assignFromFragile helper? ANd why is the second loop from iRoot+1 only?
    }

    /**
     * helper to assignFromFragile *
     */
    private void assignFromFragileHelper(int iStart, int iEnd, Node[] otherNodes) {
        for (int i = iStart; i < iEnd; i++) {
            QuasiSpeciesNode sink = (QuasiSpeciesNode)m_nodes[i];
            QuasiSpeciesNode src = (QuasiSpeciesNode)otherNodes[i];
            sink.setHeight(src.getHeight());
            sink.setParent(m_nodes[src.getParent().getNr()]);

            sink.haploAboveName = src.haploAboveName;
            sink.continuingHaploName = src.continuingHaploName;
            if (src.getAllChildNodes().size() == 0){
                ((QuasiSpeciesTip) sink).setAttachmentTimesList(((QuasiSpeciesTip) src).getAttachmentTimesList());
                ((QuasiSpeciesTip) sink).setTipTimesList(((QuasiSpeciesTip) src).getTipTimesList());
                ((QuasiSpeciesTip) sink).setTipTimesCountList(((QuasiSpeciesTip) src).getTipTimesCountList());
                ((QuasiSpeciesTip) sink).setParentHaplo(((QuasiSpeciesTip) src).getParentHaplo());
            }

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
            QuasiSpeciesTip qsNode = (QuasiSpeciesTip) node;
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
                regularTree.m_nodes[nextNodeNr]=nextNode;
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
                nextTip = new QuasiSpeciesTip();
                nextTip.setNr(nextNodeNr);
                nextTip.setID(String.valueOf(nextNodeNr));
                regularTree.m_nodes[nextNodeNr] = nextTip;
                nextNodeNr += 1;

                // Connect to parent:
                nextNode.addChild(nextTip);

                // Ensure height is set:
                nextTip.setHeight(temptiptimes[currentTipTimePosition]);
                currentTipArrayPosition++;

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
    private Tree getFlattenedHaploTree() {
        // Create new tree to modify.  Note that copy() doesn't
        // initialise the node array lists, so initArrays() must
        // be called manually.
        QuasiSpeciesTree flatTree = copy();
        flatTree.initArrays();

        for (Node node : getExternalNodes()) {
            QuasiSpeciesTip qsNode = (QuasiSpeciesTip) node;
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
    private void initFromFlatHaploTree(Tree flatHaploTree) {

        // Create quasispecies tree
        QuasiSpeciesNode[] quasiSpeciesNodes = new QuasiSpeciesNode[flatHaploTree.getNodeCount()];

        attachmentTimesList = new ArrayList<>(flatHaploTree.getExternalNodes().size());
        tipTimesList = new ArrayList<>(flatHaploTree.getExternalNodes().size());
        tipTimesCountList = new ArrayList<>(flatHaploTree.getExternalNodes().size());

        for (int i = 0; i < flatHaploTree.getExternalNodes().size(); i++){
            attachmentTimesList.add(i, null);
            tipTimesList.add(i, null);
            tipTimesCountList.add(i, null);
        }

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
                            attachmentTimesList.set(thisQuasiSpeciesNode.getNr(), tempqstimes);
                        } else if (typeObject1 instanceof double[]) {
                            attachmentTimesList.add(thisQuasiSpeciesNode.getNr(), (double[]) typeObject1);
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
                            tipTimesList.set(thisQuasiSpeciesNode.getNr(), temptiptimes);

                        } else if (typeObject2 instanceof double[]) {
                            tipTimesList.set(thisQuasiSpeciesNode.getNr(), (double[]) typeObject2);
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
                            tipTimesCountList.set(thisQuasiSpeciesNode.getNr(), temptiptimescount);

                        } else if (typeObject3 instanceof int[]) {
                            tipTimesCountList.set(thisQuasiSpeciesNode.getNr(), (int[]) typeObject3);
                        } else if (typeObject3 instanceof Double[]) {
                            int[] temptiptimescount = new int[((Double[]) typeObject3).length];
                            int count3 = 0;
                            for (Double nextEntry : (Double[]) typeObject3) {
                                int tipTimeCount = nextEntry.intValue();
                                temptiptimescount[count3] = tipTimeCount;
                                count3++;
                            }
                            tipTimesCountList.set(thisQuasiSpeciesNode.getNr(), temptiptimescount);

                        } else if (typeObject3 instanceof double[]) {
                            tipTimesCountList.set(thisQuasiSpeciesNode.getNr(), (int[]) typeObject3);
                        } else
                            throw new IllegalArgumentException("Unrecognised type metadata.");

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
        numberInternalNodes(newRoot, newRoot.getAllLeafNodes().size());

        // Assign tree topology:
        assignFromWithoutID(new QuasiSpeciesTree(newRoot));
        initArrays();

        //attachmentTimesList and tipTimesList created on a go...
        // do not use initAttachmentTimes();
        for (int i = 0; i < attachmentTimesList.size(); i++){
            double[] attachList = attachmentTimesList.get(i);
            double[] tipTimes = tipTimesList.get(i);
            int[] tipTimesCount = tipTimesCountList.get(i);
            Arrays.sort(attachList);
            // reverse the array to start with the largest value
            int totalLength = attachList.length;
            for (int j = 0; j < totalLength/2; j++){
                double tmp = attachList[j];
                attachList[j] = attachList[totalLength-1-j];
                attachList[totalLength-1-j] = tmp;
//                tmp = tipTimes[j];
//                tipTimes[j] = tipTimes[totalLength-1-j];
//                tipTimes[totalLength-1-j] = tmp;
            }
            //Manually sort tipTimes, since we in the same way need to sort the tip counts
//            Arrays.sort(tipTimes);
            double tmp;
            int tmpint;
            for (int j = 1; j < tipTimes.length; j++){
                int k = j-1;
                while (k>=0 && tipTimes[k]>tipTimes[j]){
                    tmp = tipTimes[j];
                    tipTimes[j] = tipTimes[k];
                    tipTimes[k] = tmp;
                    tmpint = tipTimesCount[j];
                    tipTimesCount[j] = tipTimesCount[k];
                    tipTimesCount[k] = tmpint;
                    k--;
                }
            }
        }
        storedAttachmentTimesList = new ArrayList<>(this.getLeafNodeCount());
        storedTipTimesList = new ArrayList<>(this.getLeafNodeCount());
        storedTipTimesCountList = new ArrayList<>(this.getLeafNodeCount());
        for (int i = 0; i<this.getLeafNodeCount(); i++){
            storedAttachmentTimesList.add(i,attachmentTimesList.get(i).clone());
            storedTipTimesList.add(i,tipTimesList.get(i).clone());
            storedTipTimesCountList.add(i,tipTimesCountList.get(i).clone());
        }

        //traverse a tree and assign nodes above and continuing haplo annotations
        for (Node tipNode : this.getExternalNodes()){
            int haplo = tipNode.getNr();
            QuasiSpeciesNode nodeBelowHaplo = (QuasiSpeciesNode) tipNode;
            double maxHaploTime = this.getAttachmentTimesList(haplo)[0];
            while (nodeBelowHaplo.getParent()!=null && nodeBelowHaplo.getParent().getHeight() < maxHaploTime){
                nodeBelowHaplo.setContinuingHaploName(haplo);
                nodeBelowHaplo = (QuasiSpeciesNode) nodeBelowHaplo.getParent();
            }
            nodeBelowHaplo.setHaploAboveName(haplo);
            nodeBelowHaplo.setContinuingHaploName(haplo);
        }

        fillParentHaplo();

        System.arraycopy(parentHaplo,0,storedParentHaplo,0,parentHaplo.length);

        startBranchCounts = countPossibleStartBranches();
        storedStartBranchCounts = new int[startBranchCounts.length];
        System.arraycopy(startBranchCounts,0,storedStartBranchCounts,0,startBranchCounts.length);

        for (Node node : this.getExternalNodes()){
            setHaploCounts(node.getNr(), attachmentTimesList.get(node.getNr()).length-1);
        }
    }

    /**
     * Initialise tree topology from Tree object with duplicate counts (but no attachment times) in trait set array
     *
     * @param uniqueHaploTree
     */
    public void initFromUniqueHaploTree(Tree uniqueHaploTree, Alignment data){

        // In unique haplo tree, there can still be duplicate sequences, if found at different points in time
        // Get the distances for the sequences:
        double[][] distanceMatrix = getSequenceDistances(data, uniqueHaploTree);

        // Build new quasi-species tree:
        ArrayList haplotypesSeen = new ArrayList<>();
        List<QuasiSpeciesNode> qsTips = new ArrayList<>();
        List<QuasiSpeciesNode> qsInternalNodes = new ArrayList<>();
        //initialize attachmentTimesList and tipTimesList
        attachmentTimesList = new ArrayList<>();
        tipTimesList = new ArrayList<>();
        tipTimesCountList = new ArrayList<>();
        setHaploCounts(haplotypeCountsInput.get(),uniqueHaploTree);
        ArrayList result = processNextNodeOfFullNewickTree(uniqueHaploTree.getRoot(),qsTips,qsInternalNodes,
                distanceMatrix,haplotypesSeen,0);

        // renumber tips to match the number of tips in the qsTree (so far matching fullTree node numbers)
        // need to match the tip times and attach time and haplo count lists!! -- this should not affect the order
        //  in these arrays, as the values were put in in the order of tips put in qsTips array
        for (int i = 0; i<qsTips.size(); i++){
            qsTips.get(i).setNr(i);
        }
        QuasiSpeciesNode newRoot = (QuasiSpeciesNode) result.get(0);
        // Number internal nodes:
        numberInternalNodes(newRoot, newRoot.getAllLeafNodes().size());

        // Assign tree topology:
        assignFromWithoutID(new QuasiSpeciesTree(newRoot));
        initArrays();

        // Make sure to properly assign the attachmentTimes
        // Set start of haplotype times as default to belong to the leaf node
        // treeNode.setHaploname(treeNode.getID());
        // done in initAndValidate!!!
        initAttachmentTimes();

        // tipTimesList created on a go...
        for (int i = 0; i < attachmentTimesList.size(); i++){
            double[] attachList = attachmentTimesList.get(i);
            double[] tipTimes = tipTimesList.get(i);
            int[] tipTimesCount = tipTimesCountList.get(i);
            Arrays.sort(attachList);
            // copy the largest bifurcation time, to indicate the haplo start time
            // done in initAttachmentTimes();
//            attachList[0]=attachList[attachList.length-1];
            // reverse the rest of the array to start with the largest value
            int totalLength = attachList.length;
            for (int j = 0; j < totalLength/2; j++){
                double tmp = attachList[j];
                attachList[j] = attachList[totalLength-1-j];
                attachList[totalLength-1-j] = tmp;
            }
            //Manually sort tipTimes, since we in the same way need to sort the tip counts
            double tmp;
            int tmpint;
            for (int j = 1; j < tipTimes.length; j++){
                int k = j-1;
                while (k>=0 && tipTimes[k]>tipTimes[j]){
                    tmp = tipTimes[j];
                    tipTimes[j] = tipTimes[k];
                    tipTimes[k] = tmp;
                    tmpint = tipTimesCount[j];
                    tipTimesCount[j] = tipTimesCount[k];
                    tipTimesCount[k] = tmpint;
                    k--;
                }
            }
        }

        storedAttachmentTimesList = new ArrayList<>(this.getLeafNodeCount());
        storedTipTimesList = new ArrayList<>(this.getLeafNodeCount());
        storedTipTimesCountList = new ArrayList<>(this.getLeafNodeCount());
        for (int i = 0; i<this.getLeafNodeCount(); i++){
            storedAttachmentTimesList.add(i,attachmentTimesList.get(i).clone());
            storedTipTimesList.add(i,tipTimesList.get(i).clone());
            storedTipTimesCountList.add(i,tipTimesCountList.get(i).clone());
        }

        //traverse a tree and assign nodes above and continuing haplo annotations
        for (Node tipNode : this.getExternalNodes()){
            int haplo = tipNode.getNr();
            QuasiSpeciesNode nodeBelowHaplo = (QuasiSpeciesNode) tipNode;
            double maxHaploTime = this.getAttachmentTimesList(haplo)[0];
            while (nodeBelowHaplo.getParent()!=null && nodeBelowHaplo.getParent().getHeight() < maxHaploTime){
                nodeBelowHaplo.setContinuingHaploName(haplo);
                nodeBelowHaplo = (QuasiSpeciesNode) nodeBelowHaplo.getParent();
            }
            nodeBelowHaplo.setHaploAboveName(haplo);
            nodeBelowHaplo.setContinuingHaploName(haplo);
        }

        fillParentHaplo();

        System.arraycopy(parentHaplo,0,storedParentHaplo,0,parentHaplo.length);

        startBranchCounts = countPossibleStartBranches();
        storedStartBranchCounts = new int[startBranchCounts.length];
        System.arraycopy(startBranchCounts,0,storedStartBranchCounts,0,startBranchCounts.length);

        for (Node node : this.getExternalNodes()){
            setHaploCounts(node.getNr(), attachmentTimesList.get(node.getNr()).length-1);
        }

    }

    /**
     * Initialise tree topology from full Tree object - have to check for duplicates in post order (children first) traversal
     *                                                  and collapse duplicates, throw error if tree does not fullfill
     *                                                  recursively monophyletic constraint of our QS model
     *
     * @param fullTree
     */
    public void initFromFullTree(Tree fullTree, Alignment data){

        // Get the distances for the sequences:
        double[][] distanceMatrix = getSequenceDistances(data, fullTree);

        // Build new quasi-species tree:
        ArrayList haplotypesSeen = new ArrayList<>();
        List<QuasiSpeciesNode> qsTips = new ArrayList<>();
        List<QuasiSpeciesNode> qsInternalNodes = new ArrayList<>();
        //initialize attachmentTimesList and tipTimesList
        attachmentTimesList = new ArrayList<>();
        tipTimesList = new ArrayList<>();
        ArrayList result = processNextNodeOfFullNewickTree(fullTree.getRoot(),qsTips,qsInternalNodes,
                                                            distanceMatrix,haplotypesSeen,0);

        // post-process tipTimesList
        processTipTimes();

        // renumber tips to match the number of tips in the qsTree (so far matching fullTree node numbers)
        // need to match the tip times and attach time and haplo count lists!! -- this should not affect the order
        //  in these arrays, as the values were put in in the order of tips put in qsTips array
        for (int i = 0; i<qsTips.size(); i++){
            qsTips.get(i).setNr(i);
        }
        QuasiSpeciesNode newRoot = (QuasiSpeciesNode) result.get(0);
        // Number internal nodes:
        numberInternalNodes(newRoot, newRoot.getAllLeafNodes().size());

        // Assign tree topology:
        assignFromWithoutID(new QuasiSpeciesTree(newRoot));
        initArrays();

        //attachmentTimesList and tipTimesList created on a go...
        // do not use initAttachmentTimes();
        for (int i = 0; i < attachmentTimesList.size(); i++){
            double[] attachList = attachmentTimesList.get(i);
            double[] tipTimes = tipTimesList.get(i);
            int[] tipTimesCount = tipTimesCountList.get(i);
            Arrays.sort(attachList);
            // copy the largest bifurcation time, to indicate the haplo start time
            attachList[0]=attachList[attachList.length-1];
            // reverse the rest of the array to start with the largest value
            int totalLength = attachList.length;
            for (int j = 0; j < totalLength/2; j++){
                double tmp = attachList[j+1];
                attachList[j+1] = attachList[totalLength-1-j];
                attachList[totalLength-1-j] = tmp;
            }
            //Manually sort tipTimes, since we in the same way need to sort the tip counts
            double tmp;
            int tmpint;
            for (int j = 1; j < tipTimes.length; j++){
                int k = j-1;
                while (k>=0 && tipTimes[k]>tipTimes[j]){
                    tmp = tipTimes[j];
                    tipTimes[j] = tipTimes[k];
                    tipTimes[k] = tmp;
                    tmpint = tipTimesCount[j];
                    tipTimesCount[j] = tipTimesCount[k];
                    tipTimesCount[k] = tmpint;
                    k--;
                }
            }
        }
        storedAttachmentTimesList = new ArrayList<>(this.getLeafNodeCount());
        storedTipTimesList = new ArrayList<>(this.getLeafNodeCount());
        storedTipTimesCountList = new ArrayList<>(this.getLeafNodeCount());
        for (int i = 0; i<this.getLeafNodeCount(); i++){
            storedAttachmentTimesList.add(i,attachmentTimesList.get(i).clone());
            storedTipTimesList.add(i,tipTimesList.get(i).clone());
            storedTipTimesCountList.add(i,tipTimesCountList.get(i).clone());
        }

        //traverse a tree and assign nodes above and continuing haplo annotations
        for (Node tipNode : this.getExternalNodes()){
            int haplo = tipNode.getNr();
            QuasiSpeciesNode nodeBelowHaplo = (QuasiSpeciesNode) tipNode;
            double maxHaploTime = this.getAttachmentTimesList(haplo)[0];
            while (nodeBelowHaplo.getParent()!=null && nodeBelowHaplo.getParent().getHeight() < maxHaploTime){
                nodeBelowHaplo.setContinuingHaploName(haplo);
                nodeBelowHaplo = (QuasiSpeciesNode) nodeBelowHaplo.getParent();
            }
            nodeBelowHaplo.setHaploAboveName(haplo);
            nodeBelowHaplo.setContinuingHaploName(haplo);
        }

        fillParentHaplo();

        System.arraycopy(parentHaplo,0,storedParentHaplo,0,parentHaplo.length);

        startBranchCounts = countPossibleStartBranches();
        storedStartBranchCounts = new int[startBranchCounts.length];
        System.arraycopy(startBranchCounts,0,storedStartBranchCounts,0,startBranchCounts.length);

        for (Node node : this.getExternalNodes()){
            setHaploCounts(node.getNr(), attachmentTimesList.get(node.getNr()).length-1);
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
     * @param nextTipNumber
     * @return
     */
    private ArrayList processNextNodeOfFullNewickTree(
            Node node, List<QuasiSpeciesNode> qsTips, List<QuasiSpeciesNode> qsInternalNodes,
            double[][] distanceMatrix, ArrayList haplotypesSeen, int nextTipNumber){

        QuasiSpeciesNode returnNode = null;
        ArrayList haplotypesAtThisNode = new ArrayList();
        // fakeHaplo the haplotype number (in the qsTree); a node in the regular that is a sequence duplicate of a previously
        //      seen haplotype is called "fakeHaplo" since all the nodes that attach to this brach should attach
        //      to the branch leading to the true tip
        int fakeHaplo = -1;
        // for leaf nodes check if the sequence has been seen at another node already
        // pass on to the parent the info on which haplo is at the tip
        if (node.isLeaf()){
            boolean skip = false;
            // check if the sequence has been seen already
            for (int i = 0; i < haplotypesSeen.size(); i++) {
                if (distanceMatrix[node.getNr()][(int) haplotypesSeen.get(i)] == 0) {
                    QuasiSpeciesNode seenNode = qsTips.get(i);
                    // check if the time of the tip is less than the uniqueHaploTree tip
                    // if not, rewrite the info on the uniqueHaploTree tip
                    if (node.getHeight() < seenNode.getHeight()) {
                        seenNode.setHeight(node.getHeight());
                        seenNode.setID(String.valueOf(node.getID()));
                    }
                    // to robustly always pick the same node as the unique haplo node even after restart from state file
                    else if (!Pattern.compile("^t").matcher(seenNode.getID()).find()
                            && Pattern.compile("^t").matcher(node.getID()).find()){
                        seenNode.setID(String.valueOf(node.getID()));
                    }
                    // if the tips do not start with t - then keep always the tip with smallest number
                    else if (seenNode.getID().compareTo(node.getID())>0){
                        seenNode.setID(String.valueOf(node.getID()));
                    }
                    // since the sequence has been seen already, assign the tip time to array
                    if (haplotypeCounts.size() != 0) {
                        boolean assigned = false;
                        // throw an error if if tip has already been found with same seq and at the same time...
                        // the user wants to use full tree? or did not correctly merge duplicate sequences?
                        for (int j = 0; j < tipTimesList.get(i).length; j++) {
                            // check if the tips with the same sequence that have been seen had also the current node's sampling time
                            if (tipTimesList.get(i)[j] == node.getHeight()) {
                                System.out.println("There are at least two tips with the same sequence and same sampling time." +
                                        "If you want to input tree with all sequences as tips, please use class " +
                                        " quasispeciestree.tree.QuasiSpeciesTreeFromFullNewick;" +
                                        " otherwise please remove duplicates and use haplotypeCounts traitset" +
                                        " to annotate the duplicate counts");
                                System.exit(0);
                            }
                            // if with different time, add to tip times and counts
                            if (!assigned && tipTimesList.get(i)[j] == 0 && tipTimesCountList.get(i)[j] == 0) {
                                tipTimesList.get(i)[j] = node.getHeight();
                                tipTimesCountList.get(i)[j] = haplotypeCounts.get(node.getNr())+1;
                                assigned = true;
                            }
                        }
                    }
                    else {
                        for (int j = 0; j < tipTimesList.get(i).length; j++) {
                            // check if the tips with the same sequence that have been seen had also the current node's sampling time
                            if (tipTimesList.get(i)[j] == 0) {
                                tipTimesList.get(i)[j] = node.getHeight();
                                break;
                            }
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
                for (int i = 0; i<distances.length; i++){
                    if (distances[i]==0)
                        newEntryLength++;
                    else
                        break;
                }
                attachmentTimesList.add(new double[newEntryLength]);
                tipTimesList.add(new double[newEntryLength]);
                tipTimesList.get(nextTipNumber)[0] = node.getHeight();
                if (haplotypeCounts.size() != 0){
                    tipTimesCountList.add(new int[newEntryLength]);
                    tipTimesCountList.get(nextTipNumber)[0] = haplotypeCounts.get(node.getNr())+1;
                }
                nextTipNumber++;
            }
        }
        else {
            ArrayList leftOut = processNextNodeOfFullNewickTree(node.getLeft(),qsTips,qsInternalNodes,
                                                                distanceMatrix,haplotypesSeen,nextTipNumber);
            nextTipNumber = (int)leftOut.get(3);
            ArrayList rightOut = processNextNodeOfFullNewickTree(node.getRight(),qsTips,qsInternalNodes,
                                                                 distanceMatrix,haplotypesSeen,nextTipNumber);
            nextTipNumber = (int)rightOut.get(3);
            QuasiSpeciesNode leftNode = (QuasiSpeciesNode) leftOut.get(0);
            QuasiSpeciesNode rightNode = (QuasiSpeciesNode) rightOut.get(0);
            ArrayList leftHaplo = (ArrayList) leftOut.get(1);
            ArrayList rightHaplo = (ArrayList) rightOut.get(1);
            int leftFakeHaplo = (int) leftOut.get(2);
            int rightFakeHaplo = (int) rightOut.get(2);
            // Case 1:  there is more than one QS identical throw exception = this is not a QS tree
            // Case 2:  both arrays have exactly one same QS this means we found a duplicate and need to record this time
            //            into attachmentTimesList
            // Case 3:  this is a real QS tree node where several haplotypes meet -- check if any of haplo has been seen already
            //            if it has, assign the internal node to the existing branch of that haplo
            ArrayList sameHaploLeftRight = new ArrayList();
            for (int i = 0; i<leftHaplo.size(); i++){
                for (int j = 0; j<rightHaplo.size(); j++){
                    if ((int) leftHaplo.get(i)==(int) rightHaplo.get(j)) {
                        sameHaploLeftRight.add(leftHaplo.get(i));
                    }
                }
            }
            //case 1
            if (sameHaploLeftRight.size()>1){
                System.out.println("The input tree is not recursively monophyletic and therefore cannot be converted to " +
                                    " a quasi-species tree. Try to input a different tree. Alternatively, input sequences only.");
                System.exit(0);
            }
            //case 2
            else if(sameHaploLeftRight.size()==1){
                // since the sequence has been seen already, assign the attachment haplo time to array
                int haplo = (int) sameHaploLeftRight.get(0);
                int tip = -1;
                for (int i = 0; i<qsTips.size(); i++){
                    if (qsTips.get(i).getNr() == haplo){
                        tip = i;
                        break;
                    }
                }
                for (int i = 0; i<attachmentTimesList.get(tip).length; i++) {
                    if (attachmentTimesList.get(tip)[i]==0) {
                        attachmentTimesList.get(tip)[i] = node.getHeight();
                        break;
                    }
                }
                haplotypesAtThisNode=leftHaplo;
                haplotypesAtThisNode.addAll(haplotypesAtThisNode.size(),rightHaplo);
                // return the first index of the recurring element
                int index1 = haplotypesAtThisNode.indexOf(haplo);
                haplotypesAtThisNode.remove(index1);
                //check what is the returnNode : will be null if both nodes are fake
                //  otherwise the higher of the two QS nodes from left and right
                //  the lower node is has already been appended to the branch leading to haplo (see case 3)
                if (leftNode==null && rightNode==null)
                    returnNode = null;
                else if (rightNode ==null)
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
                    System.out.println("There is something seriously wrong with the tree or with this algorithm.");
                    System.exit(0);
                }
                else {
                    System.out.println("What case have I forgotten?");
                    System.exit(0);
                }
            }
            //case 3
            else {
                haplotypesAtThisNode=leftHaplo;
                haplotypesAtThisNode.addAll(haplotypesAtThisNode.size(),rightHaplo);
                // look for fake haplotype (i.e. one of the nodes below this node are not qsTips
                QuasiSpeciesNode nodeToPlace = null;
                // if we have a fake node, connect the left-right haplotypes at a new node corresponding to fake haplo
                if(leftFakeHaplo!=-1 || rightFakeHaplo !=-1){
                    //find the TRUE tip belonging to the found haplo and the node that should be the parent of the new node
                    QuasiSpeciesNode checkThisNode = null;
                    QuasiSpeciesNode checkThisNodeKid = null;
                    if (leftFakeHaplo != -1){
                        checkThisNode = qsTips.get(leftFakeHaplo);
                        nodeToPlace = rightNode;
                        fakeHaplo = leftFakeHaplo;
                    }
                    if (rightFakeHaplo != -1){
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
                        if (checkThisNode.getLeft().getNr() == checkThisNodeKid.getNr())
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
        output.add(nextTipNumber);
        return output;
    }


    /**
     * Helper method used by initFromFullTree/initFromUniqueHaploTree to
     * process the tip times as array to tipTimesList and tipTimesCountList.
     */
    private void processTipTimes() {
        // post-process tipTimesList
        ArrayList tipTimesTemp = tipTimesList;
        tipTimesList = new ArrayList(tipTimesList.size());
        tipTimesCountList = new ArrayList(tipTimesList.size());
        for (int i = 0; i < tipTimesTemp.size(); i++) {
            double[] newDouble = {((double[]) tipTimesTemp.get(i))[0]};
            tipTimesList.add(i,newDouble);
            int[] newInt = {1};
            tipTimesCountList.add(i,newInt);
            double[] timesTemp = (double[]) tipTimesTemp.get(i);
            for (int j = 1; j < timesTemp.length; j++) {
                boolean haploPresentAlready = false;
                for (int k = 0; k < tipTimesList.get(i).length; k++) {
                    if (tipTimesList.get(i)[k] == timesTemp[j]) {
                        tipTimesCountList.get(i)[k]++;
                        haploPresentAlready=true;
                        break;
                    }
                }
                if (!haploPresentAlready) {
                    double[] tipTimesTempArray = tipTimesList.get(i);
                    tipTimesList.set(i,new double[tipTimesTempArray.length+1]);
                    System.arraycopy(tipTimesTempArray,0,tipTimesList.get(i),0,tipTimesTempArray.length);
                    tipTimesList.get(i)[tipTimesTempArray.length] = timesTemp[j];
                    int[] tipTimesCountTempArray = tipTimesCountList.get(i);
                    tipTimesCountList.set(i,new int[tipTimesCountTempArray.length+1]);
                    System.arraycopy(tipTimesCountTempArray,0,tipTimesCountList.get(i),0,tipTimesCountTempArray.length);
                    tipTimesCountList.get(i)[tipTimesCountTempArray.length] = 1;
                }
            }
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
    private double[][] getSequenceDistances(Alignment data, Tree tree) {
        // Get the distances for the sequences:
        Distance distance = new DifferenceCount();
        ((Distance.Base) distance).setPatterns(data);
        // make a tip sequence distance matrix
        int taxaSize = data.getTaxonCount();
        double[][] distanceMatrix = new double[taxaSize][taxaSize];
        for (int i = 0; i < taxaSize - 1; i++) {
            for (int j = i + 1; j < taxaSize; j++) {
                distanceMatrix[i][j] = distance.pairwiseDistance(data.getTaxonIndex(tree.getNode(i).getID()),
                        data.getTaxonIndex(tree.getNode(j).getID()));
                distanceMatrix[j][i] = distanceMatrix[i][j];
            }
        }
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
    private int numberInternalNodes(Node node, int nextNr) {
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
            return string.replaceAll("&", "&amp;");
        } else{
            return getFlattenedHaploTree().getRoot().toSortedNewick(new int[1], true);
        }
    }

    /////////////////////////////////////////////////
    //           StateNode implementation          //
    /////////////////////////////////////////////////
    @Override
    protected void store() {

        System.arraycopy(parentHaplo, 0, storedParentHaplo, 0, parentHaplo.length);
//        Collections.copy(storedAboveNodeHaplo, aboveNodeHaplo);
        for (int i = 0; i<storedAttachmentTimesList.size(); i++){
            System.arraycopy(attachmentTimesList.get(i),0,storedAttachmentTimesList.get(i),0,storedAttachmentTimesList.get(i).length);
            System.arraycopy(tipTimesList.get(i),0,storedTipTimesList.get(i),0,storedTipTimesList.get(i).length);
            System.arraycopy(tipTimesCountList.get(i),0,storedTipTimesCountList.get(i),0,storedTipTimesCountList.get(i).length);
        }
//        Collections.copy(storedAttachmentTimesList, attachmentTimesList);
//        Collections.copy(storedTipTimesList,tipTimesList);
        System.arraycopy(startBranchCounts, 0, storedStartBranchCounts, 0, startBranchCounts.length);

//        final ArrayList<Double[]> tmp = attachmentTimesList;
//        attachmentTimesList = storedAttachmentTimesList;
//        storedAttachmentTimesList = tmp;
//        final ArrayList<QuasiSpeciesNode> tmpparent = parentHaplo;
//        parentHaplo = storedParentHaplo;
//        storedParentHaplo = tmpparent;
//        final ArrayList<QuasiSpeciesNode> tmphaploabove = aboveNodeHaplo;
//        aboveNodeHaplo = storedAboveNodeHaplo;
//        storedAboveNodeHaplo = tmphaploabove;
//        final int[] tmpstart = startBranchCounts;
//        startBranchCounts = storedStartBranchCounts;
//        storedStartBranchCounts = tmpstart;



        // from MultiTypeTree class
        storedRoot = m_storedNodes[root.getNr()];
        int iRoot = root.getNr();

        storeNodes(0, iRoot);

        ((QuasiSpeciesNode) storedRoot).setHeight(m_nodes[iRoot].getHeight(),false);
        if (this.getLeafNodeCount()>1){
            ((QuasiSpeciesNode) storedRoot).setParent(null, false);

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
        qsStoredRoot.haploAboveName = ((QuasiSpeciesNode)m_nodes[iRoot]).haploAboveName;
        qsStoredRoot.continuingHaploName = ((QuasiSpeciesNode)m_nodes[iRoot]).continuingHaploName;

        storeNodes(iRoot+1, nodeCount);
        // from MultiTypeTree class
    }


    /**
     * helper to store *
     */

    private void storeNodes(int iStart, int iEnd) {
        // from MultiTypeTree class
        for (int i = iStart; i<iEnd; i++) {
            QuasiSpeciesNode sink = (QuasiSpeciesNode)m_storedNodes[i];
            QuasiSpeciesNode src = (QuasiSpeciesNode)m_nodes[i];
            sink.setHeight(src.getHeight(),false);
            sink.setParent(m_storedNodes[src.getParent().getNr()], false);
            if (src.getLeft()!=null) {
                sink.setLeft(m_storedNodes[src.getLeft().getNr()]);
                if (src.getRight()!=null)
                    sink.setRight(m_storedNodes[src.getRight().getNr()]);
                else
                    sink.setRight(null);
            }
        // from MultiTypeTree class
            sink.setHaploAboveName(src.getHaploAboveName());
            sink.setContinuingHaploName(src.getContinuingHaploName());

            // TODO also copy all the other values like attach times, tip times and tip counts and parent haplo
        // from MultiTypeTree class
        }
        // from MultiTypeTree class
    }


    @Override
    public void restore() {
        super.restore();
        final ArrayList<double[]> tmp = storedAttachmentTimesList;
        storedAttachmentTimesList = attachmentTimesList;
        attachmentTimesList = tmp;
        final ArrayList<double[]> tmp2 = storedTipTimesList;
        storedTipTimesList = tipTimesList;
        tipTimesList = tmp2;
        final ArrayList<int[]> tmp3 = storedTipTimesCountList;
        storedTipTimesCountList = tipTimesCountList;
        tipTimesCountList = tmp3;
        final int[] tmpparent = storedParentHaplo;
        storedParentHaplo = parentHaplo;
        parentHaplo = tmpparent;
//        final ArrayList<QuasiSpeciesNode> tmphaploabove = storedAboveNodeHaplo;
//        storedAboveNodeHaplo = aboveNodeHaplo;
//        aboveNodeHaplo = tmphaploabove;
        final int[] tmpstart = storedStartBranchCounts;
        storedStartBranchCounts = startBranchCounts;
        startBranchCounts = tmpstart;
    }
    /////////////////////////////////////////////////
    // Methods implementing the Loggable interface //
    /////////////////////////////////////////////////
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

    @Override
    public void log(int i, PrintStream printStream) {
        printStream.print("tree STATE_"+i+" = ");
        printStream.print(toString());
        printStream.print(";");


    }

    @Override
    public void close(PrintStream printStream) {
        printStream.println("End;");
    }


    /////////////////////////////////////////////////
    // Serialization and deserialization for state //
    /////////////////////////////////////////////////

    /**
     * reconstruct tree from XML fragment in the form of a DOM node *
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

}