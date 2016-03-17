package quasispeciestree.tree;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;

import java.io.PrintStream;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author Veronika Boskova created on 29/04/2015
 */

public class QuasiSpeciesTree extends Tree {

// TODO: the haplotype count can be input as Integer counts or as percentage (in case of percentage check it sums to 100)
//    public Input<TraitSet> QuasiSpeciesAttachmentTimesInput =
//            new Input<TraitSet>("HaploAttachmentTimes","Attachment times for each haplotype (excluding the one representative of each haplotype in the tree input)");
    // not necessary, when not input -- only at init, then we set all values to the time of the parent (median of the tree length)
    public Input<TraitSet> haplotypeCountsInput =
            new Input<TraitSet>("haplotypeCounts","Count of sequences for each haplotype (excluding the one representative of each haplotype in the tree input)", Input.Validate.REQUIRED);
    // necessary - we need to know which tips map to which haplotype
    public Input<RealParameter> originInput = new Input<>(
            "origin", "The time from origin to last sample (must be larger than tree height)",
            Input.Validate.REQUIRED);

    protected ArrayList<Double[]> attachmentTimesList;
    protected ArrayList<Double[]> storedAttachmentTimesList;
    protected TraitSet haplotypeCounts;
    protected String qsLabel;
    private int[] startBranchCounts;
    private int[] storedStartBranchCounts;
    private int[] parentHaplo;
    // TODO from discussion with denise on 25.01.2016
    // private instead of protected???
    protected int[] storedParentHaplo;


    public QuasiSpeciesTree() { };

    protected RealParameter origin;


    public QuasiSpeciesTree(Node rootNode) {

        if (!(rootNode instanceof QuasiSpeciesNode))
            throw new IllegalArgumentException("Attempted to instantiate "
                    + "quasi-species tree with regular root node.");

        setRoot(rootNode);
        initArrays();
    }

    public void initAndValidate(){






// init and validate from scratch in order to implement the quasi-species node -- holding the haplotype starting above
        if (m_initial.get() != null && !(this instanceof StateNodeInitialiser)) {

            if (!(m_initial.get() instanceof QuasiSpeciesTree)) {
                throw new IllegalArgumentException("Attempted to initialise "
                        + "quasi-species tree with regular tree object.");
            }

            QuasiSpeciesTree other = (QuasiSpeciesTree)m_initial.get();
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
                root = (QuasiSpeciesNode) left;
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
        haplotypeCounts = haplotypeCountsInput.get();
        qsLabel = haplotypeCounts.getTraitName();

        initAttachmentTimes();

        fillParentHaplo();

//        System.arraycopy(parentHaplo,0,storedParentHaplo,0,parentHaplo.size());

        startBranchCounts = countPossibleStartBranches();
        storedStartBranchCounts = new int[startBranchCounts.length];
        System.arraycopy(startBranchCounts,0,storedStartBranchCounts,0,startBranchCounts.length);
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
    public void initAttachmentTimes(){
        // reserves space for array list of size LeafNodeCount
        attachmentTimesList = new ArrayList<Double[]>(this.getLeafNodeCount());
        storedAttachmentTimesList = new ArrayList<Double[]>(this.getLeafNodeCount());
        // set all haploAboveNames to -1 and all continuingHaploNames to -1
        for(Node node : this.getNodesAsArray()){
            ((QuasiSpeciesNode) node).setHaploAboveName(-1);
            ((QuasiSpeciesNode) node).setContinuingHaploName(-1);
        }
        // for those nodes where haplotype arises, change the haploAboveNode to haplotype's (corresponding tip node) number
        //  and for nodes below up to the tip set continuingHaploName to the same haplotype's (tip)
        for (Node node : this.getExternalNodes()){
            // check if getNr() always returns the same >>> Node number is guaranteed not to change during an MCMC run.
            //      written in the Node class)
            // assigns to each unique haplotype an array of size haplotypeCount
            Double[] tempqstimes = new Double[(int) haplotypeCounts.getValue(node.getID()) + 1];
            // start with the star tree, and define start of haplotype at the multi-furcating node time
            for (int i=0; i<=getHaplotypeCounts((QuasiSpeciesNode)node); i++) {
//                if (this.getLeafNodeCount()>1){
                    // TODO write a function that proposes RANDOM attachment times within given interval
                    // TODO also to be used by the operators -- > is this really needed?
                    
                    tempqstimes[i]=node.getParent().getHeight()//*0.9999999
                    //if spread the haplotype at start a bit
//                    -node.getParent().getHeight()*i*(1-0.9999999);
                    -(i+1)*((node.getParent().getHeight()-node.getHeight()) /(2+getHaplotypeCounts((QuasiSpeciesNode)node)));
// testing
//                  if (node.getID()==this.getExternalNodes().get(2).getID()){
//                      tempqstimes[i]=node.getParent().getParent().getHeight()*0.9999999
//                               -node.getParent().getParent().getHeight()*i*(1-0.9999999);
//                  }
//              }
//              for (int i=0; i<=getHaplotypeCounts((QuasiSpeciesNode)node); i++) {
//                  tempqstimes[i]=node.getParent().getHeight()-0.1
//                            //if spread the haplotype at start a bit
//                           -i*0.1;
//
//                  if (node.getID()==this.getExternalNodes().get(2).getID()){
//                      tempqstimes[i]=node.getParent().getParent().getHeight()-0.1
//                                -i*0.3;
//                  }
                    // TODO this is just for orig=MRCA  = 17 example
//                } else {
//                    tempqstimes[i]=17//*0.9999999
//                            -(i+1)*(17/(2+getHaplotypeCounts((QuasiSpeciesNode)node)));
////                    System.out.println("The tree has only one haplotype. This is not accepted by the current method implementation.");
////                    System.exit(0);
//                }

            }
            ((QuasiSpeciesNode) node).setHaploAboveName(node.getNr());
            ((QuasiSpeciesNode) node).setContinuingHaploName(node.getNr());
// testing function countPossibleStartBranches()
//            if (node.getID()==this.getExternalNodes().get(2).getID()){
//                System.out.println(node.getID());
//                ((QuasiSpeciesNode) node.getParent()).setHaploName(node.getID().toString());
//            } else {
//                ((QuasiSpeciesNode) node).setHaploName(node.getID().toString());
//            }
            // attachment time list defined as "node" height, going from present (0) to past (positive height)
            attachmentTimesList.add(node.getNr(), tempqstimes);
            storedAttachmentTimesList.add(node.getNr(), tempqstimes);
// testing
//          System.out.println((int)haplotypeCounts.getValue(node.getID()));
//          for (int i=0; i<attachmentTimesList.get(node.getNr()).length; i++)
//              {System.out.println(attachmentTimesList.get(node.getNr())[i]+" ");
//          }
        }
    }

    public void setAttachmentTimesList(QuasiSpeciesNode node, Double[] tempqstimes){
        node.setAttachmentTimesList();
        attachmentTimesList.set(node.getNr(),tempqstimes);
    }

    public Double[] getAttachmentTimesList(QuasiSpeciesNode node){
        return attachmentTimesList.get(node.getNr());
    }

    public void setAttachmentTimesList(int position, Double[] tempqstimes){
        ((QuasiSpeciesNode) this.getNode(position)).setAttachmentTimesList();
        attachmentTimesList.set(position,tempqstimes);
    }

    public Double[] getAttachmentTimesList(int position){
        return attachmentTimesList.get(position);
    }

    public Double getTotalAttachmentTimes(QuasiSpeciesNode node){
        Double totalTime=0.0;
        Double[] timesList=attachmentTimesList.get(node.getNr());
        for (int i=0; i<timesList.length; i++){
            totalTime += timesList[i] ;
        }
        return totalTime;
    }

    public double getHaplotypeCounts(QuasiSpeciesNode node){
        return haplotypeCounts.getValue(node.getID());
    }

    public int getTotalAttachmentCounts(){
        int totalCount=0;
        for (Node node : this.getExternalNodes()){
            totalCount += ( haplotypeCounts.getValue(node.getID())) ;
        }
        return totalCount;
    }

    public int[] getStartBranchCounts(){
        return startBranchCounts;
    }

    public void  setStartBranchCounts(int[] startBranchCountsArray){
        ((QuasiSpeciesNode) this.getRoot()).setStartBranchCounts();
        startBranchCounts = startBranchCountsArray;
    }

    public int[] getParentHaplo() {return parentHaplo; }

    public int getParentHaplo(int position) { return parentHaplo[position]; }

    public void setParentHaplo(int[] newParentHaploArray) {
        ((QuasiSpeciesNode) this.getRoot()).setParentHaplo();
        parentHaplo = newParentHaploArray;
    }

    public void clearContinuingHaploNames() {
        for (Node node : this.getNodesAsArray()){
            ((QuasiSpeciesNode) node).setContinuingHaploName(-1);
        }
    }
    public void clearContinuingHaploNames(Node belowThisNode) {
        for (Node node : belowThisNode.getAllChildNodes()){
            ((QuasiSpeciesNode) node).setContinuingHaploName(-1);
        }
    }

    // clear and add method for parent haplo



//    /**
//     * Method to propose haplotype start time and attachment times for each of its copies
//     * for each haplotype
//     *
//     * @return Array holding for each haplotype (#this.getExternalNodes()) the times of start
//     *         and attachemtn times, held in array at position determined by node.getNr()
//     */
//    public ArrayList<Double[]> proposeHaploStartAndAttachemntTimes(QuasiSpeciesNode node, double time){
//
//
//    }







    /**
     * Method to count for a given sequence ihaploseq the possible number of attachment
     * branches of a given (same or different) haplotype
     *
     * @param node Name of the haplotype from which the haplotype sequence i is supposed to branch off
     * @param ihaploseq For which sequence (i) exactly from haplotype (node) are we estimating the possible
     *                  attachment branches?
     *                  This is important if we are trying to list all possible branches
     *                  of haplotype (node) from which the sequence can branch off, ie. deduct 1 if the
     *                  sequence age is identical to the ihaploseqage, and position of that haplotype
     *                  sequence in AttachmentTimesList is at ihaploseq.
     *                  If the haplotype i is different from haplotype at node, then set
     *                  ihaploseq=0.
     * @param haploseqage Age of a sequence for which we want to check
     *                     how many instances of haplotype (node) it can arise from
     *                     -- important if we are checking all possible branches of haplotype
     *                     (node) from which the haplotype number i can branch off, then
     *                     set haplotypeage==getAttachmentTimesList(node)[i]
     * @return Number of possible attachment branches, at position determined by node.getNr()
     *
     */
    public int countPossibleAttachmentBranches(QuasiSpeciesNode node, int ihaploseq, double haploseqage){
        // there is always one possible attachment --- the original haplotype
        int abcount = 1;
        // if node's haplotype is older than haploseqage, count how many sequences of node's haplotype are there
        // if node's haplotype is not older than haploseqage --- throw error
        Double[] AttachmentTimesTemp=getAttachmentTimesList(node);
        double nodehaplorootage = AttachmentTimesTemp[0];

        // we allow also testing of other times than the ones belonging to the node's haplotype
        if (haploseqage>nodehaplorootage){
            System.out.println("The haplotype "+ node.getID() + " of age " + nodehaplorootage
                    + " is set to start after the haplotype sequence number " + ihaploseq + " of age " + haploseqage +
                    " is set to branch off");
            System.exit(0);
        } else {
            for (int i=1; i < AttachmentTimesTemp.length; i++){
                // We do NOT allow branching from the branch with the starting starting age!!!
                if ( haploseqage < (double) AttachmentTimesTemp[i] && i != ihaploseq ){
                    abcount +=1;
                }
            }
// testing
//            System.out.println("The haplotype count " + node.getID() + " " + ihaploseq + " = " + abcount);
        }
        return abcount;
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
        int nTips=this.getLeafNodeCount();
//        // array to store the count of possible start branches for each haplotype
//        int[] startBranchCountArray = new int[nTips];
        // array to store the count of possible start branches for each internal node
        int[] startBranchCountArray = new int[this.getInternalNodeCount()];
        for (int i=0; i<this.getInternalNodeCount(); i++){
            // once determined the parent haplotype find out from how many branches
            // of the parent haplotype can it actually branch off
            if (((QuasiSpeciesNode) this.getNode(nTips + i)).getContinuingHaploName() != -1){
                startBranchCountArray[i]=countPossibleAttachmentBranches(
                        (QuasiSpeciesNode) this.getNode(((QuasiSpeciesNode) this.getNode(nTips + i)).getContinuingHaploName()),
                        0, this.getNode(nTips + i).getHeight());
            } else {
                startBranchCountArray[i]=1;
            }
// testing
//            System.out.println("Number of possible attachment branches of " + this.getNode(nTips+i).getID()
//            + " = " + startBranchCountArray[i]);
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
        // array to store the parent haplotype of each internal node, if the first haplotype, the parent haplotype is -1
        parentHaplo = new int[this.getLeafNodeCount()];
        storedParentHaplo = new int[this.getLeafNodeCount()];
        Arrays.fill(parentHaplo, -1);
        Arrays.fill(storedParentHaplo, -1);

        // starting at the root of the tree, see what is the order of the haplotype, always the haplotype below,
        // starts from the haplotype above (parent)
        //
        // if there is already haplotype arising at root-origin branch, this will have to arise from NULL (-1) haplotype anyways
        // and by default we set each entry of parentHaplo array to -1
        findParentHaplo(-1, (QuasiSpeciesNode) root, parentHaplo);

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
     * @param parentHaplo Original array holding parent haplotype for each haplotype (external node) in the tree,
     *                if not parent haplotype ... this is by default -1
//     * @param aboveNodeHaplo Original array holding node the haplotype arises above for each haplotype in the tree,
//     *                if not parent haplotype ... this is by default -1
     */
    public void findParentHaplo(int currentHaploType, QuasiSpeciesNode nextNode, int[] parentHaplo){
        // get haplotype starting at a branch directly above next node --- if no haplotype (null), check children nodes
        // check whether there is a new haplotype arising
//        int nTips=this.getLeafNodeCount();
        if (nextNode.getHaploAboveName()!=-1){
            // check whether the new haplotype is on the same or different branch than that leading to the parent haplotype
            // TODO can we do kind of binary search for String here instead of for loop???
            int newHaploType = -1;
            newHaploType = nextNode.getHaploAboveName();
            // set the haplotype from which the next haplotype arises as a parent
            parentHaplo[newHaploType] = currentHaploType;
            // set the continuingHaploName at corresponding internal nodes leading from nextNode to the respective child\
            nextNode.setContinuingHaploName(newHaploType);
            QuasiSpeciesNode settingContinuousHaplo = (QuasiSpeciesNode) this.getNode(newHaploType);
            while ( nextNode != settingContinuousHaplo ) {
                settingContinuousHaplo.setContinuingHaploName(newHaploType);
                settingContinuousHaplo = (QuasiSpeciesNode) settingContinuousHaplo.getParent();
            }

// testing
//            if (currentHaploType != -1) {
//                System.out.println("Parent haplotype of " + this.getNode(newHaploType).getID() +
//                    " is " + this.getNode(currentHaploType).getID());
//            }

            currentHaploType=newHaploType;
        }
        if (nextNode.getChildCount() > 0 ) {
            // apply the same getParentHaplo function to the child nodes, with updated current haplotype
            for (Node childNode : nextNode.getChildren()) {
                findParentHaplo(currentHaploType, (QuasiSpeciesNode) childNode, parentHaplo);
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
        m_nodes = new QuasiSpeciesNode[nodeCount];
        listNodes((QuasiSpeciesNode)root, (QuasiSpeciesNode[])m_nodes);
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
            if (node.getRight()!=null)
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
        tree.root = root.copy();
        tree.nodeCount = nodeCount;
        tree.internalNodeCount = internalNodeCount;
        tree.leafNodeCount = leafNodeCount;
        tree.qsLabel = qsLabel;
        tree.attachmentTimesList = attachmentTimesList;
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

        attachmentTimesList = qsTree.attachmentTimesList;

        QuasiSpeciesNode[] qsNodes = new QuasiSpeciesNode[qsTree.getNodeCount()];
        for (int i=0; i<qsTree.getNodeCount(); i++)
            qsNodes[i] = new QuasiSpeciesNode();

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
     *
     * @param other
     */
    @Override
    public void assignFromFragile(StateNode other) {
        QuasiSpeciesTree qsTree = (QuasiSpeciesTree) other;

//        attachmentTimesList = qsTree.attachmentTimesList;

        if (m_nodes == null) {
            initArrays();
        }
        root = (QuasiSpeciesNode) m_nodes[qsTree.root.getNr()];
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
     * Initialise tree topology from Tree object
     *
     * @param regularTree
     * @throws java.lang.Exception
     */
    public void initFromRegularTree(Tree regularTree){

        // Build new quasi-species tree:

        List<Node> activeRegularTreeNodes = new ArrayList<>();
        List<Node> nextActiveRegularTreeNodes = new ArrayList<>();
        List<QuasiSpeciesNode> activeTreeNodes = new ArrayList<>();
        List<QuasiSpeciesNode> nextActiveTreeNodes = new ArrayList<>();

        // Populate active node lists with root:
        activeRegularTreeNodes.add(regularTree.getRoot());
        QuasiSpeciesNode newRoot = new QuasiSpeciesNode();
        activeTreeNodes.add(newRoot);

        // Initialise counter used to number leaves when takeNrsFromRegularTree
        // is false:
        int nextNr = 0;

        while (!activeRegularTreeNodes.isEmpty()) {

            nextActiveRegularTreeNodes.clear();
            nextActiveTreeNodes.clear();

            for (int idx = 0; idx<activeRegularTreeNodes.size(); idx++) {
                Node regularTreeNode = activeRegularTreeNodes.get(idx);
                QuasiSpeciesNode treeNode = activeTreeNodes.get(idx);

                switch (regularTreeNode.getChildCount()) {
                    case 0:
                        // Leaf at base of branch
                        treeNode.setNr(nextNr);
                        treeNode.setID(String.valueOf(regularTreeNode.getID()));
                        // Set start of haplotype times set as default to belong to the leaf node
                        // treeNode.setHaploname(treeNode.getID());
                        // done in initAndValidate!!!

                        nextNr += 1;

                        break;

                    case 2:
                        // Non-leaf at base of branch
                        nextActiveRegularTreeNodes.add(regularTreeNode.getLeft());
                        nextActiveRegularTreeNodes.add(regularTreeNode.getRight());

                        QuasiSpeciesNode daughter = new QuasiSpeciesNode();
                        QuasiSpeciesNode son = new QuasiSpeciesNode();
                        treeNode.addChild(daughter);
                        treeNode.addChild(son);
                        nextActiveTreeNodes.add(daughter);
                        nextActiveTreeNodes.add(son);
                        // Set start of haplotype times set as default to belong to the leaf node
                        // do nothing here

                        break;
                }

                // Set node height:
                treeNode.setHeight(regularTreeNode.getHeight());
            }

            // Replace old active node lists with new:
            activeRegularTreeNodes.clear();
            activeRegularTreeNodes.addAll(nextActiveRegularTreeNodes);

            activeTreeNodes.clear();
            activeTreeNodes.addAll(nextActiveTreeNodes);

        }


        // Number internal nodes:
        numberInternalNodes(newRoot, newRoot.getAllLeafNodes().size());

        // Assign tree topology:
        assignFromWithoutID(new QuasiSpeciesTree(newRoot));
        initArrays();

    }

    /**
     * Helper method used by initFromRegularTree to assign sensible node numbers
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

//    /**
//     * Return string representation of quasi-species tree.  We use reflection
//     * here to determine whether this is being called as part of writing
//     * the state file.
//     *
//     * @return Quasi-species tree string in Newick format.
//     */
//    @Override
//    public String toString() {
//
//        // Behaves differently if writing a state file
//        StackTraceElement[] ste = Thread.currentThread().getStackTrace();
//        if (ste[2].getMethodName().equals("toXML")) {
//            // Use toShortNewick to generate Newick string without taxon labels
//            String string = getFlattenedTree(false).getRoot().toShortNewick(true);
//
//            // Sanitize ampersands if this is destined for a state file.
//            return string.replaceAll("&", "&amp;");
//        } else{
//            return getFlattenedTree(true).getRoot().toSortedNewick(new int[1], true);
//        }
//    }
    // would also need to change and import from MultiTypeTree class the getFlattenedTree function

    /////////////////////////////////////////////////
    //           StateNode implementation          //
    /////////////////////////////////////////////////
    @Override
    protected void store() {

        System.arraycopy(parentHaplo, 0, storedParentHaplo, 0, parentHaplo.length);
//        Collections.copy(storedAboveNodeHaplo, aboveNodeHaplo);
        Collections.copy(storedAttachmentTimesList, attachmentTimesList);
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
        storedRoot = (QuasiSpeciesNode) m_storedNodes[root.getNr()];
        int iRoot = root.getNr();

        storeNodes(0, iRoot);

        ((QuasiSpeciesNode) storedRoot).setHeight(m_nodes[iRoot].getHeight(),false);
        ((QuasiSpeciesNode) storedRoot).setParent(null, false);

        if (root.getLeft()!=null)
            storedRoot.setLeft(m_storedNodes[root.getLeft().getNr()]);
        else
            storedRoot.setLeft(null);
        if (root.getRight()!=null)
            storedRoot.setRight(m_storedNodes[root.getRight().getNr()]);
        else
            storedRoot.setRight(null);

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
            sink.haploAboveName = src.getHaploAboveName();
            sink.continuingHaploName = src.getContinuingHaploName();
        // from MultiTypeTree class
        }
        // from MultiTypeTree class
    }


    @Override
    public void restore() {
        super.restore();
        final ArrayList<Double[]> tmp = storedAttachmentTimesList;
        storedAttachmentTimesList = attachmentTimesList;
        attachmentTimesList = tmp;
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

//    @Override
//    public void log(int i, PrintStream printStream) {
//        printStream.print("tree STATE_"+i+" = ");
//        printStream.print(toString());
//        printStream.print(";");
//
//
//    }

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
            String sNewick = node.getTextContent().replace("&", "");

            TreeParser parser = new TreeParser();
            parser.initByName(
                    "IsLabelledNewick", false,
                    "offset", 0,
                    "adjustTipHeights", false,
                    "singlechild", true,
                    "newick", sNewick);
            //parser.m_nThreshold.setValue(1e-10, parser);
            //parser.m_nOffset.setValue(0, parser);

            initFromRegularTree(parser);

            initArrays();
        } catch (Exception ex) {
            Logger.getLogger(QuasiSpeciesTree.class.getName()).log(Level.SEVERE, null, ex);
        }
    }



}

