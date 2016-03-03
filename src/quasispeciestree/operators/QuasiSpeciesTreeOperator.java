package quasispeciestree.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesTree;
import quasispeciestree.tree.QuasiSpeciesNode;
import beast.core.parameter.RealParameter;

import java.util.ArrayList;
import java.util.List;


/**
 *  @author Veronika Boskova created on 23/07/2015
 */
@Description("This operator generates proposals for a quasi-species tree.")
public abstract class QuasiSpeciesTreeOperator extends Operator {

    public Input<QuasiSpeciesTree> quasiSpeciesTreeInput = new Input<>(
            "quasiSpeciesTree", "Quasi-species tree on which to operate.",
            Input.Validate.REQUIRED);

    public Input<RealParameter> originInput = new Input<>(
            "origin", "The time from origin to last sample (must be larger than tree height)",
            Input.Validate.REQUIRED);

    protected QuasiSpeciesTree qsTree;
    protected RealParameter origin;

    @Override
    public void initAndValidate() throws Exception {
        qsTree = quasiSpeciesTreeInput.get();
        origin = originInput.get();
    }

    /* ***********************************************************************
     * The following two methods are copied verbatim from TreeOperator.
     */

    /**
     * Obtain the sister of node "child" having parent "parent".
     *
     * @param parent the parent
     * @param child  the child that you want the sister of
     * @return the other child of the given parent.
     */
    protected Node getOtherChild(Node parent, Node child) {
        if (parent.getLeft().getNr() == child.getNr()) {
            return parent.getRight();
        } else {
            return parent.getLeft();
        }
    }

    /**
     * Replace node "child" with another node.
     *
     * @param node
     * @param child
     * @param replacement
     */
    public void replace(Node node, Node child, Node replacement) {
        node.removeChild(child);
        node.addChild(replacement);
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        replacement.makeDirty(QuasiSpeciesTree.IS_FILTHY);
    }

    /* **********************************************************************
    * The following 4 methods are adapted from MultiTypeTreeOperator.
    */

    /**
     * Disconnect edge <node,node.getParent()> by joining node's sister directly
     * to node's grandparent
     * If the haplotype in the sub-tree to be moved, i.e. sub-tree below the node,
     * does arise above the node.getParent(), choose new attachment for the haplo
     * If the haplotype in the sister sub-tree does arise above node.getParent(),
     * no need to change anything but the node.getParent().haploAboveName and node.continuingHaploName
     * and sister.haploAboveName
     *
     * @param node to detach
     * @param haplotype to detach (change haploAboveName for the corresponding node)
     * return the QuasiSpeciesNode (parent/grandparent/.../root) above which the haplotype of the subtree to be moved arises
     *                  (and for/below which the continuingHaploName have to be corrected)
     */
    public QuasiSpeciesNode disconnectBranch(QuasiSpeciesNode node, int haplotype) {

        // Check argument validity:
        QuasiSpeciesNode parent = (QuasiSpeciesNode) node.getParent();

        if (node.isRoot() || parent.isRoot())
            throw new IllegalArgumentException("Illegal argument to "
                    + "disconnectBranch().");

        QuasiSpeciesNode sister = (QuasiSpeciesNode) getOtherChild(parent, node);

        // keep track of where does the haplotype arise
        QuasiSpeciesNode nodeBelowHaploMoved = null;

        // Disconnect the haplotype chosen to be moved when disconnecting the node
        // if the haplotype to be moved is already present at node, then it must have arisen above in the tree
        if (haplotype == node.getContinuingHaploName()){
            // first check whether it arises at the node itself
            if (node.getHaploAboveName() == haplotype){
                node.setHaploAboveName(-1);
                if (!node.isLeaf())
                    node.setContinuingHaploName(-1);
                nodeBelowHaploMoved = null;
            }
            // second see whether there is continuing haplotype at the parent node
            else {
                // check at which parent it arises
                QuasiSpeciesNode nodeAbove = parent;
                while (nodeAbove.getHaploAboveName() != haplotype){
                    nodeAbove = (QuasiSpeciesNode) nodeAbove.getParent();
                }
                // set the node above which the haploMoved arises to -1, as haploMoved will be moved away
                nodeAbove.setHaploAboveName(-1);
                nodeAbove.setContinuingHaploName(-1);
                nodeBelowHaploMoved = nodeAbove;
                if (nodeBelowHaploMoved == node.getParent())
//                    nodeBelowHaploMoved = null;
                    nodeBelowHaploMoved = sister;
            }
        }
        // otherwise the haplotype to be moved must have arisen below the node (need to clear the haploAboveName)
        else{
            // start checking at the tip
            QuasiSpeciesNode nodeBelow = (QuasiSpeciesNode) qsTree.getNode(haplotype);
            while (nodeBelow.getHaploAboveName() != haplotype){
                nodeBelow = (QuasiSpeciesNode) nodeBelow.getParent();
            }
            // set the node above which the haploMoved arises to -1, as haploMoved will be moved away
            nodeBelow.setHaploAboveName(-1);
            if (!nodeBelow.isLeaf())
                nodeBelow.setContinuingHaploName(-1);
            nodeBelowHaploMoved = null;
        }

        // Add haplotype start/attachment times originally attached to parent to those attached
        // to node's sister:
        if ( parent.getHaploAboveName() != -1 ){
            if ( sister.getHaploAboveName() != -1 ){
                // if there is one haplotype arising on the parent's branch and another on sister's branch
                // there is a problem somewhere in the code!!!
                System.out.println("There must be a serious problem with the code. There are haplotypes arising at " +
                        "both parent and sister branch, while this should technically be impossible...");
                System.exit(0);
            }
            else {
                int haploAboveParentNode = parent.getHaploAboveName();
                sister.setHaploAboveName(haploAboveParentNode);
                sister.setContinuingHaploName(haploAboveParentNode);
                parent.setHaploAboveName(-1);
                parent.setContinuingHaploName(-1);
            }
        }

        // Implement topology change.
        parent.removeChild(sister);
        replace(parent.getParent(), parent, sister);
        sister.setParent(parent.getParent());


        // Ensure BEAST knows to update affected likelihoods:
        parent.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        sister.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);

        // return from which node is the haplotype moved --- we have to re-attach and rescale it!
        // and re-assign continuing haplotypes/parenthaploarray from parent node onwards
        return nodeBelowHaploMoved;
    }

    /**
     * Disconnect node from root
     * * If the haplotype in the sub-tree to be moved, i.e. sub-tree below the node,
     *   does arise at the <node,root> branch, record it and need to be reassigned
     *   starting point at the destination branch.
     *
     * @param node to detach
     * @param haplotype to detach (change haploAboveName for the corresponding node)
     * return the QuasiSpeciesNode (parent/grandparent/.../root) above which the haplotype of the subtree to be moved arises
     *                  (and for/below which the continuingHaploName have to be corrected)
     */
    public QuasiSpeciesNode disconnectBranchFromRoot(QuasiSpeciesNode node, int haplotype) {

        // Check argument validity:
        QuasiSpeciesNode parent = (QuasiSpeciesNode) node.getParent();

        if (node.isRoot() || !parent.isRoot())
            throw new IllegalArgumentException("Illegal argument to "
                    + "disconnectBranchFromRoot().");

        QuasiSpeciesNode sister = (QuasiSpeciesNode) getOtherChild(parent, node);

        // keep track of where does the haplotype arise
        QuasiSpeciesNode nodeBelowHaploMoved = null;

        // Disconnect the haplotype chosen to be moved when disconnecting the node
        // if the haplotype to be moved is already present at node, then it must have arisen above in the tree
        if (haplotype == node.getContinuingHaploName()){
            // first check whether it arises at the node itself
            if (node.getHaploAboveName() == haplotype){
                node.setHaploAboveName(-1);
                if (!node.isLeaf())
                    node.setContinuingHaploName(-1);
                nodeBelowHaploMoved = null;
            }
            // second see whether there is continuing haplotype at the parent node
            else {
                // check at which parent it arises
                QuasiSpeciesNode nodeAbove = parent;
                while (nodeAbove.getHaploAboveName() != haplotype){
                    nodeAbove = (QuasiSpeciesNode) nodeAbove.getParent();
                }
                // set the node above which the haploMoved arises to -1, as haploMoved will be moved away
                nodeAbove.setHaploAboveName(-1);
                nodeAbove.setContinuingHaploName(-1);
                nodeBelowHaploMoved = nodeAbove;
                // if nodeBelowHaploMoved is not root, then we have got a problem
                if (! nodeBelowHaploMoved.isRoot())
                    throw new IllegalArgumentException(
                            "Somehow the parent node of the haplotype to be detached"
                            + " is not the same as the root, and we are using"
                            + " disconnectBranchFromRoot() method.");
                if (nodeBelowHaploMoved == node.getParent())
//                    nodeBelowHaploMoved = null;
                    nodeBelowHaploMoved = sister;
            }
        }
        // otherwise the haplotype to be moved must have arisen below the node (need to clear the haploAboveName)
        else{
            // start checking at the tip
            QuasiSpeciesNode nodeBelow = (QuasiSpeciesNode) qsTree.getNode(haplotype);
            while (nodeBelow.getHaploAboveName() != haplotype){
                nodeBelow = (QuasiSpeciesNode) nodeBelow.getParent();
            }
            // set the node above which the haploMoved arises to -1, as haploMoved will be moved away
            nodeBelow.setHaploAboveName(-1);
            if (!nodeBelow.isLeaf())
                nodeBelow.setContinuingHaploName(-1);
            nodeBelowHaploMoved = null;
        }

        // Add haplotype start/attachment times originally attached to parent to those attached
        // to node's sister:
        if ( parent.getHaploAboveName() != -1 ){
            if ( sister.getHaploAboveName() != -1 ){
                // if there is one haplotype arising on the parent's branch and another on sister's branch
                // there is a problem somewhere in the code!!!
                System.out.println("There must be a serious problem with the code. There are haplotypes arising at " +
                        "both parent and sister branch, while this should technically be impossible...");
                System.exit(0);
            }
            else {
                int haploAboveParentNode = parent.getHaploAboveName();
                sister.setHaploAboveName(haploAboveParentNode);
                sister.setContinuingHaploName(haploAboveParentNode);
                parent.setHaploAboveName(-1);
                parent.setContinuingHaploName(-1);
            }
        }

        // Implement topology change:
        sister.setParent(null, true);
        parent.removeChild(sister);

        // Ensure BEAST knows to update affected likelihoods:
        parent.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        sister.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);

        // return from which node is the haplotype moved --- we have to re-attach and rescale it!
        // and re-assign continuing haplotypes/parenthaploarray from parent node onwards
        return nodeBelowHaploMoved;
    }

    /**
     * Creates a new branch between node and a new node at time destTime between
     * destBranchBase and its parent and sets the haploAboveName for the appropriate node.
     *
     * @param node
     * @param destBranchBase
     * @param destTime time at which the node attaches
     * @param haplotype the haplotype that has to attach to new destination
     * @param haploAttachTime the time at which the new haplotype attaches at new destination
     * return the QuasiSpeciesNode above which the moved haplotype is arising now
     *                          (and for/below which the continuingHaploName have to be corrected)
     */
    public QuasiSpeciesNode connectBranch(QuasiSpeciesNode node, QuasiSpeciesNode destBranchBase,
                                          double destTime, int haplotype, double haploAttachTime) {

        // Check argument validity:
        if (node.isRoot() || destBranchBase.isRoot())
            throw new IllegalArgumentException("Illegal argument to "
                    + "connectBranch().");

        // Obtain existing parent of node and set new time:
        QuasiSpeciesNode parent = (QuasiSpeciesNode) node.getParent();
        parent.setHeight(destTime);

        // Implement topology changes:
        replace(destBranchBase.getParent(), destBranchBase, parent);
        destBranchBase.setParent(parent);

        if (parent.getLeft() == node)
            parent.setRight(destBranchBase);
        else if (parent.getRight() == node)
            parent.setLeft(destBranchBase);

        // Ensure BEAST knows to update affected likelihoods:
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        parent.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        destBranchBase.makeDirty(QuasiSpeciesTree.IS_FILTHY);

        // if the new parent node's height is above the sister node continuingHaplo assign also to the parent
        // check whether the haploAboveNode is also above the node
        if (destBranchBase.getContinuingHaploName() != -1)
        {
            int destHaplo = destBranchBase.getContinuingHaploName();
            //check if it starts above the newly created parent node
            if (destTime < qsTree.getAttachmentTimesList((QuasiSpeciesNode) qsTree.getNode(destHaplo))[0]){
                if (destBranchBase.getHaploAboveName() != -1){
                    destBranchBase.setHaploAboveName(-1);
                    parent.setHaploAboveName(destHaplo);
                }
                parent.setContinuingHaploName(destHaplo);
            }
        }

        // set the aboveNodeHaplo for the node above which the haplotype moved arises to haplotypeNumber
        QuasiSpeciesNode nodeBelowHaploMoved = null;
        if (haploAttachTime >= destTime){
            nodeBelowHaploMoved = parent;
            while (nodeBelowHaploMoved != qsTree.getRoot() && haploAttachTime > nodeBelowHaploMoved.getParent().getHeight()){
                nodeBelowHaploMoved = (QuasiSpeciesNode) nodeBelowHaploMoved.getParent();
            }
            nodeBelowHaploMoved.setHaploAboveName(haplotype);
        }
        else{
            nodeBelowHaploMoved = (QuasiSpeciesNode) qsTree.getNode(haplotype);
            while (haploAttachTime > nodeBelowHaploMoved.getParent().getHeight()){
                nodeBelowHaploMoved = (QuasiSpeciesNode) nodeBelowHaploMoved.getParent();
            }
            nodeBelowHaploMoved.setHaploAboveName(haplotype);
            nodeBelowHaploMoved = parent;
        }

        // return the node for which the haploAboveName changed
        return nodeBelowHaploMoved;
    }

    /**
     * Set up node's parent as the new root with a height of destTime, with
     * oldRoot as node's new sister.
     *
     * @param node
     * @param oldRoot
     * @param destTime time at which the node attaches
     * @param haplotype the haplotype that has to attach to new destination
     * @param haploAttachTime the time at which the new haplotype attaches at new destination
     * return the QuasiSpeciesNode above which the moved haplotype is arising now
     *                          (and for/below which the continuingHaploName have to be corrected)
     */
    public QuasiSpeciesNode connectBranchToRoot(QuasiSpeciesNode node, QuasiSpeciesNode oldRoot,
                                    double destTime, int haplotype, double haploAttachTime) {

        // Check argument validity:
        if (node.isRoot() || !oldRoot.isRoot())
            throw new IllegalArgumentException("Illegal argument "
                    + "to connectBranchToRoot().");

        // Obtain existing parent of node and set new time:
        QuasiSpeciesNode newRoot = (QuasiSpeciesNode) node.getParent();
        newRoot.setHeight(destTime);

        // Implement topology changes:
        newRoot.setParent(null);

        if (newRoot.getLeft() == node)
            newRoot.setRight(oldRoot);
        else if (newRoot.getRight() == node)
            newRoot.setLeft(oldRoot);

        oldRoot.setParent(newRoot);

        // Ensure BEAST knows to recalculate affected likelihood:
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        newRoot.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        oldRoot.makeDirty(QuasiSpeciesTree.IS_FILTHY);

        // if the new parent node's height is above the sister node continuingHaplo assign it to the parent
        // check whether the haploAboveNode is also above the node
        if (oldRoot.getContinuingHaploName() != -1)
        {
            int destHaplo = oldRoot.getContinuingHaploName();
            //check if it starts below the newly created parent node
            if (destTime < qsTree.getAttachmentTimesList((QuasiSpeciesNode) qsTree.getNode(destHaplo))[0]){
                if (oldRoot.getHaploAboveName() != -1){
                    oldRoot.setHaploAboveName(-1);
                    newRoot.setHaploAboveName(destHaplo);
                }
                newRoot.setContinuingHaploName(destHaplo);
            }
        }

        // set the aboveNodeHaplo for the node above which the haplotype moved arises to haplotypeNumber
        QuasiSpeciesNode nodeBelowHaploMoved = null;
        if (haploAttachTime >= destTime){
            nodeBelowHaploMoved = newRoot;
            nodeBelowHaploMoved.setHaploAboveName(haplotype);
        }
        else{
            nodeBelowHaploMoved = (QuasiSpeciesNode) qsTree.getNode(haplotype);
            while (haploAttachTime > nodeBelowHaploMoved.getParent().getHeight()){
                nodeBelowHaploMoved = (QuasiSpeciesNode) nodeBelowHaploMoved.getParent();
            }
            nodeBelowHaploMoved.setHaploAboveName(haplotype);
            nodeBelowHaploMoved = newRoot;
        }


        // return the node for which the haploAboveName changed
        return nodeBelowHaploMoved;
    }

    /*
    //
    //
    //          OWN FUNCTIONS
    //
    //
    */
    /**
     * Function to find a most recent common ancestor of two nodes
     *
     * @param node1 tip node 1
     * @param node2 tip node 2
     * @param startInternalNode internal node from which to start the search, going down the tree towards to tips
     * return the QuasiSpeciesNode that is the MRCA for two haplotypes
     */
    public QuasiSpeciesNode findLastCommonAncestor(QuasiSpeciesNode node1, QuasiSpeciesNode node2, QuasiSpeciesNode startInternalNode){
        QuasiSpeciesNode returnnode = startInternalNode;

        boolean node1found = false;
        boolean node2found = false;

        for (Node tipnode : startInternalNode.getAllLeafNodes()){
            if (tipnode.getID() == node1.getID()){
                node1found = true;
            }
            if (tipnode.getID() == node2.getID()){
                node2found = true;
            }
        }
        if (node1found && node2found){
            QuasiSpeciesNode returnnode1 = findLastCommonAncestor(node1,node2, (QuasiSpeciesNode) startInternalNode.getLeft());
            QuasiSpeciesNode returnnode2 = findLastCommonAncestor(node1,node2, (QuasiSpeciesNode) startInternalNode.getRight());
            if (returnnode1 == returnnode2){
                returnnode = returnnode1;
            }
            else{
                if(returnnode1 == startInternalNode){
                    returnnode = returnnode2;
                }
                else if(returnnode2 == startInternalNode){
                    returnnode = returnnode1;
                }
                else{
                    System.out.println("Check your findLastCommonAncestor function!");
                    System.exit(0);
                }
            }
        }
//        else if ((node1found && !node2found) || (!node1found && node2found)){
        else{
            returnnode = (QuasiSpeciesNode) startInternalNode.getParent();
        }
        return (QuasiSpeciesNode) returnnode;
    }

    /**
     * Function to find a most recent common ancestor of two nodes
     */
    public void recalculateParentHaploAndCorrectContinuingHaploName(int currentHaplo, QuasiSpeciesNode nextNode){
        // recalculate parentHaplo array, re-assign continuingHaploName
        int[] newParentHaplo = new int[qsTree.getLeafNodeCount()];
        //Arrays.fill(newParentHaplo, -1);
        System.arraycopy(qsTree.getParentHaplo(), 0, newParentHaplo, 0, newParentHaplo.length);
        for(Node node : nextNode.getAllLeafNodes()){
            newParentHaplo[node.getNr()]=-1;
        }
        if (currentHaplo != -1){
            newParentHaplo[currentHaplo]=qsTree.getParentHaplo(currentHaplo);
        }

        qsTree.clearContinuingHaploNames(nextNode);
        qsTree.findParentHaplo(currentHaplo, nextNode, newParentHaplo);

        // re-set the original parentHaplo array
        qsTree.setParentHaplo(newParentHaplo);
    }

    /**
     * Function to find a maximum height until which the haplo can attach to (start from)
     *
     * @param currentHaplo the haplotype for which we search how far in the current tree it can start
     * @param startNode the node at which to start the search going up towards the origin
     * return the maximum time up to which the current haplo can be attaching
     */
    public double getMaxPossibleHaploAttachTime(int currentHaplo, QuasiSpeciesNode startNode){
        // starting from the startNode look for parent's continuingHaploName and if present
        // that node is the max possible attachment time of the haplotype
        QuasiSpeciesNode nodeToCheck = startNode;
        while (nodeToCheck.getContinuingHaploName() == -1
                || nodeToCheck.getContinuingHaploName() == currentHaplo
//                || nodeToCheck.getHaploAboveName() == currentHaplo
                ){
            if(nodeToCheck == qsTree.getRoot())
                break;
            nodeToCheck = (QuasiSpeciesNode) nodeToCheck.getParent();
        }
        if (nodeToCheck == qsTree.getRoot() &&
                (nodeToCheck.getContinuingHaploName() == currentHaplo || nodeToCheck.getContinuingHaploName() == -1))
            return origin.getValue();
        else
            return nodeToCheck.getHeight();
    }

    /**
     * Function to find a maximum height until which the haplo can attach to (start from)
     *
     * @param startNode the node at which to start the search going up towards the origin
     * @param newParentTime the time at which the new parent of startNode arises
     * return ArrayList with elements (in this order): [0] maximum node, [1] maximum time up to which there is no haplotype yet,
     *                      and if present the limiting (next) [2] haplotype
     */
    public ArrayList getMaxPossibleHaploAttachTime(QuasiSpeciesNode startNode, int haplo,double newParentTime){

        ArrayList output = new ArrayList(3);
        // starting from the startNode look for parent's continuingHaploName and if present and not the same as haplo
        // that node is the max possible attachment time of the haplotype
        QuasiSpeciesNode nodeToCheck = startNode;
        if(nodeToCheck.getContinuingHaploName() != -1 && nodeToCheck.getContinuingHaploName() != haplo &&
                qsTree.getAttachmentTimesList((QuasiSpeciesNode) qsTree.getNode(nodeToCheck.getContinuingHaploName()))[0] > newParentTime){
            // the max node is the new to be created node at the attach time!
            output.add(0,null);
            output.add(1,newParentTime);
            output.add(2,nodeToCheck.getContinuingHaploName());
            return output;
            //return newParentTime;
        }
        else {
            if (!nodeToCheck.isRoot())
                nodeToCheck = (QuasiSpeciesNode) nodeToCheck.getParent();
            while ((nodeToCheck.getContinuingHaploName() == -1 || nodeToCheck.getContinuingHaploName() == haplo)
                    && nodeToCheck != qsTree.getRoot()){
                nodeToCheck = (QuasiSpeciesNode) nodeToCheck.getParent();
            }
            if (nodeToCheck == qsTree.getRoot() &&
                    (nodeToCheck.getHaploAboveName() == -1 || nodeToCheck.getHaploAboveName() == haplo)){
                output.add(0,origin);
                output.add(1,origin.getValue());
                output.add(2,-1);
                return output;
                //return origin.getValue();
            }
            else{
                output.add(0,nodeToCheck);
                output.add(1,nodeToCheck.getHeight());
                output.add(2,nodeToCheck.getContinuingHaploName());
                return output;
                //return nodeToCheck.getHeight();
            }
        }
    }

    /**
     * Function to find a maximum node until which the haplo can attach to (start from)
     * TODO !!! WATCH OUT THIS FUNCTION RETURNS A NODE BUT ACTUALLY THERE MAY BE HAPLO ABOVE IT ALREADY, SO THE ACTUAL MAX
     * TODO POSSIBLE ATTACHMENT NODE CAN ALSO JUST BE THE NEW NODE!!!
     *
     * @param srcNode the node that attaches to destNode
     * @param destNode the node at which to start the search going up towards the origin
     * @param newParentTime the time at which the new parent of startNode arises
     * return the maximum node up to which there is no haplotype yet
     */
    public QuasiSpeciesNode getMaxPossibleHaploAttachNode(QuasiSpeciesNode srcNode, QuasiSpeciesNode destNode, double newParentTime){
        // starting from the startNode look for parent's continuingHaploName and if present
        // that node is the max possible attachment time of the haplotype
        QuasiSpeciesNode nodeToCheck = destNode;
        if(nodeToCheck.getHaploAboveName() != -1 &&
                qsTree.getAttachmentTimesList((QuasiSpeciesNode) qsTree.getNode(nodeToCheck.getHaploAboveName()))[0] > newParentTime){
            return (QuasiSpeciesNode) srcNode.getParent();
        }
        else {
            nodeToCheck = (QuasiSpeciesNode) nodeToCheck.getParent();
            while (nodeToCheck.getContinuingHaploName() == -1 && nodeToCheck != qsTree.getRoot()){
                nodeToCheck = (QuasiSpeciesNode) nodeToCheck.getParent();
            }
            return nodeToCheck;
        }
    }

    /**
     * Function to find number of haplotypes in a subtree starting at startNode that can be pulled up
     *  without the need to push any haplotypes below the common ancestor
     *
     * @param startNode the node at which to start the search going down towards the tips
     * @param parentHaplo the array of parent haplotypes for each haplotype
     * @param possibleHaplo the number of haplotypes that can be pulled up - will be re-written here
     */
    public void checkNumberOfPossibleSrcHaplo(QuasiSpeciesNode startNode, int[] parentHaplo,
                                              ArrayList<Integer> possibleHaplo){
        List<Node> children = startNode.getAllLeafNodes();
        for (Node childLeafNode : children) {
            int childLeafNodeNr = childLeafNode.getNr();
            if (parentHaplo[childLeafNodeNr]==-1){
                possibleHaplo.add(childLeafNodeNr);
            }
            else {
                boolean addnode = true;
                for (Node otherChildLeafNode : children){
                    if(otherChildLeafNode.getNr()==parentHaplo[childLeafNodeNr]){
                        addnode = false;
                        break;
                    }
                }
                if (addnode)
                    possibleHaplo.add(childLeafNodeNr);
            }
        }
    }

}
