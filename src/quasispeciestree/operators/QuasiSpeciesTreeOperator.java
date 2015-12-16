package quasispeciestree.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesTree;
import quasispeciestree.tree.QuasiSpeciesNode;
import beast.core.parameter.RealParameter;



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

    /* **********************************************************************/





    /**
     * Disconnect edge <node,node.getParent()> by joining node's sister directly
     * to node's grandmother
     * This move can only be done if the haplotype in the sub-tree to be moved,
     * i.e. sub-tree below the node, does not arise above the node.getParent().
     *
     * @param node
     */
    public void disconnectBranch(Node node) {

        // Check argument validity:
        Node parent = node.getParent();
        if (node.isRoot() || parent.isRoot())
            throw new IllegalArgumentException("Illegal argument to "
                    + "disconnectBranch().");

        Node sister = getOtherChild(parent, node);

        // Check that there are no haplotype sequences that belong to the sub-tree being
        // moved, arising at the branch leading to parent:
        String parentQSType = ((QuasiSpeciesNode) node.getParent()).getHaploName();
        if (parentQSType != null){
            for (Node childLeafNode : node.getAllLeafNodes()) {
                if (childLeafNode.getID() == parentQSType){
                    // if haplotype that is in the moved sub-tree arises on the branch
                    // to parent haplotype, throw error cannot be the case
                    System.out.println("Trying to disconnect a branch while quasi-species are evolving on it...");
                    System.exit(0);
                }
            }
        } else {
        // Add haplotype start/attachment times originally attached to parent to those attached
        // to node's sister:
            if (((QuasiSpeciesNode) sister).getHaploName() != null){
                // if there is one haplotype arising on the parent's branch and another on sister's branch
                // there is a problem somewhere in the code!!!
                System.out.println("There must be a serious problem with the code. There are haplotypes arising at " +
                        "both parent and sister branch, while this should technically be impossible...");
                System.exit(0);
            } else {
                ((QuasiSpeciesNode) sister).setHaploName(parentQSType);
            }
        }
        // Implement topology change.
        replace(parent.getParent(), parent, sister);

        // Ensure BEAST knows to update affected likelihoods:
        parent.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        sister.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);
    }

    /**
     * Disconnect node from root
     * * This move can only be done if the haplotype in the sub-tree to be moved,
     * i.e. sub-tree below the node, does not arise at the <node,root> branch.
     *
     * @param node
     */
    public void disconnectBranchFromRoot(Node node) {

        // Check argument validity:
        Node parent = node.getParent();
        if (node.isRoot() || !parent.isRoot())
            throw new IllegalArgumentException("Illegal argument to"
                    + " disconnectBranchFromRoot().");

        Node sister = getOtherChild(parent, node);

        // Check that there are no haplotype sequenes that belong to the sub-tree being
        // moved, arising at the branch leading to parent:
        String parentQSType = ((QuasiSpeciesNode) node.getParent()).getHaploName();
        if (parentQSType != null){
            for (Node childLeafNode : node.getAllLeafNodes()) {
                if (childLeafNode.getID() == parentQSType){
                    // if haplotype that is in the moved sub-tree arises on the branch
                    // to parent haplotype, throw error cannot be the case
                    System.out.println("Trying to disconnect a branch while quasi-species are evolving on it...");
                    System.exit(0);
                }
            }
        } else {
            // Add haplotype start/attachment times originally attached to parent to those attached
            // to node's sister:
            if (((QuasiSpeciesNode) sister).getHaploName() != null){
                // if there is one haplotype arising on the parent's branch and another on sister's branch
                // there is a problem somewhere in the code!!!
                System.out.println("There must be a serious problem with the code. There are haplotypes arising at " +
                        "both parent and sister branch, while this should technically be impossible...");
                System.exit(0);
            } else {
                ((QuasiSpeciesNode) sister).setHaploName(parentQSType);
            }
        }
        // Implement topology change:
        sister.setParent(null);
        parent.getChildren().remove(sister);

        // Ensure BEAST knows to update affected likelihoods:
        parent.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        sister.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);
    }

    /**
     * Creates a new branch between node and a new node at time destTime between
     * destBranchBase and its parent.
     *
     * @param node
     * @param destBranchBase
     * @param destTime
     */
    public void connectBranch(Node node, Node destBranchBase, double destTime) {

        // Check argument validity:
        if (node.isRoot() || destBranchBase.isRoot())
            throw new IllegalArgumentException("Illegal argument to "
                    + "connectBranch().");

        // Obtain existing parent of node and set new time:
        Node parent = node.getParent();
        parent.setHeight(destTime);

        // destTime must be higher than, if existant, the first haplotype arising below the node.getParent()
        for (Node childLeafNode : node.getAllLeafNodes()){
            if (destTime < qsTree.getAttachmentTimesList((QuasiSpeciesNode)childLeafNode)[0]){
                System.out.println("Cannot attach subtree to new destination. The haplotype arises above the " +
                        "destination attachment point...");
                System.exit(0);
            }
        }
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
    }

    /**
     * Set up node's parent as the new root with a height of destTime, with
     * oldRoot as node's new sister.
     *
     * @param node
     * @param oldRoot
     * @param destTime
     */
    public void connectBranchToRoot(Node node, Node oldRoot, double destTime) {

        // Check argument validity:
        if (node.isRoot() || !oldRoot.isRoot())
            throw new IllegalArgumentException("Illegal argument "
                    + "to connectBranchToRoot().");

        // destTime must be higher than, if existant, the first haplotype arising below the node.getParent()
        for (Node childLeafNode : node.getAllLeafNodes()){
            if (destTime < qsTree.getAttachmentTimesList((QuasiSpeciesNode)childLeafNode)[0]){
                System.out.println("Cannot set root of new tree to the destination time. The haplotype arises above the " +
                        "destination attachment point...");
                System.exit(0);
            }
        }
        // Obtain existing parent of node and set new time:
        Node newRoot = node.getParent();
        newRoot.setHeight(destTime);

        // Implement topology changes:
        newRoot.setParent(null);

        if (newRoot.getLeft() == node)
            newRoot.setRight(oldRoot);
        else if (newRoot.getRight() == node)
            newRoot.setLeft(oldRoot);

        oldRoot.setParent(newRoot);

        // Ensure BEAST knows to recalculate affected likelihood:
        newRoot.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        oldRoot.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        node.makeDirty(QuasiSpeciesTree.IS_FILTHY);
    }
}
