package quasispeciestree.tree;

import beast.core.Description;
import beast.evolution.tree.Node;

/**
 * @author Veronika Boskova created on 01/07/2015.
 */

@Description("A node in a quasi-species phylogenetic tree.")
public class QuasiSpeciesNode extends Node {

    /**
     * status of this node after an operation is performed on the state *
     */
    int isDirty = QuasiSpeciesTree.IS_CLEAN;

    @Override
    public void initAndValidate(){

        super.initAndValidate();

    }

//    protected String haploAboveName;
//    protected String continuingHaploName;
    protected int haploAboveName;
    protected int continuingHaploName;


    public void setAttachmentTimesList(){
        startEditing();
    }

    public void  setStartBranchCounts(){
        startEditing();
    }

    public void setParentHaplo() {
        startEditing();
    }


    /**
     * Obtain the quasi-species type/name, if any, starting on the branch above this node.
     *
     * @return quasi-species name
     */
//    public String getHaploAboveName() {return this.haploAboveName; }
    public int getHaploAboveName() {
        return this.haploAboveName;
    }

    /**
     * Sets the quasi-species starting on the branch above this node
     *
     * @param haploName New quasi-species name
     */
//    public void setHaploAboveName(String haploName) { this.haploAboveName = haploName; }
    public void setHaploAboveName(int haploName) {
        startEditing();
        this.isDirty |= QuasiSpeciesTree.IS_DIRTY;
        if (!isLeaf()) {
            ((QuasiSpeciesNode)getLeft()).isDirty |= QuasiSpeciesTree.IS_DIRTY;
            if (getRight() != null) {
                ((QuasiSpeciesNode)getRight()).isDirty |= QuasiSpeciesTree.IS_DIRTY;
            }
        }
        this.haploAboveName = haploName;
    }

    /**
     * Obtain the quasi-species type/name, of a haplotype that started earlier and continues on below this (internal) node.
     * Only non null (non -1) if this haplo arose at some previous node (so haploAboveName is non-null (non -1))
     *
     * @return quasi-species name
     */
//    public String getContinuingHaploName() { return this.continuingHaploName; }
    public int getContinuingHaploName() { return this.continuingHaploName; }

    /**
     * Sets the quasi-species haplotype as a continuing haplotype below this node
     *
     * @param haploName New quasi-species name
     */
//    public void setContinuingHaploName(String haploName) { this.continuingHaploName = haploName; }
    public void setContinuingHaploName(int haploName) {
        startEditing();
        this.isDirty |= QuasiSpeciesTree.IS_DIRTY;
        if (!isLeaf()) {
            ((QuasiSpeciesNode)getLeft()).isDirty |= QuasiSpeciesTree.IS_DIRTY;
            if (getRight() != null) {
                ((QuasiSpeciesNode)getRight()).isDirty |= QuasiSpeciesTree.IS_DIRTY;
            }
        }
        this.continuingHaploName = haploName;
    }


    /**
     * Set quasi-species tree for a copied node
     */

    public void setqsTree(QuasiSpeciesTree tree) {
        startEditing();
        this.m_tree = tree;
    }

    /**
     * @return shallow copy of node
     */
    public QuasiSpeciesNode shallowCopy() {
        QuasiSpeciesNode node = new QuasiSpeciesNode();
        node.height = height;
        node.setParent(this.getParent());
        if (getLeft()!=null) {
            node.setLeft(getLeft().copy());
            node.getLeft().setParent(node);
            if (getRight()!=null) {
                node.setRight(getRight().copy());
                node.getRight().setParent(node);
            }
        }

        node.haploAboveName = haploAboveName;
        node.continuingHaploName = continuingHaploName;

        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.ID = ID;

        return node;
    }

    /**
     * **************************
     * Methods ported from Node *
     ***************************
     */

    /**
     * methods for accessing the dirtiness state of the Node.
     * A Node is Tree.IS_DIRTY if its value (like height) has changed
     * A Node Tree.IS_if FILTHY if its parent or child has changed.
     * Otherwise the node is Tree.IS_CLEAN *
     */
    public int isDirty() {
        return isDirty;
    }

    /**
     * Sets the parent of this node
     *
     * @param parent     the node to become parent
     * @param inOperator if true, then startEditing() is called and setting the parent will make tree "filthy"
     */
    public void setParent(final Node parent, final boolean inOperator) {
        // start editing set in operator itself
        if (inOperator) startEditing();
        if (this.getParent() != parent) {
            this.setParent(parent);
            if (inOperator) this.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        }
    }

    /**
     * @return (deep) copy of node
     */
    @Override
    public QuasiSpeciesNode copy() {
        QuasiSpeciesNode node = new QuasiSpeciesNode();
        node.height = height;
        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.setParent(null);
        node.ID = ID;

        node.haploAboveName = haploAboveName;
        node.continuingHaploName = continuingHaploName;

        if (getLeft()!=null) {
            node.setLeft(getLeft().copy());
            node.getLeft().setParent(node);
            if (getRight()!=null) {
                node.setRight(getRight().copy());
                node.getRight().setParent(node);
            }
        }
        return node;
    }

    /**
     * assign values from a tree in array representation *
     * @param nodes
     * @param node
     */
    @Override
    public void assignFrom(Node[] nodes, Node node) {
        height = node.getHeight();
        labelNr = node.getNr();
        metaDataString = node.metaDataString;
        setParent(null);
        ID = node.getID();

        QuasiSpeciesNode qsNode = (QuasiSpeciesNode)node;
        haploAboveName = qsNode.haploAboveName;
        continuingHaploName = qsNode.continuingHaploName;

        if (node.getLeft()!=null) {
            setLeft(nodes[node.getLeft().getNr()]);
            getLeft().assignFrom(nodes, node.getLeft());
            getLeft().setParent(this);
            if (node.getRight()!=null) {
                setRight(nodes[node.getRight().getNr()]);
                getRight().assignFrom(nodes, node.getRight());
                getRight().setParent(this);
            }
        }
    }
    public void setHeight(final double height, boolean inOperator) {
        if (inOperator){
            this.startEditing();
        }
        this.height = height;
        this.makeDirty(QuasiSpeciesTree.IS_DIRTY);
        if (!isLeaf()) {
            getLeft().makeDirty(QuasiSpeciesTree.IS_DIRTY);
            if (getRight() != null) {
                getRight().makeDirty(QuasiSpeciesTree.IS_DIRTY);
            }
        }
    }
}
