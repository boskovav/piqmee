package quasispeciestree.tree;

import beast.core.Description;
import beast.evolution.tree.Node;

/**
 * @author Veronika Boskova created on 01/07/2015.
 */

@Description("A node in a quasi-species phylogenetic tree.")
public class QuasiSpeciesNode extends Node {

    @Override
    public void initAndValidate(){

        super.initAndValidate();

    }

    protected int haploAboveName = -1;
    protected int continuingHaploName = -1;
    private int startBranchCounts = 1;

    /**
     * Obtain the quasi-species type/name, if any, starting on the branch above this node.
     *
     * @return quasi-species name
     */
    public int getHaploAboveName() {
        return this.haploAboveName;
    }

    /**
     * Sets the quasi-species starting on the branch above this node
     *
     * @param haploName New quasi-species name
     */
    public void setHaploAboveName(int haploName) {
        startEditing();
        this.haploAboveName = haploName;
    }

    /**
     * Obtain the quasi-species type/name, of a haplotype that started earlier and continues on below this (internal) node.
     * Only non null (non -1) if this haplo arose at some previous node (so haploAboveName is non-null (non -1))
     *
     * @return quasi-species name
     */
    public int getContinuingHaploName() {
        return this.continuingHaploName;
    }

    /**
     * Sets the quasi-species haplotype as a continuing haplotype below this node
     *
     * @param haploName New quasi-species name
     */
    public void setContinuingHaploName(int haploName) {
        startEditing();
        this.continuingHaploName = haploName;
    }

    /**
     * Obtain the number of "branches" this node can attach to in the QS tree
     *
     * @return start branch count
     */
    public int getStartBranchCounts() {
        return this.startBranchCounts;
    }

    /**
     * Sets the number of branches the node can be attached to
     *
     * @param startBranchCount New quasi-species name
     */
    public void setStartBranchCounts(int startBranchCount) {
        startEditing();
        this.startBranchCounts = startBranchCount;
    }

    /**
     * Set quasi-species tree for a copied node
     */
    public void setqsTree(QuasiSpeciesTree tree) {
        startEditing();
        this.m_tree = tree;
    }

    /**
     * **************************
     * Methods ported from Node *
     ***************************
     */


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
            if (getLeft().getAllChildNodes().size() == 0)
                node.setLeft(((QuasiSpeciesTip)getLeft()).copyTip());
            else
                node.setLeft(getLeft().copy());
            node.getLeft().setParent(node);
            if (getRight()!=null) {
                if (getRight().getAllChildNodes().size() == 0)
                    node.setRight(((QuasiSpeciesTip)getRight()).copyTip());
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
            if (getLeft().getAllChildNodes().size() == 0)
                ((QuasiSpeciesTip) getLeft()).assignFromTip(node.getLeft());
            getLeft().setParent(this);
            if (node.getRight()!=null) {
                setRight(nodes[node.getRight().getNr()]);
                if (getRight().getAllChildNodes().size() == 0)
                    ((QuasiSpeciesTip) getRight()).assignFromTip(node.getRight());
                getRight().setParent(this);
            }
        }
    }

}
