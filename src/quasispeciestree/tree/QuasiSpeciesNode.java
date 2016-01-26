package quasispeciestree.tree;

import beast.core.Description;
import beast.evolution.tree.Node;

/**
 * @author Veronika Boskova created on 01/07/2015.
 */

@Description("A node in a quasi-species phylogenetic tree.")
public class QuasiSpeciesNode extends Node {


    @Override
    public void initAndValidate() throws Exception {

        super.initAndValidate();

    }

    protected String haploName;

    /**
     * Obtain the quasi-species type/name, if any, starting on the branch above this node.
     *
     * @return quasi-species name
     */
    public String getHaploName() {
        return this.haploName;
    }

    /**
     * Sets the quasi-species starting on the branch above this node
     *
     * @param haploName New quasi-species name
     */
    public void setHaploName(String haploName) {
        // start editing already set in the QuasiSpeciesHaplotypeStartSwap -- if done by the same operator several times, results in wrong stored tree object
//        startEditing();
        this.haploName = haploName;
    }

    /**
     * Set quasi-specied tree for a copied node
     */

    public void setmTree(QuasiSpeciesTree tree) {
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

        node.haploName = haploName;

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

        node.haploName = haploName;

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
        haploName = qsNode.haploName;

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

}
