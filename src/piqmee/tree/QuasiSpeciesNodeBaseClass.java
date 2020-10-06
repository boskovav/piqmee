package piqmee.tree;

import beast.core.Description;
import beast.evolution.tree.Node;

/**
 * @author Veronika Boskova created on 30/09/2020.
 */

@Description("A base class for a node in a quasi-species-type of phylogenetic tree.")
public abstract class QuasiSpeciesNodeBaseClass extends Node {

    protected int haploAboveName = -1;
    protected int continuingHaploName = -1;
    // tip times are sorted from the most recent to the most into the past
    // NOTICE: the real tip in the tree has the value of tipTimesList[0]
    protected double[] tipTimesList;
    // same order as tipTimesList count in tipTimesCountList[0] includes the real tip
    protected int[] tipTimesCountList;
    protected int parentHaplo = -1;

    @Override
    public void initAndValidate(){
        super.initAndValidate();
    }

    /**
     * ********************************************
     * Methods implementing the quasispecies node *
     **********************************************
     */

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
        this.makeDirty(QuasiSpeciesTree.IS_FILTHY);
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
        this.makeDirty(QuasiSpeciesTree.IS_FILTHY);
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
     * *******************************************
     * Methods implementing the quasispecies tip *
     *********************************************
     */

    /**
     * Obtain the total number of haplotype duplicates from the tip counts
     *
     * @return the total number of haplotype duplicates
     */
    public int getHaplotypeCountsFromTips() {
        int totalHaploCount = 0;
        for (int i = 0; i < tipTimesCountList.length; i++){
            totalHaploCount += tipTimesCountList[i];
        }
        return totalHaploCount;
    }

    /**
     * Obtain the number of "branches" this node can attach to in the QS tree
     *
     * @return parent haplo of haplotype associated with this tip
     */
    public int getParentHaplo() {
        return this.parentHaplo;
    }

    /**
     * Sets the number of branches the node can be attached to
     *
     * @param newParentHaplo New parent haplotype of the haplotype associated with this tip
     */
    public void setParentHaplo(int newParentHaplo) {
        startEditing();
        this.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        this.parentHaplo = newParentHaplo;
    }

    /**
     * **************************
     * Methods ported from Node *
     ***************************
     */

    // copy
    // assignFrom
    // scale

}
