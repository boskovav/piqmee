package quasispeciestree.tree;

import beast.core.Description;
import beast.evolution.tree.Node;

import java.util.Arrays;

/**
 * @author Veronika Boskova created on 01/07/2015.
 */

@Description("A node in a quasi-species phylogenetic tree.")
public class QuasiSpeciesNode extends Node {

    @Override
    public void initAndValidate(){

        super.initAndValidate();

        if (isLeaf()) setContinuingHaploName(this.getNr());

    }

    private int haploAboveName = -1;
    private int continuingHaploName = -1;
    private int startBranchCounts = 1;

    // attachment times are sorted from the biggest to the smallest
    private double[] attachmentTimesList;
    // tip times are sorted from the most recent to the most into the past
    // NOTICE: the real tip in the tree has the value of tipTimesList[0]
    private double[] tipTimesList;
    // same order as tipTimesList count in tipTimesCountList[0] includes the real tip
    private int[] tipTimesCountList;
    private int parentHaplo = -1;

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
        this.makeDirty(QuasiSpeciesTree.IS_FILTHY);
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
        node.startBranchCounts = startBranchCounts;

        if (attachmentTimesList != null) {
            node.attachmentTimesList = new double[attachmentTimesList.length];
            System.arraycopy(attachmentTimesList,0,node.attachmentTimesList,0,attachmentTimesList.length);
        }
        if (tipTimesList != null) {
            node.tipTimesList = new double[tipTimesList.length];
            System.arraycopy(tipTimesList,0,node.tipTimesList,0,tipTimesList.length);
            node.tipTimesCountList = new int[tipTimesCountList.length];
            System.arraycopy(tipTimesCountList,0,node.tipTimesCountList,0,tipTimesCountList.length);
        }

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
        startBranchCounts = qsNode.startBranchCounts;
        if (qsNode.attachmentTimesList != null) {
            attachmentTimesList = new double[qsNode.attachmentTimesList.length];
            System.arraycopy(qsNode.attachmentTimesList,0,attachmentTimesList,0,qsNode.attachmentTimesList.length);
        }
        if (qsNode.tipTimesList != null) {
            tipTimesList = new double[qsNode.tipTimesList.length];
            System.arraycopy(qsNode.tipTimesList,0,tipTimesList,0,qsNode.tipTimesList.length);
            tipTimesCountList = new int[qsNode.tipTimesCountList.length];
            System.arraycopy(qsNode.tipTimesCountList,0,tipTimesCountList,0,qsNode.tipTimesCountList.length);
        }
        parentHaplo = qsNode.parentHaplo;

        if (node.getLeft()!=null) {
            setLeft(nodes[node.getLeft().getNr()]);
            ((QuasiSpeciesNode) getLeft()).assignFrom(nodes, node.getLeft());
            getLeft().setParent(this);
            if (node.getRight()!=null) {
                setRight(nodes[node.getRight().getNr()]);
                ((QuasiSpeciesNode) getRight()).assignFrom(nodes, node.getRight());
                getRight().setParent(this);
            }
        }
    }


    /**
     * scale height of this node and all its descendants
     *
     * @param scale scale factor
     */
    public void scale(final double scale) {
        startEditing();
        this.makeDirty(QuasiSpeciesTree.IS_DIRTY);
        // Scale internal node heights
        if (!isLeaf() && !isFake()) {
            height *= scale;
        }
        if (!isLeaf()) {
            ((QuasiSpeciesNode) getLeft()).scale(scale);
            if (getRight() != null) {
                ((QuasiSpeciesNode) getRight()).scale(scale);
            }
            if (height < getLeft().getHeight() || height < getRight().getHeight()) {
                throw new IllegalArgumentException("Scale gives negative branch length");
            }
        }
        // Scale haplotype attachment times
        else {
            double[] tempqstimes = this.getAttachmentTimesList().clone();
            for (int i = 0; i < tempqstimes.length; i++) {
                tempqstimes[i] = tempqstimes[i] * scale;
            }
            if (tempqstimes.length > 1)
                tempqstimes[0] = tempqstimes[1];
            else
                tempqstimes[0] = this.getHeight();
            this.setAttachmentTimesList(tempqstimes);
            // Reject invalid tree scalings:
            if (scale < 1.0) {
                    // get also tip times to help define valid scalings
                    double[] temptiptimes = this.getTipTimesList();
                    int[] temptiptimescount = this.getTipTimesCountList();
                    if (tempqstimes[tempqstimes.length-1] < this.getHeight())
                        throw new IllegalArgumentException("Scale gives negative branch length");
                    // check also if branching times corresponding to samples through time are also above the respective time
                    int currentPosition = tempqstimes.length-1-(temptiptimescount[0]-1);
                    for (int i = 1; i < temptiptimes.length; i++){
                        if (tempqstimes[currentPosition] < temptiptimes[i])
                            throw new IllegalArgumentException("Scale gives negative branch length");
                        currentPosition -= temptiptimescount[i];
                    }
                    if (tempqstimes.length > 1 && tempqstimes[0] < this.getHeight()){
                        throw new IllegalStateException("problem in QuasiSpeciesTreeScale: you did not really scale the QS start apparently");
                    }
            }
        }
    }

    /////////////////////////////////////////////////
    //  Methods implementing the quasispecies tip  //
    /////////////////////////////////////////////////


    /**
     * Method to count for a given sequence ihaploseq the possible number of attachment
     * branches of a given (same or different) haplotype
     *
     * @param ihaploseq For which sequence (i) exactly from this haplotype are we estimating the possible
     *                  attachment branches?
     *                  This is important if we are trying to list all possible branches
     *                  of haplotype from which the sequence can branch off, ie. deduct 1 if the
     *                  sequence age is identical to the ihaploseqage, and position of that haplotype
     *                  sequence in AttachmentTimesList is at ihaploseq.
     *                  If the haplotype i is different from haplotype at node, then set
     *                  ihaploseq=0.
     * @param haploseqage Age of a sequence for which we want to check
     *                     how many instances of this haplotype it can arise from
     *                     -- important if we are checking all possible branches of haplotype
     *                     from which the haplotype number i can branch off, then
     *                     set haplotypeage==attachmentTimesList[i]
     * @return Number of possible attachment branches, at position determined by node.getNr()
     *
     */
    public int countPossibleAttachmentBranches(int ihaploseq, double haploseqage) {
        // there is always one possible attachment --- the original haplotype
        int abcount = 1;
        // if node's haplotype is older than haploseqage, count how many sequences of node's haplotype are there
        // if node's haplotype is not older than haploseqage --- throw error
        // we allow also testing of other times than the ones belonging to the node's haplotype
        if (haploseqage > attachmentTimesList[0]){
            throw new RuntimeException("The haplotype "+ this.getID() + " of age " + attachmentTimesList[0] +
                    " is set to start after the haplotype sequence number " + ihaploseq +
                    " of age " + haploseqage + " is set to branch off");
        } else {
            for (int i = 1; i < attachmentTimesList.length; i++){
                // We do NOT allow branching from the branch with the starting starting age!!!
                if ( haploseqage < attachmentTimesList[i] && i != ihaploseq ){
                    abcount +=1;
                }
            }
            // deduct the duplicates that stop (sampling through time) before this node attaches
            for (int i = 0; i < tipTimesList.length; i++){
                if (haploseqage <= tipTimesList[i])
                    abcount -= tipTimesCountList[i];
                    // the times are ordered from largest to smallest so stop when a smaller time than haploseqage is encountered
                else
                    break;
            }
            if(abcount==0){
                throw new RuntimeException("There is an error in counting the number of attachment points." +
                        " Please contact the developers with this error: QuasiSpeciesHaplo # 445");
            }
        }
        return abcount;
    }

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
     * Obtain the sum of branch length for the haplotype associated with this tip
     *
     * @return the sum of branch lengths defined by the attachment times list
     *          and the tip times/counts of haplotype associated with this tip
     */
    public double getTotalBranchLengths() {
        double totalTime=0.0;
        int currentTipTimePosition = 0;
        int currentTipArrayPosition = 0;
        for (int i=0; i<attachmentTimesList.length; i++){
            totalTime += attachmentTimesList[i] - tipTimesList[currentTipTimePosition];
            currentTipArrayPosition++;
            if (currentTipArrayPosition == tipTimesCountList[currentTipTimePosition]) {
                currentTipArrayPosition = 0;
                currentTipTimePosition++;
            }
        }
        return totalTime;
    }

    /**
     * Obtain the attachment times of the haplotype associated with this tip
     *
     * @return attachment times list of haplotype associated with this tip
     */
    public double[] getAttachmentTimesList() {
        return this.attachmentTimesList;
    }

    /**
     * Sets the attachment times of the haplotype associated with this tip
     *
     * @param newAttachmentTimesList New attachment times list of haplotype associated with this tip
     */
    public void setAttachmentTimesList(double[] newAttachmentTimesList) {
        startEditing();
        this.attachmentTimesList = newAttachmentTimesList;
    }

    /**
     * Sorts the attachment times of the haplotype associated with this tip in descending order
     *  and sets the first (fake) attachment time to the max attachment time
     *
     */
    public void setFirstEntryAndSortAttachTimeList() {
        Arrays.sort(attachmentTimesList);
        // copy the largest bifurcation time, to indicate the haplo start time
        attachmentTimesList[0] = attachmentTimesList[attachmentTimesList.length - 1];
        // reverse the array to start with the largest value
        sortAttachTimeList();
    }

    /**
     * Sorts the attachment times of the haplotype associated with this tip in descending order
     *
     */
    public void sortAttachTimeList() {
        Arrays.sort(attachmentTimesList);
        // reverse the array to start with the largest value
        int totalLength = attachmentTimesList.length;
        for (int j = 0; j < totalLength/2; j++){
            double tmp = attachmentTimesList[j];
            attachmentTimesList[j] = attachmentTimesList[totalLength-1-j];
            attachmentTimesList[totalLength-1-j] = tmp;
        }
    }

    /**
     * Obtain the tip times of the haplotype associated with this tip
     *
     * @return tip times list of haplotype associated with this tip
     */
    public double[] getTipTimesList() {
        return this.tipTimesList;
    }

    /**
     * Sets the tip times of the haplotype associated with this tip
     * Should only be used at initialization!!!
     *
     * @param newTipTimesList New tip times list of haplotype associated with this tip
     */
    public void setTipTimesList(double[] newTipTimesList) {
        startEditing();
        this.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        this.tipTimesList = newTipTimesList;
    }

    /**
     * Obtain the tip counts for each time of the haplotype associated with this tip
     *
     * @return tip times counts list of haplotype associated with this tip
     */
    public int[] getTipTimesCountList() {
        return this.tipTimesCountList;
    }

    /**
     * Sets the tip times of the haplotype associated with this tip
     * Should only be used at initialization!!!
     *
     * @param newTipTimesCountList New tip times count list of haplotype associated with this tip
     */
    public void setTipTimesCountList(int[] newTipTimesCountList) {
        startEditing();
        this.makeDirty(QuasiSpeciesTree.IS_FILTHY);
        this.tipTimesCountList = newTipTimesCountList;
    }

    /**
     * Sorts the tip times/counts of the haplotype associated with this tip
     *
     */
    public void sortTipTimeAndCountList() {
        //Manually sort tipTimes, since we in the same way need to sort the tip counts
        double tmp;
        int tmpint;
        for (int j = 1; j < tipTimesList.length; j++) {
            int k = j - 1;
            while (k >= 0 && tipTimesList[k] > tipTimesList[j]) {
                tmp = tipTimesList[j];
                tipTimesList[j] = tipTimesList[k];
                tipTimesList[k] = tmp;
                tmpint = tipTimesCountList[j];
                tipTimesCountList[j] = tipTimesCountList[k];
                tipTimesCountList[k] = tmpint;
                k--;
            }
        }
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
}
