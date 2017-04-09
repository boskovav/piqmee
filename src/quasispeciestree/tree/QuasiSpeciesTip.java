package quasispeciestree.tree;

import beast.core.Description;
import beast.evolution.tree.Node;

/**
 * @author Veronika Boskova created on 06/04/2017.
 */

@Description("A tip in a quasi-species phylogenetic tree.")
public class QuasiSpeciesTip extends QuasiSpeciesNode {

    @Override
    public void initAndValidate(){

        super.initAndValidate();
        setContinuingHaploName(this.getNr());

    }

    // attachment times are sorted from the biggest to the smallest
    private double[] attachmentTimesList;
    // tip times are sorted from the most recent to the most into the past
    // NOTICE: the real tip in the tree has the value of tipTimesList[0]
    private double[] tipTimesList;
    // same order as tipTimesList count in tipTimesCountList[0] includes the real tip
    private int[] tipTimesCountList;
    private int parentHaplo;


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
    protected int countPossibleAttachmentBranches(int ihaploseq, double haploseqage){
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
    public int getHaplotypeCountsFromTips(){
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
    protected double getTotalBranchLengths(){
        double totalTime=0.0;
        int currentTipTimePosition = 0;
        int currentTipArrayPosition = 0;
        for (int i=0; i<attachmentTimesList.length; i++){
            totalTime += attachmentTimesList[i]-tipTimesList[currentTipTimePosition];
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
    protected void setAttachmentTimesList(double[] newAttachmentTimesList) {
        startEditing();
        this.attachmentTimesList = newAttachmentTimesList;
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
    protected void setTipTimesList(double[] newTipTimesList) {
        startEditing();
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
    protected void setTipTimesCountList(int[] newTipTimesCountList) {
        startEditing();
        this.tipTimesCountList = newTipTimesCountList;
    }

    /**
     * Obtain the number of "branches" this node can attach to in the QS tree
     *
     * @return parent haplo of haplotype associated with this tip
     */
    protected int getParentHaplo() {
        return this.parentHaplo;
    }

    /**
     * Sets the number of branches the node can be attached to
     *
     * @param newParentHaplo New parent haplotype of the haplotype associated with this tip
     */
    protected void setParentHaplo(int newParentHaplo) {
        startEditing();
        this.parentHaplo = newParentHaplo;
    }


    /**
     * **************************
     * Methods ported from Node *
     ***************************
     */

    /**
     * @return (deep) copy of node
     */
    public QuasiSpeciesTip copyTip() {
        QuasiSpeciesTip node = new QuasiSpeciesTip();
        node.height = height;
        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.setParent(null);
        node.ID = ID;

        node.haploAboveName = haploAboveName;
        node.continuingHaploName = continuingHaploName;

        node.attachmentTimesList = attachmentTimesList;
        node.tipTimesList = tipTimesList;
        node.tipTimesCountList = tipTimesCountList;
        node.parentHaplo = parentHaplo;

        return node;
    }

    /**
     * assign values from a tree in array representation *
     * @param node
     */
    public void assignFromTip(Node node) {
        height = node.getHeight();
        labelNr = node.getNr();
        metaDataString = node.metaDataString;
        setParent(null);
        ID = node.getID();

        QuasiSpeciesTip qsNode = (QuasiSpeciesTip)node;
        haploAboveName = qsNode.haploAboveName;
        continuingHaploName = qsNode.continuingHaploName;


    }

}
