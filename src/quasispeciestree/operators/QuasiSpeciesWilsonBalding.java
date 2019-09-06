package quasispeciestree.operators;


import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.util.Randomizer;

/**
 * @author Veronika Boskova created on 15/07/2015.
 */
@Description("Implements the unweighted Wilson-Balding branch"
        +" swapping move.  This move is similar to one proposed by WILSON"
        +" and BALDING 1998 and involves removing a subtree and"
        +" re-attaching it on a new parent branch. "
        +" See <a href='http://www.genetics.org/cgi/content/full/161/3/1307/F1'>picture</a>."
        +" This version retypes each newly generated branch by drawing a"
        +" path from the migration model conditional on the types at the"
        +" branch ends.")
public class QuasiSpeciesWilsonBalding {
}


// TODO: for the normal WIlson-balding need to implement pushing down of all QS in tips of chosen subtree
// and then moving and reseeding QS on that sub-tree.
// TODO make a function in QS tree to re-seed start points of QS and their copy sequences --- can be used by above WB
// TODO: Do I need to re-calculate the array holding the number of possible attachment branches of each internal node?
// --- I think YES


/*
* TODO: need to implement check for whether the proposed branching times for each duplicate of a haplotype are
       in accordance with "appear-only-once-in-tree-history" assumption --- where to do it?  -- in operator? to be efficient
       -- randomly choose a haplotype, propose its attachment times from root to tip, if another haplotype attaches to the path (or below),
       where the QS of the first haplotype chose to attach, then allow these to only propose within attachmet node to first haplotype
       to the tip of the second haplotype, or to the tip of the haplotypes below.

       Also check in tree operators whether the QS starting/attachment times are still compatible, if not propose new QS times...
 *
 */


// import beast.util.Randomizer;
//    /**
//     * Use beast RNG to select random node from list.
//     *
//     * @param nodeList
//     * @return A randomly selected node.
//     */
//    private MultiTypeNode selectRandomNode(List<MultiTypeNode> nodeList) {
//        return nodeList.get(Randomizer.nextInt(nodeList.size()));
//    }







////public class TypedWilsonBalding extends UniformizationRetypeOperator {
//
//    public Input<Double> alphaInput = new Input<>("alpha",
//            "Root height proposal parameter", Input.Validate.REQUIRED);
//    private double alpha;
//
//    @Override
//    public void initAndValidate() throws Exception {
//        super.initAndValidate();
//
//        alpha = alphaInput.get();
//    }
//
//    @Override
//    public double proposal() {
//        // Check that operator can be applied to tree:
//        if (mtTree.getLeafNodeCount()<3)
//            throw new IllegalStateException("Tree too small for"
//                    +" TypedWilsonBalding operator.");
//
//        // Select source node:
//        Node srcNode;
//        do {
//            srcNode = mtTree.getNode(Randomizer.nextInt(mtTree.getNodeCount()));
//        } while (invalidSrcNode(srcNode));
//        Node srcNodeP = srcNode.getParent();
//        Node srcNodeS = getOtherChild(srcNodeP, srcNode);
//        double t_srcNode = srcNode.getHeight();
//        double t_srcNodeP = srcNodeP.getHeight();
//        double t_srcNodeS = srcNodeS.getHeight();
//
//        // Select destination branch node:
//        Node destNode;
//        do {
//            destNode = mtTree.getNode(Randomizer.nextInt(mtTree.getNodeCount()));
//        } while (invalidDestNode(srcNode, destNode));
//        Node destNodeP = destNode.getParent();
//        double t_destNode = destNode.getHeight();
//
//        // Handle special cases involving root:
//
//        if (destNode.isRoot()) {
//            // FORWARD ROOT MOVE
//
//            double logHR = 0.0;
//
//            // Record probability of current colouring:
//            logHR += getBranchTypeProb(srcNode);
//
//            // Record srcNode grandmother height:
//            double t_srcNodeG = srcNodeP.getParent().getHeight();
//
//            // Choose new root height:
//            double newTime = t_destNode+Randomizer.nextExponential(1.0/(alpha*t_destNode));
//
//            // Implement tree changes:
//            disconnectBranch(srcNode);
//            connectBranchToRoot(srcNode, destNode, newTime);
//            mtTree.setRoot(srcNodeP);
//
//            // Recolour root branches:
//            try {
//                logHR -= retypeRootBranches(srcNode);
//            } catch (NoValidPathException e) {
//                return Double.NEGATIVE_INFINITY;
//            }
//
//            // Return HR:
//            logHR += Math.log(alpha*t_destNode)
//                    +(1.0/alpha)*(newTime/t_destNode-1.0)
//                    -Math.log(t_srcNodeG-Math.max(t_srcNode, t_srcNodeS));
//
//            return logHR;
//        }
//
//        if (srcNodeP.isRoot()) {
//            // BACKWARD ROOT MOVE
//
//            double logHR = 0.0;
//
//            // Incorporate probability of current colouring:
//            logHR += getRootBranchTypeProb(srcNode);
//
//            // Record old srcNode parent height
//            double oldTime = t_srcNodeP;
//
//            // Choose height of new attachement point:
//            double min_newTime = Math.max(t_srcNode, t_destNode);
//            double t_destNodeP = destNodeP.getHeight();
//            double span = t_destNodeP-min_newTime;
//            double newTime = min_newTime+span*Randomizer.nextDouble();
//
//            // Implement tree changes:
//            disconnectBranchFromRoot(srcNode);
//            connectBranch(srcNode, destNode, newTime);
//            srcNodeS.setParent(null);
//            mtTree.setRoot(srcNodeS);
//
//            // Recolour new branch:
//            try {
//                logHR -= retypeBranch(srcNode);
//            } catch (NoValidPathException e) {
//                return Double.NEGATIVE_INFINITY;
//            }
//
//            // Return HR:
//            logHR += Math.log(t_destNodeP-Math.max(t_srcNode, t_destNode))
//                    -Math.log(alpha*t_srcNodeS)
//                    -(1.0/alpha)*(oldTime/t_srcNodeS-1.0);
//
//            return logHR;
//        }
//
//        // NON-ROOT MOVE
//
//        double logHR = 0.0;
//
//        // Incorporate probability of current colouring.
//        logHR += getBranchTypeProb(srcNode);
//
//        // Record srcNode grandmother height:
//        double t_srcNodeG = srcNodeP.getParent().getHeight();
//
//        // Choose height of new attachment point:
//        double min_newTime = Math.max(t_destNode, t_srcNode);
//        double t_destNodeP = destNodeP.getHeight();
//        double span = t_destNodeP-min_newTime;
//        double newTime = min_newTime+span*Randomizer.nextDouble();
//
//        // Implement tree changes:
//        disconnectBranch(srcNode);
//        connectBranch(srcNode, destNode, newTime);
//
//        // Recolour new branch:
//        try {
//            logHR -= retypeBranch(srcNode);
//        } catch (NoValidPathException e) {
//            return Double.NEGATIVE_INFINITY;
//        }
//
//        // HR contribution of topology and node height changes:
//        logHR += Math.log(t_destNodeP-Math.max(t_srcNode, t_destNode))
//                -Math.log(t_srcNodeG-Math.max(t_srcNode, t_srcNodeS));
//
//        return logHR;
//    }
//
//    /**
//     * Returns true if srcNode CANNOT be used for the CWBR move.
//     *
//     * @param srcNode
//     * @return True if srcNode invalid.
//     */
//    private boolean invalidSrcNode(Node srcNode) {
//
//        if (srcNode.isRoot())
//            return true;
//
//        Node parent = srcNode.getParent();
//
//        // This check is important for avoiding situations where it is
//        // impossible to choose a valid destNode:
//        if (parent.isRoot()) {
//
//            Node sister = getOtherChild(parent, srcNode);
//
//            if (sister.isLeaf())
//                return true;
//
//            if (srcNode.getHeight()>=sister.getHeight())
//                return true;
//        }
//
//        return false;
//    }
//
//    /**
//     * Returns true if destNode CANNOT be used for the CWBR move in conjunction
//     * with srcNode.
//     *
//     * @param srcNode
//     * @param destNode
//     * @return True if destNode invalid.
//     */
//    private boolean invalidDestNode(Node srcNode, Node destNode) {
//
//        if (destNode==srcNode
//                ||destNode==srcNode.getParent()
//                ||destNode.getParent()==srcNode.getParent())
//            return true;
//
//        Node destNodeP = destNode.getParent();
//
//        if (destNodeP!=null&&(destNodeP.getHeight()<=srcNode.getHeight()))
//            return true;
//
//        return false;
//    }
//
//    /**
//     * Retype branches with nChanges between srcNode and the root (srcNode's
//     * parent) and nChangesSister between the root and srcNode's sister.
//     *
//     * @param srcNode
//     * @return Probability of new state.
//     */
//    private double retypeRootBranches(Node srcNode) throws NoValidPathException {
//
//        double logProb = 0.0;
//
//        Node srcNodeP = srcNode.getParent();
//        Node srcNodeS = getOtherChild(srcNodeP, srcNode);
//
//        // Select new root colour:
//        ((MultiTypeNode)srcNodeP).setNodeType(Randomizer.nextInt(migModel.getNTypes()));
//
//        // Incorporate probability of new root colour:
//        logProb += Math.log(1.0/migModel.getNTypes());
//
//        // Recolour branches conditional on root type:
//        logProb += retypeBranch(srcNode);
//        logProb += retypeBranch(srcNodeS);
//
//
//        // Return probability of new colouring given boundary conditions:
//        return logProb;
//    }
//
//
//    /**
//     * Obtain joint probability of typing along branches between srcNode and
//     * the root, the sister of srcNode and the root, and the node type of the
//     * root.
//     *
//     * @param srcNode
//     * @return
//     */
//    protected double getRootBranchTypeProb(Node srcNode) {
//
//        double logProb = 0.0;
//
//        Node srcNodeP = srcNode.getParent();
//        Node srcNodeS = getOtherChild(srcNodeP, srcNode);
//
//        // Probability of node type:
//        logProb += Math.log(1.0/migModel.getNTypes());
//
//        // Probability of branch types conditional on node types:
//        logProb += getBranchTypeProb(srcNode);
//        logProb += getBranchTypeProb(srcNodeS);
//
//        return logProb;
//    }
//
//}