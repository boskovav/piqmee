package piqmee.evolution.branchratemodel;

import beast.core.*;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.tree.Tree;
import piqmee.tree.QuasiSpeciesNode;

import java.io.PrintStream;



@Description("Logger reporting categories for each branch of the tree.")
public class CategoriesLogger extends BEASTObject implements Loggable, Function {

    final public Input<GenericTreeLikelihood> likelihoodInput = new Input<>("treeLikelihood", "TreeLikelihood containing branch rate model that provides rates for a tree");
    final public Input<BranchRateModel> branchRateModelInput = new Input<>("branchratemodel", "model that provides rates for a tree", Input.Validate.XOR, likelihoodInput);
    final public Input<Tree> treeInput = new Input<>("tree", "tree for which the rates apply");

    private Tree tree = null;

    private BranchRateModel branchRateModel = null;

    static int[] categories;

    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        branchRateModel = branchRateModelInput.get();
        if (branchRateModel == null) {
            branchRateModel = likelihoodInput.get().branchRateModelInput.get();
        }

        categories = new int [tree.getNodeCount()+tree.getLeafNodeCount()];
        for (int i=0;i<tree.getNodeCount()+tree.getLeafNodeCount();i++){
            if (i<tree.getLeafNodeCount() && ((QuasiSpeciesNode) tree.getNode(i)).getAttachmentTimesList().length==1){
                categories[i]=0;
            } else if (i<tree.getNodeCount() && tree.getNode(i).isRoot()){
                categories[i]=0;
            } else if (i>=tree.getLeafNodeCount() && i<tree.getNodeCount() && ((QuasiSpeciesNode) tree.getNode(i)).getContinuingHaploName()!=-1){
                categories[i] = ((QuasiSpeciesUCRelaxedClockModel) branchRateModel).getCategories(((QuasiSpeciesNode) tree.getNode(i)).getContinuingHaploName());
            } else {
                categories[i] = ((QuasiSpeciesUCRelaxedClockModel) branchRateModel).getCategories(i);
            }
        }
    }

    @Override
    public int getDimension() {
        return categories.length;
    }

    @Override
    public double getArrayValue() {
        return categories[0];
    }

    @Override
    public double getArrayValue(final int dim) {
        if (dim > 3) {
            throw new IllegalArgumentException();
        }
        return categories[dim];
    }

    /**
     * Loggable implementation *
     */

    @Override
    public void init(final PrintStream out) {
        String id = getID();
        if (id == null) {
            id = "";
        }
        final int valueCount = getDimension();
        if (valueCount == 1) {
            out.print(getID() + "\t");
        } else {
            for (int value = 0; value < valueCount; value++) {
                out.print(getID() + (value + 1) + "\t");
            }
        }
    }

    public void log(final long sample, final PrintStream out) {
        final int valueCount = getDimension();
        for (int i=0;i<tree.getNodeCount()+tree.getLeafNodeCount();i++){
            if (i<tree.getLeafNodeCount() && ((QuasiSpeciesNode) tree.getNode(i)).getAttachmentTimesList().length==1){
                categories[i]=0;
            } else if (i<tree.getNodeCount() && tree.getNode(i).isRoot()){
                categories[i]=0;
            } else if (i>=tree.getLeafNodeCount() && i<tree.getNodeCount() && ((QuasiSpeciesNode) tree.getNode(i)).getContinuingHaploName()!=-1){
                categories[i] = ((QuasiSpeciesUCRelaxedClockModel) branchRateModel).getCategories(((QuasiSpeciesNode) tree.getNode(i)).getContinuingHaploName());
            } else {
                categories[i] = ((QuasiSpeciesUCRelaxedClockModel) branchRateModel).getCategories(i);
            }
        }
        for (int i = 0; i < valueCount; i++) {
            out.print(categories[i] + "\t");
        }
    }

    @Override
    public void close(final PrintStream out) {
        // nothing to do
    }
}
