package piqmee.util;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import piqmee.tree.QuasiSpeciesNode;
import piqmee.tree.QuasiSpeciesTree;

import java.io.PrintStream;

/**
 * @author Veronika Boskova created on 11/01/17.
 */
@Description("Allows logging the absolute height of a quasi-species tree; not only of the tree on unique haplotypes.")
public class AbsoluteTreeHeight extends CalculationNode implements Function, Loggable {

    public Input<QuasiSpeciesTree> quasiSpeciesTreeInput = new Input<>(
            "quasiSpeciesTree",
            "Quasi-species tree whose attachment times should be logged.",
            Input.Validate.REQUIRED);

    protected QuasiSpeciesTree qsTree;

    public AbsoluteTreeHeight() { };

    @Override
    public void initAndValidate() {
        qsTree = quasiSpeciesTreeInput.get();
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        double absoluteHeight = qsTree.getRoot().getHeight();
        for (Node node : qsTree.getExternalNodes()){
            double[] attachmentTimes = ((QuasiSpeciesNode) node).getAttachmentTimesList();
            if (attachmentTimes.length > 1 && attachmentTimes[1] > absoluteHeight){
                absoluteHeight = attachmentTimes[1];
            }
        }
        return absoluteHeight;
    }

    @Override
    public double getArrayValue(int iDim) {
        double absoluteHeight = qsTree.getRoot().getHeight();
        for (Node node : qsTree.getExternalNodes()){
            double[] attachmentTimes = ((QuasiSpeciesNode) node).getAttachmentTimesList();
            if (attachmentTimes.length > 1 && attachmentTimes[1] > absoluteHeight){
                absoluteHeight = attachmentTimes[1];
            }
        }
        return absoluteHeight;
    }

    @Override
    public void init(PrintStream out){
        if (getID() == null || getID().matches("\\s*")) {
            out.print(qsTree.getID() + ".absoluteHeight\t");
        } else {
            out.print(getID() + "\t");
        }
    }

    @Override
    public void log(long nSample, PrintStream out) {
        double absoluteHeight = qsTree.getRoot().getHeight();
        for (Node node : qsTree.getExternalNodes()){
            double[] attachmentTimes = ((QuasiSpeciesNode) node).getAttachmentTimesList();
            if (attachmentTimes.length > 1 && attachmentTimes[1] > absoluteHeight){
                absoluteHeight = attachmentTimes[1];
            }
        }
        out.print(absoluteHeight + "\t");
    }

    @Override
    public void close(PrintStream out) { }

}