package quasispeciestree.util;

import beast.core.*;
import beast.evolution.tree.Node;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;

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

    private QuasiSpeciesTree qsTree;

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
            Double[] attachmentTimes = qsTree.getAttachmentTimesList((QuasiSpeciesNode) node);
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
            Double[] attachmentTimes = qsTree.getAttachmentTimesList((QuasiSpeciesNode) node);
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
    public void log(int nSample, PrintStream out) {
        double absoluteHeight = qsTree.getRoot().getHeight();
        for (Node node : qsTree.getExternalNodes()){
            Double[] attachmentTimes = qsTree.getAttachmentTimesList((QuasiSpeciesNode) node);
            if (attachmentTimes.length > 1 && attachmentTimes[1] > absoluteHeight){
                absoluteHeight = attachmentTimes[1];
            }
        }
        out.print(absoluteHeight + "\t");
    }

    @Override
    public void close(PrintStream out) { }

}