package piqmee.util;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import piqmee.tree.QuasiSpeciesTree;
import piqmee.tree.QuasiSpeciesNode;

import java.io.PrintStream;

/**
 * @author Veronika Boskova created on 03/08/15.
 */
@Description("Allows logging and attachment times of the haplotypes of a quasi-species tree.")
public class AttachmentTimesAll extends CalculationNode implements Function, Loggable {

    public Input<QuasiSpeciesTree> quasiSpeciesTreeInput = new Input<>(
            "quasiSpeciesTree",
            "Quasi-species tree whose attachment times should be logged.",
            Input.Validate.REQUIRED);

    private QuasiSpeciesTree qsTree;

    public AttachmentTimesAll() { };

    @Override
    public void initAndValidate() {
        qsTree = quasiSpeciesTreeInput.get();
    }

    @Override
    public int getDimension() {
        int size = 0;
        for (Node node : qsTree.getExternalNodes()){
            size += ((QuasiSpeciesNode) node).getAttachmentTimesList().length;
        }
        return size;
    }

    @Override
    public double getArrayValue(int iDim) {
        if (iDim < getDimension()) {
            int index = iDim;
            int size;
            for (Node node : qsTree.getExternalNodes()){
                size = ((QuasiSpeciesNode) node).getAttachmentTimesList().length;
                if (index < size) {
                    return ((QuasiSpeciesNode) node).getAttachmentTimesList()[index];
                } else {
                    index -= size;
                }
            }
            return Double.NaN;
        } else
            return Double.NaN;
    }

    @Override
    public void init(PrintStream out){

        String idString = qsTree.getID();
        // print all haplo names
        for (Node node : qsTree.getExternalNodes()){
            int maxtime = ((QuasiSpeciesNode) node).getAttachmentTimesList().length;
            for (int time = 0; time < maxtime; time++) {
                String haploname = node.getID();
                out.print(idString + "_" + haploname + "_" + time + "\t");
            }
        }
    }

    @Override
    public void log(long nSample, PrintStream out) {
        // print all haplo names
        for (Node node : qsTree.getExternalNodes()) {
            double[] attachmentTimes = ((QuasiSpeciesNode) node).getAttachmentTimesList();
            for (int time = 0; time < attachmentTimes.length; time++) {
                out.print(attachmentTimes[time] + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) { }
}