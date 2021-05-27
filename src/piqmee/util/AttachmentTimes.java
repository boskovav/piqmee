package piqmee.util;

import beast.core.*;
import beast.evolution.tree.Node;
import piqmee.tree.QuasiSpeciesTree;
import piqmee.tree.QuasiSpeciesNode;

import java.io.PrintStream;

/**
 * @author Veronika Boskova created on 03/08/15.
 */
@Description("Allows logging and attachment times of the haplotypes of a quasi-species tree.")
public class AttachmentTimes extends CalculationNode implements Function, Loggable {

    public Input<QuasiSpeciesTree> quasiSpeciesTreeInput = new Input<>(
            "quasiSpeciesTree",
            "Quasi-species tree whose attachment times should be logged.",
            Input.Validate.REQUIRED);

    public Input<String> haplotypeInput = new Input<>(
            "haplotype", "Haplotype name, specifying for which haplotype" +
            " in the tree the attachement times shoudl be logged.",
            Input.Validate.REQUIRED);

    private QuasiSpeciesTree qsTree;

    private String haplotype;

    private int haploNodeNr = -1;

    private QuasiSpeciesNode haploNode;

    public AttachmentTimes() { };

    @Override
    public void initAndValidate() {
        qsTree = quasiSpeciesTreeInput.get();
        haplotype = haplotypeInput.get();
        for (Node node : qsTree.getExternalNodes()){
            if (haplotype.equals(node.getID().toString())){
                haploNodeNr = node.getNr();
            }
        }
        // check if the tip name has been found in the tree
        if (haploNodeNr < 0) {
            throw new RuntimeException("Tip "+haplotype+" could not be found. " +
                    "Maybe it was collapsed with another sequence?");
        }
    }

    @Override
    public int getDimension() {
        return ((QuasiSpeciesNode) qsTree.getNode(haploNodeNr)).getAttachmentTimesList().length;
    }

    @Override
    public double getArrayValue() {
        return ((QuasiSpeciesNode) qsTree.getNode(haploNodeNr)).getAttachmentTimesList()[0];
    }

    @Override
    public double getArrayValue(int iDim) {
        if (iDim < getDimension()) {
            return ((QuasiSpeciesNode) qsTree.getNode(haploNodeNr)).getAttachmentTimesList()[iDim];
        } else
            return Double.NaN;
    }

    @Override
    public void init(PrintStream out){

        String idString = qsTree.getID();
        haploNode = (QuasiSpeciesNode) qsTree.getNode(haploNodeNr);
        int maxtime = haploNode.getAttachmentTimesList().length;
        for (int time = 0; time < maxtime; time++) {
            String haploname = haploNode.getID();
            out.print(idString + "_" + haploname + "_" + time + "\t");
        }
    }

    @Override
    public void log(int nSample, PrintStream out) {
        double[] attachmentTimes = haploNode.getAttachmentTimesList();
        for (int time = 0; time < attachmentTimes.length; time++) {
            out.print(attachmentTimes[time] + "\t");
        }
    }

    @Override
    public void close(PrintStream out) { }
}