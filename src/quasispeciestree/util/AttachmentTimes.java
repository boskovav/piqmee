package quasispeciestree.util;

import beast.core.*;
import beast.evolution.tree.Node;
import beast.evolution.tree.QuasiSpeciesTree;
import beast.evolution.tree.QuasiSpeciesNode;

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

    private QuasiSpeciesNode haploNode;


    public AttachmentTimes() { };

    @Override
    public void initAndValidate() {
        qsTree = quasiSpeciesTreeInput.get();
        haplotype = haplotypeInput.get();
        for (Node node : qsTree.getExternalNodes()){
            if (haplotype.equals(node.getID().toString())){
                haploNode = (QuasiSpeciesNode) node;
            }
        }
    }

    @Override
    public int getDimension() {
        return qsTree.getAttachmentTimesList(haploNode).length;
    }

    @Override
    public double getArrayValue() {
        return qsTree.getAttachmentTimesList(haploNode)[0];
    }

    @Override
    public double getArrayValue(int iDim) {
        if (iDim<getDimension()) {
            return qsTree.getAttachmentTimesList(haploNode)[iDim];
        } else
            return Double.NaN;
    }

    @Override
    public void init(PrintStream out) throws Exception {

        String idString = qsTree.getID();
        int maxtime=qsTree.getAttachmentTimesList(haploNode).length;
        for (int time = 0; time < maxtime; time++) {
            String haploname = haploNode.getID();
            // TODO this is not printing
            out.print(idString + "." + haploname + "." + time + "\t");
        }
    }

    @Override
    public void log(int nSample, PrintStream out) {

        for (int time = 0; time < qsTree.getAttachmentTimesList(haploNode).length; time++) {
            out.print(qsTree.getAttachmentTimesList(haploNode)[time] + "\t");
        }
    }

    @Override
    public void close(PrintStream out) { }


}