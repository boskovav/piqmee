package piqmee.util;
import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import piqmee.tree.QuasiSpeciesNode;
import piqmee.tree.QuasiSpeciesTree;

import java.io.PrintStream;
import java.util.Arrays;

/**
 *  @author Veronika Boskova created on 13/02/20
 */

public class QuasiSpeciesNodeTreeLogger extends BEASTObject implements Loggable {

        public Input<QuasiSpeciesTree> qsTreeInput = new Input<>(
                "tree",
                "Quasi-species tree to log.",
                Input.Validate.REQUIRED);

        QuasiSpeciesTree qsTree;

        @Override
        public void initAndValidate() {
            qsTree = qsTreeInput.get();
        }

        @Override
        public void init(PrintStream out) {
            qsTree.init(out);
        }

        @Override
        public void log(long nSample, PrintStream out) {

            // Set up metadata string
            for (Node node : qsTree.getExternalNodes()) {
                QuasiSpeciesNode qsNode = (QuasiSpeciesNode) node;
                int nodeNr = node.getNr();
                double[] tempqstimes = qsNode.getAttachmentTimesList();
                qsNode.metaDataString = String.format("DuplicateBranchingTimes={%s}",Arrays.toString(tempqstimes));
            }


            out.print("tree STATE_" + nSample + " = ");
            out.print(qsTree.getRoot().toSortedNewick(new int[1], true));
            out.print(";");
        }

        @Override
        public void close(PrintStream out) {
            qsTree.close(out);
        }

}
