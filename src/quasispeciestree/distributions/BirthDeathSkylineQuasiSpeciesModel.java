package quasispeciestree.distributions;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.speciation.BirthDeathSkylineModel;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import quasispeciestree.operators.*;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;

/**
 *  @author Veronika Boskova created on 26/06/2015
 *
 *         maths: Tanja Stadler and Veronika Boskova
 *
 */
@Description("Model for calculating birth-death likelihood for quasispecies (ultrametric) trees ")

@Citation("Boskova, V., Stadler, T., Quasispecies algorithm")
public class BirthDeathSkylineQuasiSpeciesModel extends BirthDeathSkylineModel{

//    final public Input<QuasiSpeciesTree> qsTreeInput = new Input<QuasiSpeciesTree>(
//            "quasiSpeciesTree", "Quasi-Species tree over which to calculate a prior or likelihood",
//            Input.Validate.REQUIRED);

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if(SAModel || r!=null) {
            System.out.println("The sampled ancestor model has not been implemented to work with quasispecies model yet");
            System.exit(0);
        }
    }

    /*
     * Adds the number of qs duplicates to the count of tips at each of the contemporaneous sampling times ("rho" sampling time)
     * @return negative infinity if tips are found at a time when rho is zero, zero otherwise.
     */
    @Override
    public double computeN(TreeInterface tree) {

        double superReturn = super.computeN(tree);

        QuasiSpeciesTree qsTree = (QuasiSpeciesTree) tree;

        int nQS = qsTree.getTotalAttachmentCounts();
        int tipCount = tree.getLeafNodeCount();

        double[] dates = new double[nQS];

        int nCount = 0;
        for (Node node : tree.getExternalNodes()){
            double QSTimesLength = qsTree.getHaplotypeCounts((QuasiSpeciesNode)node);
            for (int j=0; j<QSTimesLength; j++){
                // TODO rewrite this to make sure all the tip sampling points are correctly assigned for qs tips
                dates[nCount+j] = node.getHeight();
            }
            nCount += QSTimesLength;
        }

        for (int k = 0; k < totalIntervals; k++) {


            for (int i = 0; i < nQS; i++) {
                if (Math.abs((times[totalIntervals - 1] - times[k]) - dates[i]) < 1e-10) {
                    if (rho[k] == 0 && psi[k] == 0) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    if (rho[k] > 0) {
                        N[k] += 1;
                    }
                }
            }
        }
        return superReturn;
    }


    /*    calculate and store Ai, Bi and p0        */
    @Override
    public Double preCalculation(TreeInterface tree) {

        QuasiSpeciesTree qsTree = (QuasiSpeciesTree) tree;

        double maxheight= tree.getRoot().getHeight();
        for (Node node : tree.getExternalNodes()){
            double[] attachmentTimes = qsTree.getAttachmentTimesList((QuasiSpeciesNode) node);
            if (attachmentTimes.length > 1 && attachmentTimes[1] > maxheight){
                maxheight = attachmentTimes[1];
            }
        }
        if (origin.get() != null && (!originIsRootEdge.get() && maxheight >= origin.get().getValue())) {
            return Double.NEGATIVE_INFINITY;
        }

        return super.preCalculation(tree);
    }


    /**
     * @param time the time
     * @param tree the tree
     * @return the number of lineages that exist at the given time in the given tree.
     */
    @Override
    public int lineageCountAtTime(double time, TreeInterface tree) {

        QuasiSpeciesTree qsTree = (QuasiSpeciesTree) tree;

        int count = super.lineageCountAtTime(time, tree);

        for (Node node : tree.getExternalNodes()){
            double[] QSTimesTemp = qsTree.getAttachmentTimesList((QuasiSpeciesNode) node);
            if (QSTimesTemp[0] > time && node.getHeight() < time){
                for (int j=1; j<QSTimesTemp.length; j++){
                    if (QSTimesTemp[j] > time) count += 1;
                }
            }
        }
        return count;
    }

    /**
     * @param time the time
     * @param tree the tree
     * @param k count the number of sampled ancestors at the given time
     * @return the number of lineages that exist at the given time in the given tree.
     */
    public int lineageCountAtTime(double time, QuasiSpeciesTree tree, int[] k) {

        QuasiSpeciesTree qsTree = (QuasiSpeciesTree) tree;

        int count = super.lineageCountAtTime(time, tree, k);

        for (Node node : tree.getExternalNodes()){
            double[] QSTimesTemp = tree.getAttachmentTimesList((QuasiSpeciesNode) node);
            if (QSTimesTemp[0] > time && node.getHeight() < time){
                for (int j=1; j<QSTimesTemp.length; j++){
                    if (QSTimesTemp[j] > time) count += 1;
                }
            }
        }
        return count;
    }

    // calculateTreeLogLikelihood also adapted from Denise's code
    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {

        QuasiSpeciesTree qsTree = (QuasiSpeciesTree) tree;

        logP = 0.;


        int nTips = tree.getLeafNodeCount();

        // nTips represent only the unique haplotypes, we use nQS for the counts of repetitive haplotypes,
        //      together nTips+nQs=N in Stadler 2013 (BDSKY) paper
        //        int nQS = 0;
        // do not use       for (int i=0;i<tree.getExternalNodes().size();i++){
        // do not use          nQS += tree.getQuasiSpeciesCounts(tree.getExternalNodes().get(i));
        // do not use       }
        //      computation separately done in computeN function and first factor of likelihood calculation below

        if (preCalculation(tree) < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // number of lineages at each time ti
        int[] n = new int[totalIntervals];

        int index = 0;
        if (times[index] < 0.)
            index = index(0.);

        double x0 = 0.;
        double temp = 0.;

        switch (conditionOn) {
            case NONE:
                temp = log_q(index, times[index], x0);
                break;
            case SURVIVAL:
                temp = p0(index, times[index], x0);
                if (temp == 1)
                    return Double.NEGATIVE_INFINITY;
                if (conditionOnRootInput.get()) {
                    temp = log_q(index, times[index], x0) - 2 * Math.log(1 - temp) - Math.log(birth[index]);
                } else {
                    temp = log_q(index, times[index], x0) - Math.log(1 - temp);
                }
                break;
            case RHO_SAMPLING:
                temp = p0hat(index, times[index], x0);
                if (temp == 1)
                    return Double.NEGATIVE_INFINITY;
                if (conditionOnRootInput.get()) {
                    temp = log_q(index, times[index], x0) - 2 * Math.log(1 - temp);
                } else {
                    temp = log_q(index, times[index], x0) - Math.log(1 - temp);
                }
                break;
            default:
                break;
        }

        logP = temp;
        if (Double.isInfinite(logP))
            return logP;

        if (printTempResults) System.out.println("first factor for origin = " + temp);


        // first product term in f[T] over all non-QS transmission times (for the tips sampled through time and at times of parameter change)
        // to start with, get array containing possible number branches the haplotype can start from
//        int[] startBranchCountArray= tree.countPossibleStartBranches();
        int[] startBranchCountArray = qsTree.getStartBranchCounts();
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {

            double x = times[totalIntervals - 1] - tree.getNode(nTips + i).getHeight();
            index = index(x);
            if (!(tree.getNode(nTips + i)).isFake()) {
                temp = Math.log(birth[index]) + log_q(index, times[index], x);
                logP += temp;
                if (printTempResults) System.out.println("1st pwd" +
                        " = " + temp + "; interval = " + i);
                if (Double.isInfinite(logP))
                    return logP;
                // term for the Quasi-Species tree likelihood calculation counting possible QS start branches
                int gamma=startBranchCountArray[i];
                temp = gamma;
                logP += Math.log(temp);
//                logP += Math.log(startBranchCountArray[i]) - Math.log(lineageCountAtTime(tree.getNode(nTips + i).getHeight(),tree));
// testing
//                System.out.println("1st pwd" +
//                        " = " + Math.log(temp) + "; QS start branches = " + tree.getNode(nTips + i).getID());
                if (printTempResults) System.out.println("1st pwd" +
                        " = " + temp + "; QS start branches = " + tree.getNode(nTips + i).getID());
                if (Double.isInfinite(logP))
                    return logP;
            }
        }
        // first product term in f[T] over all QS transmission times
        //
//        // to start with, get array containing possible number branches the haplotype can start from
//        int[] startBranchCountArray= tree.countPossibleStartBranches();
//        //
        for (Node node : tree.getExternalNodes()){
            // term for the Quasi-Species tree likelihood calculation counting possible QS start branches
//            int gamma=startBranchCountArray[node.getNr()];
//            temp = gamma;
//            logP += Math.log(temp);
//            if (printTempResults) System.out.println("1st pwd" +
//                    " = " + temp + "; QS start branches = " + node.getID());
//            if (Double.isInfinite(logP))
//                return logP;
            int nQSTemp = (int) qsTree.getHaplotypeCounts((QuasiSpeciesNode) node);
            double[] QSTimesTemp = qsTree.getAttachmentTimesList((QuasiSpeciesNode) node);
            for (int j = 1; j <= nQSTemp; j++ ){
                double x = times[totalIntervals - 1] - QSTimesTemp[j];
                index = index(x);
                // term for the Quasi-Species tree likelihood calculation counting possible attachment branches
                // THIS IS A CONSTANT FOR ANY TREE TOPOLOGY, SO OMIT THIS TERM FOR NUMERICAL STABILITY AND CALCULATION SPEED
                //int gammaj = tree.countPossibleAttachmentBranches((QuasiSpeciesNode) node, j, tree.getAttachmentTimesList((QuasiSpeciesNode) node)[j]);
                //temp = Math.log(gammaj * birth[index] * g(index, times[index], x));
                //      times 2 for left OR right -- constant factor also not taken into account
                // temp = Math.log(2 * gammaj * birth[index] * g(index, times[index], x));
//                logP += Math.log(lineageCountAtTime(QSTimesTemp[j], tree)- j) - Math.log(lineageCountAtTime(QSTimesTemp[j], tree));
                temp = Math.log(birth[index]) + log_q(index, times[index], x);
// testing
//                System.out.println(tree.getAttachmentTimesList((QuasiSpeciesNode) node)[j]);
//                System.out.println("1st pwd" +
//                        " = " + temp + "; QSinterval & QS attachment branches = " + node.getID() + " " +j);
                logP += temp;
                if (printTempResults) System.out.println("1st pwd" +
                        " = " + temp + "; QSinterval & QS attachment branches = " + node.getID() + " " +j);
                if (Double.isInfinite(logP))
                    return logP;
            }
        }

        // middle product term in f[T]
        for (int i = 0; i < nTips; i++) {
            // TODO what happens if not all tips have QS counts specified in XML --- we get warning, is this enough?
            if (!isRhoTip[i] || m_rho.get() == null) {
//            || tree.getHaplotypeCounts((QuasiSpeciesNode) tree.getNode(i))>0) {
                int nQSTemp = (int) qsTree.getHaplotypeCounts((QuasiSpeciesNode) tree.getNode(i));
                double y = times[totalIntervals - 1] - tree.getNode(i).getHeight();
                index = index(y);
// testing
//                System.out.println("2nd factor included");
//                if (!(tree.getNode(i)).isDirectAncestor()) {
//                    if (!SAModel) {
                        temp = (nQSTemp+1)*(Math.log(psi[index]) - log_q(index, times[index], y));
//                    } else {
//                        temp = Math.log(psi[index] * (r[index] + (1 - r[index]) * p0(index, times[index], y))) - log_q(index, times[index], y);
//                    }
                logP += temp;
                if (printTempResults) System.out.println("2nd PI = " + temp);
                if (psi[index] == 0 || Double.isInfinite(logP))
                    return logP;
//                } else {
//                    if (r[index] != 1) {
//                        logP += Math.log((1 - r[index])*psi[index]);
//                        if (Double.isInfinite(logP)) {
//                            return logP;
//                        }
//                    } else {
//                        //throw new Exception("There is a sampled ancestor in the tree while r parameter is 1");
//                        System.out.println("There is a sampled ancestor in the tree while r parameter is 1");
//                        System.exit(0);
//                    }
//                }
            }
        }

        // last product term in f[T], factorizing from 1 to m //
        double time;
        for (int j = 0; j < totalIntervals; j++) {
            time = j < 1 ? 0 : times[j - 1];
            int[] k = {0};
//            if (!SAModel) {
            // changed the number of lineages surviving the next parameter change period done in function lineageCountAtTime
            //  to account for the QS lineages that could be also surviving the parameter change time
            n[j] = ((j == 0) ? 0 : lineageCountAtTime(times[totalIntervals - 1] - time, tree));
//            } else {
//                n[j] = ((j == 0) ? 0 : lineageCountAtTime(times[totalIntervals - 1] - time, tree, k));
//            }
            if (n[j] > 0) {
                temp = n[j] * (log_q(j, times[j], time) + Math.log(1 - rho[j-1]));
                logP += temp;
// testing
//                System.out.println("3rd factor (nj loop) = " + temp + "; interval = " + j + "; n[j] = " + n[j]);
                if (printTempResults)
                    System.out.println("3rd factor (nj loop) = " + temp + "; interval = " + j + "; n[j] = " + n[j]);//+ "; Math.log(g(j, times[j], time)) = " + Math.log(g(j, times[j], time)));
                if (Double.isInfinite(logP))
                    return logP;

            }

//            if (SAModel && j>0 && N != null) { // term for sampled leaves and two-degree nodes at time t_i
//                logP += k[0] * (log_q(j, times[j], time) + Math.log(1-r[j])) + //here g(j,..) corresponds to q_{i+1}, r[j] to r_{i+1},
//                        (N[j-1]-k[0])*(Math.log(r[j]+ (1-r[j])*p0(j, times[j], time))); //N[j-1] to N_i, k[0] to K_i,and thus N[j-1]-k[0] to M_i
//                if (Double.isInfinite(logP)) {
//                    return logP;
//                }
//            }

            if (rho[j] > 0 && N[j] > 0) {
                temp = N[j] * Math.log(rho[j]);    // term for contemporaneous sampling
                logP += temp;
// testing
//                System.out.println("3rd factor (Nj loop) = " + temp + "; interval = " + j + "; N[j] = " + N[j]);
                if (printTempResults)
                    System.out.println("3rd factor (Nj loop) = " + temp + "; interval = " + j + "; N[j] = " + N[j]);
                if (Double.isInfinite(logP))
                    return logP;

            }
        }

//        if (SAModel) {
//            int internalNodeCount = tree.getLeafNodeCount() - ((Tree)tree).getDirectAncestorNodeCount()- 1;
//            logP +=  Math.log(2)*internalNodeCount;
//        }

        if (logP > 0)
            return Double.NEGATIVE_INFINITY;
        return logP;
    }

}
