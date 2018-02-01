package quasispeciestree.distributions;

import beast.core.Citation;
import beast.core.Description;
import beast.evolution.speciation.BirthDeathSkylineModel;

import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;

import java.util.ArrayList;

/**
 *  @author Veronika Boskova created on 26/06/2015
 *
 *         maths: Tanja Stadler and Veronika Boskova
 *
 */
@Description("Model for calculating birth-death likelihood for quasispecies (ultrametric) trees ")

@Citation("Boskova, V., Stadler, T., Quasispecies algorithm")
public class BirthDeathSkylineQuasiSpeciesModel extends BirthDeathSkylineModel{

    // Empty constructor as required:
    public BirthDeathSkylineQuasiSpeciesModel() { };

    ArrayList isRhoTip;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if(SAModel || r!=null)
            throw new IllegalArgumentException("The sampled ancestor model has not been implemented to work with quasispecies model yet");

    }

    /**
     * Factorial calculation for counting number of possible full topologies QS tree represents
     *
     * @param number the number to calculate the factorial of
     */
    public static long logfactorial(int number) {
        long result = 0;

        for (int factor = 2; factor <= number; factor++) {
            result += Math.log(factor);
        }

        return result;
    }

    /**
     * Get number of full topologies the QS tree represents between time m and m-1
     *
     * @param numberOfLineages the number of lineages at time m (closer to the present)
     * @param coalCounter number of coalescent events between times m and m-1
     */
    public double logNumerOfTrees(int numberOfLineages, int coalCounter) {
        double result = 0;

        result  = logfactorial(numberOfLineages)
                + logfactorial(numberOfLineages - 1)
                - Math.log(2) * (numberOfLineages - 1)
                - logfactorial(numberOfLineages - coalCounter)
                - logfactorial(numberOfLineages - coalCounter - 1)
                + Math.log(2) * (numberOfLineages - coalCounter - 1);

        return result;
    }



    /*
     * Adds the number of qs duplicates to the count of tips at each of the contemporaneous sampling times ("rho" sampling time)
     * @return negative infinity if tips are found at a time when rho is zero, zero otherwise.
     */
    @Override
    protected double computeN(TreeInterface tree) {

        isRhoTip = new ArrayList<Boolean>(tree.getLeafNodeCount());

        N = new int[totalIntervals];

        int tipCount = tree.getLeafNodeCount();

        double maxdate = tree.getRoot().getHeight();

        for (int i = 0; i < tipCount; i++) {

            QuasiSpeciesNode node = (QuasiSpeciesNode) tree.getNode(i);
            double[] tipTimes = node.getTipTimesList();
            int[] tipTimeCounts = node.getTipTimesCountList();

            boolean[] isRhoTipArray = new boolean[tipTimes.length];

            for (int index = 0; index < tipTimes.length; index++) {

                for (int k = 0; k < totalIntervals; k++) {

                    if (Math.abs(((times[totalIntervals - 1] - times[k]) - tipTimes[index])/maxdate) < 1e-10 ||
                            (maxdate == 0 && Math.abs((times[totalIntervals - 1] - times[k]) - tipTimes[index]) < 1e-10)) {
                        if (rho[k] == 0 && psi[k] == 0) {
                            return Double.NEGATIVE_INFINITY;
                        }
                        if (rho[k] > 0) {
                            N[k] += tipTimeCounts[index];
                            isRhoTipArray[index] = true;
                        }
                    }
                }
            }
            isRhoTip.add(i,isRhoTipArray);
        }
        return 0.;
    }


    /*    calculate and store Ai, Bi and p0        */
    @Override
    public Double preCalculation(TreeInterface tree) {

        double maxheight= tree.getRoot().getHeight();
        for (Node node : tree.getExternalNodes()){
            double[] attachmentTimes = ((QuasiSpeciesNode) node).getAttachmentTimesList();
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

        int count = 1;
        int tipCount = tree.getLeafNodeCount();
        for (int i = tipCount; i < tipCount + tree.getInternalNodeCount(); i++) {
            if (tree.getNode(i).getHeight() > time) count += 1;
        }

        for (int i = 0; i < tipCount; i++) {
            QuasiSpeciesNode node = (QuasiSpeciesNode) tree.getNode(i);
            double[] attachTimes = node.getAttachmentTimesList();
            double[] tipTimes = node.getTipTimesList();
            int[] tipTimeCounts = node.getTipTimesCountList();

            if (node.getHeight() >= time) count -= 1;
            // node.getHeight() should be the same as the tipTimes[0]
            else {
                // start at position attachTimes.length-1 and stop at 1, since position 0 is the "fake" start of the haplo
                int position = attachTimes.length-1;
                if (tipTimes[0] < time){
                    // first - 1 for if tree last tips, only 2 internal nodes
                    // second - 1 for the correct stopping point if there are two internal nodes out of 10,
                    //      corresponding to the last 3 tips go from 9 to 10
                    for (int l = position - (tipTimeCounts[0] - 1 - 1); l <= position  ; l++){
                        if (attachTimes[l] > time) count += 1;
                        else break;
                    }
                }
                position -= (tipTimeCounts[0] - 1);
                for (int j = 1; j < tipTimes.length; j++) {
                    if (tipTimes[j] < time){
                        for (int l = position - (tipTimeCounts[j] - 1); l <= position; l++){
                            if (attachTimes[l] > time) count += 1;
                            else break;
                        }
                    } else break;
                    position -= tipTimeCounts[j];
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
    @Override
    public int lineageCountAtTime(double time, TreeInterface tree, int[] k) {

        k[0]=0;
        int count = 1;
        int tipCount = tree.getLeafNodeCount();
        for (int i = tipCount; i < tipCount + tree.getInternalNodeCount(); i++) {
            if (tree.getNode(i).getHeight() > time) count += 1;
        }

        for (int i = 0; i < tipCount; i++) {
            QuasiSpeciesNode node = (QuasiSpeciesNode) tree.getNode(i);
            double[] attachTimes = node.getAttachmentTimesList();
            double[] tipTimes = node.getTipTimesList();
            int[] tipTimeCounts = node.getTipTimesCountList();

            if (node.getHeight() >= time) count -= 1;
                // node.getHeight() should be the same as the tipTimes[0]
            else {
                // start at position attachTimes.length-1 and stop at 1, since position 0 is the "fake" start of the haplo
                int position = attachTimes.length-1;
                if (tipTimes[0] < time){
                    // first - 1 for if tree last tips, only 2 internal nodes
                    // second - 1 for the correct stopping point if there are two internal nodes out of 10,
                    //      corresponding to the last 3 tips go from 9 to 10
                    for (int l = position - (tipTimeCounts[0] - 1 - 1); l <= position  ; l++){
                        if (attachTimes[l] > time) count += 1;
                        else break;
                    }
                }
                position -= (tipTimeCounts[0] - 1);
                for (int j = 1; j < tipTimes.length; j++) {
                    if (tipTimes[j] < time){
                        for (int l = position - (tipTimeCounts[j] - 1); l <= position; l++){
                            if (attachTimes[l] > time) count += 1;
                            else break;
                        }
                    } else break;
                    position -= tipTimeCounts[j];
                }
            }
            if (Math.abs(tree.getNode(i).getHeight() - time) < 1e-10) {
                count -= 1;
                if (tree.getNode(i).isDirectAncestor()) {
                    count -= 1;
                    k[0]++;
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
        // to start with, get array containing possible number of branches the true node can start from
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            QuasiSpeciesNode node = (QuasiSpeciesNode) tree.getNode(nTips + i);
            double x = times[totalIntervals - 1] - node.getHeight();
            index = index(x);
            if (!node.isFake()) {
                temp = Math.log(birth[index]) + log_q(index, times[index], x);
                logP += temp;
                if (printTempResults) System.out.println("1st pwd" +
                        " = " + temp + "; interval = " + i);
                if (Double.isInfinite(logP))
                    return logP;
                // term for the Quasi-Species tree likelihood calculation counting possible start branches (gamma)
                temp = node.getStartBranchCounts();
                logP += Math.log(temp);
                if (printTempResults) System.out.println("1st pwd" +
                        " = " + temp + "; QS start branches = " + node.getID());
                if (Double.isInfinite(logP))
                    return logP;
            }
        }
        // first product term in f[T] over all QS transmission times
        //
        for (Node node : tree.getExternalNodes()){
            int nQSTemp = qsTree.getHaplotypeCounts(node);
            double[] QSTimesTemp = ((QuasiSpeciesNode) node).getAttachmentTimesList();
            double[] QSTipTimesTemp = ((QuasiSpeciesNode) node).getTipTimesList();
            int[] QSTipCountsTemp = ((QuasiSpeciesNode) node).getTipTimesCountList();
            // start the tipTime pointer at 1 not 0, since we look at how many lineages coalesced between generation 0 and 1
            int tipTimeIndex = 1;
            int coalCounter = 0;
            int numberOfLineages = QSTipCountsTemp[0];
            double logNumberOfTrees = 0;
            for (int j = nQSTemp - 1; j > 0 ; j-- ){
                double x = times[totalIntervals - 1] - QSTimesTemp[j];
                index = index(x);
                // term for the Quasi-Species tree likelihood calculation counting possible attachment branches
                // THIS IS A CONSTANT FOR ANY TREE TOPOLOGY, SO OMIT THIS TERM FOR NUMERICAL STABILITY AND CALCULATION SPEED
                //int gammaj = tree.countPossibleAttachmentBranches((QuasiSpeciesNode) node, j, tree.getAttachmentTimesList((QuasiSpeciesNode) node)[j]);
                //    gammaj is actually 1 here
                //temp = Math.log(gammaj * birth[index] * g(index, times[index], x));
                //      times 2 for left OR right -- constant factor also not taken into account
                // temp = Math.log(2 * gammaj * birth[index] * g(index, times[index], x));
                temp = Math.log(birth[index]) + log_q(index, times[index], x);
                logP += temp;
                if (printTempResults) System.out.println("1st pwd" +
                        " = " + temp + "; QSinterval & QS attachment branches = " + node.getID() + " " + j);
                if (Double.isInfinite(logP))
                    return logP;

                // plus extra term for quasispecies model accounting for the fact that QS tree may represent several full trees
                if (tipTimeIndex >= QSTipTimesTemp.length){
                    coalCounter += 1;
                    if (j == 1 && coalCounter > 0)
                        logNumberOfTrees += logNumerOfTrees(numberOfLineages,coalCounter);
                }
                else if (QSTimesTemp[j] <= QSTipTimesTemp[tipTimeIndex]) {
                    coalCounter += 1;
                }
                else {
                    if (coalCounter > 0)
                        logNumberOfTrees += logNumerOfTrees(numberOfLineages,coalCounter);
                    numberOfLineages += QSTipCountsTemp[tipTimeIndex] - coalCounter;
                    tipTimeIndex += 1;
                    // one coal time is already above the limit of the next sampling time... so coalCounter is now set to 1
                    coalCounter = 1;
                }
            }

            logP += logNumberOfTrees;
            if (Double.isInfinite(logP))
                return logP;
        }

        // middle product term in f[T]
        for (int i = 0; i < nTips; i++) {

            QuasiSpeciesNode node = (QuasiSpeciesNode) tree.getNode(i);
            boolean[] isRhoTipArray = (boolean[]) isRhoTip.get(i);

            for (int j = 0; j < isRhoTipArray.length; j++) {
                if (!isRhoTipArray[j] || m_rho.get() == null) {
                    int nQSTemp = node.getTipTimesCountList()[j];
                    double y = times[totalIntervals - 1] - node.getTipTimesList()[j];
                    index = index(y);

//                if (!(tree.getNode(i)).isDirectAncestor()) {
//                    if (!SAModel) {
                    temp = nQSTemp * (Math.log(psi[index]) - log_q(index, times[index], y));
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
        }

        // last product term in f[T], factorizing from 1 to m //
        double time;
        for (int j = 0; j < totalIntervals; j++) {
            time = j < 1 ? 0 : times[j - 1];
//            int[] k = {0};
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

        return logP;
    }

}
