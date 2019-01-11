package piqmee.distributions;

import beast.core.Citation;
import beast.core.Description;
import beast.core.util.Log;
import beast.evolution.speciation.BirthDeathSkylineModel;

import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import piqmee.tree.QuasiSpeciesNode;
import piqmee.tree.QuasiSpeciesTree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

/**
 *  @author Veronika Boskova created on 26/06/2015
 *
 *         maths: Tanja Stadler and Veronika Boskova
 *
 */
@Description("Model for calculating birth-death likelihood for quasispecies (ultrametric) trees ")

@Citation("Boskova, V., Stadler, T., Quasispecies algorithm")
public class QuasiSpeciesBirthDeathSkylineModel extends BirthDeathSkylineModel{

    // Empty constructor as required:
    public QuasiSpeciesBirthDeathSkylineModel() { };

    ArrayList isRhoTip;
    double[] uniqueSampTimes;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        isRhoTip = new ArrayList<Boolean>(treeInput.get().getLeafNodeCount());
        for (int i = 0; i < treeInput.get().getLeafNodeCount(); i++){
            QuasiSpeciesNode node = (QuasiSpeciesNode) treeInput.get().getNode(i);
            double[] tipTimes = node.getTipTimesList();
            boolean[] isRhoTipArray = new boolean[tipTimes.length];
            isRhoTip.add(i,isRhoTipArray);
        }

        uniqueSampTimes = getArrayOfUniqueSamplingTimes(treeInput.get());

        if(SAModel || r!=null)
            throw new IllegalArgumentException("The sampled ancestor model has not been implemented to work with quasispecies model yet");

    }

    /**
     * Get number of tip that are sampled during sampling efforts
     *
     */
    public int getTotalN(){
        int totalN = 0;
        if (N != null){
            for (int i = 0; i < N.length; i++){
                totalN += N[i];
            }
        }
        return totalN;
    }

    /**
     * Sorts the array of values in descending order
     *
     */
    public void sortListDescending(double[] array) {
        Arrays.sort(array);
        // reverse the array to start with the largest value
        invertArray(array);
    }

    /**
     * Inverts the array of values
     *
     */
    public void invertArray(double[] array) {
        int totalLength = array.length;
        for (int j = 0; j < totalLength / 2; j++) {
            double tmp = array[j];
            array[j] = array[totalLength - 1 - j];
            array[totalLength - 1 - j] = tmp;
        }
    }

    /**
     * Inverts the array of values
     *
     */
    public void invertArray(int[] array) {
        int totalLength = array.length;
        for (int j = 0; j < totalLength / 2; j++) {
            int tmp = array[j];
            array[j] = array[totalLength - 1 - j];
            array[totalLength - 1 - j] = tmp;
        }
    }

    /**
     * Counts number of unique sampling times of leaves in a tree
     *
     */
    public double[] getArrayOfUniqueSamplingTimes(TreeInterface tree){

        List<Double> uniqueTimes = new ArrayList();

        int nrTips = tree.getLeafNodeCount();

        for (int i = 0; i < nrTips; i++) {
            Node node = tree.getNode(i);
            double[] samplingtimes = ((QuasiSpeciesNode) node).getTipTimesList();
            for (int time = 0; time < samplingtimes.length; time++) {
                if (!uniqueTimes.contains(samplingtimes[time])) {
                    uniqueTimes.add(samplingtimes[time]);
                }
            }
        }
        // cast from Double to double
        double[] uniqueTimesArray = new double[uniqueTimes.size()];
        for (int i = 0; i < uniqueTimes.size(); i++) {
            uniqueTimesArray[i] = uniqueTimes.get(i);
        }

        sortListDescending(uniqueTimesArray);

        return uniqueTimesArray;
    }

    /**
     * Get number of full topologies the QS tree represents
     *
     */
    public double logNumberOfQSTrees(TreeInterface tree) {

        double gamma = 0;

        QuasiSpeciesTree qsTree = (QuasiSpeciesTree) tree;

        // make 3 time arrays - 1) sampling times -- this one is made in initAndValidate already
        //                      2) true internal node times and
        //                      3) merge of these
        // store in 2 arrays a) number of bifurcations and b) total number of lineages A at time t0, t1 etc
        // at each time t0, t1 etc add the gamma contribution from A -> knowing number of bifurcations and total nr of lineages
        // add total number of QS lineages to TOTAL number of all lineages array

        // times arrays
        int nrTips = tree.getLeafNodeCount();
        int nrINTnodes = tree.getInternalNodeCount();
        double[] INT = new double[nrINTnodes];
        // fill time arrays
        for (int i = nrTips; i < nrTips+nrINTnodes; i++){
            // internal node, add its height to INT
            INT[i-nrTips] = ((tree.getNode(i)).getHeight());
        }

        // make array of sort indexes for INT --- for checking later if the node belonging to this height has QS passing through
        int[] indexes = IntStream.range(0, INT.length).boxed()
                .sorted((i, j) -> ((Double)INT[i]).compareTo(INT[j])).mapToInt(ele -> ele).toArray();
        // get in descending order
        invertArray(indexes);

        // sort time arrays in descending order
        sortListDescending(INT);
        double[] allTimes = new double[uniqueSampTimes.length + INT.length];
        System.arraycopy(uniqueSampTimes,0, allTimes,0, uniqueSampTimes.length);
        System.arraycopy(INT, 0, allTimes, uniqueSampTimes.length, INT.length);
        sortListDescending(allTimes);

        // arrays to store number of bifurcations / total number of lineages of QS and total number of all lineages
        // QS
        int[] nrqsattachments = new int[allTimes.length];
        int[] nrqslineages = new int[allTimes.length];
        // total
        //int[] nrtotallineages = new int[allTimes.length];
        //Arrays.fill(nrtotallineages,1);
        //int[] nrtotalqsattachments = new int[allTimes.length];

        // count for each haplo at each time point the n's and add to total count
        for (Node node : tree.getExternalNodes()){
            Arrays.fill(nrqsattachments,0);
            Arrays.fill(nrqslineages,0);
            double[] QSTimesTemp = ((QuasiSpeciesNode) node).getAttachmentTimesList();
            double[] QSTipTimesTemp = ((QuasiSpeciesNode) node).getTipTimesList();
            int[] QSTipTimesCountTemp = ((QuasiSpeciesNode) node).getTipTimesCountList();
            // pointers for qs arrays
            int qstimep = 1;
            int qstiptimep = QSTipTimesTemp.length - 1;
            // pointer for internal node times - to associate INT node with QS correctly
            int intp = 0;
            for (int i = 0; i < allTimes.length; i++) {
                // new bifurcations between i-1 and i, add a lineage to qs and total count
                while (qstimep < QSTimesTemp.length && QSTimesTemp[qstimep] > allTimes[i]) {
                    nrqsattachments[i] += 1;
                    qstimep++;
                }
                nrqslineages[i] += nrqsattachments[i];

                if (i+1 < allTimes.length){
                    // number of lineages at time i+1 is the same as at i ... + bifurcations (i to i+1) - sampling (i)
                    nrqslineages[i + 1] += nrqslineages[i];
                    // if the time is different from previous time - then take into account possible birth/death of lineages
                    if (allTimes[i] == allTimes[0] || allTimes[i] != allTimes[i-1]) {
                        // the time in allTimes[i] can be a time for sampling at time i, remove lineages from time i+1
                        if (qstiptimep >= 0 && QSTipTimesTemp[qstiptimep] == allTimes[i]) {
                            nrqslineages[i + 1] -= QSTipTimesCountTemp[qstiptimep];
                            qstiptimep--;
                        }
                        // the time in allTimes[i] can be the time for a real internal node
                        //              - account for all possible QS lineages it can attach to
                        // check if the internal node with this height belongs to this haplotype
                        while (intp < indexes.length && tree.getNode(indexes[intp] + nrTips).getHeight() > allTimes[i])
                            intp++;
                        if (intp < indexes.length
                                && tree.getNode(indexes[intp] + nrTips).getHeight() == allTimes[i]
                                && ((QuasiSpeciesNode) tree.getNode(indexes[intp] + nrTips)).getContinuingHaploName() == node.getNr()) {
                            // if it does, add contribution to the gamma factor for possible attachment branches
                            gamma += Math.log(nrqslineages[i] + 1);
                        }
                    }
                }

                //nrtotallineages[i] += nrqslineages[i];
                //nrtotalqsattachments[i] += nrqsattachments[i];
                // this factor is only needed for trees with tips sampled through time
                if (uniqueSampTimes.length > 1) {
                    // include all the (gammaj) factors for the possible combinations of QS lineages at each merge point
                    for (int lineage = nrqslineages[i]; lineage > nrqslineages[i] - nrqsattachments[i]; lineage--) {
                        // for qs lineages we have to have lineage + 1 since we did not count that at first split,
                        //  we created 2 lineages instead of just one, as at all later splits
                        gamma += Math.log(lineage) + Math.log(lineage + 1);
                    }
                }
            }
            // check if we checked everything in QS attachment array
            if (qstimep < qsTree.getHaplotypeCounts(node)){// && QSTimesTemp.length > 0){
                throw new RuntimeException("There is somethings wrong with accounting for attachments times of node " +
                        node.getNr() + ". Please check line 242 in QuasiSpeciesBirthDeathSkylineModel class.");
            }
        }

//        if (uniqueSampTimes.length > 1) {
//            // add factor for the possible combinations of QS lineages with all the lineages presentx in tree at time i
//            for (int i = 0; i < allTimes.length; i++) {
//                // add lineage count for when there is a birth of new lineage through true internal node
//                // This was not done above as we looped through each external node, so we had many possibilities to increase
//                //      the lineage due to the internal node. We could have opted to increase the lineage count due to the
//                //      internal node when the node's gamma contribution above was added, but it can happen that there is
//                //      no QS lineage passing through the node and we would then never have increased the lineage count.
//                for (Node node : tree.getInternalNodes()) {
//                    if ((allTimes[i] == allTimes[0] || allTimes[i] != allTimes[i-1]) && node.getHeight() == allTimes[i]) {
//                        for (int j = i + 1; j < allTimes.length; j++) {
//                            nrtotallineages[j] += 1;
//                        }
//                    }
//                }
//                //
//                for (int lineage = nrtotallineages[i]; lineage > nrtotallineages[i] - nrtotalqsattachments[i]; lineage--) {
//                    // for total lineages we have to have lineage - 1 since we did already count also the real nodes in
//                    gamma -= ((Math.log(lineage) + Math.log(lineage - 1)));
//                }
//            }
//        }

        return gamma;
    }

    /**
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

            if (node.getHeight() >= time) count -= 1;
            // node.getHeight() should be the same as the tipTimes[0]
            else {
                double[] attachTimes = node.getAttachmentTimesList();
                double[] tipTimes = node.getTipTimesList();
                int[] tipTimeCounts = node.getTipTimesCountList();
                // start at position attachTimes.length-1 and stop at 1, since position 0 is the "fake" start of the haplo
                int position = attachTimes.length-1;
                if (tipTimes[0] < time){
                    // first - 1 since if we have e.g. last three tips, we only have 2 internal nodes (3-1)=2
                    //          ... but this is still not enough, thus need a second -1
                    // second - 1 for the correct stopping point. If there are e.g. two internal nodes out of 10,
                    //          corresponding to the last 3 tips go from l = 9 to 10
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

            if (node.getHeight() >= time) count -= 1;
                // node.getHeight() should be the same as the tipTimes[0]
            else {
                double[] attachTimes = node.getAttachmentTimesList();
                double[] tipTimes = node.getTipTimesList();
                int[] tipTimeCounts = node.getTipTimesCountList();
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
                //int gamma = node.getStartBranchCounts();
                //logP += Math.log(gamma);
                // NOTE: this is all done in logNumberOfQSTrees() below
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
            //double[] QSTipTimesTemp = ((QuasiSpeciesNode) node).getTipTimesList();
            for (int j = nQSTemp - 1; j > 0; j--) {
                double x = times[totalIntervals - 1] - QSTimesTemp[j];
                index = index(x);
                temp = Math.log(birth[index]) + log_q(index, times[index], x);
                logP += temp;
                if (printTempResults)
                    System.out.println("1st pwd" + " = " + temp + "; QSinterval & QS attachment branches = " + node.getID() + " " + j);
                if (Double.isInfinite(logP))
                    return logP;
            }
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

        // factor for all possible QS trees
        logP += logNumberOfQSTrees(tree);
        if (Double.isInfinite(logP))
            return logP;

//        if (SAModel) {
//            int internalNodeCount = tree.getLeafNodeCount() - ((Tree)tree).getDirectAncestorNodeCount()- 1;
//            logP +=  Math.log(2)*internalNodeCount;
//        }

        return logP;
    }

}
