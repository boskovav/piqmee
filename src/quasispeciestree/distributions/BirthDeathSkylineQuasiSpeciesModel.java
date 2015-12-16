package quasispeciestree.distributions;

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Taxon;
import beast.evolution.tree.*;
import quasispeciestree.tree.QuasiSpeciesNode;
import quasispeciestree.tree.QuasiSpeciesTree;

import java.util.*;

/**
 *  @author Veronika Boskova created on 26/06/2015
 *
 *         maths: Tanja Stadler and Veronika Boskova
 *
 */

@Description("Model for calculating birth-death likelihood for quasispecies (ultrametric) trees ")

@Citation("Boskova, V., Stadler, T., Quasispecies algorithm \n" +
        "and \n"+
        "Stadler, T., Kuehnert, D., Bonhoeffer, S., and Drummond, A. J. (2013):\n Birth-death skyline " +
        "plot reveals temporal changes of\n epidemic spread in HIV and hepatitis C virus (HCV). PNAS 110(1): 228–33.\n" +
        "If sampled ancestors are used then please also site: Gavryushkina A, Welch D, Stadler T, Drummond AJ (2014) \n" +
        "Bayesian inference of sampled ancestor trees for epidemiology and fossil calibration. \n" +
        "PLoS Comput Biol 10(12): e1003919. doi:10.1371/journal.pcbi.1003919")


public class BirthDeathSkylineQuasiSpeciesModel extends QuasiSpeciesTreeDistribution {
    /******************************************************************************************************************/
    /* Input parameters for the BD model */
    /* same as the bdsky */
    /*
    public Input<Boolean> contemp =
            new Input<Boolean>("contemp", "Only contemporaneous sampling (i.e. all tips are from same sampling time, default false)", false);
    // do we want to implement something with serial sampling?
    */

//    public Input<RealParameter> rhoParameterInput =
//            new Input<RealParameter>("rho", "Proportion of tips sampled at present time (default = 0)"); // rho
    /*
    public Input<RealParameter> rhoParameterInput =
            new Input<RealParameter>("rho", "The proportion of lineages sampled at rho-sampling times (default 0.)");
    // for now not used, needs to be implemented otherwise (bdsky with several sampling rho)
    */

    // Set 1  (lambda, mu, psi)
//    public Input<RealParameter> birthRateParameterInput =
//            new Input<RealParameter>("birthRate", "Birth/Infection rate parameter");
//    public Input<RealParameter> deathRateParameterInput =
//            new Input<RealParameter>("deathRate", "Death rate parameter");
    /*
    public Input<RealParameter> samplingRateParameterInput =
            new Input<RealParameter>("samplingRate", "Sampling rate per individual");
    // for now not used, needs to be implemented otherwise (bdsky with several sampling rho)
    */

    // Set 2 (R0, delta, p)
//    public Input<RealParameter> R0ParameterInput =
//            new Input<RealParameter>("R0", "The basic reproductive number", Input.Validate.XOR, birthRateParameterInput);
//    public Input<RealParameter> becomeUninfectiousRateParameterInput =
//            new Input<RealParameter>("becomeUninfectiousRate", "Rate at which individuals become uninfectious, delta (through recovery or sampling)", Input.Validate.XOR, deathRateParameterInput);
    /*
    public Input<RealParameter> samplingProportionParameterInput =
            new Input<RealParameter>("samplingProportion", "The samplingProportion = samplingRate / becomeUninfectiousRate", Input.Validate.XOR, samplingRateParameterInput);
    // for now not used, needs to be implemented otherwise (bdsky with several sampling rho)
    */
    /******************************************************************************************************************/
    /* Input parameters for phylogeny */
    /* TreeInput XOR TreeIntervalsInput inherited from ancestor */
//    public Input<RealParameter> originParameterInput =
//            new Input<RealParameter>("origin", "The time from origin to last sample (must be larger than tree height)", Input.Validate.REQUIRED);
    /*
    // not sure whether this is necessary
    public Input<Boolean> originIsRootEdge =
            new Input<>("originIsRootEdge", "The origin is only the length of the root edge", false);
    */

    /******************************************************************************************************************/
    /* Random input parameters */
//    public Input<Boolean> conditionOnSurvival =
//            new Input<Boolean>("conditionOnSurvival", "if is true then condition on sampling at least one individual (psi-sampling).", true);
//    public Input<Boolean> conditionOnRhoSampling =
//            new Input<Boolean> ("conditionOnRhoSampling","if is true then condition on sampling at least one individual at present.", false);
    /******************************************************************************************************************/
    /* Shadow parameters */
//    RealParameter    rhoParameter, birthRateParameter, deathRateParameter,
//            R0Parameter, becomeUninfectiousRateParameter,
//            originParameter;
//    TraitSet QuasiSpeciesAttachmentTimes, QuasiSpeciesCounts;
    /* samplingRateParameter,samplingProportionParameter, */


    /* Model parameters */
//    protected TreeInterface tree    = null;
//    protected double [] lambda      = null,   // Infection rates
//                        mu          = null,   // Death rates
//                        psi         = null,   // Sampling rates through time
//                        QStimes     = null;   // Attachment times of the QS
//    protected double    rho         = 0,      // Proportion of tips sampled at present
//                        origin      = 0;      // Origin of the epidemic
//    protected boolean   survival    = true;  // Calculate survival probability
    //protected int       verbose     = 0;      // Not an input parameter, just for debugging (different levels)


    /******************************************************************************************************************/





    // START1 From Denise's BDSKY class

    // the interval times for the birth rate
    public Input<RealParameter> birthRateChangeTimesInput =
            new Input<RealParameter>("birthRateChangeTimes", "The times t_i specifying when birth/R rate changes occur", (RealParameter) null);

    // the interval times for the death rate
    public Input<RealParameter> deathRateChangeTimesInput =
            new Input<RealParameter>("deathRateChangeTimes", "The times t_i specifying when death/becomeUninfectious rate changes occur", (RealParameter) null);

    // the interval times for removal probability
    public Input<RealParameter> samplingRateChangeTimesInput =
            new Input<RealParameter>("samplingRateChangeTimes", "The times t_i specifying when sampling rate or sampling proportion changes occur", (RealParameter) null);

    // the interval times for sampling rate
    public Input<RealParameter> removalProbabilityChangeTimesInput =
            new Input<RealParameter>("removalProbabilityChangeTimes", "The times t_i specifying when removal probability changes occur", (RealParameter) null);

    public Input<RealParameter> intervalTimes =
            new Input<RealParameter>("intervalTimes", "The time t_i for all parameters if they are the same", (RealParameter) null);

    public Input<Boolean> birthRateChangeTimesRelativeInput =
            new Input<Boolean>("birthRateTimesRelative", "True if birth rate change times specified relative to tree height? Default false", false);

    public Input<Boolean> deathRateChangeTimesRelativeInput =
            new Input<Boolean>("deathRateTimesRelative", "True if death rate change times specified relative to tree height? Default false", false);

    public Input<Boolean> samplingRateChangeTimesRelativeInput =
            new Input<Boolean>("samplingRateTimesRelative", "True if sampling rate times specified relative to tree height? Default false", false);

    Input<Boolean> removalProbabilityChangeTimesRelativeInput =
            new Input<Boolean>("removalProbabilityTimesRelative", "True if removal probability change times specified relative to tree height? Default false", false);

    public Input<BooleanParameter> reverseTimeArraysInput =
            new Input<BooleanParameter>("reverseTimeArrays", "True if the time arrays are given in backwards time (from the present back to root). Order: 1) birth 2) death 3) sampling 4) rho 5) r. Default false." +
                    "Careful, rate array must still be given in FORWARD time (root to tips). If rhosamplingTimes given, they should be backwards and this should be true.");

    // the times for rho sampling
    public Input<RealParameter> rhoSamplingTimes =
            new Input<RealParameter>("rhoSamplingTimes", "The times t_i specifying when rho-sampling occurs", (RealParameter) null);


    public Input<RealParameter> origin =
            new Input<RealParameter>("origin", "The time from origin to last sample (must be larger than tree height)", Input.Validate.REQUIRED);

    public Input<Boolean> originIsRootEdge =
            new Input<>("originIsRootEdge", "The origin is only the length of the root edge", false);

    public Input<RealParameter> birthRate =
            new Input<RealParameter>("birthRate", "BirthRate = BirthRateVector * birthRateScalar, birthrate can change over time");
    public Input<RealParameter> deathRate =
            new Input<RealParameter>("deathRate", "The deathRate vector with birthRates between times");
    public Input<RealParameter> samplingRate =
            new Input<RealParameter>("samplingRate", "The sampling rate per individual");      // psi
    public Input<RealParameter> removalProbability =
            new Input<RealParameter>("removalProbability", "The probability of an individual to become noninfectious immediately after the sampling");

    public Input<RealParameter> m_rho =
            new Input<RealParameter>("rho", "The proportion of lineages sampled at rho-sampling times (default 0.)");
    public Input<Boolean> contemp =
            new Input<Boolean>("contemp", "Only contemporaneous sampling (i.e. all tips are from same sampling time, default false)", false);


    public Input<RealParameter> R0 =
            new Input<RealParameter>("R0", "The basic reproduction number", Input.Validate.XOR, birthRate);
    public Input<RealParameter> becomeUninfectiousRate =
            new Input<RealParameter>("becomeUninfectiousRate", "Rate at which individuals become uninfectious (throuch recovery or sampling)", Input.Validate.XOR, deathRate);
    public Input<RealParameter> samplingProportion =
            new Input<RealParameter>("samplingProportion", "The samplingProportion = samplingRate / becomeUninfectiousRate", Input.Validate.XOR, samplingRate);


    public Input<Boolean> forceRateChange =
            new Input<Boolean>("forceRateChange", "If there is more than one interval and we estimate the time of rate change, do we enforce it to be within the tree interval? Default true", true);
    public Input<Boolean> conditionOnSurvival =
            new Input<Boolean>("conditionOnSurvival", "if is true then condition on sampling at least one individual (psi-sampling).", true);
    public Input<Boolean> conditionOnRhoSampling =
            new Input<Boolean> ("conditionOnRhoSampling","if is true then condition on sampling at least one individual at present.", false);

    double t_root;
    protected double[] p0, p0hat;
    protected double[] Ai, Aihat;
    protected double[] Bi, Bihat;
    protected int[] N;   // number of leaves sampled at each time t_i

    // these four arrays are totalIntervals in length
    protected Double[] birth;
    Double[] death;
    Double[] psi;
    Double[] rho;
    Double[] r;

    // true if the node of the given index occurs at the time of a rho-sampling event
    boolean[] isRhoTip;

    /**
     * The number of change points in the birth rate
     */
    protected int birthChanges;

    /**
     * The number of change points in the death rate
     */
    int deathChanges;

    /**
     * The number of change points in the sampling rate
     */
    int samplingChanges;
    int rhoChanges;

    /**
     * The number of change point in the removal probability
     */
    int rChanges;

    /**
     * The number of times rho-sampling occurs
     */
    int rhoSamplingCount;
    Boolean constantRho;

    /**
     * Total interval count
     */
    protected int totalIntervals;

    protected List<Double> birthRateChangeTimes = new ArrayList<Double>();
    protected List<Double> deathRateChangeTimes = new ArrayList<Double>();
    protected List<Double> samplingRateChangeTimes = new ArrayList<Double>();
    protected List<Double> rhoSamplingChangeTimes = new ArrayList<Double>();
    protected List<Double> rChangeTimes = new ArrayList<Double>();

    Boolean contempData;
    //List<Interval> intervals = new ArrayList<Interval>();
    SortedSet<Double> timesSet = new TreeSet<Double>();

    protected Double[] times = new Double[]{0.};

    protected Boolean transform;
    Boolean m_forceRateChange;

    Boolean birthRateTimesRelative = false;
    Boolean deathRateTimesRelative = false;
    Boolean samplingRateTimesRelative = false;
    Boolean rTimesRelative = false;
    Boolean[] reverseTimeArrays;

//    public boolean SAModel;

    enum ConditionOn {NONE, SURVIVAL, RHO_SAMPLING};
    protected ConditionOn conditionOn=ConditionOn.SURVIVAL;

    public Boolean printTempResults;

    // END1












    public void initAndValidate() throws Exception {

        super.initAndValidate();

        // Assign shadow parameters
//        rhoParameter                    = rhoParameterInput.get();
//        birthRateParameter              = birthRateParameterInput.get();
//        deathRateParameter              = deathRateParameterInput.get();
//        R0Parameter                     = R0ParameterInput.get();
//        becomeUninfectiousRateParameter = becomeUninfectiousRateParameterInput.get();
//        originParameter                 = originParameterInput.get();


        // START2
        // start change1 treeInput to qsTreeInput
        if (!originIsRootEdge.get() && qsTreeInput.get().getRoot().getHeight() >= origin.get().getValue())
            throw new RuntimeException("Origin parameter ("+origin.get().getValue()+" ) must be larger than tree height("+qsTreeInput.get().getRoot().getHeight()+" ). Please change initial origin value!");
        // end change1

//        if (removalProbability.get() != null) SAModel = true;

        birth = null;
        death = null;
        psi = null;
        rho = null;
        r = null;
        birthRateChangeTimes.clear();
        deathRateChangeTimes.clear();
        samplingRateChangeTimes.clear();
//        if (SAModel) rChangeTimes.clear();
        totalIntervals = 0;

        m_forceRateChange = forceRateChange.get();
        birthRateTimesRelative = birthRateChangeTimesRelativeInput.get();
        deathRateTimesRelative = deathRateChangeTimesRelativeInput.get();
        samplingRateTimesRelative = samplingRateChangeTimesRelativeInput.get();
//        if (SAModel) rTimesRelative = removalProbabilityChangeTimesRelativeInput.get();

        if (reverseTimeArraysInput.get()!= null )
            reverseTimeArrays = reverseTimeArraysInput.get().getValues();
        else
            reverseTimeArrays = new Boolean[]{false, false, false, false, false};

        contempData = contemp.get();
        rhoSamplingCount = 0;
        printTempResults = false;


        if (birthRate.get() != null && deathRate.get() != null && samplingRate.get() != null) {

            transform = false;
            death = deathRate.get().getValues();
            psi = samplingRate.get().getValues();
            birth = birthRate.get().getValues();
//            if (SAModel) r = removalProbability.get().getValues();

        } else if (R0.get() != null && becomeUninfectiousRate.get() != null && samplingProportion.get() != null) {

            transform = true;

        } else {
            throw new RuntimeException("Either specify birthRate, deathRate and samplingRate OR specify R0, becomeUninfectiousRate and samplingProportion!");
        }


        if (transform) {

            if (birthChanges < 1) birthChanges = R0.get().getDimension() - 1;
            samplingChanges = samplingProportion.get().getDimension() - 1;
            deathChanges = becomeUninfectiousRate.get().getDimension() - 1;

        } else {

            if (birthChanges < 1) birthChanges = birthRate.get().getDimension() - 1;
            deathChanges = deathRate.get().getDimension() - 1;
            samplingChanges = samplingRate.get().getDimension() - 1;
        }

//        if (SAModel) rChanges = removalProbability.get().getDimension() -1;

        if (m_rho.get()!=null) {
            rho = m_rho.get().getValues();
            rhoChanges = m_rho.get().getDimension() - 1;
        }

        collectTimes();

        if (m_rho.get() != null) {

            constantRho = !(m_rho.get().getDimension() > 1);

            if (m_rho.get().getDimension() == 1 && rhoSamplingTimes.get()==null || rhoSamplingTimes.get().getDimension() < 2) {
                if (!contempData && ((samplingProportion.get() != null && samplingProportion.get().getDimension() == 1 && samplingProportion.get().getValue() == 0.) ||
                        (samplingRate.get() != null && samplingRate.get().getDimension() == 1 && samplingRate.get().getValue() == 0.))) {
                    contempData = true;
                    if (printTempResults)
                        System.out.println("Parameters were chosen for contemporaneously sampled data. Setting contemp=true.");
                }
            }

            if (contempData) {
                if (m_rho.get().getDimension() != 1)
                    throw new RuntimeException("when contemp=true, rho must have dimension 1");

                else {
                    rho = new Double[totalIntervals];
                    Arrays.fill(rho, 0.);
                    rho[totalIntervals - 1] = m_rho.get().getValue();
                    rhoSamplingCount = 1;
                }
            }

        } else {
            rho = new Double[totalIntervals];
            Arrays.fill(rho, 0.);
        }
        // start change2
        isRhoTip = new boolean[qsTreeInput.get().getLeafNodeCount()];
        // end change2

        if (conditionOnSurvival.get()) {
            conditionOn = ConditionOn.SURVIVAL;
            if (conditionOnRhoSampling.get()) {
                throw new RuntimeException("conditionOnSurvival and conditionOnRhoSampling can not be both true at the same time." +
                        "Set one of them to true and another one to false.");
            }
        } else if (conditionOnRhoSampling.get()) {
            if (!rhoSamplingConditionHolds()) {
                throw new RuntimeException("Conditioning on rho-sampling is only available for sampled ancestor analyses where r " +
                        "is set to zero and all except the last rho are zero");
            }
            conditionOn = ConditionOn.RHO_SAMPLING;
        } else {
            conditionOn = ConditionOn.NONE;
        }

        printTempResults = false;

        // END2


    }

    // START3

    /**
     * checks if r is zero, all elements of rho except the last one are
     * zero and the last one is not zero
     * @return
     */
    private boolean rhoSamplingConditionHolds() {

//        if (SAModel) {
//            for (int i=0; i<removalProbability.get().getDimension(); i++) {
//                if (removalProbability.get().getValue(i) != 0.0) {
//                    return false;
//                }
//            }
//        } else {
//        return false;
//        }
        for (int i=0; i<rho.length-1; i++) {
            if (rho[i] != 0.0) {
                return false;
            }
        }

        return (rho[rho.length-1] != 0.0);
    }

    /**
     * @return a list of intervals
     */
    public void getChangeTimes(List<Double> changeTimes, RealParameter intervalTimes, int numChanges, boolean relative, boolean reverse) {
        changeTimes.clear();

        if (printTempResults) System.out.println("relative = " + relative);

        // start change3
        double maxTime = originIsRootEdge.get()? qsTreeInput.get().getRoot().getHeight() + origin.get().getValue() :origin.get().getValue();
        // end change 3
        if (intervalTimes == null) { //equidistant

            double intervalWidth = maxTime / (numChanges + 1);

            double end;
            for (int i = 1; i <= numChanges; i++) {
                end = (intervalWidth) * i;
                changeTimes.add(end);
            }
            end = maxTime;
            changeTimes.add(end);

        } else {

            if ((!isBDSIR()) && numChanges > 0 && intervalTimes.getDimension() != numChanges + 1) {
                throw new RuntimeException("The time interval parameter should be numChanges + 1 long (" + (numChanges + 1) + ").");
            }

            int dim = intervalTimes.getDimension();

            ArrayList<Double> sortedIntervalTimes = new ArrayList<>();
            for (int i=0; i< dim; i++) {
                sortedIntervalTimes.add(intervalTimes.getValue(i));
            }
            Collections.sort(sortedIntervalTimes);

            if (!reverse && sortedIntervalTimes.get(0) != 0.0) {
                throw new RuntimeException("First time in interval times parameter should always be zero.");
            }

//            if(intervalTimes.getValue(dim-1)==maxTime) changeTimes.add(0.); //rhoSampling

            double end;
            for (int i = (reverse?0:1); i < dim; i++) {
                end = reverse ? (maxTime - sortedIntervalTimes.get(dim - i - 1)) : sortedIntervalTimes.get(i);
                if (relative) end *= maxTime;
                if (end != maxTime) changeTimes.add(end);
            }
            end = maxTime;
            changeTimes.add(end);
        }
//        }
    }

    /*
    * Counts the number of tips at each of the contemporaneous sampling times ("rho" sampling time)
    * @return negative infinity if tips are found at a time when rho is zero, zero otherwise.
    */
    // VB this function needed to be changed (@ 30.06.2015)
    // start change4
    private double computeN(QuasiSpeciesTree tree) {

        isRhoTip = new boolean[tree.getLeafNodeCount()];

        N = new int[totalIntervals];

        int nQS = 0;
//        for (int i=0;i<tree.getExternalNodes().size();i++){
//            nQS += tree.getHaplotypeCounts((QuasiSpeciesNode)tree.getExternalNodes().get(i));
//        }
        nQS = tree.getTotalAttachmentCounts();
        int tipCount = tree.getLeafNodeCount();

        double[] dates = new double[tipCount + nQS];

        for (int i = 0; i < tipCount; i++) {
            dates[i] = tree.getNode(i).getHeight();
        }

        int nCount = tipCount;
        for (Node node : tree.getExternalNodes()){
            double QSTimesLength = tree.getHaplotypeCounts((QuasiSpeciesNode)node);
            for (int j=0; j<QSTimesLength; j++){
                dates[nCount+j] = node.getHeight();
            }
            nCount += QSTimesLength;
        }

        for (int k = 0; k < totalIntervals; k++) {


            for (int i = 0; i < tipCount + nQS; i++) {
                if (Math.abs((times[totalIntervals - 1] - times[k]) - dates[i]) < 1e-10) {
                    if (rho[k] == 0 && psi[k] == 0) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    if (rho[k] > 0) {
                        N[k] += 1;
                        if (i < tipCount) isRhoTip[i] = true;
                    }
                }
            }
        }
        // end change4
        return 0.;
    }

    /**
     * Collect all the times of parameter value changes and rho-sampling events
     */
    private void collectTimes() {

        timesSet.clear();

        if (isBDSIR()) {
            birthChanges = getSIRdimension() - 1;
        }

        getChangeTimes(birthRateChangeTimes,
                birthRateChangeTimesInput.get() != null && !isSeasonalBDSIR() ? birthRateChangeTimesInput.get() : intervalTimes.get(),
                birthChanges, birthRateTimesRelative, reverseTimeArrays[0]);

        getChangeTimes(deathRateChangeTimes,
                deathRateChangeTimesInput.get() != null ? deathRateChangeTimesInput.get() : intervalTimes.get(),
                deathChanges, deathRateTimesRelative, reverseTimeArrays[1]);

        getChangeTimes(samplingRateChangeTimes,
                samplingRateChangeTimesInput.get() != null ? samplingRateChangeTimesInput.get() : intervalTimes.get(),
                samplingChanges, samplingRateTimesRelative, reverseTimeArrays[2]);

        getChangeTimes(rhoSamplingChangeTimes,
                rhoSamplingTimes.get()!=null ? rhoSamplingTimes.get() : intervalTimes.get(),
                rhoChanges, false, reverseTimeArrays[3]);

//        if (SAModel) getChangeTimes(rChangeTimes,
//                removalProbabilityChangeTimesInput.get() != null ? removalProbabilityChangeTimesInput.get() : intervalTimes.get(),
//                rChanges, rTimesRelative, reverseTimeArrays[4]);

        for (Double time : birthRateChangeTimes) {
            timesSet.add(time);
        }
        for (Double time : deathRateChangeTimes) {
            timesSet.add(time);
        }

        for (Double time : samplingRateChangeTimes) {
            timesSet.add(time);
        }

        for (Double time : rhoSamplingChangeTimes) {
            timesSet.add(time);
        }

//        if (SAModel) {
//            for (Double time : rChangeTimes) {
//                timesSet.add(time);
//            }
//        }

        if (printTempResults) System.out.println("times = " + timesSet);

        times = timesSet.toArray(new Double[timesSet.size()]);
        totalIntervals = times.length;

        if (printTempResults) System.out.println("total intervals = " + totalIntervals);

    }
    // start change5
    protected Double updateRatesAndTimes(QuasiSpeciesTree tree) {
    // end change5

        collectTimes();

        t_root = tree.getRoot().getHeight();

        if (m_forceRateChange && timesSet.last() > (originIsRootEdge.get()? t_root+ origin.get().getValue() : origin.get().getValue())) {
            return Double.NEGATIVE_INFINITY;
        }

        if (transform)
            transformParameters();
        else {

            Double[] birthRates = birthRate.get().getValues();
            Double[] deathRates = deathRate.get().getValues();
            Double[] samplingRates = samplingRate.get().getValues();
            Double[] removalProbabilities = new Double[1];
//            if (SAModel) removalProbabilities = removalProbability.get().getValues();

            birth = new Double[totalIntervals];
            death = new Double[totalIntervals];
            psi = new Double[totalIntervals];
//            if (SAModel) r =  new Double[totalIntervals];

            birth[0] = birthRates[0];

            for (int i = 0; i < totalIntervals; i++) {
                if (!isBDSIR()) birth[i] = birthRates[index(times[i], birthRateChangeTimes)];
                death[i] = deathRates[index(times[i], deathRateChangeTimes)];
                psi[i] = samplingRates[index(times[i], samplingRateChangeTimes)];
//                if (SAModel) r[i] = removalProbabilities[index(times[i], rChangeTimes)];

                if (printTempResults) {
                    if (!isBDSIR()) System.out.println("birth[" + i + "]=" + birth[i]);
                    System.out.println("death[" + i + "]=" + death[i]);
                    System.out.println("psi[" + i + "]=" + psi[i]);
//                    if (SAModel) System.out.println("r[" + i + "]=" + r[i]);
                }
            }
        }




        if (m_rho.get() != null && (m_rho.get().getDimension()==1 ||  rhoSamplingTimes.get() != null)) {

            Double[] rhos = m_rho.get().getValues();
            rho = new Double[totalIntervals];

//            rho[totalIntervals-1]=rhos[rhos.length-1];
            for (int i = 0; i < totalIntervals; i++) {

                rho[i]= //rhoSamplingChangeTimes.contains(times[i]) ? rhos[rhoSamplingChangeTimes.indexOf(times[i])] : 0.;
                        rhoChanges>0?
                                rhoSamplingChangeTimes.contains(times[i]) ? rhos[rhoSamplingChangeTimes.indexOf(times[i])] : 0.
                                : rhos[0];
            }
        }

        return 0.;
    }


    /*    calculate and store Ai, Bi and p0        */
    // start change6
    public Double preCalculation(QuasiSpeciesTree tree) {
    // end change6
        if (!originIsRootEdge.get() && tree.getRoot().getHeight() >= origin.get().getValue()) {
            return Double.NEGATIVE_INFINITY;
        }

        // updateRatesAndTimes must be called before calls to index() below
        if (updateRatesAndTimes(tree) < 0) {
            return Double.NEGATIVE_INFINITY;
        }

        if (printTempResults) System.out.println("After update rates and times");

        if (m_rho.get() != null) {
            if (contempData) {
                rho = new Double[totalIntervals];
                Arrays.fill(rho, 0.);
                rho[totalIntervals-1] = m_rho.get().getValue();
            }

        } else {
            rho = new Double[totalIntervals];
            Arrays.fill(rho, 0.0);
        }

        if (m_rho.get() != null)
            if (computeN(tree) < 0)
                return Double.NEGATIVE_INFINITY;

        int intervalCount = times.length;

        Ai = new double[intervalCount];
        Bi = new double[intervalCount];
        p0 = new double[intervalCount];

        if (conditionOn == ConditionOn.RHO_SAMPLING) {
            Aihat = new double[intervalCount];
            Bihat = new double[intervalCount];
            p0hat = new double[intervalCount];
        }

        for (int i = 0; i < intervalCount; i++) {

            Ai[i] = Ai(birth[i], death[i], psi[i]);

            if (conditionOn == ConditionOn.RHO_SAMPLING) {
                Aihat[i] = Ai(birth[i], death[i], 0.0);
            }

            if (printTempResults) System.out.println("Ai[" + i + "] = " + Ai[i] + " " + Math.log(Ai[i]));
        }

        if (printTempResults) {
            System.out.println("birth[m-1]=" + birth[totalIntervals - 1]);
            System.out.println("death[m-1]=" + death[totalIntervals - 1]);
            System.out.println("psi[m-1]=" + psi[totalIntervals - 1]);
            System.out.println("rho[m-1]=" + rho[totalIntervals - 1]);
            System.out.println("Ai[m-1]=" + Ai[totalIntervals - 1]);
        }

        Bi[totalIntervals - 1] = Bi(
                birth[totalIntervals - 1],
                death[totalIntervals - 1],
                psi[totalIntervals - 1],
                rho[totalIntervals - 1],
                Ai[totalIntervals - 1], 1.);  //  (p0[m-1] = 1)

        if (conditionOn == ConditionOn.RHO_SAMPLING) {
            Bihat[totalIntervals - 1] = Bi(
                    birth[totalIntervals - 1],
                    death[totalIntervals - 1],
                    0.0,
                    rho[totalIntervals - 1],
                    Aihat[totalIntervals - 1], 1.);  //  (p0[m-1] = 1)
        }

        if (printTempResults)
            System.out.println("Bi[m-1] = " + Bi[totalIntervals - 1] + " " + Math.log(Bi[totalIntervals - 1]));
        for (int i = totalIntervals - 2; i >= 0; i--) {

            p0[i + 1] = p0(birth[i + 1], death[i + 1], psi[i + 1], Ai[i + 1], Bi[i + 1], times[i + 1], times[i]);
            if (Math.abs(p0[i + 1] - 1) < 1e-10) {
                return Double.NEGATIVE_INFINITY;
            }
            if (conditionOn == ConditionOn.RHO_SAMPLING) {
                p0hat[i + 1] = p0(birth[i + 1], death[i + 1], 0.0, Aihat[i + 1], Bihat[i + 1], times[i + 1], times[i]);
                if (Math.abs(p0hat[i + 1] - 1) < 1e-10) {
                    return Double.NEGATIVE_INFINITY;
                }
            }
            if (printTempResults) System.out.println("p0[" + (i + 1) + "] = " + p0[i + 1]);

            Bi[i] = Bi(birth[i], death[i], psi[i], rho[i], Ai[i], p0[i + 1]);
            if (conditionOn == ConditionOn.RHO_SAMPLING) {
                Bihat[i] = Bi(birth[i], death[i], 0.0, rho[i], Aihat[i], p0hat[i + 1]);
            }

            if (printTempResults) System.out.println("Bi[" + i + "] = " + Bi[i] + " " + Math.log(Bi[i]));
        }

        if (printTempResults) {
            System.out.println("g(0, x0, 0):" + g(0, times[0], 0));
            System.out.println("g(index(1),times[index(1)],1.) :" + g(index(1), times[index(1)], 1.));
            System.out.println("g(index(2),times[index(2)],2.) :" + g(index(2), times[index(2)], 2));
            System.out.println("g(index(4),times[index(4)],4.):" + g(index(4), times[index(4)], 4));
        }

        return 0.;
    }

    public double Ai(double b, double g, double psi) {

        return Math.sqrt((b - g - psi) * (b - g - psi) + 4 * b * psi);
    }

    public double Bi(double b, double g, double psi, double r, double A, double p0) {

        return ((1 - 2 * p0 * (1 - r)) * b + g + psi) / A;
    }

    public double p0(int index, double t, double ti) {

        return p0(birth[index], death[index], psi[index], Ai[index], Bi[index], t, ti);
    }

    public double p0(double b, double g, double psi, double A, double B, double ti, double t) {

        if (printTempResults)
            System.out.println("in p0: b = " + b + "; g = " + g + "; psi = " + psi + "; A = " + A + " ; B = " + B + "; ti = " + ti + "; t = " + t);
//        return ((b + g + psi - A *((Math.exp(A*(ti - t))*(1+B)-(1-B)))/(Math.exp(A*(ti - t))*(1+B)+(1-B)) ) / (2*b));
        // formula from manuscript slightly rearranged for numerical stability
        return ((b + g + psi - A * ((1 + B) - (1 - B) * (Math.exp(A * (t - ti)))) / ((1 + B) + Math.exp(A * (t - ti)) * (1 - B))) / (2 * b));

    }

    public double p0hat(int index, double t, double ti) {

        return p0(birth[index], death[index], 0.0, Aihat[index], Bihat[index], t, ti);
    }


    public double g(int index, double ti, double t) {

//        return (Math.exp(Ai[index]*(ti - t))) / (0.25*Math.pow((Math.exp(Ai[index]*(ti - t))*(1+Bi[index])+(1-Bi[index])),2));
        // formula from manuscript slightly rearranged for numerical stability
        return (4 * Math.exp(Ai[index] * (t - ti))) / (Math.exp(Ai[index] * (t - ti)) * (1 - Bi[index]) + (1 + Bi[index])) / (Math.exp(Ai[index] * (t - ti)) * (1 - Bi[index]) + (1 + Bi[index]));
    }

    /**
     * @param t the time in question
     * @return the index of the given time in the list of times, or if the time is not in the list, the index of the
     *         next smallest time
     */
    public int index(double t, List<Double> times) {

        int epoch = Collections.binarySearch(times, t);

        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        return epoch;
    }


    /**
     * @param t the time in question
     * @return the index of the given time in the times array, or if the time is not in the array the index of the time
     *         next smallest
     */
    public int index(double t) {

        if (t >= times[totalIntervals - 1])
            return totalIntervals - 1;

        int epoch = Arrays.binarySearch(times, t);

        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        return epoch;
    }


    /**
     * @param time the time
     * @param tree the tree
     * @return the number of lineages that exist at the given time in the given tree.
     */

    // VB this function needed to be changed (@ 30.06.2015)
    // start change7
    public int lineageCountAtTime(double time, QuasiSpeciesTree tree) {
    // end change7
        int count = 1;
        int tipCount = tree.getLeafNodeCount();
        for (int i = tipCount; i < tipCount + tree.getInternalNodeCount(); i++) {
            if (tree.getNode(i).getHeight() > time) count += 1;
        }
        for (int i = 0; i < tipCount; i++) {
            if (tree.getNode(i).getHeight() >= time) count -= 1;
        }
        // start change8
        for (Node node : tree.getExternalNodes()){
            Double[] QSTimesTemp = tree.getAttachmentTimesList((QuasiSpeciesNode) node);
            if (QSTimesTemp[0] > time){
                for (int j=1; j<=QSTimesTemp.length; j++){
                    if (QSTimesTemp[j] > time) count += 1;
                }
            }
        }
        // end change8
        return count;
    }

    /**
     * @param time the time
     * @param tree the tree
     * @param k count the number of sampled ancestors at the given time
     * @return the number of lineages that exist at the given time in the given tree.
     */

    // VB this function needed to be changed (@ 30.06.2015)
    // start change9
    public int lineageCountAtTime(double time, QuasiSpeciesTree tree, int[] k) {
    // end change9
        int count = 1;
        k[0]=0;
        int tipCount = tree.getLeafNodeCount();
        for (int i = tipCount; i < tipCount + tree.getInternalNodeCount(); i++) {
            if (tree.getNode(i).getHeight() >= time) count += 1;

        }
        for (int i = 0; i < tipCount; i++) {
            if (tree.getNode(i).getHeight() > time) count -= 1;
            if (Math.abs(tree.getNode(i).getHeight() - time) < 1e-10) {
                count -= 1;
                if (tree.getNode(i).isDirectAncestor()) {
                    count -= 1;
                    k[0]++;
                }

            }
        }
        // start change10
        for (Node node : tree.getExternalNodes()){
            Double[] QSTimesTemp = tree.getAttachmentTimesList((QuasiSpeciesNode) node);
            if (QSTimesTemp[0] > time){
                for (int j=1; j<=QSTimesTemp.length; j++){
                    if (QSTimesTemp[j] > time) count += 1;
                }
            }
        }
        // end change10
        return count;
    }

    protected void transformParameters() {

        Double[] R = R0.get().getValues(); // if SAModel: R0 = lambda/delta
        Double[] b = becomeUninfectiousRate.get().getValues(); // delta = mu + psi*r
        Double[] p = samplingProportion.get().getValues(); // if SAModel: s = psi/(mu+psi*r)
        Double[] removalProbabilities = new Double[1];
//        if (SAModel) removalProbabilities = removalProbability.get().getValues();

        birth = new Double[totalIntervals];
        death = new Double[totalIntervals];
        psi = new Double[totalIntervals];
//        if (SAModel) r =  new Double[totalIntervals];

        if (isBDSIR()) birth[0] = R[0] * b[0]; // the rest will be done in BDSIR class

        for (int i = 0; i < totalIntervals; i++) {
//            if (!SAModel) {
                if (!isBDSIR()) birth[i] = R[birthChanges > 0 ? index(times[i], birthRateChangeTimes) : 0] * b[deathChanges > 0 ? index(times[i], deathRateChangeTimes) : 0];
                psi[i] = p[samplingChanges > 0 ? index(times[i], samplingRateChangeTimes) : 0] * b[deathChanges > 0 ? index(times[i], deathRateChangeTimes) : 0];
                death[i] = b[deathChanges > 0 ? index(times[i], deathRateChangeTimes) : 0] - psi[i];
//            } else {
//                birth[i] = R[birthChanges > 0 ? index(times[i], birthRateChangeTimes) : 0] * b[deathChanges > 0 ? index(times[i], deathRateChangeTimes) : 0];
//                psi[i] = p[samplingChanges > 0 ? index(times[i], samplingRateChangeTimes) : 0] * b[deathChanges > 0 ? index(times[i], deathRateChangeTimes) : 0];
//                r[i] = removalProbabilities[rChanges > 0 ? index(times[i], rChangeTimes) : 0];
//                death[i] = b[deathChanges > 0 ? index(times[i], deathRateChangeTimes) : 0] - psi[i]*r[i];
//
//
//            }


        }
    }

    // start change11
    public double calculateTreeLogLikelihood(QuasiSpeciesTree tree, Set<Taxon> exclude) {
    // end change11
        if (exclude.size() == 0) return calculateTreeLogLikelihood(tree);
        throw new RuntimeException("Not implemented!");
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

//    @Override
    public boolean canHandleTipDates() {
        return (m_rho.get() == null);
    }


    public Boolean isBDSIR() {
        return false;
    }

    public Boolean isSeasonalBDSIR() {
        return false;
    }

    public int getSIRdimension() {
        throw new RuntimeException("This is not an SIR");
    }



    // END3


// calculateTreeLogLikelihood also adapted from Denise's code
    @Override
//    public double calculateTreeLogLikelihood(TreeInterface tree, TraitSet QuasiSpeciesAttachmentTimes) {
    public double calculateTreeLogLikelihood(QuasiSpeciesTree tree) {
//first trial        System.out.println("likelihood being calculated ... ");
//first trial        System.out.println(treeInput.get());
//first trial        System.out.println(tree);
//first trial
//first trial        return 0;

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

        double x0 = 0;
        int index = 0;

        if (times[index] < 0.)
            index = index(0.);

        double temp=0;

        switch (conditionOn) {
            case NONE:
                temp = Math.log(g(index, times[index], x0));
                break;
            case SURVIVAL:
                temp = p0(index, times[index], x0);
                if (temp == 1)
                    return Double.NEGATIVE_INFINITY;
                temp = Math.log(g(index, times[index], x0) / (1 - temp));
                break;
            case RHO_SAMPLING:
                temp = p0hat(index, times[index], x0);
                if (temp == 1)
                    return Double.NEGATIVE_INFINITY;
                temp = Math.log(g(index, times[index], x0) / (1 - temp));
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
        int[] startBranchCountArray = tree.getStartBranchCounts();
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {

            double x = times[totalIntervals - 1] - tree.getNode(nTips + i).getHeight();
            index = index(x);
            if (!(tree.getNode(nTips + i)).isFake()) {
                temp = Math.log(birth[index] * g(index, times[index], x));
                logP += temp;
                if (printTempResults) System.out.println("1st pwd" +
                        " = " + temp + "; interval = " + i);
                if (Double.isInfinite(logP))
                    return logP;
                // term for the Quasi-Species tree likelihood calculation counting possible QS start branches
                int gamma=startBranchCountArray[i];
                temp = gamma;
                logP += Math.log(temp);
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
            int nQSTemp = (int) tree.getHaplotypeCounts((QuasiSpeciesNode) node);
            for (int j = 1; j <= nQSTemp; j++ ){
                Double[] QSTimesTemp = tree.getAttachmentTimesList((QuasiSpeciesNode) node);
                double x = times[totalIntervals - 1] - QSTimesTemp[j];
                index = index(x);
                // term for the Quasi-Species tree likelihood calculation counting possible attachment branches
                // THIS IS A CONSTANT FOR ANY TREE TOPOLOGY, SO OMIT THIS TERM FOR NUMERICAL STABILITY AND CALCULATION SPEED
                //int gammaj = tree.countPossibleAttachmentBranches((QuasiSpeciesNode) node, j, tree.getAttachmentTimesList((QuasiSpeciesNode) node)[j]);
                //temp = Math.log(gammaj * birth[index] * g(index, times[index], x));
                //      times 2 for left OR right -- constant factor also not taken into account
                // temp = Math.log(2 * gammaj * birth[index] * g(index, times[index], x));
                temp = Math.log(birth[index] * g(index, times[index], x));
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
            if (!isRhoTip[i] || m_rho.get() == null || tree.getHaplotypeCounts((QuasiSpeciesNode) tree.getNode(i))<1) {
                double y = times[totalIntervals - 1] - tree.getNode(i).getHeight();
                index = index(y);
// testing
//                System.out.println("2nd factor included");
                if (!(tree.getNode(i)).isDirectAncestor()) {
//                    if (!SAModel) {
                        temp = Math.log(psi[index]) - Math.log(g(index, times[index], y));
//                    } else {
//                        temp = Math.log(psi[index] * (r[index] + (1 - r[index]) * p0(index, times[index], y))) - Math.log(g(index, times[index], y));
//                    }
                    logP += temp;
                    if (printTempResults) System.out.println("2nd PI = " + temp);
                    if (psi[index] == 0 || Double.isInfinite(logP))
                        return logP;
                } else {
                    if (r[index] != 1) {
                        logP += Math.log((1 - r[index])*psi[index]);
                        if (Double.isInfinite(logP)) {
                            return logP;
                        }
                    } else {
                        //throw new Exception("There is a sampled ancestor in the tree while r parameter is 1");
                        System.out.println("There is a sampled ancestor in the tree while r parameter is 1");
                        System.exit(0);
                    }
                }
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
                temp = n[j] * (Math.log(g(j, times[j], time)) + Math.log(1 - rho[j-1]));
                logP += temp;
// testing
//                System.out.println("3rd factor (nj loop) = " + temp + "; interval = " + j + "; n[j] = " + n[j]);
                if (printTempResults)
                    System.out.println("3rd factor (nj loop) = " + temp + "; interval = " + j + "; n[j] = " + n[j]);//+ "; Math.log(g(j, times[j], time)) = " + Math.log(g(j, times[j], time)));
                if (Double.isInfinite(logP))
                    return logP;

            }

//            if (SAModel && j>0 && N != null) { // term for sampled leaves and two-degree nodes at time t_i
//                logP += k[0] * (Math.log(g(j, times[j], time)) + Math.log(1-r[j])) + //here g(j,..) corresponds to q_{i+1}, r[j] to r_{i+1},
//                        (N[j-1]-k[0])*(Math.log(r[j]+ (1-r[j])*p0(j, times[j], time))); //N[j-1] to N_i, k[0] to K_i,and thus N[j-1]-k[0] to M_i
//                if (Double.isInfinite(logP)) {
//                    return logP;
//                }
//            }

            if (rho[j] > 0 && N[j] > 0) { // term for contemporaneous sampling
                temp = N[j] * Math.log(rho[j]);
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
//            int internalNodeCount = tree.getLeafNodeCount() - ((QuasiSpeciesTree)tree).getDirectAncestorNodeCount()- 1;
//            logP +=  Math.log(2)*internalNodeCount;
//        }

        return logP;

    }

}
