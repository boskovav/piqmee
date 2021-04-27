package piqmee.likelihood;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.BeerLikelihoodCore;
import beast.evolution.likelihood.BeerLikelihoodCore4;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.likelihood.LikelihoodCore;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import piqmee.tree.QuasiSpeciesNode;
import piqmee.tree.QuasiSpeciesTree;

@Description("Calculates the probability of sequence data on a beast.piqmee.tree " +
        "given a site and substitution model using a variant of the 'peeling algorithm'. ")
public class QuasiSpeciesTreeLikelihood3 extends GenericTreeLikelihood {
	
	// TODO: this input is ignored, but required by junit test. Needs fixing.
	final public Input<Boolean> useAmbiguities = new Input<>("useAmbiguities", "flag to indicate that sites containing ambiguous states should be handled instead of ignored (the default)", false);
	
	
    public static enum Scaling {none, always, _default};
    final public Input<Scaling> scalingInput = new Input<>("scaling", "type of scaling to use, one of " + Arrays.toString(Scaling.values()) + ". If not specified, the -beagle_scaling flag is used.", Scaling._default, Scaling.values());


    /**
     * calculation engine *
     */
    protected LikelihoodCore likelihoodCore;
    public LikelihoodCore getLikelihoodCore() {return likelihoodCore;}
    protected QuasiSpeciesBeagleTreeLikelihood3 beagleLikelihood;

    /**
     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    protected SubstitutionModel substitutionModel;
    public SubstitutionModel getSubstitutionModel() {return substitutionModel;}
    
    protected SiteModel.Base siteModel;
    protected BranchRateModel.Base branchRateModel;
    int nodeCount;
    int [][] states;
    Alignment alignment;

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] branchLengths;
    protected double[] storedBranchLengths;
    
    /**
     * memory allocation for likelihoods for each of the patterns *
     */
    protected double[] patternLogLikelihoods;
    /**
     * memory allocation for the root partials *
     */
    protected double[] m_fRootPartials, rawRootPartials;
    
    /** **/
    protected double[] accumulatedLogLeafScaleFactors;
    protected double[] storedAccumulatedLogLeafScaleFactors;
    protected double[] accumulatedLeafScaleFactors;
    protected double [] scale, storedScale;

    protected int[] leafIndex;
    protected int[] storedLeafIndex;
    protected double[][][] leafLogScaleFactors;
    protected double [] logProbabilities;
    
    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    protected double[] probabilities;
    protected double[] rates;
    protected double[] storedRates;
    protected double[] tmpevectimesevals;

    protected int nStates;
    
    protected int matrixSize;

    protected boolean useScaleFactors = false;

    /**
     * flag to indicate ascertainment correction should be applied *
     */
    protected boolean useAscertainedSitePatterns = false;

    /**
     * dealing with proportion of site being invariant *
     */
    double proportionInvariant = 0;
    public double getProportionInvariant() {return proportionInvariant;}
	public void setProportionInvariant(double proportionInvariant) {this.proportionInvariant = proportionInvariant;
	}
	List<Integer> constantPattern = null;
	public List<Integer> getConstantPattern() {return constantPattern;}
    public void setConstantPattern(List<Integer> constantPattern) {this.constantPattern = constantPattern;}

    // protected double tolerance = 0;
    Node toyNode = new Node();

    @Override
    public void initAndValidate() {
        // sanity check: alignment should have same #taxa as tree
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            alignment = subset(dataInput,treeInput.get().getExternalNodes());
        }
        else
            alignment = dataInput.get();
        beagleLikelihood = null;
    	if (!scalingInput.get().equals(Scaling.always)) {
    		beagleLikelihood = new QuasiSpeciesBeagleTreeLikelihood3();
        	try {
        		beagleLikelihood.initByName(
        				"data", dataInput.get(), "tree", treeInput.get(), "siteModel", siteModelInput.get(),
        				"branchRateModel", branchRateModelInput.get(), "useAmbiguities", useAmbiguities.get(),
        				"scaling", Scaling.none); // QuasiSpeciesBeagleTreeLikelihood3 scaling does not work properly
        				// so we go back to this Java implementation if scaling is required
            	if (beagleLikelihood.getBeagle() == null) {
            		beagleLikelihood = null;
            	}
	        } catch (Exception e) {
	            // No Beagle instance was found, so we use the good old java likelihood core
	        	beagleLikelihood = null;
	    	}
		}
        // tolerance = toleranceInput.get();
        nodeCount = treeInput.get().getNodeCount();
        int leafNodeCount = treeInput.get().getLeafNodeCount();

        int nodeCount = treeInput.get().getNodeCount();
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
        	throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        siteModel = (SiteModel.Base) siteModelInput.get();
        siteModel.setDataType(alignment.getDataType());
        substitutionModel = siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        branchLengths = new double[nodeCount + leafNodeCount];
        storedBranchLengths = new double[nodeCount + leafNodeCount];

        nStates = alignment.getMaxStateCount();
        int patterns = alignment.getPatternCount();
        likelihoodCore = createLikelihoodCore(nStates);

        String className = getClass().getSimpleName();

        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));
        // print startup messages via Log.print*

        proportionInvariant = siteModel.getProportionInvariant();
        siteModel.setPropInvariantIsCategory(false);
        if (proportionInvariant > 0) {
            calcConstantPatternIndices(patterns, nStates);
        }

        initCore();

        patternLogLikelihoods = new double[patterns];
        m_fRootPartials = new double[patterns * nStates];
        rawRootPartials = new double[patterns * nStates * siteModel.getCategoryCount()];
        matrixSize = (nStates + 1) * (nStates + 1);
        probabilities = new double[(nStates + 1) * (nStates + 1)];
        Arrays.fill(probabilities, 1.0);

        if (alignment.isAscertained) {
            useAscertainedSitePatterns = true;
        }

        leafIndex = new int[leafNodeCount];        
        storedLeafIndex = new int[leafNodeCount];        
        leafLogScaleFactors = new double[2][leafNodeCount][patterns * siteModel.getCategoryCount()];
        accumulatedLogLeafScaleFactors = new double[patterns * siteModel.getCategoryCount()];
        storedAccumulatedLogLeafScaleFactors = new double[patterns * siteModel.getCategoryCount()];
        accumulatedLeafScaleFactors = new double[accumulatedLogLeafScaleFactors.length];
    	scale = new double[patterns];
    	storedScale = new double[patterns];

        logProbabilities = new double[nStates * siteModel.getCategoryCount()];

        rates = new double[nStates];
        storedRates = new double [nStates];
        tmpevectimesevals = new double[nStates * nStates];
        getNoChangeRates(rates);

    	useScaleFactors = scalingInput.get().equals(Scaling.always);
    }

    
    /**
     * Function to subset full alignment to only sequences carried at tips in the
     * unique sequence tree
     *
     * @param data  full alignment
     * @param leafs leaf nodes of the unique-sequence tree
     * @return subsetted alignment
     */
    public Alignment subset(Input<Alignment> data, List<Node> leafs){
        Alignment fullData = data.get();
        int tipCount = leafs.size();
        ArrayList sequences = new ArrayList(tipCount);

        Alignment subsetData;
        // since taxonSet is not always ordered according to the sequences, we need to reorder the alignment
        // this filtered alignment part is just for beauti to not throw its toys
        if (fullData.sequenceInput.get().size() == 0 && fullData instanceof FilteredAlignment) {
            // sort the alignment
            Alignment fullsortedAlignment =
                    new Alignment(((FilteredAlignment) fullData).alignmentInput.get().sequenceInput.get(), fullData.dataTypeInput.get());
            // select sequences for subset
            for (int i = 0; i < tipCount; i++){
                sequences.add(leafs.get(i).getNr(),fullsortedAlignment.sequenceInput.get().get(fullsortedAlignment.getTaxonIndex(leafs.get(i).getID())));
            }
            // make a new filtered alignment with this subset
            Alignment toyAlignment = new Alignment(sequences, fullData.dataTypeInput.get());
            FilteredAlignment subsetDataFiltered = new FilteredAlignment();
            if (((FilteredAlignment) fullData).constantSiteWeightsInput.get() != null) {
                subsetDataFiltered.initByName(
                        "data",        toyAlignment,
                        "filter",               ((FilteredAlignment) fullData).filterInput.get(),
                        "constantSiteWeights",  ((FilteredAlignment) fullData).constantSiteWeightsInput.get());
            } else {
                subsetDataFiltered.initByName(
                        "data",        toyAlignment,
                        "filter",               ((FilteredAlignment) fullData).filterInput.get());
            }
            // set subsetData to this new Filtered alignment
            subsetData = subsetDataFiltered;
        } else {
            // sort the alignment
            Alignment fullsortedAlignment = new Alignment(fullData.sequenceInput.get(), fullData.dataTypeInput.get());
            // select sequences for subset
            for (int i = 0; i < tipCount; i++){
                sequences.add(leafs.get(i).getNr(),fullsortedAlignment.sequenceInput.get().get(fullsortedAlignment.getTaxonIndex(leafs.get(i).getID())));
            }
            // make a new alignment with this subset
            subsetData = new Alignment(sequences, fullData.dataTypeInput.get());
        }

        return subsetData;
    }
    
    protected LikelihoodCore createLikelihoodCore(int nStates) {
        if (nStates == 4) {
            return new BeerLikelihoodCore4();
        } else {
            return new BeerLikelihoodCore(nStates);
        }
    }
    

    /**
     * Determine indices of m_fRootProbabilities that need to be updates
     * // due to sites being invariant. If none of the sites are invariant,
     * // the 'site invariant' category does not contribute anything to the
     * // root probability. If the site IS invariant for a certain character,
     * // taking ambiguities in account, there is a contribution of 1 from
     * // the 'site invariant' category.
     */
    protected void calcConstantPatternIndices(final int patterns, final int nStates) {
        constantPattern = new ArrayList<>();
        for (int i = 0; i < patterns; i++) {
            final int[] pattern = alignment.getPattern(i);
            final boolean[] isInvariant = new boolean[nStates];
            Arrays.fill(isInvariant, true);
            for (final int state : pattern) {
                final boolean[] isStateSet = alignment.getStateSet(state);
//                if (m_useAmbiguities.get() || !alignment.getDataType().isAmbiguousCode(state)) {
//                    for (int k = 0; k < nStates; k++) {
//                        isInvariant[k] &= isStateSet[k];
//                    }
//                }
            }
            for (int k = 0; k < nStates; k++) {
                if (isInvariant[k]) {
                    constantPattern.add(i * nStates + k);
                }
            }
        }
    }

    protected void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        likelihoodCore.initialize(
                nodeCount,
                alignment.getPatternCount(),
                siteModel.getCategoryCount(),
                true, false // m_useAmbiguities.get()
        );

        final int extNodeCount = nodeCount / 2 + 1;
        final int intNodeCount = nodeCount / 2;

        states = new int [extNodeCount][alignment.getPatternCount()];
//        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
//            setPartials(treeInput.get().getRoot(), alignment.getPatternCount());
//        } else {
            setStates(treeInput.get().getRoot(), alignment.getPatternCount());
//        }
        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i);
        }
    }

    /**
     * This method samples the sequences based on the tree and site model.
     */
    @Override
	public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Can't sample a fixed alignment!");
    }

    /**
     * set leaf states in likelihood core *
     */
    protected void setStates(Node node, int patternCount) {
        if (node.isLeaf()) {
            int i;
            int[] states = this.states[node.getNr()];
            int taxonIndex = getTaxonIndex(node.getID(), alignment);
            for (i = 0; i < patternCount; i++) {
                int code = alignment.getPattern(taxonIndex, i);
                int[] statesForCode = alignment.getDataType().getStatesForCode(code);
                if (statesForCode.length==1)
                    states[i] = statesForCode[0];
                else
                    states[i] = code; // Causes ambiguous states to be ignored.
            }
            likelihoodCore.setNodeStates(node.getNr(), states);

        } else {
            setStates(node.getLeft(), patternCount);
            setStates(node.getRight(), patternCount);
        }
    }

    /**
     *
     * @param taxon the taxon name as a string
     * @param data the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     *         or -1 if the taxon is not in the alignment.
     */
    private int getTaxonIndex(String taxon, Alignment data) {
        int taxonIndex = data.getTaxonIndex(taxon);
        if (taxonIndex == -1) {
        	if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
            }
            if (taxonIndex == -1) {
            	throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
            }
        }
        return taxonIndex;
	}

	/**
     * set leaf partials in likelihood core *
     */
    protected void setPartials(Node node, int patternCount) {
        if (node.isLeaf()) {
            double[] partials = new double[patternCount * nStates];
            int k = 0;
            int taxonIndex = getTaxonIndex(node.getID(), alignment);
            for (int patternIndex_ = 0; patternIndex_ < patternCount; patternIndex_++) {                
                double[] tipLikelihoods = alignment.getTipLikelihoods(taxonIndex,patternIndex_);
                if (tipLikelihoods != null) {
                	for (int state = 0; state < nStates; state++) {
                		partials[k++] = tipLikelihoods[state];
                	}
                }
                else {
                	int nStates = alignment.getPattern(taxonIndex, patternIndex_);
	                boolean[] stateSet = alignment.getStateSet(nStates);
	                for (int state = 0; state < nStates; state++) {
	                	 partials[k++] = (stateSet[state] ? 1.0 : 0.0);                
	                }
                }
            }
            likelihoodCore.setNodePartials(node.getNr(), partials);

        } else {
            setPartials(node.getLeft(), patternCount);
            setPartials(node.getRight(), patternCount);
        }
    }

    // for testing
    public double[] getRootPartials() {
        return m_fRootPartials.clone();
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    double scaleFactor = 1.0;
    double storedScaleFactor = 1.0;
    int m_nScale = 0;
    int X = 100;
    boolean startScaling = false;

    @Override
    public double calculateLogP() {
        if (beagleLikelihood != null && scaleFactor == 1.0) {
        	// only use BEAGLE implementation if not scaling
            logP = beagleLikelihood.calculateLogP();
            if (Double.isFinite(logP) || scalingInput.get().equals(Scaling.none)) {
            	return logP;
            }
            hasDirt = Tree.IS_FILTHY;
            getNoChangeRates(rates);
        }
        
        if (siteModel.isDirtyCalculation()) {
            getNoChangeRates(rates);
        }
        
        final TreeInterface tree = treeInput.get();

        try {
        	if (traverse((QuasiSpeciesNode)tree.getRoot()) != Tree.IS_CLEAN)
        		calcLogP();
        }
        catch (ArithmeticException e) {
        	return Double.NEGATIVE_INFINITY;
        }
        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
//            System.err.println("Switch off scaling");
//            m_likelihoodCore.setUseScaling(1.0);
//            m_likelihoodCore.unstore();
//            m_nHasDirt = Tree.IS_FILTHY;
//            X *= 2;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
        } else if (logP == Double.NEGATIVE_INFINITY && scaleFactor < 10 && !scalingInput.get().equals(Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
        	startScaling = true;
        	useScaleFactors = true;
            m_nScale = 0;
            scaleFactor *= 1.01;
            Log.warning.println("Turning on scaling to prevent numeric instability " + scaleFactor);
            likelihoodCore.setUseScaling(scaleFactor);
            likelihoodCore.unstore();
            System.arraycopy(storedLeafIndex, 0, leafIndex, 0, leafIndex.length);
            hasDirt = Tree.IS_FILTHY;
            traverse((QuasiSpeciesNode)tree.getRoot());
            calcLogP();
            return logP;
        }
//        System.err.println("QUASI3 :" + logP);        
        return logP;
    }

    void calcLogP() {
    	
    	accumulateLogLeafScale();
    	
    	Node root = treeInput.get().getRoot();
        final double[] proportions = siteModel.getCategoryProportions(root);
        
        likelihoodCore.getNodePartials(root.getNr(), rawRootPartials);

        integratePartials(rawRootPartials, proportions, m_fRootPartials, patternLogLikelihoods.length, proportions.length);

        if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
            proportionInvariant = siteModel.getProportionInvariant();
            // some portion of sites is invariant, so adjust root partials for this
            for (final int i : constantPattern) {
                m_fRootPartials[i] += proportionInvariant;
            }
        }

        double[] rootFrequencies = substitutionModel.getFrequencies();
        calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods, patternLogLikelihoods.length);

    	
        logP = 0.0;
        if (useAscertainedSitePatterns) {
            final double ascertainmentCorrection = alignment.getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < alignment.getPatternCount(); i++) {
                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * alignment.getPatternWeight(i);
            }
        } else {
            for (int i = 0; i < alignment.getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * alignment.getPatternWeight(i);
//            	System.out.println(patternLogLikelihoods[i] + " * "+ alignment.getPatternWeight(i) + " " + logP);
            }
        }
    }

    
    /**
     * Integrates partials across categories.
     *
     * @param inPartials  the array of partials to be integrated
     * @param proportions the proportions of sites in each category
     * @param outPartials an array into which the partials will go
     */
	protected void integratePartials(double[] inPartials, double[] proportions, double[] outPartials, 
			final int nrOfPatterns, final int nrOfMatrices) {
		
		int n = accumulatedLogLeafScaleFactors.length;    	
		boolean hasZero = false; // deals with underflows of Math.exp
		if (!useScaleFactors) {
			for (int k = 0; k < n; k++) {
				accumulatedLeafScaleFactors[k] = Math.exp(accumulatedLogLeafScaleFactors[k]);
				if (accumulatedLeafScaleFactors[k] == 0) {
					hasZero = true;
					break;
				}
			}
		}
		if (hasZero || useScaleFactors) {
			for (int k = 0; k < nrOfPatterns; k++) {
				double max = accumulatedLogLeafScaleFactors[k];
				int r = k;
		        for (int l = 1; l < nrOfMatrices; l++) {
		        	r += nrOfPatterns;
		        	max = Math.max(max, accumulatedLogLeafScaleFactors[r]);
		        }
		        scale[k] = max;
			}
			int w = 0;
	        for (int l = 0; l < nrOfMatrices; l++) {
	        	for (int k = 0; k < nrOfPatterns; k++) {
	        		accumulatedLeafScaleFactors[w] = Math.exp(accumulatedLogLeafScaleFactors[w] - scale[k]);
	        		w++;
	        	}
			}
		}
        
        
        int u = 0;
        int v = 0;
        int w = 0;
        for (int k = 0; k < nrOfPatterns; k++) {

            for (int i = 0; i < nStates; i++) {
                outPartials[u] = inPartials[v] * proportions[0]  * accumulatedLeafScaleFactors[w];
//                System.out.println(outPartials[u] +"="+ inPartials[v] +" * "+ proportions[0]  +" * "+ accumulatedLeafScaleFactors[w]);
                u++;
                v++;
            }
            w++;
        }


        for (int l = 1; l < nrOfMatrices; l++) {
            u = 0;

            for (int k = 0; k < nrOfPatterns; k++) {

                for (int i = 0; i < nStates; i++) {
                    outPartials[u] += inPartials[v] * proportions[l] * accumulatedLeafScaleFactors[w];
//                    System.out.println(outPartials[u] +"="+ inPartials[v] +" * "+ proportions[l]  +" * "+ accumulatedLeafScaleFactors[w]);
                    u++;
                    v++;
                }
                w++;
            }
//            System.out.println();
        }
//        System.out.println();
    }
    
	public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods,
			final int nrOfPatterns) {
        int v = 0;
        for (int k = 0; k < nrOfPatterns; k++) {

            double sum = 0.0;
            for (int i = 0; i < nStates; i++) {

                sum += frequencies[i] * partials[v];
                v++;
            }
            outLogLikelihoods[k] = Math.log(sum) + likelihoodCore.getLogScalingFactor(k) + scale[k];
//            System.out.println(outLogLikelihoods[k] +"=" + Math.log(sum) + " " + likelihoodCore.getLogScalingFactor(k) + " " +scale[k]);
        }
    }

	
	/**
     * Get the transition rates for the QS branches,
     * i.e. the rate of no change, i.e. the diagonal entries in the rate matrix
     *
     * @param rates the matrix where the new rates should be stored
     *
     */
    public void getNoChangeRates(double[] rates) {
        EigenDecomposition eigenDecomposition = substitutionModel.getEigenDecomposition(null);
        double[] evec = eigenDecomposition.getEigenVectors();
        double[] eval = eigenDecomposition.getEigenValues();
        double[] ievc = eigenDecomposition.getInverseEigenVectors();
        for (int i = 0; i < nStates ; i++) {
            for (int j = 0; j < nStates; j++) {
                tmpevectimesevals[i * nStates + j] = evec[i * nStates + j] * eval[j];
            }
        }
        for (int i = 0; i < nStates; i++) {
            rates[i] = 0;
            for (int j = 0; j < nStates; j++) {
                rates[i] += tmpevectimesevals[i * nStates + j] * ievc[j * nStates + i];
            }
            //            System.out.println(rates[i]);
        }
    }
    
    /* Assumes there IS a branch rate model as opposed to traverse() */
    int traverse(final QuasiSpeciesNode node){

        QuasiSpeciesTree tree = (QuasiSpeciesTree) treeInput.get();

        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = ((QuasiSpeciesNode)node).getLengthWithoutHaplo() * branchRate;
        
        if (node.isLeaf()) {
        	// TODO: VERIFY THAT WE HAVE THE RIGHT BRANCH TIME (AS LOGGED BY RELAXED CLOCK LOGGER)
            toyNode.setNr(nodeCount + nodeIndex);
            branchRate = branchRateModel.getRateForBranch(toyNode);
        	
            double totalBranchTime = ((QuasiSpeciesNode)node).getTotalBranchLengths() * branchRate;            
        	if  (update != Tree.IS_CLEAN  || totalBranchTime != branchLengths[nodeCount + nodeIndex]) {
	        	
        		branchLengths[nodeCount + nodeIndex] = totalBranchTime;
	            // likelihoodCore.setNodeMatrixForUpdate(nodeCount + nodeIndex);
	        	int k = 0;
	            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
	                final double jointBranchRate = siteModel.getRateForCategory(i, node);
	                // fill the transition probability matrix with move probabilities
	                // Arrays.fill(probabilities, 0);
	                for (int j = 0; j < nStates; j++) {
	                    logProbabilities[j + k] = totalBranchTime * jointBranchRate * rates[j];
	                    // probabilities[j + k] = Math.exp(totalBranchTime * jointBranchRate * rates[j]);
	                }
	                k += nStates;
	                // likelihoodCore.setNodeMatrix(nodeCount + nodeIndex, i, probabilities);
	            }
	            
	            // this sets node partials at top of leaf clade at partial[nodeIndex]
	            setLeafScaleForUpdate(nodeIndex);
	
	        	calculateLogLeafScale(nodeIndex, logProbabilities);
        	}
        }
        
        
        
        if (!node.isRoot() && (update != Tree.IS_CLEAN || Math.abs(branchTime - branchLengths[nodeIndex]) > 0)) {
        	branchLengths[nodeIndex] = branchTime;
            if (branchTime < 0.0) {
            	double branchRate2 = branchRateModel.getRateForBranch(node);
                throw new RuntimeException("Negative branch length: " + branchTime + " " + branchRate2);
            }

            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);

        	Node parent = node.getParent();
        	double firstBranchingTime = ((QuasiSpeciesNode)node).getFirstBranchingTime();
            for (int i = 0; i < siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), firstBranchingTime, jointBranchRate, probabilities);
                for (int j = 0; j < matrixSize; j++) {
                	probabilities[j] *= scaleFactor;
                }
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }

            update |= Tree.IS_DIRTY;
        }        


        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse((QuasiSpeciesNode) child1);

            final Node child2 = node.getRight();
            final int update2 = traverse((QuasiSpeciesNode) child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN || update != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int child1QS = ((QuasiSpeciesNode) child1).getContinuingHaploName();
                final int child1parentQS = ((QuasiSpeciesNode) child1.getParent()).getContinuingHaploName();
                final int childNum2 = child2.getNr();
                final int child2QS = ((QuasiSpeciesNode) child2).getContinuingHaploName();
                final int child2parentQS = ((QuasiSpeciesNode) child2.getParent()).getContinuingHaploName();
                if (child1parentQS != child2parentQS)
                    throw new IllegalStateException("In QuasiSpeciesTreeLikelihood - QS of parent of child 1 ne to QS of parent of child 2");

                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY)
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);

                if (siteModel.integrateAcrossCategories()) {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 632: Site categories not supported");
                }

                if (node.isRoot()) {
                    // ((QuasieSpeciesBeerLikelihoodCore)likelihoodCore).calculateOriginRootPartials(nodeIndex, child1parentQS, nodeCount, rootPartials);

                    // integrate over all possible site categories the sites can be in
//                    final double[] proportions = siteModel.getCategoryProportions(node);
//                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);xx
//
//                    if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
//                        proportionInvariant = siteModel.getProportionInvariant();
//                        // some portion of sites is invariant, so adjust root partials for this
//                        for (final int i : constantPattern) {
//                            m_fRootPartials[i] += proportionInvariant;
//                        }
//                    }
//
//                    double[] rootFrequencies = substitutionModel.getFrequencies();
//                    if (rootFrequenciesInput.get() != null) {
//                        rootFrequencies = rootFrequenciesInput.get().getFreqs();
//                    }
//                    likelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);
                }
            }
        }
        // if the tree has only one child
        else {
            if (tree.getLeafNodeCount()==1) {
            	throw new RuntimeException("Expected more than 1 leaf node in tree");
            }
        }
        return update;
    }



    protected void calculateLogLeafScale(int nodeIndex, double[] logProbabilities) {
    	double [] current = leafLogScaleFactors[leafIndex[nodeIndex]][nodeIndex];
    	int [] states = this.states[nodeIndex];
    	int patternCount = states.length;
    	int v = 0;
    	int w = 0;
    	// TODO: deal with prop invariant?
    	int categoryCount = siteModel.getCategoryCount();
    	for (int i = 0; i < categoryCount; i++) {
    		for (int j = 0; j < patternCount; j++) {
    			int state = states[j];
    			if (state >= 0 && state < nStates) {
    				// only contribute for non-ambiguous states
    				current[w++] = logProbabilities[v+state];
    			} else {
    				w++;
    			}
    		}
    		v = v + nStates;
    	}
		
	}
    
    
    protected void accumulateLogLeafScale() {
		final int n = accumulatedLogLeafScaleFactors.length;
		int leafNodeCount = nodeCount / 2 + 1;
		
		// perform delta calculation if tree is not filthy
		if (hasDirt == Tree.IS_FILTHY) {
			// recalc from scratch
			Arrays.fill(accumulatedLogLeafScaleFactors, 0.0);
			for (int j = 0; j < leafNodeCount; j++) {
				double [] x = leafLogScaleFactors[leafIndex[j]][j];
				for (int i = 0; i < n; i++) {
					accumulatedLogLeafScaleFactors[i] += x[i];
				}
			}
			
		} else {
			// calc delta
			for (int j = 0; j < leafNodeCount; j++) {
				if (leafIndex[j] != storedLeafIndex[j]) {
					double [] x = leafLogScaleFactors[leafIndex[j]][j];
					double [] oldX = leafLogScaleFactors[storedLeafIndex[j]][j];
					for (int i = 0; i < n; i++) {
						accumulatedLogLeafScaleFactors[i] += x[i] - oldX[i];
					}
				}
			}
		}

		
    }
    
	protected void setLeafScaleForUpdate(int nodeIndex) {
		leafIndex[nodeIndex] = 1 - leafIndex[nodeIndex];
	}
    
	/* return copy of pattern log likelihoods for each of the patterns in the alignment */
	public double [] getPatternLogLikelihoods() {
		if (beagleLikelihood != null && scaleFactor == 1.0) {
			return beagleLikelihood.getPatternLogLikelihoods();
		}
		return patternLogLikelihoods.clone();
	} // getPatternLogLikelihoods

    /** CalculationNode methods **/

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
        if (beagleLikelihood != null && scaleFactor == 1.0) {
            return beagleLikelihood.requiresRecalculation();
        }
        hasDirt = Tree.IS_CLEAN;
        
        if (siteModel.isDirtyCalculation() && treeInput.get().somethingIsDirty()) {
        	// if both site model and tree are dirty, it is time for a fresh recalcalation
            hasDirt = Tree.IS_FILTHY;
            return true;
        }

        if (alignment.isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            //m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        return treeInput.get().somethingIsDirty();
    }

    @Override
    public void store() {
    	storedScaleFactor = scaleFactor;
    	startScaling = false;
        if (beagleLikelihood  != null && scaleFactor == 1.0) {
        	beagleLikelihood.store();
            super.store();
            return;
        }
        if (likelihoodCore != null) {
            likelihoodCore.store();
        }
        super.store();
        if (branchLengths != null) {
        	System.arraycopy(branchLengths, 0, storedBranchLengths, 0, branchLengths.length);
        }
        System.arraycopy(accumulatedLogLeafScaleFactors, 0, storedAccumulatedLogLeafScaleFactors, 0, accumulatedLogLeafScaleFactors.length);
        System.arraycopy(leafIndex, 0, storedLeafIndex, 0, leafIndex.length);
        System.arraycopy(rates, 0, storedRates, 0, rates.length);
        System.arraycopy(scale, 0, storedScale, 0, scale.length);
    }

    @Override
    public void restore() {
    	scaleFactor = storedScaleFactor;
    	if (startScaling) {
    		startScaling = false;
    		scaleFactor = 1.0;
        	useScaleFactors = false;
            m_nScale = 0;
            Log.warning.println("Restoring: Turning off scaling");
            likelihoodCore.setUseScaling(scaleFactor);
    	}
        if (beagleLikelihood != null && scaleFactor == 1.0) {
        	beagleLikelihood.restore();
            super.restore();
            return;
        }
        if (likelihoodCore != null) {
            likelihoodCore.restore();
        }
        super.restore();
        double[] tmp;
        if (branchLengths != null) {
        	tmp = branchLengths;
        	branchLengths = storedBranchLengths;
        	storedBranchLengths = tmp;
        }
        
        tmp = accumulatedLogLeafScaleFactors; accumulatedLogLeafScaleFactors = storedAccumulatedLogLeafScaleFactors; storedAccumulatedLogLeafScaleFactors = tmp;
        tmp = rates; rates = storedRates; storedRates = tmp;
        tmp = scale; scale = storedScale; storedScale = tmp;

        int[] tmp2 = leafIndex; leafIndex = storedLeafIndex; storedLeafIndex = tmp2; 
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    @Override
	public List<String> getArguments() {
        return Collections.singletonList(alignment.getID());
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    @Override
	public List<String> getConditions() {
        return siteModel.getConditions();
    }	

} // class TreeLikelihood
