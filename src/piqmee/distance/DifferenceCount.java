package piqmee.distance;

import beast.core.Description;
import beast.evolution.alignment.distance.Distance;

/**
 *  @author Veronika Boskova created on 03/03/17
 */
@Description("Defines difference between sequences based simply on characters that differ between sequences. " +
             "Note that if characters in both sequences are unknowns '?' the distance is zero. " +
             "If collapseIdenticalSequences is true, distance 0 will be determined also for ambiguous characters that " +
             "encompass a single unique or at least one ambiguous character in the other sequence.")
public class DifferenceCount extends Distance.Base {

    /**
     * Calculate a pairwise distance
     */
    public double pairwiseDistance(int taxon1, int taxon2, boolean collapseIdenticalUptoMissingData) {
        int state1, state2;
        int[] pattern;
        double dist = 0;
        for (int i = 0; i < patterns.getPatternCount(); i++) {
            pattern = patterns.getPattern(i);
            state1 = pattern[taxon1];
            state2 = pattern[taxon2];
            // check if patterns represent a single character
            if (!dataType.isAmbiguousCode(state1) && !dataType.isAmbiguousCode(state2)
                && state1 != state2) {
                    dist += patterns.getPatternWeight(i);
            } else if (collapseIdenticalUptoMissingData){
                int[] states1 = dataType.getStatesForCode(state1);
                int[] states2 = dataType.getStatesForCode(state2);
                if (!isAtLeastOneStateIdentical(states1,states2)) {
                    // only a single character difference is enough for us to determine that sequences are not identical
                    dist += patterns.getPatternWeight(i);
                }
            }
        }
        return dist;  // / patterns.getSiteCount();
    }

    /**
     * Calculate a pairwise difference (0 = identical, 1 = not identical)
     */
    public double pairwiseDifference(int taxon1, int taxon2, boolean collapseIdenticalUptoMissingData) {
        int state1, state2;
        int[] pattern;
        for (int i = 0; i < patterns.getPatternCount(); i++) {
            pattern = patterns.getPattern(i);
            state1 = pattern[taxon1];
            state2 = pattern[taxon2];
            // check if patterns represent a single character
            if (!dataType.isAmbiguousCode(state1) && !dataType.isAmbiguousCode(state2)){
                if (state1 != state2) {
                    // only a single character difference is enough for us to determine that sequences are not identical
                    return 1;
                }
            }
            // if patterns represent multiple characters, a single match is enough to call them identical
            else if (collapseIdenticalUptoMissingData){
                int[] states1 = dataType.getStatesForCode(state1);
                int[] states2 = dataType.getStatesForCode(state2);
                if (!isAtLeastOneStateIdentical(states1,states2)) {
                    // only a single character difference is enough for us to determine that sequences are not identical
                    return 1;
                }
            }
        }
        return 0;  // / patterns.getSiteCount();
    }

    /**
     * Compares two arrays representing ambiguous sequences
     */
    public boolean isAtLeastOneStateIdentical(int[] states1, int[] states2) {
        for (int j = 0; j < states1.length; j++){
            for (int k = 0; k < states2.length; k++){
                if (states1[j] == states2[k]) {
                    return true;
                }
            }
        }
        return false;
    }

}
