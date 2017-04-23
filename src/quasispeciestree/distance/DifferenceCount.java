package quasispeciestree.distance;

import beast.core.Description;
import beast.evolution.alignment.distance.Distance;

/**
 *  @author Veronika Boskova created on 03/03/17
 */
@Description("Counts the simple number of characters that differ between sequences. " +
             "Note that unknowns are not ignored, so if both are unknowns '?' the distance is zero.")
public class DifferenceCount extends Distance.Base {

        @Override
        public double pairwiseDistance(int taxon1, int taxon2) {
            double dist = 0;
            for (int i = 0; i < patterns.getPatternCount(); i++) {
                if (patterns.getPattern(taxon1, i) != patterns.getPattern(taxon2, i)) {
                    dist += patterns.getPatternWeight(i);
                }
            }
            return dist;  // / patterns.getSiteCount();
        }
}
