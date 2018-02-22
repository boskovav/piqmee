package piqmee.util;

import beast.util.Randomizer;

/**
 * @author Veronika Boskova created on 24/10/17
 */

public class BeastRandomGenerator implements RandomGenerator {
    @Override
    public double getNext() {
        return Randomizer.nextDouble();
    }
}
