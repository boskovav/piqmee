package test.piqmee;

import piqmee.util.RandomGenerator;

/**
 * @author Veronika Boskova created on 24/10/17
 */

public class DeterministicRandomGenerator implements RandomGenerator {

    double deterministicnumber;

    public DeterministicRandomGenerator(double random) {
        this.deterministicnumber = random;
    }

    @Override
    public double getNext() {
        return deterministicnumber;
    }
}
