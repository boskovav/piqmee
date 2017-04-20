/*
* File SubstitutionModel.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
package quasispeciestree.substitutionmodel;


import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.tree.Node;
import beast.evolution.substitutionmodel.SubstitutionModel;
import quasispeciestree.tree.QuasiSpeciesNode;

@Description("Specifies substitution model from which a transition probability matrix for a given " +
        "distance can be obtained.")
public interface QuasiSpeciesSubstitutionModel extends SubstitutionModel{

    /**
     * basic implementation of a SubstitutionModel bringing together relevant super class*
     */
    @Description(value = "Base implementation of a substitution model.", isInheritable = false)
    public abstract class Base extends CalculationNode implements QuasiSpeciesSubstitutionModel {
        public Input<Frequencies> frequenciesInput =
                new Input<Frequencies>("frequencies", "substitution model equilibrium state frequencies", Input.Validate.REQUIRED);

        /**
         * shadows frequencies, or can be set by subst model *
         */
        protected Frequencies frequencies;

        /**
         * number of states *
         */
        protected int nrOfStates;

        @Override
        public void initAndValidate(){
            frequencies = frequenciesInput.get();
        }

        @Override
        public double[] getFrequencies() {
            return frequencies.getFreqs();
        }

        @Override
        public int getStateCount() {
            return nrOfStates;
        }


        @Override
        public boolean canReturnComplexDiagonalization() {
            return false;
        }

        @Override
        public double[] getRateMatrix(Node node) {
            return null;
        }

    } // class Base

    /**
     * get the complete transition probability matrix for the given distance
     * determined as (fStartTime-fEndTime)*fRate
     *
     * @param node              tree node for which to calculate the probabilities
     * @param fTime             total time of the QS
     * @param fRate             rate, includes gamma rates and branch rates
     * @param nochangematrix    an array to store the transition probability (P(no change happens at (fStartTime-fEndTime)*fRate))
     *                          So, nochangematrix must be of size n where n is number of states.
     */
    void getTransitionProbabilities(QuasiSpeciesNode node, double fTime, double fRate, double[] nochangematrix);

}
