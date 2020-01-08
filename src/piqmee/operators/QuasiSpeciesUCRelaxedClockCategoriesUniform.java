package piqmee.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.Parameter;
import beast.util.Randomizer;
import piqmee.tree.QuasiSpeciesNode;
import piqmee.tree.QuasiSpeciesTree;

@Description("Assign one parameter value to a uniformly selected value in its range - for QS UCRelaxedClock model.")

public class QuasiSpeciesUCRelaxedClockCategoriesUniform extends Operator {
    final public Input<Parameter<?>> parameterInput = new Input<>("categories", "categories of UC relaxed clock model, specified as integer parameter to sample individual values for", Input.Validate.REQUIRED, Parameter.class);
    public Input<QuasiSpeciesTree> quasiSpeciesTreeInput = new Input<>(
            "quasiSpeciesTree", "Quasi-species tree on which to operate.",
            Input.Validate.REQUIRED);

    Parameter<?> parameter;
    double lower, upper;
    int lowerIndex, upperIndex;
    protected QuasiSpeciesTree qsTree;

    @Override
    public void initAndValidate() {
        parameter = parameterInput.get();
        if (parameter instanceof IntegerParameter) {
            lowerIndex = (Integer) parameter.getLower();
            upperIndex = (Integer) parameter.getUpper();
        } else {
            throw new IllegalArgumentException("parameter should be a IntegerParameter, not " + parameter.getClass().getName());
        }
        qsTree = quasiSpeciesTreeInput.get();
    }

    @Override
    public double proposal() {
        int index = Randomizer.nextInt(parameter.getDimension());

        int leafCount = qsTree.getLeafNodeCount();
        int intNodeCount = qsTree.getInternalNodeCount();
        int nodeCount = leafCount + intNodeCount;
        while ((index >= nodeCount && ((QuasiSpeciesNode)qsTree.getNode(index-nodeCount)).getHaploAboveName() == -1)
                || (index < leafCount && ((QuasiSpeciesNode)qsTree.getNode(index)).getAttachmentTimesList().length < 2)
                || (index >= leafCount && index < nodeCount && ((QuasiSpeciesNode)qsTree.getNode(index)).getContinuingHaploName() != -1 &&
                ((QuasiSpeciesNode)qsTree.getNode(index)).getHaploAboveName() == -1 ))
            index = Randomizer.nextInt(parameter.getDimension());

        int newValue = Randomizer.nextInt(upperIndex - lowerIndex + 1) + lowerIndex; // from 0 to n-1, n must > 0,
        ((IntegerParameter) parameter).setValue(index, newValue);

        return 0.0;
    }

}
