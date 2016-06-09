package quasispeciestree.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import quasispeciestree.tree.QuasiSpeciesNode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *  @author Veronika Boskova created on 21/04/2016 finished on 24/04/2016
 */
@Description("Scale operator for quasispecies trees. Also allows additional "
        + "scalar parameters to be rescaled (either forward or inversely) "
        + "at the same time.")
public class QuasiSpeciesTreeScale extends QuasiSpeciesTreeOperator{

    public Input<List<RealParameter>> parametersInput =
        new Input<>("parameter",
            "Scale this scalar parameter by the same amount as tree.",
            new ArrayList<RealParameter>());

    public Input<List<BooleanParameter>> indicatorsInput =
        new Input<>("indicator",
            "If provided, used to specify a subset of parameter elements to scale.",
            new ArrayList<BooleanParameter>());

    public Input<List<RealParameter>> parametersInverseInput =
        new Input<>("parameterInverse",
            "Scale this scalar parameter inversely.",
            new ArrayList<RealParameter>());

    public Input<List<BooleanParameter>> indicatorsInverseInput =
        new Input<>("indicatorInverse",
            "If provided, used to specify a subset of parameter elements to scale "
            + "inversely.",
            new ArrayList<BooleanParameter>());

    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
           "Scaling is restricted to the range [1/scaleFactor, scaleFactor]");

    final public Input<Boolean> rootOnlyInput = new Input<>("rootOnly", "scale root of a tree only, ignored if tree is not specified (default false)", false);

    boolean indicatorsUsed, indicatorsInverseUsed;

    @Override
    public void initAndValidate() {

        super.initAndValidate();

        if (indicatorsInput.get().size()>0) {
            if (indicatorsInput.get().size() != parametersInput.get().size())
                throw new IllegalArgumentException("If an indicator element "
                        + "exists, the number of such elements must equal "
                        + "the number of parameter elements.");

            for (int pidx=0; pidx<parametersInput.get().size(); pidx++) {
                if (parametersInput.get().get(pidx).getDimension() !=
                        indicatorsInput.get().get(pidx).getDimension()) {
                    throw new IllegalArgumentException("The number of boolean "
                            + "values in indicator element "
                            + String.valueOf(pidx+1)
                            + " doesn't match the dimension of the "
                            + "corresponding parameter element.");
                }
            }
            indicatorsUsed = true;
        } else
            indicatorsUsed = false;

        if (indicatorsInverseInput.get().size()>0) {
            if (indicatorsInverseInput.get().size() != parametersInverseInput.get().size())
                throw new IllegalArgumentException("If an indicatorInverse element "
                        + "exists, the number of such elements must equal "
                        + "the number of parameterInverse elements.");

            for (int pidx=0; pidx<parametersInverseInput.get().size(); pidx++) {
                if (parametersInverseInput.get().get(pidx).getDimension() !=
                        indicatorsInverseInput.get().get(pidx).getDimension()) {
                    throw new IllegalArgumentException("The number of boolean "
                            + "values in indicatorInverse element "
                            + String.valueOf(pidx+1)
                            + " doesn't match the dimension of the "
                            + "corresponding parameterInverse element.");
                }
            }
            indicatorsInverseUsed = true;
        } else
            indicatorsInverseUsed = false;
    }

    @Override
    public double proposal() {

        // Choose scale factor:
        double u = Randomizer.nextDouble();
        double f = u*scaleFactorInput.get()+(1.0-u)/scaleFactorInput.get();

        // Keep track of Hastings ratio:
        double logf = Math.log(f);
        double logHastingsRatio = 0.0;

        final Node root = qsTree.getRoot();

        // scaling only the root
//        if (rootOnlyInput.get()) {
//            final Node root = qsTree.getRoot();
//            double oldHeight = root.getHeight();
//            final double newHeight = oldHeight * f;
//            logHastingsRatio -= logf;
//
//            if (newHeight < Math.max(root.getLeft().getHeight(), root.getRight().getHeight())) {
//                return Double.NEGATIVE_INFINITY;
//            }
//            root.setHeight(newHeight);
//
//            // check whether there is a haplotype arising above the root at the time of moving
//            //      if so,
//            //           1) check whether there is a second one arising "above" root, when root height changed (to lower)
//            //              if so, scale the second haplo below the newHeight of the root
//            //           2) check whether there is still one arising "above" the root, when root height changed (to higher)
//            //              if so, assign the haplo the left/right child --- depending on which subtree the haplo belongs to
//            //      else check whether there is haplotype arising "above" root, when root height changed (to lower)
//            //          if so, change the aboveNodeHaplo tag of the root & corresponding child (left or right)
//            // P(data|tree) and P(tree|eta) change only because of moved root (not because of moved haplo -- no haplo moving here)
//            // TODO but the operator move has to be reversible -- so at every move of the root height re-scale one left and one right haplo
//            int haplo = ((QuasiSpeciesNode) root).getHaploAboveName();
//
//            QuasiSpeciesNode left = (QuasiSpeciesNode)root.getLeft();
//            int haploleft = ((QuasiSpeciesNode)root.getLeft()).getHaploAboveName();
//            QuasiSpeciesNode right = (QuasiSpeciesNode)root.getRight();
//            int haploright = ((QuasiSpeciesNode)root.getRight()).getHaploAboveName();
//
//            int[] parentHaplo=qsTree.getParentHaplo();
//            ArrayList<Integer> possibleHaploLeft = new ArrayList<>();
//            ArrayList<Integer> possibleHaploRight = new ArrayList<>();
//
//            // the haplo at the root can come from the left child
//            if (haplo != -1 && haplo == left.getContinuingHaploName()){
//                if (haploleft != -1){
//                    throw new IllegalStateException("Seems like we found a haplotype at the root "
//                            + "that should come from the left subtree but yet there is another haplotype "
//                            + "arising at the left child of the root. This should not happen!");
//                }
//                else{
//                    haploleft = haplo;
//                    possibleHaploLeft.add(haplo);
//                }
//            }
//            // otherwise, choose which haplotype from the left subtree to scale
//            else if (haploleft == -1){
//                checkNumberOfPossibleSrcHaplo(left, parentHaplo, possibleHaploLeft);
//                // choose randomly from the array of haplotypes one to move up
//                haploleft = possibleHaploLeft.get(Randomizer.nextInt(possibleHaploLeft.size()));
//            }
//            // otherwise the haploleft is unique and already assigned above
//            // Incorporate probability of choosing current left haplotype to move --- only with Felsenstein
////            logHastingsRatio += Math.log(possibleHaploLeft.size());
//
//
//
//            // the haplo at the root can come from the right child
//            if (haplo == -1 && haplo == right.getContinuingHaploName()){
//                if (haploleft != -1){
//                    throw new IllegalStateException("Seems like we found a haplotype at the root "
//                            + "that should come from the right subtree but yet there is another haplotype "
//                            + "arising at the right child of the root. This should not happen!");
//                }
//                else{
//                    haploright = haplo;
//                    possibleHaploRight.add(haplo);
//                }
//            }
//            // otherwise, choose which haplotype from the left subtree to scale
//            else if (haploright == -1){
//                checkNumberOfPossibleSrcHaplo(right, parentHaplo, possibleHaploRight);
//                // choose randomly from the array of haplotypes one to move up
//                haploright = possibleHaploRight.get(Randomizer.nextInt(possibleHaploRight.size()));
//            }
//            // otherwise the haploright is unique and already assigned above
//            // Incorporate probability of choosing current right haplotype to move --- only with Felsenstein
////            logHastingsRatio += Math.log(possibleHaploRight.size());
//
//
//            double haploStartMinLeft=qsTree.getNode(haploleft).getHeight();
//            double haploStartMinRight=qsTree.getNode(haploright).getHeight();
//            double haploStartMaxLeft=0;
//            double haploStartMaxNewLeft=0;
//            double haploStartMaxRight=0;
//            double haploStartMaxNewRight=0;
//            double scalefactorLeft=0;
//            double scalefactorReft=0;
//
//            // check what was the max attachment time for left and rigth haplo
//            if (haplo==haploleft){
//                haploStartMaxLeft=origin.getValue();
//                haploStartMaxRight=oldHeight;
//            }
//            else if (haplo==haploright){
//                haploStartMaxRight=origin.getValue();
//                haploStartMaxLeft=oldHeight;
//            }
//
//            // choose which one to move till root and which one above the root
//            if (Randomizer.nextBoolean()){
//                haploStartMaxNewLeft = origin.getValue();
//                haploStartMaxNewRight = newHeight;
//            }
//            else {
//                haploStartMaxNewRight = origin.getValue();
//                haploStartMaxNewLeft = newHeight;
//            }
//
//            // reposition lefthaplo
//            // get a random number deciding where the current haplo will be moved
//            double l = Randomizer.nextDouble();
//            double tnewLeft = l*haploStartMinLeft + (1-l)*haploStartMaxLeft;
//            // reposition attachment times: attach ((time - haploStartMin) * (tHaploNew/told)) + haploStartMin
//            Double[] tempqstimesleft=qsTree.getAttachmentTimesList(haploleft).clone();
//            // get the haplotype's starting time
//            double toldQSstartleft = tempqstimesleft[0];
//            double tnewQSstartleft = 0;
//            double toldLeft=0;
//            if (tempqstimesleft.length > 1){
//                toldLeft = tempqstimesleft[1];
//                // Scale the haplotype strains
//                // scale all the other positions in the array but the 0 position (haplo start time)
//                scalefactorLeft = (tnewLeft - haploStartMinLeft)/(toldLeft - haploStartMinLeft);
//                for (int i=1; i<tempqstimesleft.length; i++) {
//                    tempqstimesleft[i] = ((tempqstimesleft[i] - haploStartMinLeft) * scalefactorLeft) + haploStartMinLeft;
//                }
//                // set the haplotype's starting time to the new time
//                double x = Randomizer.nextDouble();
//                tnewQSstartleft = x*tempqstimesleft[1] + (1-x)*haploStartMaxNewLeft;
//                tempqstimesleft[0] = tnewQSstartleft;
//                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
////                logHastingsRatio -= Math.log(haploStartMaxLeft - toldLeft);
////                logHastingsRatio += Math.log(haploStartMaxNewLeft - tempqstimesleft[1]);
//            }
//            else {
//                // set the haplotype's starting time to the new time
//                double x = Randomizer.nextDouble();
//                tnewQSstartleft = x*haploStartMinLeft + (1-x)*haploStartMaxNewLeft;
//                tempqstimesleft[0] = tnewQSstartleft;
//                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
////                logHastingsRatio -= Math.log(haploStartMaxLeft - haploStartMinLeft);
////                logHastingsRatio += Math.log(haploStartMaxNewLeft - haploStartMinLeft);
//            }
//            // rewrite the attachment times array
//            qsTree.setAttachmentTimesList(haploleft, tempqstimesleft);
//
//            // re-assign where the haplo starts at
//            QuasiSpeciesNode oldnodeBelowLeft = null;
//            QuasiSpeciesNode newnodeBelowLeft = null;
//            Node templeft = qsTree.getNode(haploleft);
//            while (oldnodeBelowLeft==null && newnodeBelowLeft==null && templeft!=root){
//                if(tnewQSstartleft>templeft.getParent().getHeight() && tnewQSstartleft<templeft.getHeight()){
//                    newnodeBelowLeft = (QuasiSpeciesNode) templeft;
//                }
//                if (((QuasiSpeciesNode)templeft).getHaploAboveName()==haploleft){
//                    oldnodeBelowLeft = (QuasiSpeciesNode) templeft;
//                }
//                templeft=templeft.getParent();
//            }
//            if (newnodeBelowLeft==null && tnewQSstartleft>root.getHeight() && tnewQSstartleft<origin.getValue())
//                newnodeBelowLeft=(QuasiSpeciesNode)root;
//            if (((QuasiSpeciesNode)root).getHaploAboveName()==haploleft)
//                oldnodeBelowLeft=(QuasiSpeciesNode)root;
//            //change the haploAboveNames for the nodes
//            if (newnodeBelowLeft!=oldnodeBelowLeft){
//                newnodeBelowLeft.setHaploAboveName(haploleft);
//                oldnodeBelowLeft.setHaploAboveName(-1);
//            }
//
//
//            // TODO this should be done after the right haplo was also repositioned
//            // Recalculate continuingHaplo and HaploAbove arrays
//            recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) root);
////            in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
//            int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
//            qsTree.setStartBranchCounts(startBranchCountsArray);
//            // TODO
//
//
//            // Incorporate probability of choosing current haplotype to move --- only with Felsenstein
////            logHastingsRatio += Math.log(possibleHaplo.size());
//            if (tempqstimesleft.length>1){
//                // Incorporate probability of current positioning of the start time of the haplotype.
//                logHastingsRatio -= Math.log(haploStartMaxLeft-haploStartMinLeft);
//                // Incorporate probability of selected positioning of the start time of the haplotype.
//                logHastingsRatio += Math.log(haploStartMaxNewLeft-haploStartMinLeft);
//                // scaling the attachment times (see Scaling Operator)
//                // assign contribution to the Hastings ratio for having different possible scales for told
//                logHastingsRatio += Math.log(haploStartMaxNewLeft/toldLeft - haploStartMinLeft/toldLeft);
//                logHastingsRatio -= Math.log(haploStartMaxNewLeft/tnewLeft - haploStartMinLeft/tnewLeft);
//                // assign contribution of each scaled attachment time
//                logHastingsRatio -= 2 * (scalefactorLeft * (haploStartMaxNewLeft - haploStartMinLeft) + haploStartMinLeft)/(haploStartMaxNewLeft);
//                logHastingsRatio += (tempqstimesleft.length-1) * Math.log(scalefactorLeft);
//            }
//            // Incorporate probability of choosing current haplotype to move back --- only with Felsenstein
//            if (left.getContinuingHaploName() != haploleft){
//                // check how many haplotypes can be pulled up
//                ArrayList<Integer> backPossibleHaplo = new ArrayList<>();
////                int[] newParentHaplo = qsTree.getParentHaplo();
////                checkNumberOfPossibleSrcHaplo((QuasiSpeciesNode) srcNode, newParentHaplo, backPossibleHaplo);
//                checkNumberOfPossibleSrcHaplo(left, qsTree.getParentHaplo(), backPossibleHaplo);
//                // account for this in the hastings ratio
////                logHastingsRatio -= Math.log(backPossibleHaplo.size());
//            }
//            // no need to be adding 0.0, so just commented out
//            // else{
//            //  logHastingsRatio -= Math.log(1);
//            // }
//
//            //reposition righthaplo
//
//                // recalcualte how many haplotypes are possible in the backward move
//
//
//
////            if (haplo != -1){
////                if (newHeight < oldHeight){
////                    // check left child
////                    if (haploleft != -1){
////                        Double[] tempqstimesleft=qsTree.getAttachmentTimesList(haploleft);
////                        if (tempqstimesleft[0]>newHeight){
////                            //return Double.NEGATIVE_INFINITY;
////                            double tminleft=left.getHeight();
////                            double tmaxleft=newHeight;
////
////                            double toldleft=0;
////                            double tnewleft=0;
////                            double scalefactorleft=0;
////
////                            // choose new time to attach
////                            if (tempqstimesleft.length >1)
////                                tnewleft = u*tminleft + (1-u)*tmaxleft;
////
////                            // reposition attachment times: attach ((time - tmin) * (tnew/told)) + tmin
////                            // get the haplotype's starting time
////                            double toldQSstartleft = tempqstimesleft[0];
////                            double tnewQSstartleft = 0;
////                            if (tempqstimesleft.length > 1){
////                                toldleft = tempqstimesleft[1];
////                                // scale all the other positions in the array but the 0 position (haplo start time)
////                                // reposition attachment times: attach ((time - tmin) * (tnew/told)) + tmin
////                                scalefactorleft = (tnewleft - tminleft)/(toldleft - tminleft);
////                                for (int i=1; i<tempqstimesleft.length; i++) {
////                                    tempqstimesleft[i] = ((tempqstimesleft[i] - tminleft) * scalefactorleft) + tminleft;
////                                }
////                                // set the haplotype's starting time to the new time
////                                double x = Randomizer.nextDouble();
////                                tnewQSstartleft = x*tempqstimesleft[1] + (1-x)*tmaxleft;
////                                tempqstimesleft[0] = tnewQSstartleft;
////                                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//////                logHastingsRatio -= Math.log(tmaxleft - toldleft);
//////                logHastingsRatio += Math.log(tmaxleft - tempqstimesleft[1]);
////                            }
////                            else {
////                                // set the haplotype's starting time to the new time
////                                double x = Randomizer.nextDouble();
////                                tnewQSstartleft = x*tminleft + (1-x)*tmaxleft;
////                                tempqstimesleft[0] = tnewQSstartleft;
////                                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
////                                //              logHastingsRatio -= Math.log(tmaxleft - tminleft);
////                                //              logHastingsRatio += Math.log(tmaxleft - tminleft);
////                                // since they are identical, no need to change Hastings ratio
////                            }
////                            // rewrite the attachment times array
////                            qsTree.setAttachmentTimesList(left, tempqstimesleft);
////
////                            // Incorporate the probability of scaling all the attachment times
////                            if (tempqstimesleft.length > 1){
////                                // assign contribution to the Hastings ratio for having different possible scales for told
////                                logHastingsRatio += Math.log(tmaxleft/toldleft - tminleft/toldleft);
////                                logHastingsRatio -= Math.log(tmaxleft/tnewleft - tminleft/tnewleft);
////                                // assign contribution of each scaled attachment time
////                                logHastingsRatio -= 2 * (Math.log (scalefactorleft * (toldleft - tminleft) + tminleft) - Math.log(toldleft));
////                                logHastingsRatio += (tempqstimesleft.length-1) * Math.log(scalefactorleft);
////                            }
////
////                            // change internal node haploName - both remove and add if necessary
////                            Node intnode = left;
////                            Node intnodeparent;
////                            if (tnewQSstartleft > toldQSstartleft) {
////                                intnodeparent = intnode.getParent();
////                                while (intnodeparent!=null && intnodeparent.getHeight() < tnewQSstartleft){
////                                    intnode = intnodeparent;
////                                    intnodeparent = intnodeparent.getParent();
////                                }
////                            }
////                            else {
////                                intnode = qsTree.getNode(haploleft);
////                                intnodeparent = intnode.getParent();
////                                while (intnodeparent != null && intnodeparent.getHeight() < tnewQSstartleft){
////                                    intnode = intnodeparent;
////                                    intnodeparent = intnodeparent.getParent();
////                                }
////                            }
////                            if (intnode != left){
////                                // change the internal nodes' haploName, where necessary
////                                left.setHaploAboveName(-1);
////                                ((QuasiSpeciesNode) intnode).setHaploAboveName(haploleft);
////                            }
////                        }
////                    }
////                    // check right child
////                    if (haploright != -1){
////                        Double[] tempqstimesright=qsTree.getAttachmentTimesList(haploright);
////                        if (tempqstimesright[0]>newHeight){
////                            //return Double.NEGATIVE_INFINITY;
////                            double tminright=right.getHeight();
////                            double tmaxright=newHeight;
////
////                            double toldright=0;
////                            double tnewright=0;
////                            double scalefactorright=0;
////
////                            // choose new time to attach
////                            if (tempqstimesright.length >1)
////                                tnewright = u*tminright + (1-u)*tmaxright;
////
////                            // reposition attachment times: attach ((time - tmin) * (tnew/told)) + tmin
////                            // get the haplotype's starting time
////                            double toldQSstartright = tempqstimesright[0];
////                            double tnewQSstartright = 0;
////                            if (tempqstimesright.length > 1){
////                                toldright = tempqstimesright[1];
////                                // scale all the other positions in the array but the 0 position (haplo start time)
////                                // reposition attachment times: attach ((time - tmin) * (tnew/told)) + tmin
////                                scalefactorright = (tnewright - tminright)/(toldright - tminright);
////                                for (int i=1; i<tempqstimesright.length; i++) {
////                                    tempqstimesright[i] = ((tempqstimesright[i] - tminright) * scalefactorright) + tminright;
////                                }
////                                // set the haplotype's starting time to the new time
////                                double x = Randomizer.nextDouble();
////                                tnewQSstartright = x*tempqstimesright[1] + (1-x)*tmaxright;
////                                tempqstimesright[0] = tnewQSstartright;
////                                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//////                logHastingsRatio -= Math.log(tmaxright - toldright);
//////                logHastingsRatio += Math.log(tmaxright - tempqstimesright[1]);
////                            }
////                            else {
////                                // set the haplotype's starting time to the new time
////                                double x = Randomizer.nextDouble();
////                                tnewQSstartright = x*tminright + (1-x)*tmaxright;
////                                tempqstimesright[0] = tnewQSstartright;
////                                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
////                                //              logHastingsRatio -= Math.log(tmaxright - tminright);
////                                //              logHastingsRatio += Math.log(tmaxright - tminright);
////                                // since they are identical, no need to change Hastings ratio
////                            }
////                            // rewrite the attachment times array
////                            qsTree.setAttachmentTimesList(right, tempqstimesright);
////
////                            // Incorporate the probability of scaling all the attachment times
////                            if (tempqstimesright.length > 1){
////                                // assign contribution to the Hastings ratio for having different possible scales for told
////                                logHastingsRatio += Math.log(tmaxright/toldright - tminright/toldright);
////                                logHastingsRatio -= Math.log(tmaxright/tnewright - tminright/tnewright);
////                                // assign contribution of each scaled attachment time
////                                logHastingsRatio -= 2 * (Math.log (scalefactorright * (toldright - tminright) + tminright) - Math.log(toldright));
////                                logHastingsRatio += (tempqstimesright.length-1) * Math.log(scalefactorright);
////                            }
////
////                            // change internal node haploName - both remove and add if necessary
////                            Node intnode = right;
////                            Node intnodeparent;
////                            if (tnewQSstartright > toldQSstartright) {
////                                intnodeparent = intnode.getParent();
////                                while (intnodeparent!=null && intnodeparent.getHeight() < tnewQSstartright){
////                                    intnode = intnodeparent;
////                                    intnodeparent = intnodeparent.getParent();
////                                }
////                            }
////                            else {
////                                intnode = qsTree.getNode(haploright);
////                                intnodeparent = intnode.getParent();
////                                while (intnodeparent != null && intnodeparent.getHeight() < tnewQSstartright){
////                                    intnode = intnodeparent;
////                                    intnodeparent = intnodeparent.getParent();
////                                }
////                            }
////                            if (intnode != right){
////                                // change the internal nodes' haploName, where necessary
////                                right.setHaploAboveName(-1);
////                                ((QuasiSpeciesNode) intnode).setHaploAboveName(haploright);
////                            }
////                        }
////                    }
////
////                    // recalculate parentHaplo array, re-assign continuingHaploName
////                    recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) root);
////
////                    // in any case (changed or not the aboveNodeHaplo/parentHaplo array) recalculate countPossibleStartBranches
////                    int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
////                    qsTree.setStartBranchCounts(startBranchCountsArray);
////                }
////                else if (newHeight > qsTree.getAttachmentTimesList(haplo)[0]){
////                    ((QuasiSpeciesNode) root).setHaploAboveName(-1);
////                    ((QuasiSpeciesNode) root).setContinuingHaploName(-1);
////                    if (left.getContinuingHaploName() == haplo)
////                        left.setHaploAboveName(haplo);
////                    else if (right.getContinuingHaploName() == haplo)
////                        right.setHaploAboveName(haplo);
////                }
////            }
////            else {
////                if (newHeight < oldHeight){
////                    // check if only one haplotype above the root
////                    if (haploleft != -1 && haploright != -1){
////                        Double[] tempqstimesleft=qsTree.getAttachmentTimesList(haploleft);
////                        Double[] tempqstimesright=qsTree.getAttachmentTimesList(haploright);
////                        if (tempqstimesleft[0]>newHeight && tempqstimesright[0]>newHeight)
////                            return Double.NEGATIVE_INFINITY;
////                    }
////                    // check left child
////                    else if (haploleft != -1){
////                        Double[] tempqstimesleft=qsTree.getAttachmentTimesList(haploleft);
////                        if (tempqstimesleft[0]>newHeight){
////                            ((QuasiSpeciesNode) root).setHaploAboveName(haploleft);
////                            ((QuasiSpeciesNode) root).setContinuingHaploName(haploleft);
////                            left.setHaploAboveName(-1);
////                        }
////                    }
////                    // check right child
////                    if (haploright != -1){
////                        Double[] tempqstimesright=qsTree.getAttachmentTimesList(haploright);
////                        if (tempqstimesright[0]>newHeight){
////                            ((QuasiSpeciesNode) root).setHaploAboveName(haploright);
////                            ((QuasiSpeciesNode) root).setContinuingHaploName(haploright);
////                            right.setHaploAboveName(-1);
////                        }
////                    }
////                }
////            }
//            // recalculate parentHaplo array
//            recalculateParentHaploAndCorrectContinuingHaploName(-1, (QuasiSpeciesNode) root);
//            // recalculate countPossibleStartBranches (at least for root this may have changed)
////            int[] startBranchCountsArray = qsTree.countPossibleStartBranches();
//            qsTree.setStartBranchCounts(startBranchCountsArray);
//
//            // Return Hastings ratio:
//            return logHastingsRatio;
//









        // scaling the entire tree
//        } else {
            logHastingsRatio -= 2 * logf;

            // Scale internal node heights
            for (Node node : qsTree.getInternalNodes()) {
                node.setHeight(node.getHeight()*f);
                logHastingsRatio += logf;
            }

            // Scale haplotype attachment times
            for (Node leaf : qsTree.getExternalNodes()) {
                Double[] tempqstimes=qsTree.getAttachmentTimesList((QuasiSpeciesNode)leaf).clone();
                for (int i=0; i<tempqstimes.length; i++) {
                    tempqstimes[i] = tempqstimes[i] * f;
                }
                // reject the move as soon as the first QS attachment time exceeds the origin height
                if (tempqstimes.length>1 && tempqstimes[1]>origin.getValue())
                    return Double.NEGATIVE_INFINITY;
                // select a new QS start time for QS above the root --- always
                if (((QuasiSpeciesNode) root).getContinuingHaploName()==leaf.getNr()){
                    double x = Randomizer.nextDouble();
                    if (tempqstimes.length>1){
                        if (root.getHeight()>tempqstimes[1]){
                            tempqstimes[0] = x*root.getHeight() + (1-x)*origin.getValue();
                        // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                      logHastingsRatio -= Math.log(origin.getValue() - root.getHeight()/f);
//                      logHastingsRatio += Math.log(origin.getValue() - root.getHeight());
                        }
                        else {
                            tempqstimes[0] = x*tempqstimes[1] + (1-x)*origin.getValue();
                        // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                      logHastingsRatio -= Math.log(origin.getValue() - tempqstimes[1]/f);
//                      logHastingsRatio += Math.log(origin.getValue() - tempqstimes[1]);
                        }
                    }
                    else {
                        tempqstimes[0] = x*root.getHeight() + (1-x)*origin.getValue();
                        // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                      logHastingsRatio -= Math.log(origin.getValue() - root.getHeight()/f);
//                      logHastingsRatio += Math.log(origin.getValue() - root.getHeight());
                    }
                }

                // make sure the QS haplo start does not interfere with the probability of acceptance of the move
                if (((QuasiSpeciesNode) leaf).getHaploAboveName()==leaf.getNr()){
                    double y = Randomizer.nextDouble();
                    if (tempqstimes.length>1){
                        tempqstimes[0] = y*tempqstimes[1] + (1-y)*leaf.getParent().getHeight();
                        // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                      logHastingsRatio -= Math.log(leaf.getParent().getHeight() - tempqstimes[1]/f);
//                      logHastingsRatio += Math.log(leaf.getParent().getHeight() - tempqstimes[1]);
                    }
                    else {
                        tempqstimes[0] = y*leaf.getHeight() + (1-y)*leaf.getParent().getHeight();
                        // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                      logHastingsRatio -= Math.log(leaf.getParent().getHeight() - leaf.getHeight()/f);
//                      logHastingsRatio += Math.log(leaf.getParent().getHeight() - leaf.getHeight());
                    }
                }

                qsTree.setAttachmentTimesList(leaf.getNr(), tempqstimes);
                // assign contribution of the QS start to the Hastings ratio --- only with Felsenstein
//                logHastingsRatio += logf;
                // contribution from haplotype duplicate attach time scaling
                if (tempqstimes.length > 1)
                    logHastingsRatio += (tempqstimes.length-1) * logf;
            }
//        }

        // Scale parameters:
        for (int pidx=0; pidx<parametersInput.get().size(); pidx++) {
            RealParameter param = parametersInput.get().get(pidx);
            for (int i=0; i<param.getDimension(); i++) {
                if (!indicatorsUsed ||
                        indicatorsInput.get().get(pidx).getValue(i)) {
                    double oldValue = param.getValue(i);
                    double newValue = oldValue*f;
                    if (newValue < param.getLower() || newValue > param.getUpper())
                        return Double.NEGATIVE_INFINITY;

                    param.setValue(i, newValue);
                    logHastingsRatio += logf;
                }
            }
        }

        // Scale parameters inversely:
        for (int pidx=0; pidx<parametersInverseInput.get().size(); pidx++) {
            RealParameter param = parametersInverseInput.get().get(pidx);
            for (int i=0; i<param.getDimension(); i++) {
                if (!indicatorsInverseUsed ||
                        indicatorsInverseInput.get().get(pidx).getValue(i)) {
                    double oldValue = param.getValue(i);
                    double newValue = oldValue/f;
                    if (newValue < param.getLower() || newValue > param.getUpper())
                        return Double.NEGATIVE_INFINITY;

                    param.setValue(i, newValue);
                    logHastingsRatio -= logf;
                }
            }
        }

        // Reject invalid tree scalings:
        if (f<1.0) {
            for (Node leaf : qsTree.getExternalNodes()) {
                Double[] tempqstimes=qsTree.getAttachmentTimesList(leaf.getNr());
                if (tempqstimes.length>1 && tempqstimes[0]<leaf.getHeight()){
                    System.out.println("problem in hereeeeee you did not really scale the QS start apparently");
                    System.exit(0);
                }
                if (leaf.getParent().getHeight()<leaf.getHeight() || tempqstimes[tempqstimes.length-1]<leaf.getHeight())
                    return Double.NEGATIVE_INFINITY;
            }
        }else {
            if (root.getHeight()>origin.getValue())
                return Double.NEGATIVE_INFINITY;
        }

        // Return Hastings ratio:
        return logHastingsRatio;
    }

}
