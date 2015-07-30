package quasispeciestree.operators;

/**
 * @author Veronika Boskova created on 15/07/2015.
 */
public class QuasiSpeciesWilsonBalding {
}



/*
* TODO: need to implement check for whether the proposed branching times for each duplicate of a haplotype are
       in accordance with "appear-only-once-in-tree-history" assumption --- where to do it?  -- in operator? to be efficient
       -- randomly choose a haplotype, propose its attachment times from root to tip, if another haplotype attaches to the path (or below),
       where the QS of the first haplotype chose to attach, then allow these to only propose within attachmet node to first haplotype
       to the tip of the second haplotype, or to the tip of the haplotypes below.

       Also check in tree operators whether the QS starting/attachment times are still compatible, if not propose new QS times...
 *
 */


// import beast.util.Randomizer;
//    /**
//     * Use beast RNG to select random node from list.
//     *
//     * @param nodeList
//     * @return A randomly selected node.
//     */
//    private MultiTypeNode selectRandomNode(List<MultiTypeNode> nodeList) {
//        return nodeList.get(Randomizer.nextInt(nodeList.size()));
//    }