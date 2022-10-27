package test.piqmee;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import piqmee.tree.QuasiSpeciesTree;
import piqmee.tree.QuasiSpeciesTreeFromNewick;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Veronika Boskova created on 21/02/2018
 */
public class QuasiSpeciesTestCase {

    static public Alignment getAlignment(String[] data) {
        List<Sequence> seqList = new ArrayList<Sequence>();

        for (int i = 0; i < data.length; i++) {
            String taxonID = "t" + (i);
            seqList.add(new Sequence(taxonID, data[i]));
        }

        Alignment alignment = new Alignment(seqList, "nucleotide");
        return alignment;
    }

    static public QuasiSpeciesTree setTreeFromFullNewick(String inputTree, String[] data) {
        int Nleaves = data.length;
        StringBuilder traitSB = new StringBuilder();
        List<Sequence> seqList = new ArrayList<Sequence>();

        for (int i = 0; i < Nleaves; i++) {
            String taxonID = "t" + (i);
            seqList.add(new Sequence(taxonID, data[i]));

            if (i > 0)
                traitSB.append(",");
            traitSB.append(taxonID).append("=").append(i+1);
        }

        Alignment alignment = new Alignment(seqList, "nucleotide");
        TaxonSet taxonSet = new TaxonSet(alignment);

        QuasiSpeciesTree tree = new QuasiSpeciesTreeFromNewick();

        tree.setInputValue("newick", inputTree);
        tree.setInputValue("adjustTipHeights", "false");
        tree.setInputValue("taxonset", taxonSet);
        tree.setInputValue("data", alignment);
        tree.initAndValidate();

        return tree;
    }

    static public QuasiSpeciesTree setTreeFromNewick(String inputTree, String[] data) {
        int Nleaves = data.length;
        StringBuilder traitSB = new StringBuilder();
        List<Sequence> seqList = new ArrayList<Sequence>();

        for (int i = 0; i < Nleaves; i++) {
            String taxonID = "t" + (i);
            seqList.add(new Sequence(taxonID, data[i]));

            if (i > 0)
                traitSB.append(",");
            traitSB.append(taxonID).append("=").append(i+1);
        }

        Alignment alignment = new Alignment(seqList, "nucleotide");
        TaxonSet taxonSet = new TaxonSet(alignment);
        TraitSet haploCounts = new TraitSet();

        haploCounts.initByName(
                "traitname", "qscounts",
                "taxa", taxonSet,
                "value", traitSB.toString());

        QuasiSpeciesTree tree = new QuasiSpeciesTreeFromNewick();

        tree.setInputValue("newick", inputTree);
        tree.setInputValue("adjustTipHeights", "false");
        tree.setInputValue("taxonset", taxonSet);
        tree.setInputValue("haplotypeCounts",haploCounts);
        tree.setInputValue("data",alignment);
        tree.initAndValidate();

        return tree;
    }
}
