package beast.app.piqmee.beauti;

import beast.app.beauti.BeautiDoc;
import beast.app.beauti.BeautiPanelConfig;
import beast.app.beauti.TipDatesInputEditor;
import beast.app.draw.InputEditor;
import beast.core.BEASTInterface;
import beast.core.Input;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import beast.app.beauti.GuessPatternDialog;

import java.awt.event.ActionEvent;
import java.util.List;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.AbstractTableModel;
import beast.evolution.tree.Tree;
import piqmee.tree.QuasiSpeciesTree;


/**
 * BEAUti input editor for QuasiSpeciesTree haplotype copies count trait.
 *
 *  @author Veronika Boskova
 *  */
public class HaplotypeCountsTraitSetInputEditor extends TipDatesInputEditor {

    HaplotypeCountsTraitTableModel tableModel;
    TraitSet traitSet;
    TaxonSet taxonSet;

    public HaplotypeCountsTraitSetInputEditor(BeautiDoc doc) {
        super(doc);
    }

    @Override
    public Class<?> type() {
        return QuasiSpeciesTree.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {

        if (plugin.getClass().equals("beast.app.beauti.BeautiConfig") && !((BeautiPanelConfig) plugin).getName().equals("Sequence Counts")){
            super.init(input, plugin, itemNr, bExpandOption, bAddButtons);
        } else {

            traitSet = ((QuasiSpeciesTree) input.get()).getHaplotypeCountsTrait();
            Tree tree = null;
            if (itemNr >= 0) {
                tree = (Tree) ((List<?>) input.get()).get(itemNr);
            } else {
                tree = (Tree) input.get();
            }
            traitSet.setID("haplotypeCounts.t:" + BeautiDoc.parsePartition(tree.getID()));
            taxonSet = traitSet.taxaInput.get();
            tableModel = new HaplotypeCountsTraitTableModel(traitSet);
            JTable table = new JTable(tableModel);
            JButton guessButton = new JButton("Guess");
            guessButton.addActionListener((ActionEvent e) -> {
                GuessPatternDialog dlg = new GuessPatternDialog(null, ".*(\\d\\d\\d\\d).*");

                String traitString = "";
                switch (dlg.showDialog("Guess counts")) {
                    case canceled:
                        return;
                    case trait:
                        traitString = dlg.getTrait();
                        break;
                    case pattern:
                        StringBuilder traitStringBuilder = new StringBuilder();
                        for (String taxonName : taxonSet.asStringList()) {
                            String matchString = dlg.match(taxonName);
                            if (matchString == null)
                                return;

                            if (traitStringBuilder.length() > 0)
                                traitStringBuilder.append(",");

                            traitStringBuilder.append(taxonName).append("=").append(matchString);
                        }
                        traitString = traitStringBuilder.toString();
                        break;
                }
                traitSet.traitsInput.setValue(traitString, traitSet);
                try {
                    traitSet.initAndValidate();
                } catch (Exception ex) {
                    System.err.println("Error setting type trait.");
                }
                refreshPanel();
            });

            JButton clearButton = new JButton("Clear");
            clearButton.addActionListener((ActionEvent e) -> {
                StringBuilder traitStringBuilder = new StringBuilder();
                for (String taxonName : taxonSet.asStringList()) {
                    if (traitStringBuilder.length() > 0)
                        traitStringBuilder.append(",");
                    traitStringBuilder.append(taxonName).append("=0");
                }
                traitSet.traitsInput.setValue(traitStringBuilder.toString(), traitSet);
                try {
                    traitSet.initAndValidate();
                } catch (Exception ex) {
                    System.err.println("Error clearing haplotype count trait.");
                }
                refreshPanel();
            });

            Box boxVert = Box.createVerticalBox();

            Box boxHoriz = Box.createHorizontalBox();
            boxHoriz.add(Box.createHorizontalGlue());
            boxHoriz.add(guessButton);
            boxHoriz.add(clearButton);
            boxVert.add(boxHoriz);
            boxVert.add(new JScrollPane(table));

            add(boxVert);
        }
    }

    class HaplotypeCountsTraitTableModel extends AbstractTableModel {

        TraitSet haplotypeCountsTraitSet;

        public HaplotypeCountsTraitTableModel(TraitSet haplotypeCountsTraitSet) {
            this.haplotypeCountsTraitSet = haplotypeCountsTraitSet;
        }

        @Override
        public int getRowCount() {
            return haplotypeCountsTraitSet.taxaInput.get().getTaxonCount();
        }

        @Override
        public int getColumnCount() {
            return 2;
        }

        @Override
        public Object getValueAt(int rowIndex, int columnIndex) {
            if (columnIndex<0 || columnIndex>=getRowCount())
                return null;

            switch(columnIndex) {
                case 0:
                    // Taxon name:
                    return haplotypeCountsTraitSet.taxaInput.get().getTaxonId(rowIndex);
                case 1:
                    // Type:
                    return haplotypeCountsTraitSet.getStringValue(rowIndex);
                default:
                    return null;
            }
        }

        @Override
        public boolean isCellEditable(int rowIndex, int columnIndex) {
            return columnIndex == 1;
        }

        @Override
        public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
            String taxon = taxonSet.getTaxonId(rowIndex);
            String traitString = traitSet.traitsInput.get();
            int startIdx = traitString.indexOf(taxon + "=");
            int endIdx = traitString.indexOf(",", startIdx);

            if (Double.valueOf((String) aValue) < 1) {
                System.err.println("Error setting haplotype counts trait value. Make sure to input counts, not percentages");
            }

            String newTraitString = traitString.substring(0, startIdx);
            newTraitString += taxon + "=" + (String)aValue;
            if (endIdx>=0)
                newTraitString += traitString.substring(endIdx);

            traitSet.traitsInput.setValue(newTraitString, traitSet);
            try {
                traitSet.initAndValidate();
            } catch (Exception ex) {
                System.err.println("Error setting haplotype counts trait value.");
            }

            fireTableCellUpdated(rowIndex, columnIndex);
        }

        @Override
        public String getColumnName(int column) {
            switch(column) {
                case 0:
                    return "Name";
                case 1:
                    return "Count";
                default:
                    return null;
            }
        }
    }
}