/*
 * Copyright (C) 2015 Tim Vaughan (tgvaughan@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package piqmee.app.beauti;


import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TraitSet;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BeautiPanelConfig;
import beastfx.app.inputeditor.GuessPatternDialog;
import beastfx.app.inputeditor.TipDatesInputEditor;
import beastfx.app.util.Alert;
import beastfx.app.util.FXUtils;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.control.Button;
import javafx.scene.control.SelectionMode;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import piqmee.tree.QuasiSpeciesTree;

/**
 * Adapted from BEAUti input editor for MultiTypeTree/Mascot type traits.
 *
 * @author Tim Vaughan (tgvaughan@gmail.com)
 * @author Walter Xie
 * @author Jordan Douglas
 */
public class HaplotypeCountsInputEditor extends TipDatesInputEditor {

    TraitSet traitSet;
    TaxonSet taxonSet;
    TableView<HaploCount> table;
    private ObservableList<HaploCount> counts;

    public HaplotypeCountsInputEditor(BeautiDoc doc) {
        super(doc);

    }

    @Override
    public Class<?> type() {
        return QuasiSpeciesTree.class;
    }

    @Override
    public void init(Input<?> input, BEASTInterface plugin, int itemNr, ExpandOption bExpandOption, boolean bAddButtons) {
        if (plugin.getClass().getName().equals("beastfx.app.inputeditor.BeautiPanelConfig") && 
        	!((BeautiPanelConfig) plugin).getName().equals("Sequence Counts")){
        	// behave like a normal TipDatesInputEditor
            super.init(input, plugin, itemNr, bExpandOption, bAddButtons);
        } else {

            traitSet = ((QuasiSpeciesTree) input.get()).getHaplotypeCountsTrait();

	    	counts = FXCollections.observableArrayList();
	    	
	        taxonSet = traitSet.taxaInput.get();
	
	        table = new TableView<>();
	        table.setEditable(true);
	        table.setPrefWidth(1024);
	        table.getSelectionModel().setSelectionMode(SelectionMode.MULTIPLE);
	        table.setItems(counts);
	        
	        TableColumn<HaploCount, String> column1 = new TableColumn<>("Name");
	        column1.setCellValueFactory(new PropertyValueFactory<>("name"));
	        column1.setPrefWidth(500);
	        column1.setEditable(false);
	
	        TableColumn<HaploCount, String> column2 =new TableColumn<>("Count");
	        column2.setPrefWidth(500);
	        column2.setEditable(true);
	        column2.setCellValueFactory(new PropertyValueFactory<>("count"));
	        column2.setCellFactory(TextFieldTableCell.forTableColumn());
	        
	        column2.setOnEditCommit(
	                t -> t.getTableView().getItems().
	                        get(t.getTablePosition().getRow()).
	                        setCount(t.getNewValue())
	        );
	
	        
	        table.getColumns().add(column1);
	        table.getColumns().add(column2);
	
	        traitSetToCounts();
	
	        Button guessButton = new Button("Guess");
	        guessButton.setOnAction(e -> {
	            GuessPatternDialog dlg = new GuessPatternDialog(null,
	                ".*(\\d\\d\\d\\d).*");
	            
	            String traitString = "";
	            switch(dlg.showDialog("Guess counts")) {
	                case canceled:
	                    return;
	                case trait:
	                    traitString = dlg.getTrait(); //TODO replace to what?
	                    break;
	                case pattern:
	                    StringBuilder traitStringBuilder = new StringBuilder();
	                    for (String taxonName : taxonSet.asStringList()) {
	                        String matchString = dlg.match(taxonName);
	                        if (matchString == null)
	                            return;
	                        
	                        if (traitStringBuilder.length()>0)
	                            traitStringBuilder.append(",");
	                        
	                        traitStringBuilder.append(taxonName)
	                            .append("=")
	                            .append(matchString);
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
	
	            traitSetToCounts();
	
	        });
	
	        Button clearButton = new Button("Clear");
	        clearButton.setOnAction(e -> {
	//            StringBuilder traitStringBuilder = new StringBuilder();
	//            for (String taxonName : taxonSet.asStringList()) {
	//                if (traitStringBuilder.length()>0)
	//                    traitStringBuilder.append(",");
	//                traitStringBuilder.append(taxonName).append("=0");
	//            }
	//            traitSet.traitsInput.setValue(traitStringBuilder.toString(), traitSet);
	//            try {
	//                traitSet.initAndValidate();
	//            } catch (Exception ex) {
	//                System.err.println("Error clearing type trait.");
	//            }
	
	            saveTraitSet(true);
	
	            traitSetToCounts();
	
	        });
	
	        VBox boxVert = FXUtils.newVBox();
	
	        HBox boxHoriz = FXUtils.newHBox();
	        boxHoriz.getChildren().add(guessButton);
	        boxHoriz.getChildren().add(clearButton);
	
	        boxVert.getChildren().add(boxHoriz);
	        boxVert.getChildren().add(table);
	        
	        this.pane = FXUtils.newHBox();
	        getChildren().add(pane);
	
	        pane.getChildren().add(boxVert);
        }
    }

    /**
     * Save {@link HaploCount} to {@link TraitSet#traitsInput}.
     * @param toClear  true to set trait to 0, false to set to count
     */
    private void saveTraitSet(boolean toClear) {
        StringBuilder traitStringBuilder = new StringBuilder();
        for (HaploCount c : counts) {
            if (traitStringBuilder.length()>0) traitStringBuilder.append(",");
            
            String count = c.getCount();
            String taxon = c.getName();
            traitStringBuilder.append(taxon).append("=");
            if (toClear)
                traitStringBuilder.append("0");
            else
                traitStringBuilder.append(count);
        }
        traitSet.traitsInput.setValue(traitStringBuilder.toString(), traitSet);
        try {
            traitSet.initAndValidate();
        } catch (Exception ex) {
            String action = toClear ? "clearing" : "setting";
            System.err.println("Error " + action + " type trait.");
            Alert.showMessageDialog(this,
                    "Error " + action + " type trait : " + traitSet.getTraitName() + " !",
                    "Type Trait Editor", Alert.ERROR_MESSAGE);
        }
    }

    /**
     * Import {@link TraitSet} values into {@link HaploCount}.
     */
    private void traitSetToCounts() {
        counts.clear();
        // add data
        for (String taxonName : taxonSet.asStringList()) {
            String loc = traitSet.getStringValue(taxonName);
            counts.add(new HaploCount(taxonName, loc));
        }
        
        table.refresh();
        // refreshPanel();
    }

    /**
     * Table model
     */
    public class HaploCount {
        String name;
        String count;

        public HaploCount(String name, String count) {
            this.name = name;
            this.count = count;
        }

        public String getName() {
            return name;
        }

        public void setName(String name) {
            this.name = name;
        }

        public String getCount() {
            return count;
        }

        public void setCount(String location) {
            this.count = location;
            saveTraitSet(false);
        }
    }

}
