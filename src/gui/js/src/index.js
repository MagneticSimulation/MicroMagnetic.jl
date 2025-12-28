import { CellSelector } from './CellSelector.js';
import GUIManager from './GUIManager.js';
import Visualization from './Visualization.js';

// Export to global scope
window.Visualization = Visualization ;
window.GUIManager = GUIManager;

// Initialize GUI when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    console.log('Initializing MicroMagneticGUI...');
    
    // Create and initialize visualization
    const container = document.getElementById('visualization-container');
    if (container) {
        // Create visualization instance
        const visualization = new Visualization();
        visualization.init(document.getElementById('visualization-container'));
        
        // Create GUI manager instance
        const guiManager = new GUIManager(visualization);
        guiManager.initWebSocketClient();
        
        // Create relax task manager using GUIManager
        const relaxTaskManager = guiManager.initRelaxTaskManager('#cells-container');
        
        // Initialize CellSelector
        const cellSelector = new CellSelector(relaxTaskManager);
        
        // Update cell selector options when cell type changes
        const cellSelectorElement = document.getElementById('cell-selector');
        if (cellSelectorElement) {
            cellSelectorElement.addEventListener('change', () => {
                cellSelector.updateOptions();
            });
        }
        
        // Update cell selector options when cell is selected
        const cellsContainerElement = document.getElementById('cells-container');
        if (cellsContainerElement) {
            cellsContainerElement.addEventListener('click', (event) => {
                if (event.target.closest('.cell')) {
                    cellSelector.updateOptions();
                }
            });
        }
        
        // Add selection change listener to update cell selector options
        relaxTaskManager.addSelectionChangeListener((selectedCell) => {
            console.log('Selection changed, updating cell selector options:', selectedCell.name);
            cellSelector.updateOptions();
        });
        
        // Update cell selector options initially
        cellSelector.updateOptions();
        
        // Add event listener for add-cell-btn
        const addCellButton = document.getElementById('add-cell-btn');
        if (addCellButton) {
            addCellButton.addEventListener('click', () => {
                console.log('Add new cell button clicked');
                const selectedCell = relaxTaskManager.getSelectedCell();
                if (selectedCell) {
                    // Insert new cell below the selected cell
                    const newCell = new Cell('');
                    relaxTaskManager.insertCellBelowSelected(newCell);
                } else {
                    // If no cell is selected, just add it to the end
                    const newCell = new Cell('');
                    relaxTaskManager.addCell(newCell);
                }
                cellSelector.updateOptions();
            });
        }
        
        console.log('MicroMagneticGUI initialized successfully!');
    } else {
        console.error('Visualization container not found');
    }
});