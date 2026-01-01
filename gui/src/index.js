import { CellSelector } from './CellSelector.js';
import GUIManager from './GUIManager.js';
import Visualization from './Visualization.js';
import { Cell } from './Cell.js';

// Export to global scope
window.Visualization = Visualization ;
window.GUIManager = GUIManager;

// Initialize GUI when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    console.log('Initializing MicroMagneticGUI...');
    
    // Create and initialize visualization
    const containerId = 'visualization-container';
    if (document.getElementById(containerId)) {
        // Create visualization instance with container ID
        const visualization = new Visualization(containerId);
        
        // Create GUI manager instance
        const guiManager = new GUIManager(visualization);
        guiManager.initWebSocketClient();
        
        // Create relax task manager using GUIManager
        const relaxTaskManager = guiManager.initRelaxTaskManager('#cells-container');
        
        // Initialize CellSelector
        const cellSelector = new CellSelector(relaxTaskManager);
        
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