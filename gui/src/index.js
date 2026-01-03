import { CellSelector } from './CellSelector.js';
import GUIManager from './GUIManager.js';
import Visualization from './Visualization.js';
import { Cell } from './Cell.js';
import { examples } from './tasks.js';

// Export to global scope
window.Visualization = Visualization ;
window.GUIManager = GUIManager;

function filterExamplesByTaskType(taskType) {
    const filteredExamples = [];
    for (const [key, example] of Object.entries(examples)) {
        if (example.task === taskType) {
            filteredExamples.push({ value: key, label: example.title });
        }
    }
    return filteredExamples;
}

function updateExampleSelector() {
    const taskTypeSelect = document.getElementById('task-type');
    const exampleSelector = document.getElementById('example-selector');
    
    if (!taskTypeSelect || !exampleSelector) {
        console.error('Task type or example selector not found');
        return;
    }
    
    const selectedTaskType = taskTypeSelect.value;
    const filteredExamples = filterExamplesByTaskType(selectedTaskType);
    
    exampleSelector.innerHTML = '<option value="" selected>-- Select Example --</option>';
    
    filteredExamples.forEach(example => {
        const option = document.createElement('option');
        option.value = example.value;
        option.textContent = example.label;
        exampleSelector.appendChild(option);
    });
}

// Initialize GUI when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    console.log('Initializing MicroMagneticGUI...');
    
    updateExampleSelector();

    const taskTypeSelect = document.getElementById('task-type');
    if (taskTypeSelect) {
        taskTypeSelect.addEventListener('change', updateExampleSelector);
    }
    
    // Create and initialize visualization
    const containerId = 'visualization-container';
    if (document.getElementById(containerId)) {
        // Create visualization instance with container ID
        const visualization = new Visualization(containerId);
        
        // Create GUI manager instance
        const guiManager = new GUIManager(visualization);
        guiManager.initWebSocketClient();
        
        // Create relax task manager using GUIManager (default)
        const relaxTaskManager = guiManager.initRelaxTaskManager('#cells-container');
        
        // Initialize CellSelector
        const cellSelector = new CellSelector(relaxTaskManager);
        
        // Add event listener for add-cell-btn
        const addCellButton = document.getElementById('add-cell-btn');
        if (addCellButton) {
            addCellButton.addEventListener('click', () => {
                console.log('Add new cell button clicked');
                const selectedCell = guiManager.taskManager.getSelectedCell();
                if (selectedCell) {
                    // Insert new cell below the selected cell
                    const newCell = new Cell('');
                    guiManager.taskManager.insertCellBelowSelected(newCell);
                } else {
                    // If no cell is selected, just add it to the end
                    const newCell = new Cell('');
                    guiManager.taskManager.addCell(newCell);
                }
                cellSelector.updateOptions();
            });
        }
        
        // Add event listener for example selection to update CellSelector
        const exampleSelector = document.getElementById('example-selector');
        if (exampleSelector) {
            exampleSelector.addEventListener('change', () => {
                // Update CellSelector with the current taskManager
                if (guiManager.taskManager) {
                    cellSelector.updateCellManager(guiManager.taskManager);
                }
            });
        }
        
        console.log('MicroMagneticGUI initialized successfully!');
    } else {
        console.error('Visualization container not found');
    }
});