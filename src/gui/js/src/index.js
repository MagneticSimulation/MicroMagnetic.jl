import * as THREE from 'three';
import { GUI } from 'lil-gui';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CellManager } from './CellManager.js';
import { CellSelectorManager } from './CellSelectorManager.js';
import { createRelaxTaskManager, initTaskTemplates } from './tasks.js';
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
        
        // Add event listener for run-code button
        const runButton = document.getElementById('run-code');
        if (runButton) {
            runButton.addEventListener('click', () => {
                console.log('Run button clicked');
                
                if (guiManager.codeEditor) {
                    const code = guiManager.codeEditor.getValue();
                    guiManager.runCode(code);
                }
            });
        }
        
        // Add event listener for run-simulation-btn button
        const runSimulationButton = document.getElementById('run-simulation-btn');
        if (runSimulationButton) {
            runSimulationButton.addEventListener('click', () => {
                console.log('Run simulation button clicked');
                
                if (guiManager.codeEditor) {
                    const code = guiManager.codeEditor.getValue();
                    guiManager.runCode(code);
                }
            });
        }
        
        // Add event listeners for test buttons
        const testConnectionButton = document.getElementById('test-connection-btn');
        if (testConnectionButton) {
            testConnectionButton.addEventListener('click', () => {
                console.log('Test connection button clicked');
                
                // Update status message
                const statusElement = document.getElementById('status-message');
                if (statusElement) {
                    statusElement.textContent = 'Testing connection to Julia server...';
                }
                
                // Send test command
                guiManager.sendCommand('test_connection', {});
            });
        }
        
        const sendCommandButton = document.getElementById('send-command-btn');
        if (sendCommandButton) {
            sendCommandButton.addEventListener('click', () => {
                console.log('Send test command button clicked');
                
                // Update status message
                const statusElement = document.getElementById('status-message');
                if (statusElement) {
                    statusElement.textContent = 'Sending test command to Julia server...';
                }
                
                // Send test command with parameters
                guiManager.sendCommand('test_command', {
                    message: 'Hello from JavaScript!',
                    timestamp: Date.now(),
                    numbers: [1, 2, 3, 4, 5]
                });
            });
        }
        
        const getStatusButton = document.getElementById('get-status-btn');
        if (getStatusButton) {
            getStatusButton.addEventListener('click', () => {
                console.log('Get Julia status button clicked');
                
                // Update status message
                const statusElement = document.getElementById('status-message');
                if (statusElement) {
                    statusElement.textContent = 'Getting Julia server status...';
                }
                
                // Send status command
                guiManager.sendCommand('get_status', {});
            });
        }
        
        const stopServerButton = document.getElementById('stop-server-btn');
        if (stopServerButton) {
            stopServerButton.addEventListener('click', () => {
                console.log('Close connection button clicked');
                
                // Update status message
                const statusElement = document.getElementById('status-message');
                if (statusElement) {
                    statusElement.textContent = 'Closing connection...';
                }
                
                // Send stop command
                guiManager.sendCommand('stop_server', {});
            });
        }
        
        // Initialize task templates
        initTaskTemplates();
        
        // Create relax task manager and populate cells-container with relax task steps
        const relaxTaskManager = createRelaxTaskManager('#cells-container');
        
        // Initialize CellSelectorManager
        const cellSelectorManager = new CellSelectorManager(relaxTaskManager);
        
        // Add selection change listener to update cell selector options
        relaxTaskManager.addSelectionChangeListener((selectedCell) => {
            console.log('Selection changed, updating cell selector options:', selectedCell.name);
            cellSelectorManager.updateOptions();
        });
        
        // Update cell selector options initially
        cellSelectorManager.updateOptions();
        
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
                cellSelectorManager.updateOptions();
            });
        }
        
        console.log('MicroMagneticGUI initialized successfully!');
    } else {
        console.error('Visualization container not found');
    }
});