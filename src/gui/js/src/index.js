import * as THREE from 'three';
import { GUI } from 'lil-gui';
import WebSocketClient from './client.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CellManager } from './CellManager.js';
import { CellSelectorManager } from './CellSelectorManager.js';
import { createRelaxTaskManager, initTaskTemplates } from './tasks.js';

/**
 * MagneticVisualization class for visualizing magnetization distributions
 */
class MagneticVisualization {
    constructor() {
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.controls = null;
        this.cube = null;
        this.animationId = null;
        this.arrows = [];
        this.gridSize = [10, 10, 10];
        this.dimensions = [10, 10, 10];
        this.arrowGroup = new THREE.Group();
        this.arrowPositions = null;
        this.webSocketClient = null;
        this.wsClient = null;
        this.isConnected = false;
        this.connectionAttempts = 0;
        this.maxConnectionAttempts = 5;
    }

    /**
     * Initialize visualization scene
     * @param {Object} container - DOM element containing renderer
     * @param {Object} data - Initialization data
     */
    init(container, data = null) {
        // Create scene
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0xf0f0f0);

        // Create camera
        const width = container.clientWidth;
        const height = container.clientHeight;
        this.camera = new THREE.PerspectiveCamera(75, width / height, 0.1, 1000);
        this.camera.position.z = 5;

        // Create renderer
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(width, height);
        container.innerHTML = '';
        container.appendChild(this.renderer.domElement);

        // Add helpers
        this.scene.add(new THREE.GridHelper(10, 10));
        this.scene.add(new THREE.AxesHelper(5));

        // Create cube
        const geometry = new THREE.BoxGeometry();
        const material = new THREE.MeshStandardMaterial({ 
            color: 0x0077ff,
            metalness: 0.3,
            roughness: 0.4,
            transparent: true,
            opacity: 0.3
        });
        this.cube = new THREE.Mesh(geometry, material);
        this.scene.add(this.cube);

        // Add lights
        this.scene.add(new THREE.AmbientLight(0xffffff, 0.5));
        
        const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
        directionalLight.position.set(5, 5, 5);
        this.scene.add(directionalLight);

        // Add arrow group
        this.scene.add(this.arrowGroup);

        // Add orbit controls
        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;

        // Initialize Julia communicator
        this.initWebSocketClient();

        // Initialize magnetization if data provided
        if (data) {
            this.updateMagnetization(data);
        }

        // Start animation
        this.animate();

        console.log('MagneticVisualization initialized successfully!');
    }

    /**
     * Initialize WebSocket client for communicating with Julia server
     */
    initWebSocketClient() {
        try {
            console.log('Initializing WebSocket client');
            
            // Create WebSocket client instance
            this.webSocketClient = new WebSocketClient({
                autoConnect: false // Disable auto-connect to allow manual control
            });
            
            // Initialize WebSocket UI elements
            this.initWebSocketUI();
            
            // Set up event listeners
            this.webSocketClient.on('connect', () => {
                console.log('WebSocket connected to Julia server');
                this.updateWebSocketStatus(true);
                
                const statusElement = document.getElementById('status-message');
                if (statusElement) {
                    statusElement.textContent = 'Connected to Julia server';
                    statusElement.className = 'status-connected';
                }
            });
            
            this.webSocketClient.on('disconnect', () => {
                console.log('WebSocket disconnected from Julia server');
                this.updateWebSocketStatus(false);
                
                const statusElement = document.getElementById('status-message');
                if (statusElement) {
                    statusElement.textContent = 'Disconnected from Julia server';
                    statusElement.className = 'status-disconnected';
                }
            });
            
            this.webSocketClient.on('error', (error) => {
                console.error('WebSocket error:', error);
                
                const statusElement = document.getElementById('status-message');
                if (statusElement) {
                    statusElement.textContent = `Error: ${error.message}`;
                    statusElement.className = 'status-error';
                }
            });
            
            // Set up message handlers
            this.webSocketClient.on('simulation_status', (data) => {
                this.handleJuliaMessage('simulation_status', data);
            });
            
            this.webSocketClient.on('magnetization_data', (data) => {
                this.handleJuliaMessage('magnetization_data', data);
            });
            
            this.webSocketClient.on('command_response', (data) => {
                this.handleJuliaMessage('command_response', data);
            });
            
            this.webSocketClient.on('error', (data) => {
                this.handleJuliaMessage('error', data);
            });
            
        } catch (error) {
            console.error('Failed to initialize WebSocket client:', error);
            
            const statusElement = document.getElementById('status-message');
            if (statusElement) {
                statusElement.textContent = `Error: ${error.message}`;
                statusElement.className = 'status-error';
            }
        }
    }

    /**
     * Initialize WebSocket UI elements and event listeners
     */
    initWebSocketUI() {
        // Get UI elements
        const statusIndicator = document.getElementById('ws-status-indicator');
        const statusText = document.getElementById('ws-status-text');
        const connectBtn = document.getElementById('ws-connect-btn');
        
        if (!statusIndicator || !statusText || !connectBtn) {
            console.error('WebSocket UI elements not found');
            return;
        }
        
        // Initialize status
        this.updateWebSocketStatus(false);
        
        // Add connect button event listener
        connectBtn.addEventListener('click', async () => {
            try {
                if (this.webSocketClient && this.webSocketClient.isConnected()) {
                    // Disconnect
                    this.webSocketClient.disconnect();
                } else {
                    // Connect
                    await this.webSocketClient.connect();
                }
            } catch (error) {
                console.error('WebSocket connection error:', error);
                this.updateWebSocketStatus(false);
                
                const statusElement = document.getElementById('status-message');
                if (statusElement) {
                    statusElement.textContent = `Connection error: ${error.message}`;
                    statusElement.className = 'status-error';
                }
            }
        });
        
        // Initialize simulation type and task type selections
        this.initSimulationTypeSelection();
        this.initTaskTypeSelection();
    }

    /**
     * Update WebSocket status display
     * @param {boolean} isConnected - Connection status
     */
    updateWebSocketStatus(isConnected) {
        const statusIndicator = document.getElementById('ws-status-indicator');
        const statusText = document.getElementById('ws-status-text');
        const connectBtn = document.getElementById('ws-connect-btn');
        
        if (!statusIndicator || !statusText || !connectBtn) {
            return;
        }
        
        // Update status indicator
        if (isConnected) {
            statusIndicator.className = 'status-indicator connected';
            statusText.textContent = 'Connected';
            connectBtn.className = 'connect-btn connected';
            connectBtn.textContent = 'Disconnect';
        } else {
            statusIndicator.className = 'status-indicator disconnected';
            statusText.textContent = 'Disconnected';
            connectBtn.className = 'connect-btn disconnected';
            connectBtn.textContent = 'Connect';
        }
    }

    /**
     * Initialize simulation type selection
     */
    initSimulationTypeSelection() {
        const selectElement = document.getElementById('simulation-type');
        
        selectElement.addEventListener('change', (event) => {
            const selectedType = event.target.value;
            console.log('Selected simulation type:', selectedType);
            // Store selected simulation type
            this.simulationType = selectedType;
        });
        
        // Set default simulation type
        this.simulationType = 'fd';
    }

    /**
     * Initialize task type selection
     */
    initTaskTypeSelection() {
        const selectElement = document.getElementById('task-type');
        
        selectElement.addEventListener('change', (event) => {
            const selectedType = event.target.value;
            console.log('Selected task type:', selectedType);
            // Store selected task type
            this.taskType = selectedType;
        });
        
        // Set default task type
        this.taskType = 'relax';
    }

    /**
     * Handle message from Julia server
     * @param {string} type - Message type
     * @param {Object} data - Message data
     */
    handleJuliaMessage(type, data) {
        console.log('Received message from Julia:', { type, data });
        
        switch (type) {
            case 'simulation_status':
                this.updateSimulationStatus(data);
                break;
            case 'magnetization_data':
                this.updateMagnetization(data);
                break;
            case 'command_response':
                this.handleCommandResponse(data);
                break;
            case 'error':
                this.handleError(data);
                break;
            default:
                console.warn('Unknown message type:', type);
        }
    }

    /**
     * Update simulation status display
     * @param {Object} status - Simulation status
     */
    updateSimulationStatus(status) {
        console.log('Simulation status:', status);
        
        // Update status message
        const statusElement = document.getElementById('status-message');
        if (statusElement) {
            if (status.has_sim) {
                statusElement.textContent = `Status: ${status.status}\nRunning: ${status.is_running}`;
            } else {
                statusElement.textContent = 'Status: No simulation created';
            }
        }
    }

    /**
     * Handle command response from Julia server
     * @param {Object} response - Command response
     */
    handleCommandResponse(response) {
        console.log('Command response:', response);
        
        // Update status message
        const statusElement = document.getElementById('status-message');
        if (statusElement) {
            if (response.success) {
                statusElement.textContent = `Command "${response.command}" executed successfully\nResult: ${JSON.stringify(response.result)}`;
            } else {
                statusElement.textContent = `Error executing command "${response.command}":\n${response.error}`;
            }
        }
    }

    /**
     * Handle error from Julia server
     * @param {Object} error - Error message
     */
    handleError(error) {
        console.error('Julia error:', error);
        
        // Update status message
        const statusElement = document.getElementById('status-message');
        if (statusElement) {
            statusElement.textContent = `Error: ${error.message}`;
        }
    }

    /**
     * Send command to Julia server
     * @param {string} command - Command to execute
     * @param {Object} params - Command parameters
     */
    sendCommand(command, params = {}) {
        console.log(`Sending command to Julia: ${command}`, params);
        
        if (this.webSocketClient) {
            this.webSocketClient.sendCommand(command, params).then(response => {
                this.handleCommandResponse(response);
            });
        }
    }
    /**
     * Animation loop
     */
    animate() {
        this.animationId = requestAnimationFrame(() => this.animate());

        if (this.renderer && this.scene && this.camera) {
            // Update controls
            if (this.controls) {
                this.controls.update();
            }
            
            this.renderer.render(this.scene, this.camera);
        }
    }

    /**
     * Update magnetization data
     * @param {Object} data - Magnetization data
     * @param {Object} options - Visualization options
     */
    updateMagnetization(data, options = {}) {
        console.log('Updating magnetization data:', data);
        console.log('Options:', options);
    
        // Update grid and dimensions
        let shouldRecalculatePositions = false;
        if (data.cells) {
            if (!this.gridSize || JSON.stringify(this.gridSize) !== JSON.stringify(data.cells)) {
                this.gridSize = data.cells;
                shouldRecalculatePositions = true;
            }
        }
        if (data.dimensions) {
            if (!this.dimensions || JSON.stringify(this.dimensions) !== JSON.stringify(data.dimensions)) {
                this.dimensions = data.dimensions;
                shouldRecalculatePositions = true;
            }
        }
    
        // Calculate cell size
        const cellSize = [
            this.dimensions[0] / this.gridSize[0],
            this.dimensions[1] / this.gridSize[1],
            this.dimensions[2] / this.gridSize[2]
        ];
    
        // Calculate arrow scale
        const baseScale = Math.min(...cellSize) * 0.8;
        const arrowScale = baseScale * (options.arrowScaleFactor || 1.0);
    
        // Filter cells by selection
        const selection = options.selection || { type: 'full' };
    
        // Collect arrow data
        const arrowData = [];
        const nx = this.gridSize[0];
        const ny = this.gridSize[1];
        const nz = this.gridSize[2];
    
        // Recalculate positions if needed
        if (shouldRecalculatePositions || !this.arrowPositions) {
            console.log('Recalculating arrow positions');
            this.arrowPositions = [];
            
            for (let i = 0; i < nx; i++) {
                for (let j = 0; j < ny; j++) {
                    for (let k = 0; k < nz; k++) {
                        if (!this.isInSelection(i, j, k, selection)) continue;

                        const index = i * ny * nz + j * nz + k;
                        const magnetization = data.magnetization[index];

                        if (magnetization) {
                            // Calculate arrow position (cell center)
                            const x = (i - (nx - 1) / 2) * cellSize[0];
                            const y = (j - (ny - 1) / 2) * cellSize[1];
                            const z = (k - (nz - 1) / 2) * cellSize[2];
                            
                            this.arrowPositions.push([x, y, z]);
                            
                            arrowData.push({
                                position: [x, y, z],
                                direction: magnetization
                            });
                        }
                    }
                }
            }
        } else {
            // Use stored positions, update directions only
            console.log('Using stored arrow positions, updating directions');
            let positionIndex = 0;
            
            for (let i = 0; i < nx; i++) {
                for (let j = 0; j < ny; j++) {
                    for (let k = 0; k < nz; k++) {
                        if (!this.isInSelection(i, j, k, selection)) continue;

                        const index = i * ny * nz + j * nz + k;
                        const magnetization = data.magnetization[index];

                        if (magnetization && positionIndex < this.arrowPositions.length) {
                            arrowData.push({
                                position: this.arrowPositions[positionIndex],
                                direction: magnetization
                            });
                            positionIndex++;
                        }
                    }
                }
            }
        }
    
        console.log(`Processing ${arrowData.length} arrows`);
    
        // Handle arrow count changes
        const currentArrowCount = this.arrows.length > 0 ? this.arrows[0].count : 0;
        
        if (currentArrowCount !== arrowData.length) {
            console.log(`Arrow count changed: clearing ${currentArrowCount} arrows, creating ${arrowData.length} new ones`);
            this.clearArrows();
            if (arrowData.length > 0) {
                this.createArrowInstances(arrowData, arrowScale);
            }
        } else if (arrowData.length > 0) {
            console.log(`Updating ${arrowData.length} existing arrow instances`);
            this.updateArrowInstances(arrowData, arrowScale);
        }
    }

    /**
     * Create arrow instances using InstancedMesh
     * @param {Array} arrowData - Arrow data
     * @param {number} arrowScale - Scale factor
     */
    createArrowInstances(arrowData, arrowScale) {
        // Create geometries
        const coneGeometry = new THREE.ConeGeometry(0.05, 0.2, 32);
        coneGeometry.translate(0, -0.2, 0);
        
        const cylinderGeometry = new THREE.CylinderGeometry(0.01, 0.01, 0.2, 32);
        cylinderGeometry.translate(0, -0.2, 0);
    
        // Create material
        const material = new THREE.MeshStandardMaterial({ 
            color: 0x0077ff,
            metalness: 0.3,
            roughness: 0.4
        });
    
        // Create instanced meshes
        const coneMesh = new THREE.InstancedMesh(coneGeometry, material, arrowData.length);
        const cylinderMesh = new THREE.InstancedMesh(cylinderGeometry, material, arrowData.length);
    
        // Set up matrices and colors
        const coneMatrix = new THREE.Matrix4();
        const cylinderMatrix = new THREE.Matrix4();
        const color = new THREE.Color();
        const up = new THREE.Vector3(0, 1, 0);
    
        for (let i = 0; i < arrowData.length; i++) {
            const data = arrowData[i];
            const position = new THREE.Vector3(...data.position);
            const direction = new THREE.Vector3(...data.direction).normalize();

            // Calculate rotation
            const quaternion = new THREE.Quaternion();
            quaternion.setFromUnitVectors(up, direction);

            // Calculate positioning
            const arrowDirection = direction.clone();
            const totalArrowLength = 0.4 * arrowScale;
            const offset = arrowDirection.clone().multiplyScalar(totalArrowLength * 0.5);
            
            // Position cylinder
            const cylinderPosition = position.clone().sub(offset)
                .add(arrowDirection.clone().multiplyScalar(0.2 * arrowScale));
            cylinderMatrix.compose(cylinderPosition, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            cylinderMesh.setMatrixAt(i, cylinderMatrix);

            // Position cone
            const conePosition = position.clone().sub(offset)
                .add(arrowDirection.clone().multiplyScalar(0.4 * arrowScale));
            coneMatrix.compose(conePosition, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            coneMesh.setMatrixAt(i, coneMatrix);

            // Set color based on direction
            color.setRGB(Math.abs(direction.x), Math.abs(direction.y), Math.abs(direction.z));
            coneMesh.setColorAt(i, color);
            cylinderMesh.setColorAt(i, color);
        }
    
        // Add to scene
        this.arrowGroup.add(coneMesh);
        this.arrowGroup.add(cylinderMesh);
    
        // Store references
        this.arrows.push(coneMesh, cylinderMesh);
    }

    /**
     * Update existing arrow instances
     * @param {Array} arrowData - Arrow data
     * @param {number} arrowScale - Scale factor
     */
    updateArrowInstances(arrowData, arrowScale) {
        if (this.arrows.length < 2) {
            console.error('Not enough arrow meshes');
            return;
        }
    
        const coneMesh = this.arrows[0];
        const cylinderMesh = this.arrows[1];
        const coneMatrix = new THREE.Matrix4();
        const cylinderMatrix = new THREE.Matrix4();
        const color = new THREE.Color();
        const up = new THREE.Vector3(0, 1, 0);
    
        for (let i = 0; i < arrowData.length; i++) {
            const data = arrowData[i];
            const position = new THREE.Vector3(...data.position);
            const direction = new THREE.Vector3(...data.direction).normalize();

            // Calculate rotation
            const quaternion = new THREE.Quaternion();
            quaternion.setFromUnitVectors(up, direction);

            // Calculate positioning
            const arrowDirection = direction.clone();
            const totalArrowLength = 0.6 * arrowScale;
            const offset = arrowDirection.clone().multiplyScalar(totalArrowLength * 0.5);
            
            // Update cylinder
            const cylinderPosition = position.clone().sub(offset)
                .add(arrowDirection.clone().multiplyScalar(0.2 * arrowScale));
            cylinderMatrix.compose(cylinderPosition, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            cylinderMesh.setMatrixAt(i, cylinderMatrix);

            // Update cone
            const conePosition = position.clone().sub(offset)
                .add(arrowDirection.clone().multiplyScalar(0.4 * arrowScale));
            coneMatrix.compose(conePosition, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            coneMesh.setMatrixAt(i, coneMatrix);

            // Update color
            color.setRGB(Math.abs(direction.x), Math.abs(direction.y), Math.abs(direction.z));
            coneMesh.setColorAt(i, color);
            cylinderMesh.setColorAt(i, color);
        }
    
        // Mark for update
        coneMesh.instanceMatrix.needsUpdate = true;
        coneMesh.instanceColor.needsUpdate = true;
        cylinderMesh.instanceMatrix.needsUpdate = true;
        cylinderMesh.instanceColor.needsUpdate = true;
    }

    /**
     * Check if cell is in selection area
     * @param {number} i - X index
     * @param {number} j - Y index
     * @param {number} k - Z index
     * @param {Object} selection - Selection config
     * @returns {boolean} In selection
     */
    isInSelection(i, j, k, selection) {
        switch (selection.type) {
            case 'full':
                return true;
            case 'slice':
                if (selection.axis === 'x' && i === selection.position) return true;
                if (selection.axis === 'y' && j === selection.position) return true;
                if (selection.axis === 'z' && k === selection.position) return true;
                return false;
            case 'box':
                return i >= selection.min[0] && i < selection.max[0] &&
                       j >= selection.min[1] && j < selection.max[1] &&
                       k >= selection.min[2] && k < selection.max[2];
            default:
                return true;
        }
    }

    /**
     * Clear all arrows
     */
    clearArrows() {
        for (const arrow of this.arrows) {
            this.arrowGroup.remove(arrow);
            arrow.geometry.dispose();
            arrow.material.dispose();
        }
        this.arrows = [];
    }

    /**
     * Resize renderer
     * @param {number} width - New width
     * @param {number} height - New height
     */
    resize(width, height) {
        if (this.camera) {
            this.camera.aspect = width / height;
            this.camera.updateProjectionMatrix();
        }

        if (this.renderer) {
            this.renderer.setSize(width, height);
        }
    }

    /**
     * Get visualization status
     * @returns {Object} Status info
     */
    getStatus() {
        return {
            initialized: !!this.scene,
            dimensions: this.dimensions,
            cells: this.gridSize,
            arrowCount: this.arrows.length > 0 ? this.arrows[0].count : 0
        };
    }

    /**
     * Clean up resources
     */
    /**
     * Run Julia code
     * @param {string} code - Julia code to execute
     */
    runCode(code) {
        console.log('Running Julia code:', code);
        
        // Update status message
        const statusElement = document.getElementById('status-message');
        if (statusElement) {
            statusElement.textContent = 'Running Julia code...';
        }
        
        // Send code to Julia server
        this.sendCommand('run_code', { code });
    }

    /**
     * Dispose resources
     */
    dispose() {
        if (this.animationId) {
            cancelAnimationFrame(this.animationId);
        }

        this.clearArrows();

        if (this.renderer) {
            this.renderer.dispose();
        }

        if (this.webSocketClient) {
            this.webSocketClient.disconnect();
        }

        console.log('MagneticVisualization disposed');
    }
}

// Export to global scope
window.MagneticVisualization = MagneticVisualization;

// Initialize GUI when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    console.log('Initializing MicroMagneticGUI...');
    
    // Create and initialize visualization
    const container = document.getElementById('visualization-container');
    if (container) {
        const visualization = new MagneticVisualization();
        visualization.init(container);
        
        // Add event listener for run-code button
        const runButton = document.getElementById('run-code');
        if (runButton) {
            runButton.addEventListener('click', () => {
                console.log('Run button clicked');
                
                if (visualization.codeEditor) {
                    const code = visualization.codeEditor.getValue();
                    visualization.runCode(code);
                }
            });
        }
        
        // Add event listener for run-simulation-btn button
        const runSimulationButton = document.getElementById('run-simulation-btn');
        if (runSimulationButton) {
            runSimulationButton.addEventListener('click', () => {
                console.log('Run simulation button clicked');
                
                if (visualization.codeEditor) {
                    const code = visualization.codeEditor.getValue();
                    visualization.runCode(code);
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
                visualization.sendCommand('test_connection', {});
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
                visualization.sendCommand('test_command', {
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
                visualization.sendCommand('get_status', {});
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
                visualization.sendCommand('stop_server', {});
            });
        }
        
        // Initialize task templates
        initTaskTemplates();
        
        // Create relax task manager and populate cells-container with relax task steps
        const relaxTaskManager = createRelaxTaskManager('#cells-container');
        
        // Initialize CellSelectorManager
        const cellSelectorManager = new CellSelectorManager(relaxTaskManager);
        
        // Update cell selector options initially
        cellSelectorManager.updateOptions();
        
        // Add event listener for add-cell-btn
        const addCellButton = document.getElementById('add-cell-btn');
        if (addCellButton) {
            addCellButton.addEventListener('click', () => {
                console.log('Add new cell button clicked');
                const newCell = new Cell('');
                relaxTaskManager.addCell(newCell);
                cellSelectorManager.updateOptions();
            });
        }
        
        // Add event listener for update-cell-btn
        const updateCellButton = document.getElementById('update-cell-btn');
        if (updateCellButton) {
            updateCellButton.addEventListener('click', () => {
                console.log('Update cell button clicked');
                const selectedCell = relaxTaskManager.getSelectedCell();
                if (selectedCell) {
                    // Here you can implement the update logic
                    console.log('Updating cell:', selectedCell.id);
                    // For example, you could refresh the cell's content or update its properties
                    // Remove the refresh() call as it doesn't exist
                } else {
                    console.log('No cell selected for update');
                }
            });
        }
        
        console.log('MicroMagneticGUI initialized successfully!');
    } else {
        console.error('Visualization container not found');
    }
});