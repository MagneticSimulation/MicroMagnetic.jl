import WebSocketClient from './client.js';
import { createExampleTaskManager, examples } from './tasks.js';

/**
 * GUIManager class for handling user interactions and WebSocket communication
 */
class GUIManager {
    constructor(visualization) {
        this.visualization = visualization;
        this.wsClient = null;
        this.isConnected = false;
        this.connectionAttempts = 0;
        this.maxConnectionAttempts = 5;
        this.simulationType = 'fd';
        this.taskType = 'relax';
        this.example = '';
        this.taskManager = null;
    }

    /**
     * Initialize WebSocket client for communicating with Julia server
     */
    initWebSocketClient() {
        try {
            console.log('Initializing WebSocket client');
            
            this.wsClient = new WebSocketClient();
            
            // Initialize WebSocket UI elements
            this.initWebSocketUI();
            
            // Set up event listeners
            this.setupEventListeners();
            
        } catch (error) {
            console.error('Failed to initialize WebSocket client:', error);
            this.updateExecutionUI('error', "WebSocket error!", error.toString());
        }
    }

    /**
     * Set up WebSocket event listeners
     */
    setupEventListeners() {
        // WebSocket connection events
        this.wsClient.on('connect', () => {
            console.log('WebSocket connected');
            this.isConnected = true;
            this.updateWebSocketStatus(true);
        });
        
        this.wsClient.on('disconnect', () => {
            console.log('WebSocket disconnected');
            this.isConnected = false;
            this.updateWebSocketStatus(false);
        });
        
        this.wsClient.on('error', (error) => {
            console.error('WebSocket error:', error);
            this.isConnected = false;
            this.updateWebSocketStatus(false);
            this.updateExecutionUI('error', "Connection error!", error.message);
        });
        
        // Data events
        this.wsClient.on('m_data', this.handleMagnetizationData.bind(this));
        this.wsClient.on('run_code_response', this.handleCommandResponse.bind(this));
    }

    /**
     * Handle magnetization data message
     * @param {Object} data - Magnetization data
     */
    handleMagnetizationData(data) {
        console.log('Received magnetization data:', data);
        
        if (this.visualization && data.m_data) {
            this.visualization.updateMagnetization(data.m_data);
        }
    }

    /**
     * Initialize relax task manager
     * @param {string} containerSelector - CSS selector for the container element
     */
    initRelaxTaskManager(containerSelector) {
        this.taskManager = createExampleTaskManager(containerSelector, 'std4', this);
        
        // Pass wsClient to taskManager if needed
        if (this.taskManager && this.wsClient) {
            this.taskManager.wsClient = this.wsClient;
        }
        
        return this.taskManager;
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
                if (this.wsClient && this.wsClient.isConnected()) {
                    // Disconnect
                    this.wsClient.disconnect();
                } else {
                    // Connect
                    await this.wsClient.connect();
                    this.updateExecutionUI('success', "Connection successful!", "");
                }
            } catch (error) {
                console.error('WebSocket connection error:', error);
                this.updateWebSocketStatus(false);
                this.updateExecutionUI('error', "Connection error!", error.message);
            }
        });
        
        // Initialize simulation type and task type selections
        this.initSimulationTypeSelection();
        this.initTaskTypeSelection();
        this.initExampleSelection();
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
     * Initialize example selection
     */
    initExampleSelection() {
        const selectElement = document.getElementById('example-selector');
        
        selectElement.addEventListener('change', (event) => {
            const selectedExample = event.target.value;
            console.log('Selected example:', selectedExample);
            // Store selected example
            this.example = selectedExample;
            
            // If a valid example is selected, create a task manager for it
            if (selectedExample && examples[selectedExample]) {
                this.taskManager = createExampleTaskManager('#cells-container', selectedExample, this);
                
                // Pass wsClient to taskManager if needed
                if (this.taskManager && this.wsClient) {
                    this.taskManager.wsClient = this.wsClient;
                }
            }
        });
        
        // Set default example
        this.example = '';
    }

    /**
     * Handle command response from Julia server
     * @param {Object} response - Command response
     */
    handleCommandResponse(response) {
        console.log('Command response:', response);
        
        // Update execution status UI
        if (this.updateExecutionUI) {
            if (response.success) {
                this.updateExecutionUI('success', 'Execution completed', response.stdout || 'No output');
            } else {
                this.updateExecutionUI('error', 'Execution failed', response.error || 'Unknown error occurred');
            }
        }

        // Handle different response types with validation
        if (response.type === 'fd_mesh_data') {
           if (this.visualization && response.fd_mesh_data) {
               this.visualization.displayFDMesh(response.fd_mesh_data);
           }
        } else if (response.type === 'Ms_data') {
            if (this.visualization && response.Ms_data) {
                this.visualization.displayCustomSurface(response.Ms_data);
            }
        } else if (response.type === 'm_data') {
            if (this.visualization && response.m_data) {
                this.visualization.updateMagnetization(response.m_data);
            }
        }
    }

    /**
     * Run Julia code
     * @param {string} code - Julia code to execute
     * @returns {Promise} - Promise that resolves when code execution message is sent
     */
    runCode(code, desc = '') {        
        this.updateExecutionUI('running', 'Executing code...', '');
        // Send code to Julia server and return the promise
        return this.wsClient.sendCommand('run_code', { code, desc });
    }

    /**
     * Update execution UI (status bar and unified output)
     * @param {string} status - Status type: 'running', 'success', 'error'
     * @param {string} statusMessage - Status bar message
     * @param {string} outputMessage - Output area message
     */
    updateExecutionUI(status, statusMessage, outputMessage) {
        // Update status bar
        const executionStatusElement = document.getElementById('execution-status');
        if (executionStatusElement) {
            executionStatusElement.className = 'status-value';
            executionStatusElement.classList.add(status);
            executionStatusElement.textContent = statusMessage;
        }
        
        // Update unified output only if not in running status
        if (status !== 'running') {
            const unifiedOutput = document.getElementById('unified-output');
            if (unifiedOutput) {
                const outputClasses = {
                    success: 'output-success',
                    error: 'output-error'
                };
                
                const outputClass = outputClasses[status] || '';
                const messageElement = `<pre>${outputMessage}</pre>`;
                
                unifiedOutput.innerHTML = `
                    <div class="${outputClass}">
                        ${messageElement}
                    </div>
                `;
            }
        }
    }

    /**
     * Dispose resources
     */
    dispose() {
        if (this.wsClient) {
            this.wsClient.disconnect();
        }

        console.log('GUIManager disposed');
    }
}

export default GUIManager;