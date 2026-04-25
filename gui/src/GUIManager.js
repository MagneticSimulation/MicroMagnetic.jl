import WebSocketClient from './client.js';


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
        this.simStatePanel = null;
        this.codeEditorPanel = null;
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
        this.wsClient.on('visualization_update', this.handleVisualizationUpdate.bind(this));
        this.wsClient.on('sim_state_update', this.handleSimStateUpdate.bind(this));
        this.wsClient.on('run_code_response', this.handleCommandResponse.bind(this));
    }

    /**
     * Handle visualization update message
     * @param {Object} data - Visualization data (mesh, spin, Ms, etc.)
     */
    handleVisualizationUpdate(data) {
        console.log('Visualization update:', data);
        if (this.visualization) {
            this.visualization.updateFromVisData(data);
        }
    }

    /**
     * Handle simulation state update message
     * @param {Object} data - Simulation state data
     */
    handleSimStateUpdate(data) {
        console.log('Sim state update:', data);
        if (this.simStatePanel) {
            this.simStatePanel.update(data);
        }
    }

    /**
     * Handle magnetization data message (legacy)
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
                this.updateExecutionUI('error', 'Execution failed', response.stderr || 'Unknown error occurred');
            }
        }
        
        // Update CodeEditorPanel output if available
        if (this.codeEditorPanel) {
            const outputText = response.stdout || response.stderr || '';
            this.codeEditorPanel.showOutput(outputText, !response.success);
        }
        
        // Note: Visualization updates are now handled by the visualization_update event
        // Mesh and magnetization data are no longer sent through run_code_response
    }

    /**
     * Run Julia code
     * @param {string} code - Julia code to execute
     * @returns {Promise} - Promise that resolves when code execution message is sent
     */
    runCode(code) {
        this.updateExecutionUI('running', 'Executing code...', '');
        // Send code to Julia server and return the promise
        return this.wsClient.sendCommand('run_code', { code });
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
     * Set CodeEditorPanel reference
     */
    setCodeEditorPanel(panel) {
        this.codeEditorPanel = panel;
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