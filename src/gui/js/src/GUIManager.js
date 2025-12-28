import WebSocketClient from './client.js';

/**
 * GUIManager class for handling user interactions and WebSocket communication
 */
class GUIManager {
    constructor(visualization) {
        this.visualization = visualization;
        this.webSocketClient = null;
        this.wsClient = null;
        this.isConnected = false;
        this.connectionAttempts = 0;
        this.maxConnectionAttempts = 5;
        this.simulationType = 'fd';
        this.taskType = 'relax';
    }

    /**
     * Initialize WebSocket client for communicating with Julia server
     */
    initWebSocketClient() {
        try {
            console.log('Initializing WebSocket client');
            
            this.webSocketClient = new WebSocketClient();
            
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
                if (this.visualization) {
                    this.visualization.updateMagnetization(data);
                }
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
        if (this.webSocketClient) {
            this.webSocketClient.disconnect();
        }

        console.log('GUIManager disposed');
    }
}

export default GUIManager;