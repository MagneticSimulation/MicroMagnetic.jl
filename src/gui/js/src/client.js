// WebSocket client for MicroMagnetic.jl server
// Supports multiple users, reconnection, and command execution

class WebSocketClient {
    /**
     * Create a new WebSocket client
     * @param {Object} options - Configuration options
     * @param {string} [options.serverUrl] - WebSocket server URL (auto-detected by default)
     * @param {number} [options.heartbeatInterval=30000] - Heartbeat interval in ms
     * @param {boolean} [options.autoConnect=true] - Automatically connect on creation
     * @param {string} [options.sessionId] - Session ID for client identification
     */
    constructor(options = {}) {
        // 自动根据当前页面协议选择 WebSocket 协议
        const protocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:';
        this.config = {
            serverUrl: options.serverUrl || `${protocol}//${window.location.host}`,
            heartbeatInterval: options.heartbeatInterval || 30000,
            sessionId: options.sessionId || this.generateSessionId()
        };

        // Connection state
        this.ws = null;
        this.clientId = null;
        this.connected = false;
        this.heartbeatTimer = null;
        
        // Event listeners
        this.listeners = new Map();
        
        // Command callbacks
        this.commandCallbacks = new Map();
        this.commandIdCounter = 1;
        
        // Logging
        this.logLevel = options.logLevel || 'info';
        
        // Auto-connect if enabled       
        if (this.config.autoConnect) {
            this.connect();
        }
    }
    
    /**
     * Generate a session ID
     * @returns {string} Generated session ID
     */
    generateSessionId() {
        return 'session_' + Math.random().toString(36).substring(2, 11);
    }
    
    /**
     * Log a message with timestamp
     * @param {string} level - Log level
     * @param {string} message - Message to log
     * @param {...any} args - Additional arguments
     */
    log(level, message, ...args) {
        const timestamp = new Date().toISOString();
        const logMessage = `[${timestamp}] [${level.toUpperCase()}] ${message}`;
        
        // Filter by log level
        const levels = ['debug', 'info', 'warn', 'error'];
        const currentLevel = levels.indexOf(this.logLevel);
        const messageLevel = levels.indexOf(level);
        
        if (messageLevel >= currentLevel) {
            if (level === 'error') {
                console.error(logMessage, ...args);
            } else if (level === 'warn') {
                console.warn(logMessage, ...args);
            } else {
                console.log(logMessage, ...args);
            }
        }
    }
    
    /**
     * Connect to the WebSocket server
     * @returns {Promise} Promise that resolves when connected
     */
    connect() {
        return new Promise((resolve, reject) => {
            if (this.connected) {
                resolve();
                return;
            }
            
            // Add session ID to the connection URL
            const url = new URL(this.config.serverUrl);
            url.searchParams.append('session_id', this.config.sessionId);
            
            this.log('info', `Connecting to ${url.toString()}`);
            
            try {
                this.ws = new WebSocket(url.toString());
                
                this.ws.onopen = (event) => {
                    this.log('info', 'WebSocket connection established');
                    this.connected = true;
                    this.emit('connect', event);
                    this.startHeartbeat();
                    resolve(event);
                };                
                this.ws.onmessage = (event) => {
                    this.handleMessage(event.data);
                };
                
                this.ws.onerror = (event) => {
                    this.log('error', 'WebSocket error', event);
                    this.emit('error', event);
                    reject(event);
                };
                
                this.ws.onclose = (event) => {
                    this.log('info', `WebSocket connection closed: ${event.code} ${event.reason}`);
                    this.connected = false;
                    this.clientId = null;
                    this.stopHeartbeat();
                    this.emit('disconnect', event);
                };
                
            } catch (error) {
                this.log('error', 'Failed to create WebSocket connection', error);
                reject(error);
            }
        });
    }
    
    /**
     * Attempt to reconnect to the server
     */
    attemptReconnect() {
        // Removed automatic reconnection logic
        this.log('info', 'Reconnection attempt disabled');
    }
    
    /**
     * Start heartbeat mechanism
     */
    startHeartbeat() {
        if (this.heartbeatTimer) {
            clearInterval(this.heartbeatTimer);
        }
        
        this.heartbeatTimer = setInterval(() => {
            if (this.connected) {
                this.sendHeartbeat();
            }
        }, this.config.heartbeatInterval);
    }
    
    /**
     * Stop heartbeat mechanism
     */
    stopHeartbeat() {
        if (this.heartbeatTimer) {
            clearInterval(this.heartbeatTimer);
            this.heartbeatTimer = null;
        }
    }
    
    /**
     * Send heartbeat to server
     */
    sendHeartbeat() {
        const heartbeatMsg = {
            type: 'heartbeat',
            timestamp: Date.now()
        };
        
        this.send(heartbeatMsg).catch(error => {
            this.log('warn', 'Heartbeat failed', error);
        });
    }
    
    /**
     * Handle incoming WebSocket message
     * @param {string} data - Raw message data
     */
    handleMessage(data) {
        try {
            const message = JSON.parse(data);
            this.log('debug', 'Received message', message);
            
            // Handle different message types based on server implementation
            const type = message.type;
            
            if (type === 'connected') {
                // Handle connection confirmation
                this.clientId = message.data.client_id;
                this.log('info', `Connected with client ID: ${this.clientId}`);
                this.emit('connection', message.data);
            } else if (type === 'command_response') {
                this.handleCommandResponse(message);
            } else if (type === 'heartbeat_response') {
                this.emit('heartbeat', message.data);
            } else if (type === 'error') {
                this.log('error', `Server error: ${message.data.message || 'Unknown error'}`);
                this.emit('error', message.data);
            } else {
                // Emit all other message types
                this.emit(type, message.data);
            }
            
        } catch (error) {
            this.log('error', 'Failed to parse message', error, data);
        }
    }
    
    /**
     * Handle command response from server
     * @param {Object} message - Command response message
     */
    handleCommandResponse(message) {
        const data = message.data;
        const commandId = data.id;
        const callback = this.commandCallbacks.get(commandId);
        
        if (callback) {
            if (data.success) {
                callback.resolve(data.result);
            } else {
                callback.reject(new Error(data.error || 'Command failed'));
            }
            this.commandCallbacks.delete(commandId);
        } else {
            this.log('warn', `No callback found for command ID: ${commandId}`);
        }
        
        this.emit('command_response', data);
    }
    
    /**
     * Send a raw message to the server
     * @param {Object} message - Message to send
     * @returns {Promise} Promise that resolves when message is sent
     */
    send(message) {
        return new Promise((resolve, reject) => {
            if (!this.connected || !this.ws) {
                reject(new Error('Not connected to server'));
                return;
            }
            
            try {
                const jsonStr = JSON.stringify(message);
                this.ws.send(jsonStr);
                this.log('debug', 'Sent message', message);
                resolve();
            } catch (error) {
                this.log('error', 'Failed to send message', error);
                reject(error);
            }
        });
    }
    
    /**
     * Send a command to the server
     * @param {string} command - Command name
     * @param {Object} data - Command data
     * @returns {Promise} Promise that resolves with command result
     */
    sendCommand(command, data = {}) {
        return new Promise((resolve, reject) => {
            const commandId = this.commandIdCounter++;
            
            const message = {
                type: 'command',
                data: {
                    id: commandId.toString(),
                    command: command,
                    ...data
                }
            };
            
            // Store callback
            this.commandCallbacks.set(commandId.toString(), { resolve, reject });
            
            // Send message
            this.send(message).catch(reject);
        });
    }
    
    /**
     * Send a custom message to the server
     * @param {string} type - Message type
     * @param {Object} data - Message data
     * @returns {Promise} Promise that resolves when message is sent
     */
    sendMessage(type, data = {}) {
        const message = {
            type: type,
            data: data
        };
        
        return this.send(message);
    }
    
    /**
     * Test connection to server
     * @returns {Promise} Promise that resolves with connection status
     */
    testConnection() {
        return this.sendMessage('ping').then(() => ({ status: 'connected' }));
    }
    
    /**
     * Get server status
     * @returns {Promise} Promise that resolves with server status
     */
    getStatus() {
        return this.sendMessage('get_status');
    }
    
    /**
     * Start a simulation
     * @param {Object} simulationData - Simulation parameters
     * @returns {Promise} Promise that resolves with simulation info
     */
    startSimulation(simulationData = {}) {
        return this.sendMessage('simulation_start', simulationData);
    }
    
    /**
     * Stop a simulation
     * @param {Object} simulationData - Simulation data
     * @returns {Promise} Promise that resolves when simulation is stopped
     */
    stopSimulation(simulationData = {}) {
        return this.sendMessage('simulation_stop', simulationData);
    }
    
    /**
     * Get simulation data
     * @param {Object} params - Parameters for data retrieval
     * @returns {Promise} Promise that resolves with simulation data
     */
    getData(params = {}) {
        return this.sendMessage('get_data', params);
    }
    
    /**
     * Disconnect from server
     */
    disconnect() {
        this.log('info', 'Disconnecting from server');
        
        this.stopHeartbeat();
        
        if (this.ws) {
            this.ws.close(1000, 'Client disconnect');
        }
        
        this.connected = false;
        this.clientId = null;
        
        // Clear command callbacks
        this.commandCallbacks.forEach((callback) => {
            callback.reject(new Error('Disconnected'));
        });
        this.commandCallbacks.clear();
        
        this.emit('disconnected');
    }
    
    /**
     * Add event listener
     * @param {string} event - Event name
     * @param {Function} callback - Event callback
     */
    on(event, callback) {
        if (!this.listeners.has(event)) {
            this.listeners.set(event, []);
        }
        this.listeners.get(event).push(callback);
    }
    
    /**
     * Remove event listener
     * @param {string} event - Event name
     * @param {Function} callback - Event callback to remove
     */
    off(event, callback) {
        if (this.listeners.has(event)) {
            const callbacks = this.listeners.get(event);
            const index = callbacks.indexOf(callback);
            if (index !== -1) {
                callbacks.splice(index, 1);
            }
        }
    }
    
    /**
     * Emit an event
     * @param {string} event - Event name
     * @param {...any} args - Event arguments
     */
    emit(event, ...args) {
        if (this.listeners.has(event)) {
            this.listeners.get(event).forEach(callback => {
                try {
                    callback(...args);
                } catch (error) {
                    this.log('error', `Error in event listener for '${event}'`, error);
                }
            });
        }
    }
    
    /**
     * Check if connected to server
     * @returns {boolean} Connection status
     */
    isConnected() {
        return this.connected;
    }
    
    /**
     * Get client ID
     * @returns {string|null} Client ID
     */
    getClientId() {
        return this.clientId;
    }
    
    /**
     * Get session ID
     * @returns {string} Session ID
     */
    getSessionId() {
        return this.config.sessionId;
    }
}

// Export for different module systems
if (typeof module !== 'undefined' && module.exports) {
    // CommonJS/Node.js
    module.exports = WebSocketClient;
} else if (typeof define === 'function' && define.amd) {
    // AMD
    define([], function() {
        return WebSocketClient;
    });
} else {
    // Browser global
    window.WebSocketClient = WebSocketClient;
}