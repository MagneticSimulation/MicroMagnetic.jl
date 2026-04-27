class CodeEditorPanel {
    constructor(containerId, guiManager) {
        this.container = document.getElementById(containerId);
        this.guiManager = guiManager;
        this.editor = null;
        this.history = [];
        this.historyVisible = false;
        
        if (!this.container) {
            console.error(`CodeEditorPanel container with ID '${containerId}' not found.`);
            return;
        }
        
        this.buildUI();
        this.initEditor();
        this.bindEvents();
    }
    
    buildUI() {
        this.container.innerHTML = `
            <div class="code-editor-panel">
                <div class="editor-header">
                    <h4>Code Editor</h4>
                    <div class="editor-controls">
                    <button id="run-button" class="run-btn">Run</button>
                    <button id="clear-editor" class="toggle-btn">Clear Editor</button>
                    <button id="toggle-history" class="toggle-btn">History</button>
                </div>
                </div>
                
                <div class="editor-content">
                    <div id="editor-container"></div>
                    
                    <div id="history-panel" class="history-panel hidden">
                        <h5>History</h5>
                        <div id="history-list" class="history-list"></div>
                    </div>
                </div>
            </div>
        `;
        
        this.runButton = this.container.querySelector('#run-button');
        this.clearEditorButton = this.container.querySelector('#clear-editor');
        this.toggleHistoryButton = this.container.querySelector('#toggle-history');
        this.editorContainer = this.container.querySelector('#editor-container');
        this.historyPanel = this.container.querySelector('#history-panel');
        this.historyList = this.container.querySelector('#history-list');
    }
    
    initEditor() {
        // Check if CodeMirror is available
        if (typeof CodeMirror === 'undefined') {
            console.error('CodeMirror is not loaded. Please include the CodeMirror library.');
            this.editorContainer.innerHTML = '<div class="editor-error">CodeMirror library not found.</div>';
            return;
        }
        
        // Create CodeMirror instance
        this.editor = CodeMirror(this.editorContainer, {
            mode: 'julia',
            theme: 'default',
            lineNumbers: true,
            lineWrapping: false,
            viewportMargin: Infinity,
            indentUnit: 4,
            smartIndent: true,
            value: '',
            placeholder: 'Enter Julia code here...',
            extraKeys: {
                'Ctrl-Enter': () => this.runCode(),
                'Shift-Enter': () => this.runCurrentLine()
            }
        });
        
        // Set initial size
        this.editor.setSize('100%', '250px');
        
        // Focus the editor
        setTimeout(() => {
            this.editor.focus();
        }, 100);
    }
    
    bindEvents() {
        // Run button click
        if (this.runButton) {
            this.runButton.addEventListener('click', () => this.runCode());
        }
        
        // Clear editor button
        if (this.clearEditorButton) {
            this.clearEditorButton.addEventListener('click', () => this.clearEditor());
        }
        
        // Toggle history visibility
        if (this.toggleHistoryButton) {
            this.toggleHistoryButton.addEventListener('click', () => this.toggleHistory());
        }
    }
    
    runCode() {
        const code = this.editor.getValue().trim();
        if (!code) return;
        
        // Add to history
        this.addToHistory(code);
        
        // Run the code using GUIManager
        if (this.guiManager && this.guiManager.runCode) {
            this.guiManager.runCode(code);
        } else {
            console.error('Error: GUIManager not available');
        }
    }
    
    runCurrentLine() {
        const cursor = this.editor.getDoc().getCursor();
        const line = this.editor.getLine(cursor.line);
        const code = line.trim();
        
        if (code) {
            // Add to history
            this.addToHistory(code);
            
            // Run the line
            if (this.guiManager && this.guiManager.runCode) {
                this.guiManager.runCode(code, 'relax');
            }
        }
    }
    

    
    toggleHistory() {
        if (!this.historyPanel || !this.toggleHistoryButton) return;
        
        this.historyVisible = !this.historyVisible;
        this.historyPanel.classList.toggle('hidden', !this.historyVisible);
        this.toggleHistoryButton.textContent = this.historyVisible ? 'Hide History' : 'History';
        
        // Update history list when showing
        if (this.historyVisible) {
            this.renderHistory();
        }
    }
    
    addToHistory(code) {
        // Add to history array (avoid duplicates)
        if (this.history.length === 0 || this.history[this.history.length - 1] !== code) {
            this.history.push(code);
            
            // Limit history to 50 items
            if (this.history.length > 50) {
                this.history.shift();
            }
        }
        
        // Update history panel if visible
        if (this.historyVisible) {
            this.renderHistory();
        }
    }
    
    renderHistory() {
        if (!this.historyList) return;
        
        if (this.history.length === 0) {
            this.historyList.innerHTML = '<div class="history-empty">No history available</div>';
            return;
        }
        
        this.historyList.innerHTML = this.history.map((code, index) => {
            const truncatedCode = code.length > 80 ? code.substring(0, 80) + '...' : code;
            return `
                <div class="history-item" data-index="${index}">
                    <div class="history-index">${index + 1}</div>
                    <div class="history-code">${this.escapeHtml(truncatedCode)}</div>
                </div>
            `;
        }).join('');
        
        // Bind click events to history items
        this.historyList.querySelectorAll('.history-item').forEach(item => {
            item.addEventListener('click', () => {
                const index = parseInt(item.dataset.index);
                const code = this.history[index];
                this.editor.setValue(code);
                this.editor.focus();
            });
        });
    }
    
    clearHistory() {
        this.history = [];
        this.renderHistory();
    }
    
    insertAtCursor(text) {
        if (!this.editor) return;
        
        const doc = this.editor.getDoc();
        const cursor = doc.getCursor();
        doc.replaceRange(text, cursor);
        this.editor.focus();
    }

    clearEditor() {
        if (!this.editor) return;
        this.editor.setValue("");
        this.editor.focus();
    }

    escapeHtml(str) {
        const div = document.createElement('div');
        div.appendChild(document.createTextNode(str));
        return div.innerHTML;
    }
}

export default CodeEditorPanel;