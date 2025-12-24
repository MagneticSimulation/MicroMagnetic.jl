/**
 * Cell class - Simple code cell using CodeMirror
 */

class Cell {
    /**
     * Constructor
     * @param {string} content - Cell content
     * @param {string} description - Cell description
     */
    constructor(content = '', description = 'Code Cell') {
        this.id = `cell-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`;
        this.content = content;
        this.description = description;
        this.element = null;
        this.editor = null;
        this.outputElement = null;
        this.output = '';
        this.collapsed = false; // Default expanded state
    }

    /**
     * Render cell
     * @returns {HTMLElement} Cell DOM element
     */
    render() {
        // Create cell container
        this.element = document.createElement('div');
        this.element.className = 'cell';
        this.element.dataset.cellId = this.id;

        // Create cell header
        const header = document.createElement('div');
        header.className = 'cell-header';
        
        // Add click event for collapse/expand
        header.addEventListener('click', () => {
            this.toggleCollapse();
        });

        // Create description
        const description = document.createElement('div');
        description.className = 'cell-description';
        description.textContent = this.description;

        // Create action buttons
        const actions = document.createElement('div');
        actions.className = 'cell-actions';

        // Run button
        const runBtn = document.createElement('button');
        runBtn.className = 'cell-btn run-btn';
        runBtn.textContent = '▶';
        runBtn.title = 'Run';
        runBtn.addEventListener('click', (e) => {
            e.stopPropagation(); // Prevent bubble, avoid triggering header click event
            this.run();
        });

        // Delete button
        const deleteBtn = document.createElement('button');
        deleteBtn.className = 'cell-btn delete-btn';
        deleteBtn.textContent = '✕';
        deleteBtn.title = 'Delete';
        deleteBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            this.delete();
        });

        actions.appendChild(runBtn);
        actions.appendChild(deleteBtn);
        header.appendChild(description);
        header.appendChild(actions);

        // Create content area
        const content = document.createElement('div');
        content.className = 'cell-content';

        // Create editor container
        const editorContainer = document.createElement('div');
        editorContainer.className = 'cell-editor';

        // Create textarea for CodeMirror
        const textarea = document.createElement('textarea');
        textarea.value = this.content;
        editorContainer.appendChild(textarea);

        content.appendChild(editorContainer);

        // Create output area
        this.outputElement = document.createElement('div');
        this.outputElement.className = 'cell-output';
        this.outputElement.style.display = 'none';

        this.element.appendChild(header);
        this.element.appendChild(content);
        this.element.appendChild(this.outputElement);

        // Initialize CodeMirror editor
        this.initEditor(textarea);

        return this.element;
    }

    /**
     * Initialize CodeMirror editor
     */
    initEditor(textarea) {
        this.editor = CodeMirror.fromTextArea(textarea, {
            mode: 'julia',
            theme: 'material',
            lineNumbers: false,
            lineWrapping: true,
            viewportMargin: Infinity,
            indentUnit: 4,
            smartIndent: true,
            extraKeys: {
                'Shift-Enter': () => this.run(),
                'Ctrl-Enter': () => {
                    this.run();
                    // Can add logic to create new cell here
                }
            }
        });

        // Listen for content changes
        this.editor.on('change', () => {
            this.content = this.editor.getValue();
        });

        // Set appropriate font size
        this.editor.setSize('100%', 'auto');
    }

    /**
     * Run cell code
     */
    run() {
        const code = this.editor.getValue();
        this.outputElement.style.display = 'block';
        
        try {
            // Can integrate Julia execution logic here
            // Temporarily use simulated execution
            this.outputElement.innerHTML = `
                <div class="output-success">
                    <strong>Output:</strong>
                    <pre>Executed code (${code.length} characters)</pre>
                </div>
            `;
            this.outputElement.className = 'cell-output success';
        } catch (error) {
            this.outputElement.innerHTML = `
                <div class="output-error">
                    <strong>Error:</strong>
                    <pre>${error.message}</pre>
                </div>
            `;
            this.outputElement.className = 'cell-output error';
        }
    }

    /**
     * Delete cell
     */
    delete() {
        if (this.element && this.element.parentNode) {
            this.element.parentNode.removeChild(this.element);
        }
    }

    /**
     * Toggle collapse/expand state
     */
    toggleCollapse() {
        this.collapsed = !this.collapsed;
        const content = this.element.querySelector('.cell-content');
        const output = this.outputElement;
        
        if (this.collapsed) {
            content.style.display = 'none';
            output.style.display = 'none';
        } else {
            content.style.display = 'block';
            if (this.output) {
                output.style.display = 'block';
            }
        }
    }

    /**
     * Get cell content
     */
    getValue() {
        return this.editor ? this.editor.getValue() : this.content;
    }

    /**
     * Set cell content
     */
    setValue(content) {
        this.content = content;
        if (this.editor) {
            this.editor.setValue(content);
        }
    }
}

/**
 * CellManager class - Manage multiple cells
 */

class CellManager {
    /**
     * Constructor
     * @param {string} containerSelector - Container selector
     */
    constructor(containerSelector) {
        this.container = document.querySelector(containerSelector);
        this.cells = [];
        this.title = 'Code Cells';
        this.collapsed = false; // Default expanded state
        
        if (this.container) {
            this.init();
        }
    }

    /**
     * Initialize CellManager
     */
    init() {
        // Create header
        this.createHeader();
        
        // Create default cells
        this.createDefaultCells();
        
        // Set add cell button event
        const addBtn = document.getElementById('add-cell-btn');
        if (addBtn) {
            addBtn.addEventListener('click', () => this.addCell());
        }
    }

    /**
     * Create header
     */
    createHeader() {
        // Save original container content
        const originalContent = this.container.innerHTML;
        
        // Clear container
        this.container.innerHTML = '';
        
        // Create header
        const header = document.createElement('div');
        header.className = 'cell-manager-header';
        
        // Add click event for collapse/expand
        header.addEventListener('click', () => {
            this.toggleCollapse();
        });

        // Create title text
        const titleText = document.createElement('div');
        titleText.className = 'cell-manager-title';
        titleText.textContent = this.title;
        
        // Create action buttons container
        const headerActions = document.createElement('div');
        headerActions.className = 'cell-manager-actions';
        
        // Run all button
        const runAllBtn = document.createElement('button');
        runAllBtn.className = 'cell-manager-btn run-all-btn';
        runAllBtn.textContent = 'Run All';
        runAllBtn.title = 'Run All Cells';
        runAllBtn.addEventListener('click', (e) => {
            e.stopPropagation(); // Prevent bubble, avoid triggering header click event
            this.runAll();
        });

        // Clear all outputs button
        const clearAllBtn = document.createElement('button');
        clearAllBtn.className = 'cell-manager-btn clear-all-btn';
        clearAllBtn.textContent = 'Clear All';
        clearAllBtn.title = 'Clear All Outputs';
        clearAllBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            this.clearAllOutputs();
        });

        headerActions.appendChild(runAllBtn);
        headerActions.appendChild(clearAllBtn);
        
        // Add title and buttons to header
        header.appendChild(titleText);
        header.appendChild(headerActions);
        
        // Create content container
        this.contentArea = document.createElement('div');
        this.contentArea.className = 'cell-manager-content';
        this.contentArea.innerHTML = originalContent;
        
        // Add header and content to container
        this.container.appendChild(header);
        this.container.appendChild(this.contentArea);
    }

    /**
     * Create default cells
     */
    createDefaultCells() {
        // Add a default cell
        const defaultCell = new Cell(
            '# Welcome to MicroMagnetic.jl\n# This is a code cell. Press Shift+Enter to run.',
            'Default Cell'
        );
        this.addCell(defaultCell);
    }

    /**
     * Add cell
     * @param {Cell} cell - Cell instance
     */
    addCell(cell) {
        if (!cell) {
            // Create new cell
            cell = new Cell(
                '# New code cell\nprintln("Hello from new cell!")',
                'New Cell'
            );
        }
        
        const cellElement = cell.render();
        
        if (this.contentArea) {
            this.contentArea.appendChild(cellElement);
            this.cells.push(cell);
        }
    }

    /**
     * Get all cells
     */
    getAllCells() {
        return this.cells;
    }

    /**
     * Run all cells
     */
    runAll() {
        this.cells.forEach(cell => {
            if (cell.run) {
                cell.run();
            }
        });
    }

    /**
     * Clear all outputs
     */
    clearAllOutputs() {
        this.cells.forEach(cell => {
            if (cell.outputElement) {
                cell.outputElement.style.display = 'none';
                cell.outputElement.innerHTML = '';
            }
        });
    }

    /**
     * Toggle collapse/expand state
     */
    toggleCollapse() {
        this.collapsed = !this.collapsed;
        
        if (this.contentArea) {
            if (this.collapsed) {
                this.contentArea.style.display = 'none';
            } else {
                this.contentArea.style.display = 'block';
            }
        }
    }
}

// Export classes - Support ES modules, CommonJS and browser globals
export { Cell, CellManager };

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { Cell, CellManager };
} else {
    window.Cell = Cell;
    window.CellManager = CellManager;
}