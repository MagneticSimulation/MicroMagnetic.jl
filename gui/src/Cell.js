/**
 * Cell class - Simple code cell using CodeMirror
 */

class Cell {
    /**
     * Constructor
     * @param {string} content - Cell content
     * @param {string} description - Cell description
     * @param {CellManager} cellManager - Associated CellManager instance
     * @param {string} id - Optional cell ID
     * @param {string} cellType - Optional short description
     * @param {boolean} isInteraction - Whether this cell is an interaction
     * @param {boolean} allowClose - Whether this cell can be closed/deleted
     */
    constructor(content = '', description = 'Code Cell', cellManager = null, id = null, cellType = '', isInteraction = false, allowClose = true) {
        this.id = id || this.generateId();
        this.content = content;
        this.defaultContent = content;
        this.description = description;
        this.cellType = cellType;
        this.isInteraction = isInteraction;
        this.cellManager = cellManager;
        this.element = null;
        this.editor = null;
        this.outputElement = null;
        this.output = '';
        this.collapsed = true; // Default expanded state
        this.selected = false; // Track selection state
        this.allowClose = allowClose;
        this.runBtn = null; // Reference to run button element
    }

    /**
     * Generate a unique cell ID
     * @returns {string} Unique cell ID
     */
    generateId() {
        // Use timestamp + 8 random base36 characters
        const timestamp = Date.now().toString(36);
        const random = Math.random().toString(36).substring(2, 10); // 8 characters
        return `cell-${timestamp}-${random}`;
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
        
        header.addEventListener('click', (e) => {
            e.stopPropagation(); 
            this.toggleCollapse();
            if (!this.selected) {
                this.select();
            }
        });

        // Create description
        const description = document.createElement('div');
        description.className = 'cell-description';
        description.textContent = this.description;

        // Create action buttons
        const actions = document.createElement('div');
        actions.className = 'cell-actions';

        this.runBtn = document.createElement('button');
        this.runBtn.className = 'cell-btn run-btn';
        this.runBtn.textContent = '▶';
        this.runBtn.title = 'Run';
        this.runBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            this.run();
        });

        // Reset button
        const resetBtn = document.createElement('button');
        resetBtn.className = 'cell-btn reset-btn';
        resetBtn.textContent = '↻';
        resetBtn.title = 'Reset to default code';
        resetBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            this.reset();
        });

        actions.appendChild(this.runBtn);
        actions.appendChild(resetBtn);

        if (this.allowClose) {
            const deleteBtn = document.createElement('button');
            deleteBtn.className = 'cell-btn delete-btn';
            deleteBtn.textContent = '✕';
            deleteBtn.title = 'Delete';
            deleteBtn.addEventListener('click', (e) => {
                e.stopPropagation();
                this.delete();
            });
            actions.appendChild(deleteBtn);
        }

        header.appendChild(description);
        header.appendChild(actions);

        // Create content area
        const content = document.createElement('div');
        content.className = 'cell-content';

        // Create editor container
        const editorContainer = document.createElement('div');
        editorContainer.className = 'cell-editor';

        content.appendChild(editorContainer);

        // Create output area
        this.outputElement = document.createElement('div');
        this.outputElement.className = 'cell-output';
        this.outputElement.style.display = 'none';

        this.element.appendChild(header);
        this.element.appendChild(content);
        this.element.appendChild(this.outputElement);

        // Add click event for cell selection
        // Simplified version - only add to the cell element
        this.element.addEventListener('click', (e) => {
            if (!this.selected) {
                this.select();
            }
        });

        // Initialize CodeMirror editor directly on the container
        this.initEditor(editorContainer);

        // Apply initial collapse state
        if (this.collapsed) {
            const content = this.element.querySelector('.cell-content');
            const output = this.outputElement;
            content.style.display = 'none';
            output.style.display = 'none';
        }

        return this.element;
    }

    /**
     * Initialize CodeMirror editor
     */
    initEditor(container) {
        this.editor = CodeMirror(container, {
            mode: 'julia',
            theme: 'default',
            lineNumbers: false,
            lineWrapping: true,
            viewportMargin: Infinity,
            indentUnit: 4,
            smartIndent: true,
            value: this.content, // Set initial content directly
            extraKeys: {
                'Shift-Enter': () => this.run(),
                'Ctrl-Enter': () => this.run()
            },
            // Add these options to ensure proper initialization
            //autoFocus: true,
            //readOnly: false,
            placeholder: 'Enter code here...'
        });

        // Force editor to render content immediately
        this.editor.refresh();
        
        // Ensure the editor has the correct dimensions
        this.editor.setSize('100%', '100%');
        
        // Explicitly set the value again to ensure it's displayed
        this.editor.setValue(this.content);
        
        // Force another refresh to make sure everything is rendered
        this.editor.refresh();

        // Listen for content changes
        this.editor.on('change', () => {
            this.content = this.editor.getValue();
            // Enable run button when content changes
            if (this.runBtn) {
                this.runBtn.disabled = false;
            }
        });

        // Set appropriate font size and height
        this.editor.setSize('100%', 'auto');
        this.editor.refresh();
    }

    /**
     * Reset cell to default content
     */
    reset() {
        this.setValue(this.defaultContent);
        // Clear output when resetting
        this.outputElement.style.display = 'none';
        this.outputElement.innerHTML = '';
        this.output = '';
        // Enable run button when resetting
        if (this.runBtn) {
            this.runBtn.disabled = false;
        }
    }

    /**
     * Run cell code
     */
    run() {
        // Disable run button when running
        if (this.runBtn) {
            this.runBtn.disabled = true;
        }
        this.cellManager.runCell(this);
    }

    /**
     * Delete cell
     */
    delete() {
        // First notify CellManager to remove this cell from its management
        if (this.cellManager) {
            this.cellManager.removeCell(this);
        }
        
        // Then remove the element from DOM
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
            
            // Refresh CodeMirror when expanding
            setTimeout(() => {
                if (this.editor) {
                    this.editor.refresh();
                    // Ensure the editor has the correct dimensions
                    this.editor.setSize('100%', 'auto');
                    this.editor.refresh();
                }
            }, 10); // Small delay to ensure DOM is updated
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

    /**
     * Select this cell
     */
    select() {
            console.log('Cell.select() called for cell:', this.id);
            if (this.cellManager) {
                console.log('Calling cellManager.selectCellById:', this.id);
                this.cellManager.selectCellById(this.id);
            } else {
                console.log('Cell has no cellManager reference');
            }
        
    }
}

export { Cell };

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { Cell };
} else {
    window.Cell = Cell;
}