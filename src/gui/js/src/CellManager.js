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
        this.selectedCell = null; // Track selected cell
        this.selectionChangeListeners = []; // Listeners for selection changes
        
        if (this.container) {
            this.container.classList.add('cell-manager');
            this.init();
        }
    }

    /**
     * Initialize CellManager
     */
    init() {
        // Create header
        this.createHeader();
        
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
        
        // Add title to header
        header.appendChild(titleText);
        
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
                'New Cell',
                this // Pass CellManager instance
            );
        } else {
            // Ensure existing cell has reference to CellManager
            cell.cellManager = this;
        }
        
        const cellElement = cell.render();
        
        if (this.contentArea) {
            this.contentArea.appendChild(cellElement);
            this.cells.push(cell);
        }
    }

    /**
     * Insert a new cell below the currently selected cell
     * @param {Cell} cell - Cell instance to insert
     */
    insertCellBelowSelected(cell) {
        if (!cell) {
            // Create new cell
            cell = new Cell(
                '# New code cell\nprintln("Hello from new cell!")',
                'New Cell',
                this // Pass CellManager instance
            );
        } else {
            // Ensure existing cell has reference to CellManager
            cell.cellManager = this;
        }
        
        const cellElement = cell.render();
        
        if (this.contentArea && this.selectedCell) {
            // Find the index of the selected cell
            const selectedIndex = this.cells.findIndex(c => c.id === this.selectedCell.id);
            
            if (selectedIndex !== -1) {
                // Insert the new cell after the selected cell in the array
                this.cells.splice(selectedIndex + 1, 0, cell);
                
                // Insert the new cell element after the selected cell element in the DOM
                const selectedElement = this.contentArea.querySelector(`[data-cell-id="${this.selectedCell.id}"]`);
                if (selectedElement) {
                    selectedElement.parentNode.insertBefore(cellElement, selectedElement.nextSibling);
                }
            }
        } else if (this.contentArea) {
            // If no cell is selected, just add it to the end
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
     * Get a cell by its ID
     * @param {string} id - The ID of the cell to get
     * @returns {Cell|null} - The cell with the specified ID, or null if not found
     */
    getCellById(id) {
        return this.cells.find(cell => cell.id === id);
    }

    /**
     * Get selected cell
     */
    getSelectedCell() {
        return this.selectedCell;
    }

    /**
     * Add a listener for selection changes
     * @param {Function} listener - The listener function to call when selection changes
     */
    addSelectionChangeListener(listener) {
        if (typeof listener === 'function') {
            this.selectionChangeListeners.push(listener);
        }
    }

    /**
     * Remove a selection change listener
     * @param {Function} listener - The listener function to remove
     */
    removeSelectionChangeListener(listener) {
        this.selectionChangeListeners = this.selectionChangeListeners.filter(l => l !== listener);
    }

    /**
     * Notify all selection change listeners
     */
    notifySelectionChange() {
        this.selectionChangeListeners.forEach(listener => {
            try {
                listener(this.selectedCell);
            } catch (error) {
                console.error('Error in selection change listener:', error);
            }
        });
    }

    /**
     * Select a cell by its ID
     * @param {string} id - The ID of the cell to select
     */
    selectCellById(id) {
        // Update selected cell reference
        this.selectedCell = this.cells.find(cell => cell && cell.id === id);
        
        // Update UI to reflect selection
        this.updateSelectionUI();
        
        // Notify listeners of selection change
        this.notifySelectionChange();
    }

    /**
     * Update UI to reflect current selection state
     */
    updateSelectionUI() {
        // Update all cells to reflect their selection state
        this.cells.forEach((cell) => {
            console.log('Updating cell:', cell.description);
            if (cell && cell.element) {
                const isSelected = cell === this.selectedCell;
                cell.selected = isSelected; // Update cell's selected property
                cell.element.classList.toggle('selected', isSelected);
            }
        });
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
    /**
     * Remove a cell from the manager
     * @param {Cell} cell - The cell to remove
     */
    removeCell(cell) {
        if (cell) {
            // Remove from cells array
            this.cells = this.cells.filter(c => c && c.id !== cell.id);
            
            // If the removed cell was selected, clear selection
            if (this.selectedCell === cell) {
                this.selectedCell = null;
            }
            
            // Update selection UI
            this.updateSelectionUI();
        }
    }
}


// 导出 CellManager 类（ES 模块）
export { CellManager };

// 导出 CellManager 类（CommonJS）
if (typeof module !== 'undefined' && module.exports) {
    module.exports = { CellManager };
}

// 如果在浏览器环境中，挂载到全局对象
if (typeof window !== 'undefined') {
    window.CellManager = CellManager;
}