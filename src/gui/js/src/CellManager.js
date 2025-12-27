import { CellSelectorManager } from './CellSelectorManager.js';

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