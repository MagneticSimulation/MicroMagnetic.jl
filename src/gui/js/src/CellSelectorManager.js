/**
 * CellSelectorManager class - Manages cell selector functionality
 */
class CellSelectorManager {
    /**
     * Constructor
     * @param {CellManager} cellManager - Associated CellManager instance
     */
    constructor(cellManager) {
        this.cellManager = cellManager;
        this.selectorElement = document.getElementById('cell-selector');
        this.initEventListeners();
    }

    /**
     * Initialize event listeners
     */
    initEventListeners() {
        if (this.selectorElement) {
            this.selectorElement.addEventListener('change', this.handleSelectorChange.bind(this));
        }
    }

    /**
     * Handle selector change event
     */
    handleSelectorChange() {
        if (!this.selectorElement) return;
        
        const selectedId = this.selectorElement.value;
        console.log('Cell selector changed to:', selectedId);
        
        if (selectedId) {
            console.log('Calling cellManager.selectCellById with ID:', selectedId);
            this.cellManager.selectCellById(selectedId);
        }
    }

    /**
     * Update selector options based on the selected cell
     * For demo purposes, only show two options: default and the selected cell's name
     */
    updateOptions() {
        if (!this.selectorElement) return;
        
        this.selectorElement.innerHTML = '';
        
        // Add default option
        const defaultOption = document.createElement('option');
        defaultOption.value = '';
        defaultOption.textContent = '-- Select Cell --';
        this.selectorElement.appendChild(defaultOption);
        
        // Get the currently selected cell
        const selectedCell = this.cellManager.getSelectedCell();
        
        // If a cell is selected, add an option for it
        if (selectedCell) {
            const cellOption = document.createElement('option');
            cellOption.value = selectedCell.id;
            // Use cell name if available, otherwise use description or ID
            cellOption.textContent = selectedCell.name || selectedCell.description || `Cell ${selectedCell.id}`;
            this.selectorElement.appendChild(cellOption);
        }
        
        // Update selected state
        this.updateSelectedOption();
    }

    /**
     * Update selected option
     */
    updateSelectedOption() {
        if (!this.selectorElement) return;
        
        const selectedCell = this.cellManager.getSelectedCell();
        if (selectedCell) {
            this.selectorElement.value = selectedCell.id;
        } else {
            this.selectorElement.value = '';
        }
    }

    /**
     * Add specific attributes to cell for selector
     * @param {Cell} cell - Cell to add attributes to
     */
    addCellSpecificAttributes(cell) {
        // Add specific data attributes based on cell content or type
        if (cell.content.includes('relax')) {
            cell.element.setAttribute('data-task-type', 'relax');
        } else if (cell.content.includes('llg')) {
            cell.element.setAttribute('data-task-type', 'llg');
        }
    }
}

// Export CellSelectorManager class (ES modules)
export { CellSelectorManager };

// Export CellSelectorManager class (CommonJS)
if (typeof module !== 'undefined' && module.exports) {
    module.exports = CellSelectorManager;
}

// If in browser environment, mount to global object
if (typeof window !== 'undefined') {
    window.CellSelectorManager = CellSelectorManager;
}