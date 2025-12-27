// CellSelectorManager.js

// Cell templates dictionary - moved outside the class for better scalability
const cellTemplates = {
    fdmesh_open: {
        content: `using MicroMagnetic

mesh = FDMesh(nx=10, ny=10, nz=10, dx=1e-9, dy=1e-9, dz=1e-9)`,
        cellType: "mesh"
    },
    fdmesh_pbc: {
        content: `using MicroMagnetic

mesh = FDMesh(nx=10, ny=10, nz=10, dx=1e-9, dy=1e-9, dz=1e-9, pbc="xy")`,
        cellType: "mesh"
    }
};

class CellSelectorManager {
    constructor(cellManager) {
        this.cellManager = cellManager;
        this.selectorElement = document.getElementById('cell-selector');
        
        // Add event listener for selector changes
        this.selectorElement.addEventListener('change', () => this.handleSelectorChange());
        
        // Initialize with default options
        this.updateOptions();
    }
    
    updateOptions() {
        const selector = this.selectorElement;
        const selectedCell = this.cellManager.getSelectedCell();
        
        // Clear existing options
        selector.innerHTML = '';
        
        // Add default option
        const defaultOption = document.createElement('option');
        defaultOption.value = '';
        defaultOption.textContent = '-- Select Cell --';
        selector.appendChild(defaultOption);
        
        // If there's a selected cell, add appropriate options
        if (selectedCell) {
            if (selectedCell.cellType === 'mesh') {
                // For mesh cells, add special options
                const emptyOption = document.createElement('option');
                emptyOption.value = 'empty_cell';
                emptyOption.textContent = 'Empty Cell';
                selector.appendChild(emptyOption);
                
                const fdmeshOpenOption = document.createElement('option');
                fdmeshOpenOption.value = 'fdmesh_open';
                fdmeshOpenOption.textContent = 'FDMesh (Open)';
                selector.appendChild(fdmeshOpenOption);
                
                const fdmeshPbcOption = document.createElement('option');
                fdmeshPbcOption.value = 'fdmesh_pbc';
                fdmeshPbcOption.textContent = 'FDMesh (PBC)';
                selector.appendChild(fdmeshPbcOption);
            } else {
                // For regular cells, just show the cell name or ID
                const cellOption = document.createElement('option');
                cellOption.value = selectedCell.id;
                cellOption.textContent = selectedCell.name || selectedCell.description || selectedCell.id;
                selector.appendChild(cellOption);
            }
        }
        
        // Update button states based on selected option
        this.updateButtonStates();
    }
    
    handleSelectorChange() {
        const selectedValue = this.selectorElement.value;
        console.log('Cell selector changed:', selectedValue);
        
        // Handle different types of selections
        if (selectedValue === 'empty_cell') {
            // Handle empty cell selection
            console.log('Empty cell selected');
        } else if (selectedValue === 'fdmesh_open' || selectedValue === 'fdmesh_pbc') {
            // Handle mesh template selection - call updateCellContent
            console.log('Mesh template selected:', selectedValue);
            this.updateCellContent(selectedValue);
        } else if (selectedValue) {
            // Handle regular cell selection
            console.log('Cell selected:', selectedValue);
        }
        
        // Update button states based on selected option
        this.updateButtonStates();
    }
    
    updateButtonStates() {
        const selectedValue = this.selectorElement.value;
        const addButton = document.getElementById('add-cell-btn');
        const updateButton = document.getElementById('update-cell-btn');
        
        // Update button states based on selected value
        if (selectedValue === 'empty_cell') {
            // Allow adding cells when empty cell is selected
            if (addButton) addButton.disabled = false;
            if (updateButton) updateButton.disabled = true;
        } else if (selectedValue === 'fdmesh_open' || selectedValue === 'fdmesh_pbc') {
            // Disable add button when mesh templates are selected
            if (addButton) addButton.disabled = true;
            if (updateButton) updateButton.disabled = false;
        } else {
            // Default state
            if (addButton) addButton.disabled = false;
            if (updateButton) updateButton.disabled = false;
        }
    }
    
    updateCellContent(templateKey) {
        const selectedCell = this.cellManager.getSelectedCell();
        if (!selectedCell) {
            console.error('No cell selected to update');
            return;
        }
        
        // Get template from the cellTemplates dictionary
        const template = cellTemplates[templateKey];
        if (!template) {
            console.error('Template not found:', templateKey);
            return;
        }
        
        // Update cell content
        selectedCell.content = template.content;
        selectedCell.defaultContent = template.content;
        //selectedCell.description = template.description;
        selectedCell.cellType = template.cellType;
        
        // Update the cell's editor content
        if (selectedCell.editor) {
            selectedCell.editor.setValue(template.content);
        }
        
        console.log('Cell content updated with template:', templateKey);
    }
    
    updateSelectedOption() {
        // Set the default selected option to empty_cell if available
        const selector = this.selectorElement;
        for (let i = 0; i < selector.options.length; i++) {
            if (selector.options[i].value === 'empty_cell') {
                selector.selectedIndex = i;
                break;
            }
        }
        
        // Update button states after setting the selected option
        this.updateButtonStates();
    }
}

// Export the CellSelectorManager class and cellTemplates dictionary
export { CellSelectorManager, cellTemplates };