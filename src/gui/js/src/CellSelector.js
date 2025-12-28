// CellSelector.js

// Cell templates dictionary - moved outside the class for better scalability
const cellTemplates = {
    fdmesh_open: {
        content: `mesh = FDMesh(nx=10, ny=10, nz=10, dx=1e-9, dy=1e-9, dz=1e-9)`,
        cellType: "mesh"
    },
    fdmesh_pbc: {
        content: `mesh = FDMesh(nx=10, ny=10, nz=10, dx=1e-9, dy=1e-9, dz=1e-9, pbc="xy")`,
        cellType: "mesh"
    },
    ms_constant: {
        content: `set_Ms(sim, 8e5) #A/m`,
        cellType: "ms"
    },
    ms_spatial: {
        content: `function circular_Ms(x, y, z)
    if x^2 + y^2 <= (50nm)^2
        return 8e5
    end
    return 0.0
end
set_Ms(sim, circular_Ms)`,
        cellType: "ms"
    },
    ms_grid: {
        content: `function circular_Ms(i,j,k,dx,dy,dz)
    if (i-50.5)^2 + (j-50.5)^2 <= 50^2
        return 8e5
    end
    return 0.0
end
set_Ms(sim, circular_Ms)`,
        cellType: "ms"
    },
    sim_sd: {
        content: `sim = Sim(mesh; driver="SD", name="std")`,
        cellType: "sim"
    },
    sim_llg: {
        content: `sim = Sim(mesh; driver="LLG", name="std", integrator="DormandPrince")
sim.driver.alpha = 0.05
sim.driver.gamma = 2.21e5
sim.driver.integrator.tol = 1e-6`,
        cellType: "sim"
    },
    m0_constant: {
        content: `init_m0(sim, (1, 0.25, 0.1))`,
        cellType: "m0"
    },
    m0_grid: {
        content: `function init_fun(i, j, k, dx, dy, dz)
    x = i - 10
    y = j - 10
    r = (x^2 + y^2)^0.5
    if r < 2
        return (0, 0, 1)
    end
    return (-y / r, x / r, 0)
end
init_m0(sim, init_fun)`,
        cellType: "m0"
    },
    m0_spatial: {
        content: `function init_fun(x, y, z)
    r = (x^2 + y^2)^0.5
    if r < 10e-9
        return (0, 0, 1)
    end
    return (-y / r, x / r, 0)
end
init_m0(sim, init_fun)`,
        cellType: "m0"
    }
};

class CellSelector {
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
        defaultOption.value = 'empty_cell';
        defaultOption.textContent = 'Empty Cell';
        selector.appendChild(defaultOption);
        
        // Add template options based on cell type
        if (selectedCell && selectedCell.cellType) {
            Object.entries(cellTemplates).forEach(([key, template]) => {
                if (template.cellType === selectedCell.cellType) {
                    const option = document.createElement('option');
                    option.value = key;
                    option.textContent = this.capitalizeFirstLetter(key.replace(/_/g, ' '));
                    selector.appendChild(option);
                }
            });
        }
        
        // Update selected option and button states
        this.updateSelectedOption();
    }
    
    capitalizeFirstLetter(string) {
        return string.charAt(0).toUpperCase() + string.slice(1);
    }
    
    handleSelectorChange() {
        const selectedValue = this.selectorElement.value;
        console.log('Selected template:', selectedValue);
        
        // Update cell content if a template is selected
        if (selectedValue && selectedValue !== 'empty_cell') {
            this.updateCellContent(selectedValue);
        }
        
        // Update button states based on the selected value
        this.updateButtonStates();
    }
    
    updateButtonStates() {
        const selectedValue = this.selectorElement.value;
        const addButton = document.getElementById('add-cell-btn');
        
        // Update button states based on selected value
        if (selectedValue === 'empty_cell') {
            // Allow adding cells when empty cell is selected
            if (addButton) addButton.disabled = false;
        } else {
            if (addButton) addButton.disabled = true;
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

// Export the CellSelector class and cellTemplates dictionary
export { CellSelector, cellTemplates };