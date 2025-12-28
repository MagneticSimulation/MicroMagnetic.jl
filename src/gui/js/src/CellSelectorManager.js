// CellSelectorManager.js

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
        defaultOption.value = 'empty_cell';
        defaultOption.textContent = 'Empty Cell';
        selector.appendChild(defaultOption);
        
        // If there's a selected cell, add appropriate options
        if (selectedCell) {
            if (selectedCell.cellType === 'mesh') {                
                const fdmeshOpenOption = document.createElement('option');
                fdmeshOpenOption.value = 'fdmesh_open';
                fdmeshOpenOption.textContent = 'FDMesh (Open)';
                selector.appendChild(fdmeshOpenOption);
                
                const fdmeshPbcOption = document.createElement('option');
                fdmeshPbcOption.value = 'fdmesh_pbc';
                fdmeshPbcOption.textContent = 'FDMesh (PBC)';
                selector.appendChild(fdmeshPbcOption);
            } else if (selectedCell.cellType === 'ms') {                
                const msConstantOption = document.createElement('option');
                msConstantOption.value = 'ms_constant';
                msConstantOption.textContent = 'Constant Ms';
                selector.appendChild(msConstantOption);
                
                const msSpatialOption = document.createElement('option');
                msSpatialOption.value = 'ms_spatial';
                msSpatialOption.textContent = 'Spatial Ms (x,y,z)';
                selector.appendChild(msSpatialOption);
                
                const msGridOption = document.createElement('option');
                msGridOption.value = 'ms_grid';
                msGridOption.textContent = 'Grid Ms (i,j,k)';
                selector.appendChild(msGridOption);
            } else if (selectedCell.cellType === 'sim') {                
                const simSdOption = document.createElement('option');
                simSdOption.value = 'sim_sd';
                simSdOption.textContent = 'Sim (SD Driver)';
                selector.appendChild(simSdOption);
                
                const simLlgOption = document.createElement('option');
                simLlgOption.value = 'sim_llg';
                simLlgOption.textContent = 'Sim (LLG Driver)';
                selector.appendChild(simLlgOption);
            } else if (selectedCell.cellType === 'm0') {                
                const m0ConstantOption = document.createElement('option');
                m0ConstantOption.value = 'm0_constant';
                m0ConstantOption.textContent = 'Constant m0';
                selector.appendChild(m0ConstantOption);
                
                const m0GridOption = document.createElement('option');
                m0GridOption.value = 'm0_grid';
                m0GridOption.textContent = 'Grid m0 (i,j,k)';
                selector.appendChild(m0GridOption);
                
                const m0SpatialOption = document.createElement('option');
                m0SpatialOption.value = 'm0_spatial';
                m0SpatialOption.textContent = 'Spatial m0 (x,y,z)';
                selector.appendChild(m0SpatialOption);
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
        } else if (selectedValue === 'ms_constant' || selectedValue === 'ms_spatial' || selectedValue === 'ms_grid') {
            // Handle Ms template selection - call updateCellContent
            console.log('Ms template selected:', selectedValue);
            this.updateCellContent(selectedValue);
        } else if (selectedValue === 'sim_sd' || selectedValue === 'sim_llg') {
            // Handle Sim template selection - call updateCellContent
            console.log('Sim template selected:', selectedValue);
            this.updateCellContent(selectedValue);
        } else if (selectedValue === 'm0_constant' || selectedValue === 'm0_grid' || selectedValue === 'm0_spatial') {
            // Handle m0 template selection - call updateCellContent
            console.log('m0 template selected:', selectedValue);
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

// Export the CellSelectorManager class and cellTemplates dictionary
export { CellSelectorManager, cellTemplates };