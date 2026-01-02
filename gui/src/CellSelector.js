// CellSelector.js
import { FDMeshPanel } from './FDMeshPanel.js';

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
        this.fdMeshPanel = null;
        this.showPanelBtn = null;
        this.lastSelectedCell = null; 
        
        // Initialize features
        this.initializeCellSelector();
        this.initializeFDMeshPanel();
        this.setupCellListeners();
        
        console.log('CellSelector initialized');
    }
    
    /**
     * Initialize Cell selector functionality
     */
    initializeCellSelector() {
        if (!this.selectorElement) {
            console.warn('Cell selector element not found');
            return;
        }
        
        // Add event listener for selector changes
        this.selectorElement.addEventListener('change', () => this.handleSelectorChange());
        
        // Add event listener for show panel button
        this.showPanelBtn = document.getElementById('show-panel-btn');
        if (this.showPanelBtn) {
            this.showPanelBtn.addEventListener('click', () => {
                this.showPanel();
            });
        }
        
        // Initialize with default options
        this.updateOptions();
    }
    
    /**
     * Initialize FDMesh panel
     */
    initializeFDMeshPanel() {
        this.fdMeshPanel = new FDMeshPanel(this.cellManager);
        this.fdMeshPanel.init();
    }
    
    /**
     * Set up cell listeners
     */
    setupCellListeners() {
        // Add selection change listener to update cell selector options
        this.cellManager.addSelectionChangeListener((cell) => {
            this.updateOptions();
            this.onCellSelected(cell);
        });
    }
    
    /**
     * Update the cellManager reference
     * @param {CellManager} newCellManager - New CellManager instance
     */
    updateCellManager(newCellManager) {
        this.cellManager = newCellManager;
        this.setupCellListeners();
        this.updateOptions();
    }
    
    /**
     * Update options in the cell selector dropdown
     */
    updateOptions() {
        const selector = this.selectorElement;
        if (!selector) return;
        
        // Clear existing options
        selector.innerHTML = '';
        
        // Add default option
        const defaultOption = document.createElement('option');
        defaultOption.value = '';
        defaultOption.textContent = '-- Select Cell --';
        selector.appendChild(defaultOption);
        
        // Add template options based on cell type
        const selectedCell = this.cellManager.getSelectedCell();
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
            return;
        }
        
        const template = cellTemplates[templateKey];
        if (!template) {
            return;
        }
        
        // Update cell content with template
        selectedCell.content = template.content;
        selectedCell.cellType = template.cellType;
        
        if (selectedCell.editor) {
            selectedCell.editor.setValue(template.content);
        }
        
        // Refresh the cell selector options and update selected option
        this.updateOptions();
        this.updateSelectedOption();
        
        this.lastSelectedCell = selectedCell;
        
        // Refresh panel data if it's already visible
        if (this.fdMeshPanel && this.fdMeshPanel.container) {
            this.fdMeshPanel.refreshFromCellCode(selectedCell);
        }
        
        console.log('Cell content updated with template:', templateKey);
    }
    
    updateSelectedOption() {
        // Set the default selected option to the first option (which is "-- Select Cell --")
        if (this.selectorElement && this.selectorElement.options.length > 0) {
            this.selectorElement.selectedIndex = 0;
        }
        
        // Update button states after setting the selected option
        this.updateButtonStates();
    }
    
    /**
     * Show panel for the current selected cell
     */
    showPanel() {
        const selectedCell = this.cellManager.getSelectedCell() || this.lastSelectedCell;
        if (!selectedCell) {
            console.warn('No cell selected to show panel for');
            return;
        }
        
        // Show appropriate panel based on cell type
        if (selectedCell.cellType === 'mesh') {
            if (this.fdMeshPanel) {
                this.fdMeshPanel.showPanel(selectedCell);
            }
        } else {
            console.log('No panel available for cell type:', selectedCell.cellType);
        }
    }
    
    /**
     * Handle cell selection changes
     */
    onCellSelected(cell) {
        this.lastSelectedCell = cell;
        
        // Refresh panel data if it's already visible
        if (cell && this.fdMeshPanel) {
            this.fdMeshPanel.refreshFromCellCode(cell);
        }
    }
}

export { CellSelector, cellTemplates };