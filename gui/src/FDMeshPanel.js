class FDMeshPanel {
    constructor(cellManager) {
        this.cellManager = cellManager;
        this.container = null;
        this.elements = {};
    }


    init() {
        this.createPanel();
        this.bindEvents();
    }


    createPanel() {
        if (document.getElementById('fdmesh-panel')) return;

        const panel = document.createElement('div');
        panel.className = 'control-section';
        panel.id = 'fdmesh-panel';
        panel.innerHTML = `
            <div id="fdmesh-header">
                <h3 id="fdmesh-title">FDMesh Configuration</h3>
                <button id="fdmesh-close">&times;</button>
            </div>
            <div id="fdmesh-content">
                <div class="fdmesh-section">
                    <small class="fdmesh-label">Grid Dimensions (nx, ny, nz)</small>
                    <div class="fdmesh-input-grid">
                        <input type="number" id="fdmesh-nx" class="fdmesh-input" value="10" min="1">
                        <input type="number" id="fdmesh-ny" class="fdmesh-input" value="10" min="1">
                        <input type="number" id="fdmesh-nz" class="fdmesh-input" value="10" min="1">
                    </div>
                </div>
                <div class="fdmesh-section">
                    <small class="fdmesh-label">Cell Sizes (dx, dy, dz) (nm)</small>
                    <div class="fdmesh-input-grid">
                        <input type="number" id="fdmesh-dx" class="fdmesh-input" value="2.5" min="0.1">
                        <input type="number" id="fdmesh-dy" class="fdmesh-input" value="2.5" min="0.1">
                        <input type="number" id="fdmesh-dz" class="fdmesh-input" value="2.5" min="0.1">
                    </div>
                </div>
                <div class="fdmesh-section">
                    <div class="fdmesh-pbc-row">
                        <small class="fdmesh-label">Periodic Boundary Conditions</small>
                        <div class="fdmesh-checkbox-group">
                            <label class="fdmesh-checkbox-label">
                                <input type="checkbox" id="fdmesh-pbc-x"> X
                            </label>
                            <label class="fdmesh-checkbox-label">
                                <input type="checkbox" id="fdmesh-pbc-y"> Y
                            </label>
                            <label class="fdmesh-checkbox-label">
                                <input type="checkbox" id="fdmesh-pbc-z"> Z
                            </label>
                        </div>
                    </div>
                </div>
            </div>
        `;

        document.body.appendChild(panel);
        
        this.container = panel;
        
        this.elements = {};
        
        const elementKeys = ['nx', 'ny', 'nz', 'dx', 'dy', 'dz', 'pbcX', 'pbcY', 'pbcZ', 'closeBtn'];
        const specialIds = {
            closeBtn: 'fdmesh-close',
            pbcX: 'fdmesh-pbc-x',
            pbcY: 'fdmesh-pbc-y',
            pbcZ: 'fdmesh-pbc-z'
        };
        
        elementKeys.forEach(key => {
            const id = specialIds[key] || `fdmesh-${key.toLowerCase()}`;
            this.elements[key] = document.getElementById(id);
        });
        
        this.bindCloseEvents();
        this.hidePanel();
    }
    
    bindCloseEvents() {
        if (this.elements.closeBtn) {
            this.elements.closeBtn.addEventListener('click', () => this.hidePanel());
        }
        
        document.addEventListener('keydown', (e) => {
            if (e.key === 'Escape' && this.container.classList.contains('expanded')) {
                this.hidePanel();
            }
        });
    }
    
    showPanel(selectedCell = null) {
        if (!selectedCell) return;

        // Refresh panel with selected cell data
        this.refreshFromCellCode(selectedCell);

        // Show the panel
        this.container.classList.add('expanded');
        this.isVisible = true;

        // Position panel relative to selected cell
        const cellRect = selectedCell.element.getBoundingClientRect();
        const panelRect = this.container.getBoundingClientRect();
        
        // Position panel to the right of the cell, aligned with top
        // Adjust position to ensure it stays within viewport
        let left = cellRect.right + 10;
        let top = cellRect.top;
        
        // Check if panel would go off right edge of viewport
        if (left + panelRect.width > window.innerWidth) {
            left = cellRect.left - panelRect.width - 10;
        }
        
        // Check if panel would go off bottom edge of viewport
        if (top + panelRect.height > window.innerHeight) {
            top = Math.max(10, window.innerHeight - panelRect.height - 10);
        }
        
        // Apply new position
        this.container.style.position = 'fixed';
        this.container.style.top = `${top}px`;
        this.container.style.left = `${left}px`;
        this.container.style.right = 'auto';
        this.container.style.bottom = 'auto';
        this.container.style.transform = 'none';
    }
    
    hidePanel() {
        if (this.container) {
            this.container.classList.remove('expanded');
        }
    }


    bindEvents() {
        const inputs = ['nx', 'ny', 'nz', 'dx', 'dy', 'dz', 'pbcX', 'pbcY', 'pbcZ'];        
        inputs.forEach(key => {
            const input = this.elements[key];
            if (input) {
                input.addEventListener('change', () => this.applyConfiguration());
                if (!key.startsWith('pbc')) {
                    input.addEventListener('input', () => this.applyConfiguration());
                }
            }
        });
    }

    /**
     * Apply configuration - generate code and add to cell
     */
    applyConfiguration() {
        const code = this.generateCode();
        
        const selectedCell = this.cellManager.getSelectedCell();
        if (selectedCell) {
            selectedCell.content = code;
            selectedCell.defaultContent = code;
            
            if (selectedCell.editor) {
                selectedCell.editor.setValue(code);
            }
        }
    }


    /**
     * Generate Julia code for FDMesh
     * @returns {string} - Julia code
     */
    generateCode() {
        const e = this.elements;
        const nx = parseInt(e.nx.value) || 10;
        const ny = parseInt(e.ny.value) || 10;
        const nz = parseInt(e.nz.value) || 10;
        const dx = parseFloat(e.dx.value) || 1;
        const dy = parseFloat(e.dy.value) || 1;
        const dz = parseFloat(e.dz.value) || 1;
        
        const pbc = [];
        if (e.pbcX.checked) pbc.push('x');
        if (e.pbcY.checked) pbc.push('y');
        if (e.pbcZ.checked) pbc.push('z');    
        const pbcStr = pbc.join('') || '';
        
        let code = `mesh = FDMesh(nx=${nx}, ny=${ny}, nz=${nz}, dx=${dx}e-9, dy=${dy}e-9, dz=${dz}e-9`;
        if (pbcStr.length > 0) {
            code += `, pbc="${pbcStr}"`;
        }
        code += ');';
        
        return code;
    }
    
    refreshFromCellCode(cell) {
        if (!cell || cell.cellType !== 'mesh') return;

        const e = this.elements;
        
        const code = cell.content;
        
        const nxMatch = code.match(/nx=(\d+)/i);
        const nyMatch = code.match(/ny=(\d+)/i);
        const nzMatch = code.match(/nz=(\d+)/i);
        
        if (nxMatch) e.nx.value = nxMatch[1];
        if (nyMatch) e.ny.value = nyMatch[1];
        if (nzMatch) e.nz.value = nzMatch[1];
        
        const dxMatch = code.match(/dx=([0-9.]+(?:e[-+]?[0-9]+)?)/i);
        const dyMatch = code.match(/dy=([0-9.]+(?:e[-+]?[0-9]+)?)/i);
        const dzMatch = code.match(/dz=([0-9.]+(?:e[-+]?[0-9]+)?)/i);
        
        if (dxMatch) e.dx.value = parseFloat(dxMatch[1]) * 1e9; 
        if (dyMatch) e.dy.value = parseFloat(dyMatch[1]) * 1e9; 
        if (dzMatch) e.dz.value = parseFloat(dzMatch[1]) * 1e9;
        
        const pbcMatch = code.match(/pbc="([xyz]*)"/i);
        const pbcStr = pbcMatch ? pbcMatch[1] : '';
        
        e.pbcX.checked = pbcStr.includes('x');
        e.pbcY.checked = pbcStr.includes('y');
        e.pbcZ.checked = pbcStr.includes('z');
    }
}

export { FDMeshPanel };

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { FDMeshPanel };
} else {
    window.FDMeshPanel = FDMeshPanel;
}