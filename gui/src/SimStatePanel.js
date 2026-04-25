import TreeView from './TreeView.js';

class SimStatePanel {
    constructor(containerId) {
        this.container = document.getElementById(containerId);
        this.defaultHtml = '<div class="sim-state-empty">No simulation loaded.</div>';
        
        // Only initialize if container exists
        if (this.container) {
            this.container.innerHTML = this.defaultHtml;
            this.treeView = new TreeView(this.container);
        } else {
            console.error(`SimStatePanel container with ID '${containerId}' not found.`);
        }
    }

    update(state) {
        // Only update if container exists
        if (!this.container) {
            console.error('SimStatePanel container not initialized.');
            return;
        }
        
        if (!state || Object.keys(state).length === 0) {
            this.container.innerHTML = this.defaultHtml;
            return;
        }
        
        // Use TreeView to render the state data
        this.treeView.render(state);
    }
}

export default SimStatePanel;