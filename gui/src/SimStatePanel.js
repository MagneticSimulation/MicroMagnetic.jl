class SimStatePanel {
    constructor(containerId) {
        this.container = document.getElementById(containerId);
        this.defaultHtml = '<div class="state-empty">No simulation loaded.</div>';
        
        // Only set innerHTML if container exists
        if (this.container) {
            this.container.innerHTML = this.defaultHtml;
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
        
        let html = '<div class="state-list">';
        if (state.mesh) {
            html += `<div class="state-row"><span class="state-key">Mesh:</span> ${state.mesh}</div>`;
        }
        if (state.sim_name) {
            html += `<div class="state-row"><span class="state-key">Simulation:</span> ${state.sim_name}</div>`;
        }
        if (state.driver) {
            html += `<div class="state-row"><span class="state-key">Driver:</span> ${state.driver}</div>`;
        }
        if (state.Ms_max !== undefined) {
            const msVal = Number(state.Ms_max).toExponential(2);
            html += `<div class="state-row"><span class="state-key">Ms max:</span> ${msVal} A/m</div>`;
        }
        if (state.m0) {
            html += `<div class="state-row"><span class="state-key">Initial M:</span> ${state.m0}</div>`;
        }
        if (state.interactions && state.interactions.length > 0) {
            html += '<div class="state-row"><span class="state-key">Interactions:</span><ul class="interaction-list">';
            state.interactions.forEach(inter => {
                html += `<li>${inter}</li>`;
            });
            html += '</ul></div>';
        }
        html += '</div>';
        
        this.container.innerHTML = html;
    }
}

export default SimStatePanel;