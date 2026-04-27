/**
 * PlotPanel - Floating plot window using Plotly.js
 */
class PlotPanel {
    constructor() {
        this.container = null;
        this.titleBar = null;
        this.content = null;
        this.plotDiv = null;
        this.isCollapsed = false;
        this.isDragging = false;
        this.dragOffset = { x: 0, y: 0 };
        this.currentData = null;
        this.mirrorLoop = false;
        this._visible = false;

        this.createDOM();
        this.setupEventListeners();
    }

    get visible() {
        return this._visible;
    }

    set visible(value) {
        this._visible = value;
        if (value) {
            this.show();
        } else {
            this.hide();
        }
        if (this._guiController) {
            this._guiController.updateDisplay();
        }
    }

    createDOM() {
        this.container = document.createElement('div');
        this.container.style.cssText = `
            position: fixed; bottom: 20px; right: 20px;
            width: 400px; height: 300px; z-index: 1000;
            display: none; flex-direction: column;
            background: #1a1a2e; border-radius: 8px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.4); overflow: hidden;
        `;

        this.titleBar = document.createElement('div');
        this.titleBar.style.cssText = `
            display: flex; justify-content: space-between; align-items: center;
            padding: 8px 12px; background: #16213e; color: #fff; cursor: move;
        `;

        const title = document.createElement('span');
        title.textContent = 'Plot';

        this.collapseBtn = document.createElement('button');
        this.collapseBtn.textContent = '−';
        this.collapseBtn.style.cssText = 'background:none; border:none; color:#fff; font-size:18px; cursor:pointer;';

        this.closeBtn = document.createElement('button');
        this.closeBtn.textContent = '×';
        this.closeBtn.style.cssText = 'background:none; border:none; color:#fff; font-size:20px; cursor:pointer;';

        this.titleBar.appendChild(title);
        this.titleBar.appendChild(this.collapseBtn);
        this.titleBar.appendChild(this.closeBtn);

        this.content = document.createElement('div');
        this.content.style.cssText = 'flex:1; overflow:hidden;';

        this.plotDiv = document.createElement('div');
        this.plotDiv.style.cssText = 'width:100%; height:100%;';

        this.content.appendChild(this.plotDiv);
        this.container.appendChild(this.titleBar);
        this.container.appendChild(this.content);
        document.body.appendChild(this.container);
    }

    setupEventListeners() {
        this.closeBtn.addEventListener('click', () => this.hide());
        this.titleBar.addEventListener('dblclick', () => this.toggleCollapse());
        this.collapseBtn.addEventListener('click', () => this.toggleCollapse());

        this.titleBar.addEventListener('mousedown', (e) => {
            if (e.target === this.closeBtn || e.target === this.collapseBtn) return;
            this.isDragging = true;
            this.dragOffset.x = e.clientX - this.container.offsetLeft;
            this.dragOffset.y = e.clientY - this.container.offsetTop;
        });

        document.addEventListener('mousemove', (e) => {
            if (!this.isDragging) return;
            this.container.style.left = `${e.clientX - this.dragOffset.x}px`;
            this.container.style.top = `${e.clientY - this.dragOffset.y}px`;
            this.container.style.right = 'auto';
            this.container.style.bottom = 'auto';
        });

        document.addEventListener('mouseup', () => this.isDragging = false);
    }

    toggleCollapse() {
        this.isCollapsed = !this.isCollapsed;
        this.content.style.display = this.isCollapsed ? 'none' : 'flex';
        this.collapseBtn.textContent = this.isCollapsed ? '+' : '−';
        this.container.style.height = this.isCollapsed ? 'auto' : '350px';
    }

    show() {
        this._visible = true;
        this.container.style.display = 'flex';
    }

    hide() {
        this._visible = false;
        this.container.style.display = 'none';
    }

    setMirrorLoop(enabled) {
        this.mirrorLoop = enabled;
        if (this.currentData) this.update(this.currentData);
    }

    clear() {
        this.currentData = null;
        this._visible = false;
        this.hide();
    }

    update(data) {
        if (!data) { this.hide(); return; }
        this.show();
        this.currentData = data;
        if (data.title) this.titleBar.querySelector('span').textContent = data.title;

        const mirrorLoop = data.mirror_loop !== undefined ? data.mirror_loop : this.mirrorLoop;
        // Create a copy of traces to avoid modifying original data
        let traces = data.traces ? JSON.parse(JSON.stringify(data.traces)) : [{ x: data.x || [], y: data.y || [], mode: data.mode || 'lines+markers', name: data.name || 'trace' }];

        if (mirrorLoop && traces.length > 0) {
            const orig = traces[0];
            traces.push({
                x: orig.x.map(v => -v), y: orig.y.map(v => -v),
                mode: orig.mode || 'lines+markers', showlegend: false
            });
        }

        const xTitle = typeof data.layout?.xaxis?.title === 'string' ? data.layout.xaxis.title : (data.layout?.xaxis?.title?.text || 'X');
        const yTitle = typeof data.layout?.yaxis?.title === 'string' ? data.layout.yaxis.title : (data.layout?.yaxis?.title?.text || 'Y');

        const layout = {
            paper_bgcolor: 'transparent', plot_bgcolor: '#1a1a2e',
            font: { color: '#e0e0e0', family: 'Open Sans, sans-serif', size: 12 },
            margin: { t: 40, r: 20, b: 50, l: 60 },
            showlegend: false,
            xaxis: { title: { text: xTitle, font: { color: '#64b5f6', size: 13 } },
                tickfont: { color: '#b0b0b0', size: 11 }, gridcolor: 'rgba(255,255,255,0.1)', linecolor: 'rgba(255,255,255,0.3)' },
            yaxis: { title: { text: yTitle, font: { color: '#81c784', size: 13 } },
                tickfont: { color: '#b0b0b0', size: 11 }, gridcolor: 'rgba(255,255,255,0.1)', linecolor: 'rgba(255,255,255,0.3)' }
        };

        const config = {
            responsive: true, displayModeBar: true, displaylogo: false,
            modeBarButtonsToRemove: ['select2d', 'lasso2d'],
            modeBarButtons: [['zoom2d', 'resetScale2d'], ['zoomIn2d', 'zoomOut2d'], ['toImage']]
        };

        try { Plotly.react(this.plotDiv, traces, layout, config); }
        catch (e) { console.error('PlotPanel update error:', e); }
    }
}

export default PlotPanel;