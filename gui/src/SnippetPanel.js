class SnippetPanel {
    constructor(containerId, codeEditor) {
        this.container = document.getElementById(containerId);
        this.codeEditor = codeEditor;
        this.snippets = [];
        this.examples = [];
        this.snippetsActiveTag = 'all';
        this.examplesActiveTag = 'all';
        this.mode = 'snippets'; // 'snippets' or 'examples'
        this.buildUI();
        this.loadSnippets();
    }

    buildUI() {
        this.container.innerHTML = `
            <div class="snippet-panel">
                <div class="snippet-header">
                    <div class="snippet-tabs">
                        <button class="snippet-tab active" data-mode="snippets">Snippets</button>
                        <button class="snippet-tab" data-mode="examples">Examples</button>
                    </div>
                    <input type="text" id="snippet-search" placeholder="Search..." />
                </div>
                <div id="snippet-tags" class="snippet-tags"></div>
                <div id="snippet-list" class="snippet-list"></div>
            </div>
            <!-- Example Modal -->
            <div id="example-modal" class="example-modal">
                <div class="example-modal-content">
                    <div class="example-modal-header">
                        <h4 id="example-modal-title"></h4>
                        <span class="example-modal-close">&times;</span>
                    </div>
                    <div id="example-modal-body" class="example-modal-body"></div>
                    <div class="example-modal-footer">
                        <button id="example-insert-all">Insert All</button>
                    </div>
                </div>
            </div>
        `;
        this.searchInput = this.container.querySelector('#snippet-search');
        this.tagsContainer = this.container.querySelector('#snippet-tags');
        this.listContainer = this.container.querySelector('#snippet-list');
        this.searchInput.addEventListener('input', () => this.filterAndRender());
        
        // Tab switching event listeners
        this.container.querySelectorAll('.snippet-tab').forEach(tab => {
            tab.addEventListener('click', () => {
                this.switchMode(tab.dataset.mode);
            });
        });
        
        // Modal close event
        this.modal = this.container.querySelector('#example-modal');
        this.modalClose = this.container.querySelector('.example-modal-close');
        this.modalInsertAll = this.container.querySelector('#example-insert-all');
        
        this.modalClose.addEventListener('click', () => {
            this.modal.style.display = 'none';
        });
        
        // Click outside modal to close
        window.addEventListener('click', (e) => {
            if (e.target === this.modal) {
                this.modal.style.display = 'none';
            }
        });
        
        // Insert all button event
        this.modalInsertAll.addEventListener('click', () => {
            this.insertAllExampleSteps();
        });
    }

    async loadSnippets() {
        try {
            const response = await fetch('/api/snippets');
            if (!response.ok) throw new Error('Failed to load snippets');
            this.snippets = await response.json();
            if (this.mode === 'snippets') {
                this.renderTags();
                this.filterAndRender();
            }
        } catch (err) {
            if (this.mode === 'snippets') {
                this.listContainer.innerHTML = '<div class="snippet-error">Could not load snippets.</div>';
            }
            console.error(err);
        }
    }

    async loadExamples() {
        try {
            const response = await fetch('/api/examples');
            if (!response.ok) throw new Error('Failed to load examples');
            this.examples = await response.json();
            this.renderTags();
            this.filterAndRender();
        } catch (err) {
            this.listContainer.innerHTML = '<div class="snippet-error">Could not load examples.</div>';
            console.error(err);
        }
    }

    switchMode(mode) {
        this.mode = mode;
        
        // Update tab UI
        this.container.querySelectorAll('.snippet-tab').forEach(tab => {
            tab.classList.toggle('active', tab.dataset.mode === mode);
        });
        
        // Load data if needed
        if (mode === 'examples' && this.examples.length === 0) {
            this.loadExamples();
        } else {
            this.renderTags();
            this.filterAndRender();
        }
    }

    renderTags() {
        const items = this.mode === 'snippets' ? this.snippets : this.examples;
        const activeTag = this.mode === 'snippets' ? this.snippetsActiveTag : this.examplesActiveTag;
        const tagsSet = new Set();
        items.forEach(item => item.tags.forEach(tag => tagsSet.add(tag)));
        const tags = ['all', ...Array.from(tagsSet).sort()];
        
        this.tagsContainer.innerHTML = tags.map(tag =>
            `<span class="snippet-tag ${activeTag===tag?'active':''}" data-tag="${tag}">${tag}</span>`
        ).join('');
        
        this.tagsContainer.querySelectorAll('.snippet-tag').forEach(el => {
            el.addEventListener('click', () => {
                if (this.mode === 'snippets') {
                    this.snippetsActiveTag = el.dataset.tag;
                } else {
                    this.examplesActiveTag = el.dataset.tag;
                }
                this.filterAndRender();
            });
        });
    }

    filterAndRender() {
        const query = this.searchInput.value.toLowerCase();
        const items = this.mode === 'snippets' ? this.snippets : this.examples;
        const activeTag = this.mode === 'snippets' ? this.snippetsActiveTag : this.examplesActiveTag;
        
        const filtered = items.filter(item => {
            if (activeTag !== 'all' && !item.tags.includes(activeTag)) return false;
            if (query && !item.title.toLowerCase().includes(query) && !item.description.toLowerCase().includes(query)) return false;
            return true;
        });
        
        if (filtered.length === 0) {
            const emptyText = this.mode === 'snippets' ? 'No snippets found.' : 'No examples found.';
            this.listContainer.innerHTML = `<div class="snippet-empty">${emptyText}</div>`;
            return;
        }
        
        if (this.mode === 'snippets') {
            // Render snippets
            this.listContainer.innerHTML = filtered.map(s => `
                <div class="snippet-item" data-snippet-id="${s.title}">
                    <div class="snippet-title">${s.title}</div>
                    <div class="snippet-desc">${s.description}</div>
                    <div class="snippet-code-preview"><code>${this.escapeHtml(s.code.substring(0, 80))}${s.code.length>80?'...':''}</code></div>
                </div>
            `).join('');

            // Bind click events for snippets
            this.listContainer.querySelectorAll('.snippet-item').forEach(el => {
                el.addEventListener('click', () => {
                    const title = el.dataset.snippetId;
                    const snippet = this.snippets.find(s => s.title === title);
                    this.handleSnippetClick(snippet);
                });
            });
        } else {
            // Render examples
            this.listContainer.innerHTML = filtered.map(example => `
                <div class="snippet-item" data-example-id="${example.title}">
                    <div class="snippet-title">${example.title}</div>
                    <div class="snippet-desc">${example.description}</div>
                    <div class="snippet-code-preview">
                        <code>${example.steps.length} steps: ${this.escapeHtml(example.steps[0].code.substring(0, 60))}${example.steps[0].code.length>60?'...':''}</code>
                    </div>
                    <div class="snippet-steps-count">${example.steps.length} steps</div>
                </div>
            `).join('');
            
            // Bind click events for examples
            this.listContainer.querySelectorAll('.snippet-item').forEach(el => {
                el.addEventListener('click', () => {
                    const title = el.dataset.exampleId;
                    const example = this.examples.find(e => e.title === title);
                    this.handleExampleClick(example);
                });
            });
        }
        
        // Sync tag highlighting
        this.tagsContainer.querySelectorAll('.snippet-tag').forEach(el => {
            el.classList.toggle('active', el.dataset.tag === activeTag);
        });
    }

    async handleSnippetClick(snippet) {
        let finalCode = snippet.code;
        if (snippet.arguments && snippet.arguments.length > 0) {
            const result = await this.showArgsDialog(snippet);
            if (!result) return;
            const { action, values } = result;
            for (const [key, val] of Object.entries(values)) {
                finalCode = finalCode.replace(new RegExp(`\\$\\{${key}\\}|\\$${key}`, 'g'), val);
            }
            if (action === 'replace') {
                if (this.codeEditor && this.codeEditor.clearEditor) {
                    this.codeEditor.clearEditor();
                }
            } else {
                finalCode = "\n" + finalCode;
            }
        } else {
            finalCode = "\n" + finalCode;
        }

        if (this.codeEditor && this.codeEditor.insertAtCursor) {
            this.codeEditor.insertAtCursor(finalCode);
        } else {
            console.warn('CodeEditorPanel insertAtCursor not available');
        }
    }

    showArgsDialog(snippet) {
        return new Promise((resolve) => {
            const overlay = document.createElement('div');
            overlay.className = 'snippet-args-overlay';
            let formHtml = `<div class="snippet-args-dialog"><h4>${snippet.title}</h4>`;
            snippet.arguments.forEach(arg => {
                formHtml += `
                    <label>${arg.name}
                        <input type="text" id="arg-${arg.name}" value="${arg.default ?? ''}" />
                    </label>`;
            });
            formHtml += `
                <div class="args-buttons">
                    <button id="args-cancel">Cancel</button>
                    <button id="args-replace">Replace</button>
                    <button id="args-ok">Insert</button>
                </div></div>`;
            overlay.innerHTML = formHtml;
            document.body.appendChild(overlay);
            const cancel = overlay.querySelector('#args-cancel');
            const ok = overlay.querySelector('#args-ok');
            const replace = overlay.querySelector('#args-replace');
            cancel.onclick = () => { document.body.removeChild(overlay); resolve(null); };
            ok.onclick = () => {
                const values = {};
                snippet.arguments.forEach(arg => {
                    values[arg.name] = overlay.querySelector(`#arg-${arg.name}`).value;
                });
                document.body.removeChild(overlay);
                resolve({ action: 'insert', values });
            };
            replace.onclick = () => {
                const values = {};
                snippet.arguments.forEach(arg => {
                    values[arg.name] = overlay.querySelector(`#arg-${arg.name}`).value;
                });
                document.body.removeChild(overlay);
                resolve({ action: 'replace', values });
            };
        });
    }
    
    handleExampleClick(example) {
        // Store current example for modal
        this.currentExample = example;
        
        // Update modal title
        const modalTitle = this.container.querySelector('#example-modal-title');
        modalTitle.textContent = example.title;
        
        // Update modal body with steps
        const modalBody = this.container.querySelector('#example-modal-body');
        modalBody.innerHTML = example.steps.map((step, index) => 
            `<div class="example-step">
                <div class="example-step-header">
                    <span class="example-step-number">${index + 1}</span>
                    <h6 class="example-step-title">${step.title}</h6>
                </div>
                <pre class="example-step-code">${this.escapeHtml(step.code)}</pre>
            </div>`
        ).join('');
        
        // Show modal
        this.modal.style.display = 'block';
    }
    
    insertAllExampleSteps() {
        if (!this.currentExample || !this.codeEditor || !this.codeEditor.insertAtCursor) {
            return;
        }

        if (this.codeEditor.clearEditor) {
            this.codeEditor.clearEditor();
        }

        const allCode = this.currentExample.steps.map(step => step.code).join('\n\n');

        this.codeEditor.insertAtCursor(allCode);

        this.modal.style.display = 'none';
    }

    escapeHtml(str) {
        const div = document.createElement('div');
        div.appendChild(document.createTextNode(str));
        return div.innerHTML;
    }
}

export default SnippetPanel;