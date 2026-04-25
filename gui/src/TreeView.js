class TreeView {
    constructor(container) {
        this.container = container;
        this.initStyles();
    }

    // Inject CSS styles for the tree view
    initStyles() {
        // Check if styles are already injected
        if (document.getElementById('tree-view-styles')) {
            return;
        }

        const style = document.createElement('style');
        style.id = 'tree-view-styles';
        style.textContent = `
            .tree-view {
                font-family: monospace;
                font-size: 12px;
                color: #333;
                margin: 0;
                padding: 0;
            }
            .tree-list {
                list-style: none;
                margin: 0;
                padding: 0;
            }
            .tree-node {
                list-style: none;
                margin: 0;
                padding: 2px 0;
                position: relative;
            }
            .tree-toggle {
                cursor: pointer;
                user-select: none;
                margin-right: 4px;
                display: inline-block;
                width: 20px;
                text-align: center;
                color: #666;
                font-size: 20px;
                line-height: 1;
            }
            .tree-toggle:hover {
                color: #007bff;
            }
            .tree-key {
                font-weight: 600;
                color: #881280;
            }
            .tree-value {
                color: #1a1aa6;
            }
            .tree-closure {
                padding-left: 1.2em; /* Align with the key name */
                font-weight: normal;
                color: #666;
            }
            .tree-children {
                margin-left: 1.2em; /* Indent child nodes */
                padding-left: 0.8em;
            }
            .tree-bracket {
                color: #666;
            }
            .tree-index {
                color: #007bff;
                font-weight: 600;
            }
        `;
        document.head.appendChild(style);
    }

    // Render data as a tree
    render(data) {
        this.container.innerHTML = '';
        
        // Handle root data
        const rootList = document.createElement('ul');
        rootList.className = 'tree-list tree-view';
        
        if (typeof data === 'object' && data !== null) {
            // If root is an object, iterate through its keys
            const keys = Object.keys(data);
            keys.forEach(key => {
                const rootNode = this.buildNode(data[key], key, true, true); // true for isRoot
                rootList.appendChild(rootNode);
            });
        } else {
            // If root is a simple value
            const rootNode = this.buildNode(data, null, false, true);
            rootList.appendChild(rootNode);
        }
        
        this.container.appendChild(rootList);
    }

    // Build tree node recursively
    buildNode(value, key = '', isRoot = false, isLast = false) {
        const li = document.createElement('li');
        li.className = 'tree-node';

        // Handle null values
        if (value === null) {
            const keySpan = this.createKeySpan(key);
            li.appendChild(keySpan);
            const valueSpan = document.createElement('span');
            valueSpan.className = 'tree-value';
            valueSpan.textContent = 'null';
            li.appendChild(valueSpan);
            return li;
        }

        // Handle simple values (string, number, boolean)
        if (typeof value !== 'object') {
            const keySpan = this.createKeySpan(key);
            li.appendChild(keySpan);
            const valueSpan = document.createElement('span');
            valueSpan.className = 'tree-value';
            valueSpan.textContent = typeof value === 'string' ? `"${value}"` : value;
            li.appendChild(valueSpan);
            return li;
        }

        // Handle objects and arrays
        const isArray = Array.isArray(value);
        
        // Create toggle button
        const toggle = document.createElement('span');
        toggle.className = 'tree-toggle';
        toggle.textContent = isRoot ? '▾' : '▸'; // Root nodes are expanded by default
        li.appendChild(toggle);

        // Create key span
        const keySpan = this.createKeySpan(key);
        li.appendChild(keySpan);

        // Create opening bracket
        const openingBracket = document.createElement('span');
        openingBracket.className = 'tree-bracket';
        openingBracket.textContent = isArray ? '[' : '{';
        li.appendChild(openingBracket);

        // Create child nodes container
        const childrenUl = document.createElement('ul');
        childrenUl.className = 'tree-list tree-children';
        
        // Determine keys to iterate through
        const keys = isArray ? Array.from({length: value.length}, (_, i) => i) : Object.keys(value);
        
        // Add child nodes
        keys.forEach((k, idx) => {
            const childValue = isArray ? value[k] : value[k];
            const isChildLast = idx === keys.length - 1;
            const childNode = this.buildNode(childValue, k, false, isChildLast);
            childrenUl.appendChild(childNode);
        });
        
        li.appendChild(childrenUl);

        // Create closing bracket line
        const closeLi = document.createElement('li');
        closeLi.className = 'tree-closure';
        closeLi.textContent = isArray ? ']' : '}';
        li.appendChild(closeLi);

        // Set initial visibility state
        if (isRoot) {
            // Root nodes are fully expanded
            childrenUl.style.display = '';
            closeLi.style.display = '';
        } else {
            // Child nodes are collapsed by default
            childrenUl.style.display = 'none';
            closeLi.style.display = 'none';
        }

        // Add toggle functionality
        toggle.addEventListener('click', () => {
            const isVisible = childrenUl.style.display !== 'none';
            childrenUl.style.display = isVisible ? 'none' : '';
            closeLi.style.display = isVisible ? 'none' : '';
            toggle.textContent = isVisible ? '▸' : '▾';
        });

        return li;
    }

    // Create key span with proper styling
    createKeySpan(key) {
        if (!key && key !== 0) {
            const span = document.createElement('span');
            span.textContent = '';
            return span;
        }
        
        const span = document.createElement('span');
        if (typeof key === 'number') {
            span.className = 'tree-index';
            span.textContent = `${key}: `;
        } else {
            span.className = 'tree-key';
            span.textContent = `${key}: `;
        }
        return span;
    }
}

export default TreeView;