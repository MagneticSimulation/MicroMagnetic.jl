import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

/**
 * MagneticVisualization class for visualizing magnetization distributions
 */
class MagneticVisualization {
    constructor() {
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.controls = null;
        this.cube = null;
        this.animationId = null;
        this.arrows = [];
        this.gridSize = [10, 10, 10];
        this.dimensions = [10, 10, 10];
        this.arrowGroup = new THREE.Group();
        this.arrowPositions = null;
        this.codeEditor = null;
    }

    /**
     * Initialize visualization scene
     * @param {Object} container - DOM element containing renderer
     * @param {Object} data - Initialization data
     */
    init(container, data = null) {
        // Create scene
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0xf0f0f0);

        // Create camera
        const width = container.clientWidth;
        const height = container.clientHeight;
        this.camera = new THREE.PerspectiveCamera(75, width / height, 0.1, 1000);
        this.camera.position.z = 5;

        // Create renderer
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(width, height);
        container.innerHTML = '';
        container.appendChild(this.renderer.domElement);

        // Add helpers
        this.scene.add(new THREE.GridHelper(10, 10));
        this.scene.add(new THREE.AxesHelper(5));

        // Create cube
        const geometry = new THREE.BoxGeometry();
        const material = new THREE.MeshStandardMaterial({ 
            color: 0x0077ff,
            metalness: 0.3,
            roughness: 0.4,
            transparent: true,
            opacity: 0.3
        });
        this.cube = new THREE.Mesh(geometry, material);
        this.scene.add(this.cube);

        // Add lights
        this.scene.add(new THREE.AmbientLight(0xffffff, 0.5));
        
        const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
        directionalLight.position.set(5, 5, 5);
        this.scene.add(directionalLight);

        // Add arrow group
        this.scene.add(this.arrowGroup);

        // Add orbit controls
        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;

        // Initialize code editor
        this.initCodeEditor();

        // Initialize magnetization if data provided
        if (data) {
            this.updateMagnetization(data);
        }

        // Start animation
        this.animate();

        console.log('MagneticVisualization initialized successfully!');
    }

    /**
     * Initialize code editor
     */
    initCodeEditor() {
        const editorElement = document.getElementById('code-editor');
        if (editorElement) {
            this.codeEditor = ace.edit(editorElement);
            
            // Set editor options
            this.codeEditor.setOptions({
                mode: 'ace/mode/julia',
                theme: 'ace/theme/github',
                fontSize: '14px',
                tabSize: 4,
                useSoftTabs: true,
                enableBasicAutocompletion: true,
                enableSnippets: true,
                enableLiveAutocompletion: true,
                wrapBehavioursEnabled: true,
                wrap: true,
                showPrintMargin: false,
                showGutter: true,
                highlightActiveLine: true,
                autoScrollEditorIntoView: true
            });
            
            // Set default Julia code
            this.codeEditor.setValue(`# MicroMagnetic.jl Simulation
using MicroMagnetic

# Create mesh
mesh = BoxMesh([10, 10, 10], [2, 2, 2])

# Set material parameters
set_material_parameters(mesh, 
    Ms=8e5,
    A=1e-11,
    K=0.0,
    u=[0, 0, 1]
)

# Set initial magnetization
set_magnetization(mesh, [1, 0, 0])

# Create driver
driver = LLGDriver(alpha=0.1)

# Run simulation
run_simulation(mesh, driver, 1e-9, 100)

# Visualize results
visualize(mesh)
`, 1); // 1 means set cursor at the beginning
            
            // Add event listener for code changes
            this.codeEditor.on('change', () => {
                const code = this.codeEditor.getValue();
                console.log('Code changed:', code);
            });
            
            console.log('Code editor initialized successfully!');
        } else {
            console.warn('Code editor container not found');
        }
    }

    /**
     * Animation loop
     */
    animate() {
        this.animationId = requestAnimationFrame(() => this.animate());

        if (this.renderer && this.scene && this.camera) {
            // Update controls
            if (this.controls) {
                this.controls.update();
            }
            
            this.renderer.render(this.scene, this.camera);
        }
    }

    /**
     * Update magnetization data
     * @param {Object} data - Magnetization data
     * @param {Object} options - Visualization options
     */
    updateMagnetization(data, options = {}) {
        console.log('Updating magnetization data:', data);
        console.log('Options:', options);
    
        // Update grid and dimensions
        let shouldRecalculatePositions = false;
        if (data.cells) {
            if (!this.gridSize || JSON.stringify(this.gridSize) !== JSON.stringify(data.cells)) {
                this.gridSize = data.cells;
                shouldRecalculatePositions = true;
            }
        }
        if (data.dimensions) {
            if (!this.dimensions || JSON.stringify(this.dimensions) !== JSON.stringify(data.dimensions)) {
                this.dimensions = data.dimensions;
                shouldRecalculatePositions = true;
            }
        }
    
        // Calculate cell size
        const cellSize = [
            this.dimensions[0] / this.gridSize[0],
            this.dimensions[1] / this.gridSize[1],
            this.dimensions[2] / this.gridSize[2]
        ];
    
        // Calculate arrow scale
        const baseScale = Math.min(...cellSize) * 0.8;
        const arrowScale = baseScale * (options.arrowScaleFactor || 1.0);
    
        // Filter cells by selection
        const selection = options.selection || { type: 'full' };
    
        // Collect arrow data
        const arrowData = [];
        const nx = this.gridSize[0];
        const ny = this.gridSize[1];
        const nz = this.gridSize[2];
    
        // Recalculate positions if needed
        if (shouldRecalculatePositions || !this.arrowPositions) {
            console.log('Recalculating arrow positions');
            this.arrowPositions = [];
            
            for (let i = 0; i < nx; i++) {
                for (let j = 0; j < ny; j++) {
                    for (let k = 0; k < nz; k++) {
                        if (!this.isInSelection(i, j, k, selection)) continue;

                        const index = i * ny * nz + j * nz + k;
                        const magnetization = data.magnetization[index];

                        if (magnetization) {
                            // Calculate arrow position (cell center)
                            const x = (i - (nx - 1) / 2) * cellSize[0];
                            const y = (j - (ny - 1) / 2) * cellSize[1];
                            const z = (k - (nz - 1) / 2) * cellSize[2];
                            
                            this.arrowPositions.push([x, y, z]);
                            
                            arrowData.push({
                                position: [x, y, z],
                                direction: magnetization
                            });
                        }
                    }
                }
            }
        } else {
            // Use stored positions, update directions only
            console.log('Using stored arrow positions, updating directions');
            let positionIndex = 0;
            
            for (let i = 0; i < nx; i++) {
                for (let j = 0; j < ny; j++) {
                    for (let k = 0; k < nz; k++) {
                        if (!this.isInSelection(i, j, k, selection)) continue;

                        const index = i * ny * nz + j * nz + k;
                        const magnetization = data.magnetization[index];

                        if (magnetization && positionIndex < this.arrowPositions.length) {
                            arrowData.push({
                                position: this.arrowPositions[positionIndex],
                                direction: magnetization
                            });
                            positionIndex++;
                        }
                    }
                }
            }
        }
    
        console.log(`Processing ${arrowData.length} arrows`);
    
        // Handle arrow count changes
        const currentArrowCount = this.arrows.length > 0 ? this.arrows[0].count : 0;
        
        if (currentArrowCount !== arrowData.length) {
            console.log(`Arrow count changed: clearing ${currentArrowCount} arrows, creating ${arrowData.length} new ones`);
            this.clearArrows();
            if (arrowData.length > 0) {
                this.createArrowInstances(arrowData, arrowScale);
            }
        } else if (arrowData.length > 0) {
            console.log(`Updating ${arrowData.length} existing arrow instances`);
            this.updateArrowInstances(arrowData, arrowScale);
        }
    }

    /**
     * Create arrow instances using InstancedMesh
     * @param {Array} arrowData - Arrow data
     * @param {number} arrowScale - Scale factor
     */
    createArrowInstances(arrowData, arrowScale) {
        // Create geometries
        const coneGeometry = new THREE.ConeGeometry(0.05, 0.2, 32);
        coneGeometry.translate(0, -0.2, 0);
        
        const cylinderGeometry = new THREE.CylinderGeometry(0.01, 0.01, 0.2, 32);
        cylinderGeometry.translate(0, -0.2, 0);
    
        // Create material
        const material = new THREE.MeshStandardMaterial({ 
            color: 0x0077ff,
            metalness: 0.3,
            roughness: 0.4
        });
    
        // Create instanced meshes
        const coneMesh = new THREE.InstancedMesh(coneGeometry, material, arrowData.length);
        const cylinderMesh = new THREE.InstancedMesh(cylinderGeometry, material, arrowData.length);
    
        // Set up matrices and colors
        const coneMatrix = new THREE.Matrix4();
        const cylinderMatrix = new THREE.Matrix4();
        const color = new THREE.Color();
        const up = new THREE.Vector3(0, 1, 0);
    
        for (let i = 0; i < arrowData.length; i++) {
            const data = arrowData[i];
            const position = new THREE.Vector3(...data.position);
            const direction = new THREE.Vector3(...data.direction).normalize();

            // Calculate rotation
            const quaternion = new THREE.Quaternion();
            quaternion.setFromUnitVectors(up, direction);

            // Calculate positioning
            const arrowDirection = direction.clone();
            const totalArrowLength = 0.4 * arrowScale;
            const offset = arrowDirection.clone().multiplyScalar(totalArrowLength * 0.5);
            
            // Position cylinder
            const cylinderPosition = position.clone().sub(offset)
                .add(arrowDirection.clone().multiplyScalar(0.2 * arrowScale));
            cylinderMatrix.compose(cylinderPosition, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            cylinderMesh.setMatrixAt(i, cylinderMatrix);

            // Position cone
            const conePosition = position.clone().sub(offset)
                .add(arrowDirection.clone().multiplyScalar(0.4 * arrowScale));
            coneMatrix.compose(conePosition, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            coneMesh.setMatrixAt(i, coneMatrix);

            // Set color based on direction
            color.setRGB(Math.abs(direction.x), Math.abs(direction.y), Math.abs(direction.z));
            coneMesh.setColorAt(i, color);
            cylinderMesh.setColorAt(i, color);
        }
    
        // Add to scene
        this.arrowGroup.add(coneMesh);
        this.arrowGroup.add(cylinderMesh);
    
        // Store references
        this.arrows.push(coneMesh, cylinderMesh);
    }

    /**
     * Update existing arrow instances
     * @param {Array} arrowData - Arrow data
     * @param {number} arrowScale - Scale factor
     */
    updateArrowInstances(arrowData, arrowScale) {
        if (this.arrows.length < 2) {
            console.error('Not enough arrow meshes');
            return;
        }
    
        const coneMesh = this.arrows[0];
        const cylinderMesh = this.arrows[1];
        const coneMatrix = new THREE.Matrix4();
        const cylinderMatrix = new THREE.Matrix4();
        const color = new THREE.Color();
        const up = new THREE.Vector3(0, 1, 0);
    
        for (let i = 0; i < arrowData.length; i++) {
            const data = arrowData[i];
            const position = new THREE.Vector3(...data.position);
            const direction = new THREE.Vector3(...data.direction).normalize();

            // Calculate rotation
            const quaternion = new THREE.Quaternion();
            quaternion.setFromUnitVectors(up, direction);

            // Calculate positioning
            const arrowDirection = direction.clone();
            const totalArrowLength = 0.6 * arrowScale;
            const offset = arrowDirection.clone().multiplyScalar(totalArrowLength * 0.5);
            
            // Update cylinder
            const cylinderPosition = position.clone().sub(offset)
                .add(arrowDirection.clone().multiplyScalar(0.2 * arrowScale));
            cylinderMatrix.compose(cylinderPosition, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            cylinderMesh.setMatrixAt(i, cylinderMatrix);

            // Update cone
            const conePosition = position.clone().sub(offset)
                .add(arrowDirection.clone().multiplyScalar(0.4 * arrowScale));
            coneMatrix.compose(conePosition, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            coneMesh.setMatrixAt(i, coneMatrix);

            // Update color
            color.setRGB(Math.abs(direction.x), Math.abs(direction.y), Math.abs(direction.z));
            coneMesh.setColorAt(i, color);
            cylinderMesh.setColorAt(i, color);
        }
    
        // Mark for update
        coneMesh.instanceMatrix.needsUpdate = true;
        coneMesh.instanceColor.needsUpdate = true;
        cylinderMesh.instanceMatrix.needsUpdate = true;
        cylinderMesh.instanceColor.needsUpdate = true;
    }

    /**
     * Check if cell is in selection area
     * @param {number} i - X index
     * @param {number} j - Y index
     * @param {number} k - Z index
     * @param {Object} selection - Selection config
     * @returns {boolean} In selection
     */
    isInSelection(i, j, k, selection) {
        switch (selection.type) {
            case 'full':
                return true;
            case 'slice':
                if (selection.axis === 'x' && i === selection.position) return true;
                if (selection.axis === 'y' && j === selection.position) return true;
                if (selection.axis === 'z' && k === selection.position) return true;
                return false;
            case 'box':
                return i >= selection.min[0] && i < selection.max[0] &&
                       j >= selection.min[1] && j < selection.max[1] &&
                       k >= selection.min[2] && k < selection.max[2];
            default:
                return true;
        }
    }

    /**
     * Clear all arrows
     */
    clearArrows() {
        for (const arrow of this.arrows) {
            this.arrowGroup.remove(arrow);
            arrow.geometry.dispose();
            arrow.material.dispose();
        }
        this.arrows = [];
    }

    /**
     * Resize renderer
     * @param {number} width - New width
     * @param {number} height - New height
     */
    resize(width, height) {
        if (this.camera) {
            this.camera.aspect = width / height;
            this.camera.updateProjectionMatrix();
        }

        if (this.renderer) {
            this.renderer.setSize(width, height);
        }
    }

    /**
     * Get visualization status
     * @returns {Object} Status info
     */
    getStatus() {
        return {
            initialized: !!this.scene,
            dimensions: this.dimensions,
            cells: this.gridSize,
            arrowCount: this.arrows.length > 0 ? this.arrows[0].count : 0
        };
    }

    /**
     * Clean up resources
     */
    dispose() {
        if (this.animationId) {
            cancelAnimationFrame(this.animationId);
        }

        this.clearArrows();

        if (this.renderer) {
            this.renderer.dispose();
        }

        console.log('MagneticVisualization disposed');
    }
}

// Export to global scope
window.MagneticVisualization = MagneticVisualization;

// Initialize GUI when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    console.log('Initializing MicroMagneticGUI...');
    
    // Create and initialize visualization
    const container = document.getElementById('visualization-container');
    if (container) {
        const visualization = new MagneticVisualization();
        visualization.init(container);
        
        // Add sample data for testing
        const sampleData = {
            cells: [10, 10, 10],
            dimensions: [10, 10, 10],
            magnetization: []
        };
        
        // Generate sample magnetization data
        for (let i = 0; i < 10; i++) {
            for (let j = 0; j < 10; j++) {
                for (let k = 0; k < 10; k++) {
                    // Create some interesting pattern
                    const x = Math.sin(i * 0.5) * Math.cos(k * 0.3);
                    const y = Math.cos(j * 0.5) * Math.sin(i * 0.3);
                    const z = Math.sin(k * 0.5) * Math.cos(j * 0.3);
                    
                    // Normalize
                    const norm = Math.sqrt(x*x + y*y + z*z);
                    sampleData.magnetization.push([x/norm, y/norm, z/norm]);
                }
            }
        }
        
        // Update visualization with sample data
        visualization.updateMagnetization(sampleData);
        
        console.log('MicroMagneticGUI initialized successfully!');
    } else {
        console.error('Visualization container not found');
    }
});