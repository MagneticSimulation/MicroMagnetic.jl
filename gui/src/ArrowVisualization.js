import * as THREE from 'three';

/**
 * Handles arrow visualization for magnetization
 */
class ArrowVisualization {
    constructor(scene) {
        this.scene = scene;
        this.arrowGroup = new THREE.Group();
        this.scene.add(this.arrowGroup);
        
        this.arrows = [];
        this.arrowPositions = null;
        this.currentSpinData = null;
        
        // Arrow visualization settings
        this.settings = {
            arrowSize: 1.0,
            arrowSampling: 'cartesian',
            sampleDensity: { nx: 10, ny: 10, nz: 10 }
        };
    }

    /**
     * Update magnetization data with arrows
     * @param {Array} data - Magnetization data
     * @param {Object} gridSize - Grid dimensions [nx, ny, nz]
     * @param {Object} dimensions - Physical dimensions [width, height, depth]
     * @param {Object} options - Visualization options
     */
    updateMagnetization(data, gridSize, dimensions, options = {}) {
        console.log('ArrowVisualizer: updateMagnetization called with data length:', data.length);

        const [nx, ny, nz] = gridSize;
        
        if (data.length !== nx * ny * nz * 3) {
            console.error('ArrowVisualizer: Data length mismatch');
            return;
        }
        
        // Store current spin data
        this.currentSpinData = data;
        this.gridSize = gridSize;
        this.dimensions = dimensions;
        
        // Calculate cell size
        const cellSize = [
            dimensions[0] / gridSize[0],
            dimensions[1] / gridSize[1],
            dimensions[2] / gridSize[2]
        ];
        
        // Calculate arrow scale
        const baseScale = Math.min(...cellSize) * 0.8;
        const arrowScale = baseScale * (options.arrowScaleFactor || 1.0) * this.settings.arrowSize;
        
        // Filter cells by selection
        const selection = options.selection || { type: 'full' };
        
        // Calculate arrow positions
        if (!this.arrowPositions || this.arrowPositions.length !== nx * ny * nz) {
            this.calculateArrowPositions(gridSize, dimensions, selection);
        }
        
        // Collect arrow data
        const arrowData = this.collectArrowData(data, gridSize, selection);
        
        console.log(`ArrowVisualizer: Processing ${arrowData.length} arrows`);
        
        // Create or update arrows
        const currentArrowCount = this.arrows.length > 0 ? this.arrows[0].count : 0;
        
        if (currentArrowCount !== arrowData.length) {
            console.log(`ArrowVisualizer: Creating ${arrowData.length} arrows`);
            this.clearArrows();
            if (arrowData.length > 0) {
                this.createArrowInstances(arrowData, arrowScale);
            }
        } else if (arrowData.length > 0) {
            console.log(`ArrowVisualizer: Updating ${arrowData.length} arrow directions`);
            this.updateArrowInstances(arrowData, arrowScale);
        }
    }

    /**
     * Calculate arrow positions based on grid and selection
     */
    calculateArrowPositions(gridSize, dimensions, selection) {
        const [nx, ny, nz] = gridSize;
        const cellSize = [
            dimensions[0] / gridSize[0],
            dimensions[1] / gridSize[1],
            dimensions[2] / gridSize[2]
        ];
        
        this.arrowPositions = [];
        
        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny; j++) {
                for (let k = 0; k < nz; k++) {
                    if (!this.isInSelection(i, j, k, selection)) continue;

                    const x = (i - (nx - 1) / 2) * cellSize[0];
                    const y = (j - (ny - 1) / 2) * cellSize[1];
                    const z = (k - (nz - 1) / 2) * cellSize[2];
                    
                    this.arrowPositions.push([x, y, z]);
                }
            }
        }
    }

    /**
     * Collect arrow data from magnetization data
     */
    collectArrowData(data, gridSize, selection) {
        const [nx, ny, nz] = gridSize;
        const arrowData = [];
        let positionIndex = 0;
        
        for (let i = 0; i < nx; i++) {
            for (let j = 0; j < ny; j++) {
                for (let k = 0; k < nz; k++) {
                    if (!this.isInSelection(i, j, k, selection)) continue;

                    const index = i + j * nx + k * nx * ny;
                    const magIndex = index * 3;

                    const mx = data[magIndex];
                    const my = data[magIndex + 1];
                    const mz = data[magIndex + 2];

                    if (positionIndex < this.arrowPositions.length) {
                        arrowData.push({
                            position: this.arrowPositions[positionIndex],
                            direction: [mx, my, mz]
                        });
                        positionIndex++;
                    }
                }
            }
        }
        
        return arrowData;
    }

    /**
     * Create arrow instances using InstancedMesh
     */
    createArrowInstances(arrowData, arrowScale) {
        // Create geometries
        const coneGeometry = new THREE.ConeGeometry(0.05, 0.2, 32);
        coneGeometry.translate(0, -0.2, 0);
        
        const cylinderGeometry = new THREE.CylinderGeometry(0.01, 0.01, 0.2, 32);
        cylinderGeometry.translate(0, -0.2, 0);
        
        // Create shared material
        const material = new THREE.MeshStandardMaterial({ 
            color: 0x0077ff,
            metalness: 0.3,
            roughness: 0.4
        });
        
        // Create instanced meshes
        const coneMesh = new THREE.InstancedMesh(coneGeometry, material, arrowData.length);
        const cylinderMesh = new THREE.InstancedMesh(cylinderGeometry, material, arrowData.length);
        
        // Set up transformation matrices
        const coneMatrix = new THREE.Matrix4();
        const cylinderMatrix = new THREE.Matrix4();
        const color = new THREE.Color();
        const up = new THREE.Vector3(0, 1, 0);
        const arrowLength = 0.6 * arrowScale;
        const offset = arrowLength * 0.5;
        
        for (let i = 0; i < arrowData.length; i++) {
            const data = arrowData[i];
            const position = new THREE.Vector3(...data.position);
            const direction = new THREE.Vector3(...data.direction).normalize();
            
            // Calculate rotation
            const quaternion = new THREE.Quaternion();
            quaternion.setFromUnitVectors(up, direction);
            
            // Calculate positions
            const arrowDirection = direction.clone();
            const offsetVector = arrowDirection.clone().multiplyScalar(offset);
            
            // Position cylinder
            const cylinderPosition = position.clone().sub(offsetVector)
                .add(arrowDirection.clone().multiplyScalar(0.2 * arrowScale));
            cylinderMatrix.compose(cylinderPosition, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            cylinderMesh.setMatrixAt(i, cylinderMatrix);
            
            // Position cone
            const conePosition = position.clone().sub(offsetVector)
                .add(arrowDirection.clone().multiplyScalar(0.4 * arrowScale));
            coneMatrix.compose(conePosition, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            coneMesh.setMatrixAt(i, coneMatrix);
            
            // Set color
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
     */
    updateArrowInstances(arrowData, arrowScale) {
        if (this.arrows.length < 2) {
            console.error('ArrowVisualizer: Not enough arrow meshes');
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
            
            // Calculate positions
            const arrowDirection = direction.clone();
            const totalArrowLength = 0.6 * arrowScale;
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
     * Set arrow size
     */
    setArrowSize(size) {
        this.settings.arrowSize = size;
        if (this.currentSpinData) {
            this.updateMagnetization(this.currentSpinData, this.gridSize, this.dimensions);
        }
    }

    /**
     * Set visibility
     */
    setVisible(visible) {
        this.arrowGroup.visible = visible;
    }

    /**
     * Dispose resources
     */
    dispose() {
        this.clearArrows();
        this.scene.remove(this.arrowGroup);
    }
}

export default ArrowVisualization;
