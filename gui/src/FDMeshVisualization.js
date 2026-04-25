import * as THREE from 'three';
import { getColor, normalizeValue } from './colormaps.js';

/**
 * Handles FD mesh visualization with cell-based coloring
 */
class FDMeshVisualization {
    constructor(scene) {
        this.scene = scene;
        
        // Create a single group for all cells
        this.cellGroup = new THREE.Group();
        this.scene.add(this.cellGroup);
        
        this.scaleFactor = 0.02;
        this.gridSize = [10, 10, 10];
        this.dimensions = [10, 10, 10];
        
        // Cell properties
        this.cellInstances = null;
        this.colormap = 'viridis';
        this.threshold = 1e-5;
        
        // Cache for current data
        this.currentMeshData = null;
        this.currentMsData = null;
        this.currentRegions = null;
    }

    /**
     * Update FD Mesh with cell-based coloring
     * @param {Object} meshData - Mesh data containing nx, ny, nz, dx, dy, dz
     * @param {Array} msData - Ms data array (optional, used for coloring)
     * @param {Object} options - Additional options
     */
    updateMesh(meshData, msData = null, options = {}) {
        // Clear existing cells
        this.clearCells();
        
        // Extract mesh parameters
        const { nx, ny, nz, dx, dy, dz, regions } = meshData;
        const scale = this.scaleFactor;
        
        // Calculate dimensions with scale factor
        const scaledDx = dx * scale;
        const scaledDy = dy * scale;
        const scaledDz = dz * scale;
        
        const width = nx * scaledDx;
        const height = ny * scaledDy;
        const depth = nz * scaledDz;
        
        // Update grid info
        this.gridSize = [nx, ny, nz];
        this.dimensions = [width, height, depth];
        
        // Calculate cell size
        const cellSize = Math.min(scaledDx, scaledDy, scaledDz) * 0.9;
        
        // Store current data for later updates
        this.currentMeshData = meshData;
        this.currentMsData = msData;
        this.currentRegions = regions;
        
        // Determine which data to use for coloring
        const useMsData = msData !== null;
        const useRegions = !useMsData && regions !== undefined;
        
        // Prepare array to hold cell positions and colors
        const cellData = [];
        
        // Helper function to get cell index
        const getCellIndex = (i, j, k) => i + j * nx + k * nx * ny;
        
        // Normalize Ms data if available
        let minMs = 0;
        let maxMs = 1;
        if (useMsData) {
            if (msData.length !== nx * ny * nz) {
                console.warn('Ms data length does not match grid size');
                return;
            }
            minMs = Math.min(...msData);
            maxMs = Math.max(...msData);
        }
        
        // Generate colors for regions if available
        const regionColors = this.generateRegionColors();
        
        // Iterate through all cells
        for (let k = 0; k < nz; k++) {
            for (let j = 0; j < ny; j++) {
                for (let i = 0; i < nx; i++) {
                    const idx = getCellIndex(i, j, k);
                    
                    // Skip empty cells based on threshold
                    if (useMsData && msData[idx] < this.threshold) {
                        continue;
                    }
                    
                    // Calculate cell center position
                    const x = (i - nx / 2 + 0.5) * scaledDx;
                    const y = (j - ny / 2 + 0.5) * scaledDy;
                    const z = (k - nz / 2 + 0.5) * scaledDz;
                    
                    // Determine cell color
                    let color;
                    if (useMsData) {
                        // Use colormap for Ms data
                        const normalizedValue = normalizeValue(msData[idx], minMs, maxMs);
                        color = getColor(normalizedValue, this.colormap);
                    } else if (useRegions) {
                        // Use discrete colors for regions
                        const regionId = regions[idx] || 0;
                        color = regionColors[regionId % regionColors.length];
                    } else {
                        // Default color if no data
                        color = { r: 0.7, g: 0.7, b: 0.7 }; // Light gray
                    }
                    
                    // Add cell data
                    cellData.push({
                        position: { x, y, z },
                        color
                    });
                }
            }
        }
        
        if (cellData.length === 0) {
            console.log('No cells to display above threshold');
            return;
        }
        
        console.log(`Found ${cellData.length} cells to display out of ${nx * ny * nz} total cells`);
        
        // Create geometry and material for instanced mesh
        const geometry = new THREE.BoxGeometry(cellSize, cellSize, cellSize);
        
        // Use a material that supports vertex colors for instanced coloring
        const material = new THREE.MeshLambertMaterial({
            vertexColors: true,
            transparent: true,
            opacity: 0.8,
            side: THREE.DoubleSide
        });
        
        // Create instanced mesh
        this.cellInstances = new THREE.InstancedMesh(geometry, material, cellData.length);
        
        // Create color attribute buffer
        const colors = new Float32Array(cellData.length * 3);
        
        // Set matrices and colors for each instance
        const matrix = new THREE.Matrix4();
        for (let i = 0; i < cellData.length; i++) {
            const cell = cellData[i];
            
            // Set position matrix
            matrix.makeTranslation(cell.position.x, cell.position.y, cell.position.z);
            this.cellInstances.setMatrixAt(i, matrix);
            
            // Set color
            colors[i * 3] = cell.color.r;
            colors[i * 3 + 1] = cell.color.g;
            colors[i * 3 + 2] = cell.color.b;
        }
        
        // Apply colors to geometry
        this.cellInstances.geometry.setAttribute('color', new THREE.InstancedBufferAttribute(colors, 3));
        
        // Update instance matrices
        this.cellInstances.instanceMatrix.needsUpdate = true;
        
        // Add to scene
        this.cellGroup.add(this.cellInstances);
        
        return this.cellInstances;
    }

    /**
     * Generate colors for different regions
     * @returns {Array} Array of {r, g, b} color objects
     */
    generateRegionColors() {
        return [
            { r: 0.8, g: 0.2, b: 0.2 }, // Red
            { r: 0.2, g: 0.8, b: 0.2 }, // Green
            { r: 0.2, g: 0.2, b: 0.8 }, // Blue
            { r: 0.8, g: 0.8, b: 0.2 }, // Yellow
            { r: 0.8, g: 0.2, b: 0.8 }, // Magenta
            { r: 0.2, g: 0.8, b: 0.8 }, // Cyan
            { r: 0.6, g: 0.4, b: 0.2 }, // Orange
            { r: 0.4, g: 0.2, b: 0.6 }  // Purple
        ];
    }

    /**
     * Adjust camera to fit the mesh
     * @param {THREE.Camera} camera - Camera to adjust
     * @param {THREE.OrbitControls} controls - Orbit controls
     */
    adjustCameraToFitMesh(camera, controls) {
        if (!this.cellGroup || this.cellGroup.children.length === 0) return;
        
        // Calculate bounding box of all cells
        const boundingBox = new THREE.Box3();
        boundingBox.setFromObject(this.cellGroup);
        
        // Calculate size
        const size = boundingBox.getSize(new THREE.Vector3());
        
        // Calculate distance needed
        const maxDim = Math.max(size.x, size.y, size.z) * 1.5;
        const fov = camera.fov * (Math.PI / 180);
        const distance = maxDim / (2 * Math.tan(fov / 2));
        
        // Set camera position
        camera.position.set(-distance / 10, distance / 2, distance);
        
        // Look at center
        const center = boundingBox.getCenter(new THREE.Vector3());
        camera.lookAt(center);
        
        if (controls) {
            controls.target.copy(center);
            controls.update();
        }
    }

    /**
     * Clear all cells from the scene
     */
    clearCells() {
        while (this.cellGroup.children.length > 0) {
            const child = this.cellGroup.children[0];
            if (child.geometry) child.geometry.dispose();
            if (child.material) child.material.dispose();
            this.cellGroup.remove(child);
        }
        this.cellInstances = null;
    }

    /**
     * Set the colormap for Ms data visualization
     * @param {string} name - Name of the colormap to use
     */
    setColormap(name) {
        this.colormap = name;
        
        // Redraw with current data if available
        if (this.currentMeshData) {
            this.updateMesh(this.currentMeshData, this.currentMsData);
        }
    }

    /**
     * Set the threshold for displaying cells
     * @param {number} value - Threshold value
     */
    setThreshold(value) {
        this.threshold = value;
        
        // Redraw with current data if available
        if (this.currentMeshData) {
            this.updateMesh(this.currentMeshData, this.currentMsData);
        }
    }

    /**
     * Set visibility of the entire mesh
     * @param {boolean} visible - Whether the mesh should be visible
     */
    setVisible(visible) {
        this.cellGroup.visible = visible;
    }
    
    /**
     * Set the scale factor for the mesh
     * @param {number} scale - New scale factor
     */
    setScaleFactor(scale) {
        this.scaleFactor = scale;
        if (this.currentMeshData) {
            // Redraw the mesh with new scale
            this.updateMesh(this.currentMeshData, this.currentMsData);
        }
    }

    /**
     * Get grid information
     * @returns {Object} Grid information including size and dimensions
     */
    getGridInfo() {
        return {
            gridSize: this.gridSize,
            dimensions: this.dimensions
        };
    }

    /**
     * Dispose all resources
     */
    dispose() {
        this.clearCells();
        this.scene.remove(this.cellGroup);
        
        // Clear cached data
        this.currentMeshData = null;
        this.currentMsData = null;
        this.currentRegions = null;
    }
}

export default FDMeshVisualization;