import * as THREE from 'three';
import { getColor, normalizeValue } from './colormaps.js';

/**
 * Handles FD mesh visualization with surface-based optimization for large systems.
 *
 * Supports three display modes:
 * - 'outline': Display mesh outline as wireframe (light yellow, default for performance)
 * - 'cells': Display individual cells as filled boxes (memory intensive)
 * - 'box': Display a transparent box surface (light gray)
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
        this.wireframeInstances = null;
        this.boxHelper = null;
        this.boxSurface = null;
        this.colormap = 'viridis';
        this.threshold = 1e-5;

        // Display type: 'outline' | 'cells' | 'box'
        this.displayType = 'outline';

        // Cache for current data
        this.currentMeshData = null;
        this.currentMsData = null;
    }

    /**
     * Set display type
     * @param {string} type - Display type: 'cells', 'outline', or 'box'
     */
    setDisplayType(type) {
        if (this.displayType === type) return;
        this.displayType = type;
        if (this.currentMeshData) {
            this.updateMesh(this.currentMeshData, this.currentMsData);
        }
    }

    /**
     * Update FD Mesh with surface-based optimization
     *
     * @param {Object} meshData - Mesh data containing nx, ny, nz, dx, dy, dz
     * @param {Array} msData - Ms data array (optional, used for coloring)
     * @param {Object} options - Additional options
     */
    updateMesh(meshData, msData = null, options = {}) {
        // Clear existing cells
        this.clearCells();

        // Extract mesh parameters
        const { nx, ny, nz, dx, dy, dz } = meshData;
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

        // Determine if Ms data is available
        const useMsData = msData !== null;

        // Helper function to get cell index
        const getCellIndex = (i, j, k) => i + j * nx + k * nx * ny;

        // Helper function to check bounds
        const inBounds = (i, j, k) =>
            i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz;

        // Helper function to check if cell is empty
        const isCellEmpty = (i, j, k) => {
            if (!inBounds(i, j, k)) return true;
            const idx = getCellIndex(i, j, k);
            if (useMsData) {
                return msData[idx] < this.threshold;
            }
            return false;
        };

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

        // Default color
        const defaultColor = { r: 0.7, g: 0.7, b: 0.7 };

        // Handle different display types
        if (this.displayType === 'outline') {
            this.displayOutline(width, height, depth);
            return;
        } else if (this.displayType === 'box') {
            this.displayTransparentBox(width, height, depth);
            return;
        }

        // 'cells' mode: Count surface cells (cells with at least one empty neighbor)
        let surfaceCellCount = 0;
        for (let k = 0; k < nz; k++) {
            for (let j = 0; j < ny; j++) {
                for (let i = 0; i < nx; i++) {
                    const idx = getCellIndex(i, j, k);
                    const ms = useMsData ? msData[idx] : 1.0;

                    // Skip empty cells (Ms < threshold)
                    if (useMsData && ms < this.threshold) {
                        continue;
                    }

                    // Check if cell is on surface (has at least one empty neighbor)
                    const isSurface = (
                        isCellEmpty(i + 1, j, k) || // right
                        isCellEmpty(i - 1, j, k) || // left
                        isCellEmpty(i, j + 1, k) || // top
                        isCellEmpty(i, j - 1, k) || // bottom
                        isCellEmpty(i, j, k + 1) || // front
                        isCellEmpty(i, j, k - 1)    // back
                    );

                    if (isSurface) {
                        surfaceCellCount++;
                    }
                }
            }
        }

        if (surfaceCellCount === 0) {
            console.log('No surface cells to display');
            return;
        }

        console.log(`Found ${surfaceCellCount} surface cells out of ${nx * ny * nz} total cells`);

        // Create geometry
        const geometry = new THREE.BoxGeometry(cellSize, cellSize, cellSize);

        // Create material for filled cells
        const material = new THREE.MeshLambertMaterial({
            vertexColors: true,
            transparent: true,
            opacity: 0.8,
            side: THREE.DoubleSide
        });

        // Create instanced mesh
        this.cellInstances = new THREE.InstancedMesh(geometry, material, surfaceCellCount);

        // Create color attribute buffer
        const colors = new Float32Array(surfaceCellCount * 3);

        // Set matrices and colors for each instance
        const matrix = new THREE.Matrix4();
        let instanceIndex = 0;

        for (let k = 0; k < nz; k++) {
            for (let j = 0; j < ny; j++) {
                for (let i = 0; i < nx; i++) {
                    const idx = getCellIndex(i, j, k);
                    const ms = useMsData ? msData[idx] : 1.0;

                    // Skip empty cells
                    if (useMsData && ms < this.threshold) {
                        continue;
                    }

                    // Check if cell is on surface
                    const isSurface = (
                        isCellEmpty(i + 1, j, k) ||
                        isCellEmpty(i - 1, j, k) ||
                        isCellEmpty(i, j + 1, k) ||
                        isCellEmpty(i, j - 1, k) ||
                        isCellEmpty(i, j, k + 1) ||
                        isCellEmpty(i, j, k - 1)
                    );

                    if (!isSurface) {
                        continue;
                    }

                    // Calculate cell center position
                    const x = (i - nx / 2 + 0.5) * scaledDx;
                    const y = (j - ny / 2 + 0.5) * scaledDy;
                    const z = (k - nz / 2 + 0.5) * scaledDz;

                    // Set position matrix
                    matrix.makeTranslation(x, y, z);
                    this.cellInstances.setMatrixAt(instanceIndex, matrix);

                    // Set color based on Ms value
                    let color;
                    if (useMsData) {
                        const normalizedValue = normalizeValue(ms, minMs, maxMs);
                        color = getColor(normalizedValue, this.colormap);
                    } else {
                        color = defaultColor;
                    }

                    colors[instanceIndex * 3] = color.r;
                    colors[instanceIndex * 3 + 1] = color.g;
                    colors[instanceIndex * 3 + 2] = color.b;

                    instanceIndex++;
                }
            }
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
     * Display mesh outline as wireframe (light yellow)
     */
    displayOutline(width, height, depth) {
        const geometry = new THREE.BoxGeometry(width, height, depth);
        const edges = new THREE.EdgesGeometry(geometry);
        const material = new THREE.LineBasicMaterial({ color: 0xffffaa, linewidth: 2 });
        this.wireframeInstances = new THREE.LineSegments(edges, material);
        this.cellGroup.add(this.wireframeInstances);
        console.log(`Displayed outline wireframe: ${width}x${height}x${depth}`);
    }

    /**
     * Display a transparent box surface
     */
    displayTransparentBox(width, height, depth) {
        const geometry = new THREE.BoxGeometry(width, height, depth);
        const material = new THREE.MeshBasicMaterial({
            color: 0xcccccc,
            transparent: true,
            opacity: 0.3,
            side: THREE.DoubleSide,
            depthWrite: false
        });
        this.boxSurface = new THREE.Mesh(geometry, material);
        this.cellGroup.add(this.boxSurface);
        console.log(`Displayed transparent box surface: ${width}x${height}x${depth}`);
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
        this.wireframeInstances = null;
        this.boxHelper = null;
        this.boxSurface = null;
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
    }
}

export default FDMeshVisualization;