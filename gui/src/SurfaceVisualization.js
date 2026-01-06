import * as THREE from 'three';
import { MarchingCubes } from 'three/addons/objects/MarchingCubes.js';
import { getColorArrayRGB, getColormap } from './colormaps.js';

/**
 * Handles both surface and isosurface visualization for magnetization components
 */
class SurfaceVisualization {
    constructor(scene) {
        this.scene = scene;
        this.surfaceGroup = new THREE.Group();
        this.scene.add(this.surfaceGroup);

        this.gridSize = null;
        this.dimensions = null;
        this.spin = null;

        // Current visualization settings
        this.currentType = 'surface'; // 'surface' or 'isosurface'
        this.currentComponent = 'mx'; // 'mx', 'my', 'mz'
        this.isoValue = 0.5; // For isosurface
        this.position = 1; // For surface index
        this.surfaceAxis = 'z'; // 'x', 'y', 'z'
        this.currentColormap = 'viridis'; // Colormap name

        // Visual elements
        this.currentMesh = null;

        // Colormap cache
        this.colormapCache = new Map();
    }

    /**
     * Set grid information
     */
    setGridInfo(gridSize, dimensions) {
        this.gridSize = gridSize;
        this.dimensions = dimensions;
    }

    /**
     * Update visualization
     */
    updateVisualization(data, options = {}) {
        this.spin = data;
        
        const type = options.type || 'surface';

        this.currentType = type;
        this.currentComponent = options.component || 'mx';
        this.isoValue = options.isoValue || 0.0;
        this.position = options.position || 1;
        this.surfaceAxis = options.direction || 'z';

        this.clearSurface();

        if (!this.spin || !this.gridSize || !this.dimensions) {
            console.warn('SurfaceVisualization: Missing data or grid info');
            return;
        }

        if (type === 'surface') {
            this.displaySurface();
        } else if (type === 'isosurface') {
            this.displayIsosurface();
        }
    }

    /**
     * Display 2D surface slice
     */
    displaySurface() {
        const [nx, ny, nz] = this.gridSize;
        const [dimX, dimY, dimZ] = this.dimensions;

        const componentIndex = this.getComponentIndex(this.currentComponent);
        const sliceIndex = Math.min(Math.max(this.position, 0),
            this.surfaceAxis === 'x' ? nx - 1 :
                this.surfaceAxis === 'y' ? ny - 1 : nz - 1);

        // Create or reuse canvas for texture
        if (!this.canvas) {
            this.canvas = document.createElement('canvas');
            this.ctx = this.canvas.getContext('2d');
        }

        // Determine image dimensions based on axis
        let width, height, getDataIndex;

        switch (this.surfaceAxis) {
            case 'x':
                width = ny;
                height = nz;
                getDataIndex = (i, j) => (j * nx * ny + i * nx + sliceIndex) * 3 + componentIndex;
                break;
            case 'y':
                width = nx;
                height = nz;
                getDataIndex = (i, j) => (j * nx * ny + sliceIndex * nx + i) * 3 + componentIndex;
                break;
            case 'z':
            default:
                width = nx;
                height = ny;
                getDataIndex = (i, j) => (sliceIndex * nx * ny + j * nx + i) * 3 + componentIndex;
                break;
        }

        // Set canvas size
        this.canvas.width = width;
        this.canvas.height = height;

        // Get image data
        const imageData = this.ctx.createImageData(width, height);
        const data = imageData.data;

        // Extract slice data and find min/max for normalization
        // Data format: reshape(spin, 3, nx, ny, nz) -> [component, x, y, z]
        let minVal = Infinity;
        let maxVal = -Infinity;
        const sliceData = new Float32Array(width * height);

        // Extract slice based on axis
        for (let i = 0; i < width; i++) {
            for (let j = 0; j < height; j++) {
                const dataIndex = getDataIndex(i, j);
                const value = parseFloat(this.spin[dataIndex]);
                const idx = i + j * width;
                sliceData[idx] = value;

                if (!isNaN(value) && isFinite(value)) {
                    minVal = Math.min(minVal, value);
                    maxVal = Math.max(maxVal, value);
                }
            }
        }

        // Check if we have valid data
        if (!this.spin || this.spin.length === 0) {
            console.warn('SurfaceVisualization: No valid spin data');
            return;
        }

        // Handle case where all values are the same or no valid values found
        if (Math.abs(maxVal - minVal) < 1e-10 || maxVal === -Infinity) {
            minVal = -1;
            maxVal = 1;
        }

        // Create image with colormap
        // Get RGB array directly
        const colormap = getColormap(this.currentColormap);
        const colormapArray = getColorArrayRGB(colormap, 256);

        for (let i = 0; i < width; i++) {
            for (let j = 0; j < height; j++) {
                const idx = i + j * width;
                const value = sliceData[idx];

                // Normalize to [0, 1]
                const normalized = (value - minVal) / (maxVal - minVal);

                // Map to colormap index (0-255)
                const cmapIdx = Math.min(255, Math.max(0, Math.floor(normalized * 255)));

                // Get color from colormap (flat array: r,g,b,r,g,b,...)
                const colorIdx = cmapIdx * 3;
                const r = Math.floor(colormapArray[colorIdx] * 255);
                const g = Math.floor(colormapArray[colorIdx + 1] * 255);
                const b = Math.floor(colormapArray[colorIdx + 2] * 255);

                // Set pixel in image data (flip Y for correct orientation)
                const pixelIdx = (i + (height - 1 - j) * width) * 4;
                data[pixelIdx] = r;     // R
                data[pixelIdx + 1] = g; // G
                data[pixelIdx + 2] = b; // B
                data[pixelIdx + 3] = 255; // A
            }
        }

        // Put image data to canvas
        this.ctx.putImageData(imageData, 0, 0);

        // Create or update texture
        if (!this.texture) {
            this.texture = new THREE.CanvasTexture(this.canvas);
            this.texture.minFilter = THREE.LinearFilter;
            this.texture.magFilter = THREE.LinearFilter;
            this.texture.wrapS = THREE.ClampToEdgeWrapping;
            this.texture.wrapT = THREE.ClampToEdgeWrapping;
        } else {
            this.texture.needsUpdate = true;
        }

        // Create plane geometry with correct aspect ratio
        let planeWidth, planeHeight;
        switch (this.surfaceAxis) {
            case 'x':
                planeWidth = dimY;
                planeHeight = dimZ;
                break;
            case 'y':
                planeWidth = dimX;
                planeHeight = dimZ;
                break;
            case 'z':
            default:
                planeWidth = dimX;
                planeHeight = dimY;
                break;
        }

        const geometry = new THREE.PlaneGeometry(planeWidth, planeHeight);
        const material = new THREE.MeshBasicMaterial({
            map: this.texture,
            side: THREE.DoubleSide,
            transparent: true,
            opacity: 0.9
        });

        this.currentMesh = new THREE.Mesh(geometry, material);

        // Position the surface
        if (this.surfaceAxis === 'x') {
            const cellSize = dimX / nx;
            this.currentMesh.position.x = (this.position - (nx - 1) / 2) * cellSize;
            this.currentMesh.rotation.z = Math.PI / 2;
            this.currentMesh.rotation.y = Math.PI / 2;
        } else if (this.surfaceAxis === 'y') {
            const cellSize = dimY / ny;
            this.currentMesh.position.y = (this.position - (ny - 1) / 2) * cellSize;
            this.currentMesh.rotation.x = Math.PI / 2; 
        } else {
            const cellSize = dimZ / nz;
            this.currentMesh.position.z = (this.position - (nz - 1) / 2) * cellSize;
        }

        this.surfaceGroup.add(this.currentMesh);

    }

    displayIsosurface() {
        if (!this.spin || this.spin.length === 0) {
            console.warn('SurfaceVisualization: No valid spin data for isosurface');
            return;
        }

        const [nx, ny, nz] = this.gridSize;
        const [dimX, dimY, dimZ] = this.dimensions;
        const totalPoints = nx * ny * nz;

        // Extract data
        const componentIndex = this.getComponentIndex(this.currentComponent);
        const field = new Float32Array(totalPoints);

        // Copy data directly
        for (let idx = 0; idx < totalPoints; idx++) {
            field[idx] = this.spin[idx * 3 + componentIndex];
        }
        // Clear previous mesh
        this.clearSurface();

        // Create MarchingCubes
        const resolution = Math.max(nx, ny, nz);

        const material = new THREE.MeshPhongMaterial({
            color: 0xff0000,
            specular: 0x111111,
            shininess: 30,
            side: THREE.DoubleSide,
            transparent: true,
            opacity: 0.7
        });

        const marchingCubes = new MarchingCubes(
            resolution,
            material,
            false,   // enableUvs
            false,   // enableColors
            100000   // maxPolyCount
        );

        // Set isolation value
        marchingCubes.isolation = this.isoValue;
        // Reset and fill data
        marchingCubes.reset();

        // Fill MarchingCubes with setCell
        this.fillMarchingCubesWithSetCell(marchingCubes, field, nx, ny, nz, resolution);

        // Generate isosurface
        marchingCubes.update();

        // Check results
        if (marchingCubes.geometry) {
            const positions = marchingCubes.geometry.attributes.position;
            console.log(`Generated ${positions.count} vertices`);

            if (positions.count === 0) {
                console.warn('No vertices generated. Trying different isolation value.');
                this.tryDifferentIsolation(marchingCubes, field, nx, ny, nz, resolution);
                return;
            }
        }

        // Scale and position
        marchingCubes.scale.set(dimX / 2, dimY / 2, dimZ / 2);

        // Add to scene
        this.currentMesh = marchingCubes;
        this.surfaceGroup.add(this.currentMesh);

        console.log('Surface created successfully');
    }

    /**
     * Fill MarchingCubes using setCell method
     */
    fillMarchingCubesWithSetCell(marchingCubes, field, nx, ny, nz, resolution) {
        console.log('Filling MarchingCubes with setCell...');

        let cellsFilled = 0;

        if (resolution === nx && resolution === ny && resolution === nz) {
            // Grid sizes match, fill directly
            for (let x = 0; x < nx; x++) {
                for (let y = 0; y < ny; y++) {
                    for (let z = 0; z < nz; z++) {
                        const idx = x + nx * (y + ny * z);
                        const value = field[idx];
                        marchingCubes.setCell(x, y, z, value);
                        cellsFilled++;
                    }
                }
            }
        } else {
            // Grid sizes don't match, need interpolation
            for (let x = 0; x < resolution; x++) {
                for (let y = 0; y < resolution; y++) {
                    for (let z = 0; z < resolution; z++) {
                        // Calculate position in original data
                        const srcX = Math.floor(x * nx / resolution);
                        const srcY = Math.floor(y * ny / resolution);
                        const srcZ = Math.floor(z * nz / resolution);

                        const idx = srcX + nx * (srcY + ny * srcZ);
                        const value = field[idx];
                        marchingCubes.setCell(x, y, z, value);
                        cellsFilled++;
                    }
                }
            }
        }

        console.log(`Filled ${cellsFilled} cells`);
    }


    /**
     * Get component index from component name
     */
    getComponentIndex(component) {
        switch (component) {
            case 'mx': return 0;
            case 'my': return 1;
            case 'mz': return 2;
            default: return 0;
        }
    }


    /**
     * Clear current surface/isosurface
     */
    clearSurface() {
        if (this.currentMesh) {
            // Dispose geometry and material
            this.currentMesh.geometry.dispose();
            this.currentMesh.material.dispose();

            // Dispose texture to prevent memory leaks
            if (this.texture) {
                this.texture.dispose();
                this.texture = null;
            }

            this.surfaceGroup.remove(this.currentMesh);
            this.currentMesh = null;
        }
    }

    /**
     * Set visibility
     */
    setVisible(visible) {
        this.surfaceGroup.visible = visible;
    }

    /**
     * Dispose resources
     */
    dispose() {
        this.clearSurface();
        this.scene.remove(this.surfaceGroup);
    }
}

export default SurfaceVisualization;