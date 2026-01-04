import * as THREE from 'three';

class VolumeVisualization {
    constructor(scene) {
        this.scene = scene;
        this.volumeGroup = new THREE.Group();
        this.scene.add(this.volumeGroup);
        
        this.gridSize = [10, 10, 10];
        this.dimensions = [10, 10, 10];
    }

    /**
     * Display only surface cubes (cells with at least one empty neighbor)
     */
    displayVolume(volumeData, gridSize, dimensions) {
        this.clearVolume();
        
        this.gridSize = gridSize;
        this.dimensions = dimensions;
        
        const [nx, ny, nz] = gridSize;
        
        // Calculate cell sizes
        const dx = dimensions[0] / gridSize[0];
        const dy = dimensions[1] / gridSize[1];
        const dz = dimensions[2] / gridSize[2];
        
        // Helper function to check if cell is empty
        const getCellIndex = (i, j, k) => i + j * nx + k * nx * ny;
        const isCellEmpty = (i, j, k) => {
            // Check bounds
            if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) {
                return true;
            }
            const index = getCellIndex(i, j, k);
            return volumeData[index] < 1e-5;
        };
        
        // Check if cell is on surface (has at least one empty neighbor)
        const isSurfaceCell = (i, j, k) => {
            const index = getCellIndex(i, j, k);
            
            // Skip empty cells
            if (volumeData[index] < 1e-5) {
                return false;
            }
            
            // Check all 6 neighbors
            return (
                isCellEmpty(i + 1, j, k) || // right
                isCellEmpty(i - 1, j, k) || // left
                isCellEmpty(i, j + 1, k) || // top
                isCellEmpty(i, j - 1, k) || // bottom
                isCellEmpty(i, j, k + 1) || // front
                isCellEmpty(i, j, k - 1)    // back
            );
        };
        
        // Count surface cells
        let surfaceCellCount = 0;
        for (let k = 0; k < nz; k++) {
            for (let j = 0; j < ny; j++) {
                for (let i = 0; i < nx; i++) {
                    if (isSurfaceCell(i, j, k)) {
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
        
        // Use instanced mesh for surface cells only
        const cubeSize = Math.min(dx, dy, dz) * 0.9;
        const geometry = new THREE.BoxGeometry(cubeSize, cubeSize, cubeSize);
        
        // Create material with soft color
        const material = new THREE.MeshLambertMaterial({
            color: 0x88ccff, // Very light blue
            transparent: true,
            opacity: 0.7,
            side: THREE.DoubleSide
        });
        
        const instancedMesh = new THREE.InstancedMesh(geometry, material, surfaceCellCount);
        
        const matrix = new THREE.Matrix4();
        let instanceIndex = 0;
        
        for (let k = 0; k < nz; k++) {
            for (let j = 0; j < ny; j++) {
                for (let i = 0; i < nx; i++) {
                    if (!isSurfaceCell(i, j, k)) {
                        continue;
                    }
                    
                    // Calculate cell center position
                    const x = (i - nx / 2 + 0.5) * dx;
                    const y = (j - ny / 2 + 0.5) * dy;
                    const z = (k - nz / 2 + 0.5) * dz;
                    
                    // Set position
                    matrix.makeTranslation(x, y, z);
                    instancedMesh.setMatrixAt(instanceIndex, matrix);
                    
                    instanceIndex++;
                }
            }
        }
        
        instancedMesh.instanceMatrix.needsUpdate = true;
        this.volumeGroup.add(instancedMesh);
        
        console.log(`Displayed ${surfaceCellCount} surface cells`);
    }

    /**
     * Clear all displayed volume cells
     */
    clearVolume() {
        while (this.volumeGroup.children.length > 0) {
            const child = this.volumeGroup.children[0];
            if (child.geometry) child.geometry.dispose();
            if (child.material) child.material.dispose();
            this.volumeGroup.remove(child);
        }
    }

    setVisible(visible) {
        this.volumeGroup.visible = visible;
    }

    dispose() {
        this.clearVolume();
        this.scene.remove(this.volumeGroup);
    }
}

export default VolumeVisualization;