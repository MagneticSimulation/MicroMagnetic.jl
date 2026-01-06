import * as THREE from 'three';

/**
 * Handles FD mesh visualization
 */
class FDMeshVisualization {
    constructor(scene) {
        this.scene = scene;
        this.meshGroup = new THREE.Group();
        this.scene.add(this.meshGroup);
        
        this.scaleFactor = 0.02;
        this.gridSize = [10, 10, 10];
        this.dimensions = [10, 10, 10];
        this.currentMesh = null;
        this.edgeLines = null;
    }

    /**
     * Display FD Mesh
     * @param {Object} meshData - Mesh data containing nx, ny, nz, dx, dy, dz
     */
    displayFDMesh(meshData) {
        this.clearMesh();
        
        const { nx, ny, nz, dx, dy, dz } = meshData;
        const scale = this.scaleFactor;
        const width = nx * dx * scale;
        const height = ny * dy * scale;
        const depth = nz * dz * scale;

        this.gridSize = [nx, ny, nz];
        this.dimensions = [width, height, depth];
        const boxGeometry = new THREE.BoxGeometry(width, height, depth, 1, 1, 1);
        const boxMaterial = new THREE.MeshPhongMaterial({ 
            color: 0xffffff, 
            transparent: true,
            opacity: 0.5, 
        });
        
        const mesh = new THREE.Mesh(boxGeometry, boxMaterial);
        
        // Enable shadows
        mesh.castShadow = true;
        mesh.receiveShadow = true;
        
        // Adjust vertex positions
        const position = mesh.geometry.attributes.position;
        for (let i = 0; i < position.count; i++) {
            position.setXYZ(i,
                position.getX(i) * -1,
                position.getY(i) * -1,
                position.getZ(i) * -1,
            );
        }
        
        // Update geometry
        mesh.geometry.attributes.position.needsUpdate = true;
        this.currentMesh = mesh;
        
        // Add mesh to group
        this.meshGroup.add(mesh);
        
        // Add yellow edges
        const edgesGeometry = new THREE.EdgesGeometry(boxGeometry);
        this.edgeLines = new THREE.LineSegments(
            edgesGeometry,
            new THREE.LineBasicMaterial({ color: 0xffff00 })
        );
        this.meshGroup.add(this.edgeLines);
        
        return mesh;
    }

    /**
     * Adjust camera to fit the mesh
     * @param {THREE.Camera} camera - Camera to adjust
     * @param {THREE.OrbitControls} controls - Orbit controls
     */
    adjustCameraToFitMesh(camera, controls) {
        if (!this.currentMesh) return;
        
        // Calculate bounding box
        const boundingBox = new THREE.Box3().setFromObject(this.currentMesh);
        
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
     * Clear mesh
     */
    clearMesh() {
        while (this.meshGroup.children.length > 0) {
            const child = this.meshGroup.children[0];
            if (child.geometry) child.geometry.dispose();
            if (child.material) child.material.dispose();
            this.meshGroup.remove(child);
        }
        this.currentMesh = null;
        this.edgeLines = null;
    }

    /**
     * Set visibility
     */
    setVisible(visible) {
        this.meshGroup.visible = visible;
    }

    /**
     * Get grid information
     */
    getGridInfo() {
        return {
            gridSize: this.gridSize,
            dimensions: this.dimensions
        };
    }

    /**
     * Dispose resources
     */
    dispose() {
        this.clearMesh();
        this.scene.remove(this.meshGroup);
    }
}

export default FDMeshVisualization;