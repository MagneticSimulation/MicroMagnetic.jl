import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

/**
 * Visualization class for visualizing magnetization distributions
 */
class Visualization {
    constructor(containerId) {
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
        this.meshGroup = new THREE.Group(); 
        this.gridHelper = null;
        this.scaleFactor = 0.02;

        this.container = document.getElementById(containerId);
        if (this.container) {
            this.init(this.container);
        } else {
            console.error(`Container with ID '${containerId}' not found`);
        }
    }

    /**
     * Initialize visualization scene
     * @param {Object} container - DOM element containing renderer
     * @param {Object} data - Initialization data
     */
    init(container, data = null) {
        // Create scene
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x000000);

        // Create camera
        const width = container.clientWidth;
        const height = container.clientHeight;
        this.camera = new THREE.PerspectiveCamera(75, width / height, 0.1, 1000);
        this.camera.position.set(0, 10, 10);
        this.camera.lookAt(0, 0, 0);

        this.scene.rotation.x = -Math.PI / 2; 

        // Create renderer
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(width, height);
        container.innerHTML = '';
        container.appendChild(this.renderer.domElement);
        
        this.gridHelper = new THREE.GridHelper(100, 100, 0x666666, 0x555555);
        this.gridHelper.rotation.x = Math.PI / 2;
        this.scene.add(this.gridHelper);
    
        this.scene.add(new THREE.AmbientLight(0xffffff, 0.8));
        const light = new THREE.HemisphereLight(0xffffff, 0x222222, 1.5);
        this.scene.add(light);
        // Add arrow group
        this.scene.add(this.arrowGroup);
        this.scene.add(this.meshGroup);

        // Add orbit controls
        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;

        // Initialize with data if provided
        if (data) {
            if (data.type === 'mesh') {
                this.displayFDMesh(data);
            } else if (data.type === 'magnetization') {
                this.updateMagnetization(data);
            }
        }

        // Start animation
        this.animate();

        // Handle window resize
        window.addEventListener('resize', () => this.resize());

        console.log('Visualization initialized successfully!');
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
     * Display FD Mesh
     * @param {Object} meshData - Mesh data containing nx, ny, nz, dx, dy, dz
     */
    displayFDMesh(meshData) {
        this.clearArrows();
        this.clearMesh();
        
        const { nx, ny, nz, dx, dy, dz} = meshData;
        const scale = this.scaleFactor;
        const width = nx * dx * scale;
        const height = ny * dy * scale;
        const depth = nz * dz * scale;
                
        const boxGeometry = new THREE.BoxGeometry(width, height, depth, 1, 1, 1);
        const boxMaterial = new THREE.MeshPhongMaterial({ 
            color: 0xffffff, 
            transparent: true,
            opacity: 0.5, 
        });
        const mesh = new THREE.Mesh(boxGeometry, boxMaterial);
        
        // Enable shadows for better visual effect
        mesh.castShadow = true;
        mesh.receiveShadow = true;
        
        const position = mesh.geometry.attributes.position;
        for (let i = 0; i < position.count; i++) {
            position.setXYZ(i,
                position.getX(i) * -1,
                position.getY(i) * -1,
                position.getZ(i) * -1,
            );
        }
        // Update the geometry after modifying vertices
        mesh.geometry.attributes.position.needsUpdate = true;
        //mesh.geometry.computeVertexNormals();

        // Add mesh to group
        this.meshGroup.add(mesh);
        
        this.gridHelper.position.set(0, 0, -depth/2);
        // Adjust camera to fit the mesh
        this.adjustCameraToFitMesh(mesh);
    }

    /**
     * Create custom surface visualization based on surfaceData
     * @param {number} nx - Number of cells in x direction
     * @param {number} ny - Number of cells in y direction
     * @param {number} nz - Number of cells in z direction
     * @param {number} dx - Cell size in x direction
     * @param {number} dy - Cell size in y direction
     * @param {number} dz - Cell size in z direction
     * @param {Array} surfaceData - Array of cell status information (length = nx*ny*nz). 
     *                              Each element is Ms value (0 if cell should be excluded)
     */
    createCustomSurface(nx, ny, nz, dx, dy, dz, surfaceData) {
        // Create a material for surface faces
        const surfaceMaterial = new THREE.MeshBasicMaterial({ 
            color: 0xe74c3c,
            transparent: true,
            opacity: 0.7, 
            side: THREE.DoubleSide
        });
        
        // Create a material for surface edges
        const edgeMaterial = new THREE.LineBasicMaterial({ 
            color: 0xc0392b, 
            transparent: false,
            linewidth: 1.5 
        });
        
        // Prepare vertices and faces arrays for custom surface geometry
        const vertices = [];
        const faces = [];
        let vertexIndex = 0;
        
        // Helper function to get cell index from coordinates
        const getCellIndex = (i, j, k) => i + j * nx + k * nx * ny;
        
        // Helper function to check if a cell has Ms=0
        const isCellEmpty = (i, j, k) => {
            // Check if cell is out of bounds
            if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) {
                return true; // Treat out of bounds as empty
            }
            // Check if cell has Ms=0
            const index = getCellIndex(i, j, k);
            return surfaceData[index] === 0;
        };
        
        // Helper function to add a face to the geometry
        const addFace = (v1, v2, v3, v4) => {
            // Add vertices to array
            vertices.push(v1.x, v1.y, v1.z);
            vertices.push(v2.x, v2.y, v2.z);
            vertices.push(v3.x, v3.y, v3.z);
            vertices.push(v4.x, v4.y, v4.z);
            
            // Add faces (two triangles per quad)
            faces.push(vertexIndex, vertexIndex + 1, vertexIndex + 2);
            faces.push(vertexIndex, vertexIndex + 2, vertexIndex + 3);
            
            // Increment vertex index
            vertexIndex += 4;
        };
        
        // Calculate cell half sizes for positioning
        const halfDx = dx / 2;
        const halfDy = dy / 2;
        const halfDz = dz / 2;
        
        // Traverse all cells
        for (let k = 0; k < nz; k++) {
            for (let j = 0; j < ny; j++) {
                for (let i = 0; i < nx; i++) {
                    const cellIndex = getCellIndex(i, j, k);
                    
                    // Skip cells with Ms=0
                    if (surfaceData[cellIndex] === 0) {
                        continue;
                    }
                    
                    // Calculate cell center position
                    const x = (i - nx / 2 + 0.5) * dx;
                    const y = (j - ny / 2 + 0.5) * dy;
                    const z = (k - nz / 2 + 0.5) * dz;
                    
                    // Define cell vertices
                    const v1 = new THREE.Vector3(x - halfDx, y - halfDy, z - halfDz);
                    const v2 = new THREE.Vector3(x + halfDx, y - halfDy, z - halfDz);
                    const v3 = new THREE.Vector3(x + halfDx, y + halfDy, z - halfDz);
                    const v4 = new THREE.Vector3(x - halfDx, y + halfDy, z - halfDz);
                    const v5 = new THREE.Vector3(x - halfDx, y - halfDy, z + halfDz);
                    const v6 = new THREE.Vector3(x + halfDx, y - halfDy, z + halfDz);
                    const v7 = new THREE.Vector3(x + halfDx, y + halfDy, z + halfDz);
                    const v8 = new THREE.Vector3(x - halfDx, y + halfDy, z + halfDz);
                    
                    // Check and add faces for each side of the cell
                    
                    // Front face (z+)
                    if (isCellEmpty(i, j, k + 1)) {
                        addFace(v5, v6, v7, v8);
                    }
                    
                    // Back face (z-)
                    if (isCellEmpty(i, j, k - 1)) {
                        addFace(v1, v2, v3, v4);
                    }
                    
                    // Right face (x+)
                    if (isCellEmpty(i + 1, j, k)) {
                        addFace(v2, v6, v7, v3);
                    }
                    
                    // Left face (x-)
                    if (isCellEmpty(i - 1, j, k)) {
                        addFace(v5, v1, v4, v8);
                    }
                    
                    // Top face (y+)
                    if (isCellEmpty(i, j + 1, k)) {
                        addFace(v4, v3, v7, v8);
                    }
                    
                    // Bottom face (y-)
                    if (isCellEmpty(i, j - 1, k)) {
                        addFace(v5, v6, v2, v1);
                    }
                }
            }
        }
        
        // Create geometry from collected vertices and faces
        const customGeometry = new THREE.BufferGeometry();
        
        // Only create surface if we have vertices
        if (vertices.length > 0) {
            // Create vertex buffer
            const vertexBuffer = new Float32Array(vertices);
            customGeometry.setAttribute('position', new THREE.BufferAttribute(vertexBuffer, 3));
            
            // Create face buffer
            const faceBuffer = new Uint32Array(faces);
            customGeometry.setIndex(new THREE.BufferAttribute(faceBuffer, 1));
            
            // Compute normals for proper lighting (if needed)
            customGeometry.computeVertexNormals();
            
            // Create and add the surface mesh
            const customMesh = new THREE.Mesh(customGeometry, surfaceMaterial);
            
            this.meshGroup.add(customMesh);
            
            // Create wireframe for custom surface
            const customWireframe = new THREE.EdgesGeometry(customGeometry);
            const customEdges = new THREE.LineSegments(customWireframe, edgeMaterial);
            
            this.meshGroup.add(customEdges);
        } else {
            console.log('No custom surface faces to render. All cells may be empty.');
        }
    }

    /**
     * Adjust camera to fit the mesh
     * @param {THREE.Mesh} mesh - Mesh object to fit
     */
    adjustCameraToFitMesh(mesh) {
        // Calculate bounding box of the mesh
        const boundingBox = new THREE.Box3().setFromObject(mesh);
        
        // Calculate the size of the bounding box
        const size = boundingBox.getSize(new THREE.Vector3());
        
        // Calculate the distance needed to fit the entire mesh
        const maxDim = Math.max(size.x, size.y, size.z)*1.5;
        const fov = this.camera.fov * (Math.PI / 180);
        const distance = maxDim / (2 * Math.tan(fov / 2));
        
        // Set camera position
        this.camera.position.set(-distance/10, distance/2, distance);
        
        // Look at the center of the mesh
        const center = boundingBox.getCenter(new THREE.Vector3());
        this.camera.lookAt(center);
        this.controls.target.copy(center);
        this.controls.update();
    }

    /**
     * Check if cell is in selection area
     * @param {number i} - X index
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
     * Clear all mesh elements
     */
    clearMesh() {
        while (this.meshGroup.children.length > 0) {
            const child = this.meshGroup.children[0];
            if (child.geometry) child.geometry.dispose();
            if (child.material) child.material.dispose();
            this.meshGroup.remove(child);
        }
    }

    /**
     * Resize renderer
     */
    resize() {
        if (!this.container || !this.camera || !this.renderer) return;
        
        const width = this.container.clientWidth;
        const height = this.container.clientHeight;
        
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
            arrowCount: this.arrows.length > 0 ? this.arrows[0].count : 0,
            meshElementCount: this.meshGroup.children.length
        };
    }

    /**
     * Dispose resources
     */
    dispose() {
        if (this.animationId) {
            cancelAnimationFrame(this.animationId);
        }

        this.clearArrows();
        this.clearMesh();

        if (this.renderer) {
            this.renderer.dispose();
        }

        console.log('Visualization disposed');
    }
}

export default Visualization;