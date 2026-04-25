import * as THREE from 'three';

/**
 * Handles Finite Element Mesh visualization
 */
class FEMeshVisualization {
    constructor(scene) {
        this.scene = scene;
        this.meshGroup = new THREE.Group();
        this.scene.add(this.meshGroup);
        
        this.scaleFactor = 0.02;
        this.currentMesh = null;
        this.nodes = null;
        this.edges = null;
        this.coordinates = [];
        this.cellVerts = [];
        this.regionIds = [];
        this.surfaceIds = [];
    }

    /**
     * Clear the current mesh from the scene
     */
    clearMesh() {
        if (this.nodes) {
            this.scene.remove(this.nodes);
            this.nodes.geometry.dispose();
            this.nodes.material.dispose();
            this.nodes = null;
        }
        
        if (this.edges) {
            this.scene.remove(this.edges);
            this.edges.geometry.dispose();
            this.edges.material.dispose();
            this.edges = null;
        }
        
        if (this.currentMesh) {
            this.scene.remove(this.currentMesh);
            this.currentMesh.geometry.dispose();
            this.currentMesh.material.dispose();
            this.currentMesh = null;
        }
        
        // Clear stored data
        this.coordinates = [];
        this.cellVerts = [];
        this.regionIds = [];
        this.surfaceIds = [];
    }

    /**
     * Display FE Mesh
     * @param {Object} meshData - Mesh data containing coordinates, cell_verts, region_ids
     */
    displayFEMesh(meshData) {
        this.clearMesh();
        
        // Extract mesh data
        const { coordinates, cell_verts, region_ids } = meshData;
        
        if (!coordinates || !cell_verts) {
            console.error('Missing required mesh data: coordinates and cell_verts are required');
            return;
        }
        
        this.coordinates = coordinates;
        this.cellVerts = cell_verts;
        this.regionIds = region_ids || [];
        
        // Scale the coordinates
        const scaledCoords = this.scaleCoordinates(coordinates);
        
        // Create geometry for tetrahedrons
        const geometry = new THREE.BufferGeometry();
        const positions = [];
        const colors = [];
        const indices = [];
        
        // Generate unique materials for different regions
        const regionColors = this.generateRegionColors();
        
        // Process each tetrahedron
        for (let cellIdx = 0; cellIdx < cell_verts.length; cellIdx++) {
            const verts = cell_verts[cellIdx];
            const regionId = region_ids[cellIdx] || 0;
            const color = regionColors[regionId % regionColors.length];
            
            // Add vertices for this tetrahedron
            const baseIndex = positions.length / 3;
            
            for (let i = 0; i < 4; i++) {
                const vertIdx = verts[i];
                const [x, y, z] = scaledCoords[vertIdx];
                
                // Add position
                positions.push(x, y, z);
                
                // Add color
                colors.push(color.r, color.g, color.b);
            }
            
            // Define the tetrahedron faces (4 triangular faces)
            const tetraIndices = [
                0, 1, 2,
                0, 2, 3,
                0, 3, 1,
                1, 3, 2
            ];
            
            // Add indices for this tetrahedron
            for (const idx of tetraIndices) {
                indices.push(baseIndex + idx);
            }
        }
        
        // Set geometry attributes
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
        geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
        geometry.setIndex(indices);
        geometry.computeVertexNormals();
        
        // Create material
        const material = new THREE.MeshPhongMaterial({
            vertexColors: true,
            transparent: true,
            opacity: 0.8,
            side: THREE.DoubleSide,
            flatShading: false
        });
        
        // Create mesh
        this.currentMesh = new THREE.Mesh(geometry, material);
        this.currentMesh.castShadow = true;
        this.currentMesh.receiveShadow = true;
        this.scene.add(this.currentMesh);
        
        // Create nodes
        this.createNodes(scaledCoords);
        
        // Create edges
        this.createEdges(scaledCoords, cell_verts);
    }
    
    /**
     * Scale coordinates by the scale factor
     */
    scaleCoordinates(coordinates) {
        const scaledCoords = [];
        for (const [x, y, z] of coordinates) {
            scaledCoords.push([
                x * this.scaleFactor,
                y * this.scaleFactor,
                z * this.scaleFactor
            ]);
        }
        return scaledCoords;
    }
    
    /**
     * Create node points visualization
     */
    createNodes(scaledCoords) {
        const geometry = new THREE.BufferGeometry();
        const positions = [];
        
        for (const [x, y, z] of scaledCoords) {
            positions.push(x, y, z);
        }
        
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
        
        const material = new THREE.PointsMaterial({
            color: 0x000000,
            size: 0.2 * this.scaleFactor,
            transparent: true,
            opacity: 0.8
        });
        
        this.nodes = new THREE.Points(geometry, material);
        this.scene.add(this.nodes);
    }
    
    /**
     * Create edges visualization
     */
    createEdges(scaledCoords, cellVerts) {
        // Create a set of unique edges
        const edgesSet = new Set();
        
        for (const verts of cellVerts) {
            // Define all 6 edges of a tetrahedron
            const edges = [
                [verts[0], verts[1]],
                [verts[0], verts[2]],
                [verts[0], verts[3]],
                [verts[1], verts[2]],
                [verts[1], verts[3]],
                [verts[2], verts[3]]
            ];
            
            // Add edges to the set (sorted to avoid duplicates)
            for (let [v1, v2] of edges) {
                if (v1 > v2) {
                    [v1, v2] = [v2, v1];
                }
                edgesSet.add(`${v1},${v2}`);
            }
        }
        
        // Create edge geometry
        const geometry = new THREE.BufferGeometry();
        const positions = [];
        
        for (const edgeStr of edgesSet) {
            const [v1, v2] = edgeStr.split(',').map(Number);
            const [x1, y1, z1] = scaledCoords[v1];
            const [x2, y2, z2] = scaledCoords[v2];
            
            positions.push(x1, y1, z1);
            positions.push(x2, y2, z2);
        }
        
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
        
        const material = new THREE.LineBasicMaterial({
            color: 0x333333,
            transparent: true,
            opacity: 0.5,
            linewidth: 1
        });
        
        this.edges = new THREE.LineSegments(geometry, material);
        this.scene.add(this.edges);
    }
    
    /**
     * Generate colors for different regions
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
     * Toggle nodes visibility
     */
    toggleNodes(visible) {
        if (this.nodes) {
            this.nodes.visible = visible;
        }
    }
    
    /**
     * Toggle edges visibility
     */
    toggleEdges(visible) {
        if (this.edges) {
            this.edges.visible = visible;
        }
    }
    
    /**
     * Toggle mesh faces visibility
     */
    toggleFaces(visible) {
        if (this.currentMesh) {
            this.currentMesh.visible = visible;
        }
    }
    
    /**
     * Set the scale factor for the mesh
     */
    setScaleFactor(scale) {
        this.scaleFactor = scale;
        if (this.coordinates.length > 0 && this.cellVerts.length > 0) {
            // Redraw the mesh with new scale
            this.displayFEMesh({
                coordinates: this.coordinates,
                cell_verts: this.cellVerts,
                region_ids: this.regionIds
            });
        }
    }
    
    /**
     * Set visibility of all mesh components
     * @param {boolean} visible - Whether the mesh should be visible
     */
    setVisible(visible) {
        this.toggleNodes(visible);
        this.toggleEdges(visible);
        this.toggleFaces(visible);
    }
    
    /**
     * Get grid information for other components
     */
    getGridInfo() {
        return {
            type: 'fem',
            coordinates: this.coordinates,
            cellVerts: this.cellVerts,
            regionIds: this.regionIds,
            scaleFactor: this.scaleFactor
        };
    }
}

export default FEMeshVisualization;
