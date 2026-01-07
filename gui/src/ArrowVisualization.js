import * as THREE from 'three';
import { getColor } from './colormaps.js';

class ArrowVisualization {
    constructor(scene) {
        this.scene = scene;
        this.arrowGroup = new THREE.Group();
        this.scene.add(this.arrowGroup);
        
        this.arrows = [];
        // Store all arrow positions calculated by calculateArrowPositions
        this.arrowPositions = null;
        this.data = null;
        this.gridSize = null;
        this.dimensions = null;
        this.cellSize = null;
        
        // Store initial positions for optimization in updateArrowInstances
        this.initialCylinderPositions = null;
        this.initialConePositions = null;
        this.initialArrowScale = null;
        // Color mapping properties
        this.component = 'mx'; // 'mx', 'my', 'mz'
        this.colormap = 'viridis'; // Colormap name
        
        this.settings = {
            arrowSize: 1.0,
            sampling: 'cartesian',
            sampleNx: 10,
            sampleNy: 10,
            sampleNz: 5,
            radius: 10,
            ringNum: 4,
            ringStep: 5,
        };
    }

    setGridInfo(gridSize, dimensions) {
        this.gridSize = gridSize;
        this.dimensions = dimensions;
        this.cellSize = [
            dimensions[0] / gridSize[0],
            dimensions[1] / gridSize[1],
            dimensions[2] / gridSize[2]
        ];
        this.settings.radius = Math.min(dimensions[1], dimensions[0]) / 2.0;
        this.settings.radius -= Math.max(this.cellSize[0], this.cellSize[1]) / 2;
    }

    updateMagnetization(data, updatePosition=false, options = {}) {
        if (!data || !this.gridSize || !this.dimensions) {
            return;
        }
                
        const [nx, ny, nz] = this.gridSize;
        
        if (data.length !== nx * ny * nz * 3) {
            return;
        }
        
        this.data = data;
        this.settings.sampling = options.sampling;
        if (options.sampling === 'cartesian') {
            this.settings.sampleNx = options.sampleNx;
            this.settings.sampleNy = options.sampleNy;
        } else {
            this.settings.ringNum = options.sampleNx;
             this.settings.ringStep = options.sampleNy;
        }   
        this.settings.sampleNz = options.sampleNz;
        this.settings.arrowSize = options.arrowSize;
        this.component = options.component;
        this.colormap = options.colormap;
        
        // Recalculate arrow positions if sampling settings changed
        if ( updatePosition || !this.arrowPositions) {
            this.calculateArrowPositions();
            this.initialCylinderPositions = null;
            this.initialConePositions = null;
            this.initialArrowScale = null;
        }
        const arrowData = this.collectArrowData();
        if (arrowData.length === 0) {
            return;
        }
        
        const baseScale = Math.min(...this.cellSize) * 3;
        const arrowScale = baseScale * this.settings.arrowSize;

        const currentArrowCount = this.arrows.length > 0 ? this.arrows[0].count : 0;
        
        if (updatePosition || currentArrowCount !== arrowData.length) {
            this.clearArrows();
            this.createArrowInstances(arrowData, arrowScale);
        } else {
            this.updateArrowInstances(arrowData, arrowScale);
        }
    }

    calculateArrowPositions() {
        const [nx, ny, nz] = this.gridSize;
        const [dimX, dimY, dimZ] = this.dimensions;
        const cellSize = [dimX / nx, dimY / ny, dimZ / nz];

        this.arrowPositions = [];
        const { sampling, sampleNz } = this.settings;
        const stepZ = (nz - 1) / Math.max(1, sampleNz - 1);

        if (sampling === 'cylindrical') {
            const { radius, ringNum, ringStep } = this.settings;
            for (let k = 0; k < sampleNz; k++) {
                       
                const gridK = (sampleNz === 1) ? (nz - 1) / 2 : k * stepZ;
                const z = (gridK - (nz - 1) / 2) * cellSize[2];
                
                this.arrowPositions.push([0, 0, z]);
                
                const dR = radius / ringNum;
                for (let ring = 1; ring <= ringNum; ring++) {
                    const r = dR * ring;
                    const num = ringStep * ring;
                    const dTheta = (2 * Math.PI) / num;

                    for (let j = 0; j < num; j++) {
                        const theta = dTheta * j;
                        this.arrowPositions.push([r * Math.cos(theta), r * Math.sin(theta), z]);
                    }
                }
            }
        } else {
            const { sampleNx, sampleNy } = this.settings;
            const stepX = (nx - 1) / Math.max(1, sampleNx - 1);
            const stepY = (ny - 1) / Math.max(1, sampleNy - 1);

            for (let i = 0; i < sampleNx; i++) {
                for (let j = 0; j < sampleNy; j++) {
                    for (let k = 0; k < sampleNz; k++) {
                        const gridI = (sampleNx === 1) ? (nx - 1) / 2 : i * stepX;
                        const gridJ = (sampleNy === 1) ? (ny - 1) / 2 : j * stepY;
                        const gridK = (sampleNz === 1) ? (nz - 1) / 2 : k * stepZ;
                    
                        const x = (gridI - (nx - 1) / 2) * cellSize[0];
                        const y = (gridJ - (ny - 1) / 2) * cellSize[1];
                        const z = (gridK - (nz - 1) / 2) * cellSize[2];
                    
                        this.arrowPositions.push([x, y, z]);
                    }
                }
            }
        }
    }

    trilinearInterpolation(x, y, z) {
        // Store frequently accessed values in local variables
        const gridSize = this.gridSize;
        const data = this.data;
        
        const nx = gridSize[0];
        const ny = gridSize[1];
        const nz = gridSize[2];
        
        const nx1 = nx - 1;
        const ny1 = ny - 1;
        const nz1 = nz - 1;
        
        // Inline clamp function to reduce function calls
        const xClamped = Math.max(0, Math.min(nx1 - 1e-10, x));
        const yClamped = Math.max(0, Math.min(ny1 - 1e-10, y));
        const zClamped = Math.max(0, Math.min(nz1 - 1e-10, z));
        
        const x0 = Math.floor(xClamped);
        const x1 = Math.min(x0 + 1, nx1);
        const y0 = Math.floor(yClamped);
        const y1 = Math.min(y0 + 1, ny1);
        const z0 = Math.floor(zClamped);
        const z1 = Math.min(z0 + 1, nz1);
        
        const dx = xClamped - x0;
        const dy = yClamped - y0;
        const dz = zClamped - z0;
        
        const nxy = nx * ny;
        
        // Calculate indices once
        const idx000 = (z0 * nxy + y0 * nx + x0) * 3;
        const idx001 = (z1 * nxy + y0 * nx + x0) * 3;
        const idx010 = (z0 * nxy + y1 * nx + x0) * 3;
        const idx011 = (z1 * nxy + y1 * nx + x0) * 3;
        const idx100 = (z0 * nxy + y0 * nx + x1) * 3;
        const idx101 = (z1 * nxy + y0 * nx + x1) * 3;
        const idx110 = (z0 * nxy + y1 * nx + x1) * 3;
        const idx111 = (z1 * nxy + y1 * nx + x1) * 3;
        
        // Pre-allocate result array
        const result = [0, 0, 0];
        
        // Unroll the loop for better performance
        for (let comp = 0; comp < 3; comp++) {
            // Get values for this component
            const v000 = data[idx000 + comp];
            const v001 = data[idx001 + comp];
            const v010 = data[idx010 + comp];
            const v011 = data[idx011 + comp];
            const v100 = data[idx100 + comp];
            const v101 = data[idx101 + comp];
            const v110 = data[idx110 + comp];
            const v111 = data[idx111 + comp];
            
            // Interpolate along z-axis
            const c00 = v000 + (v001 - v000) * dz;
            const c01 = v010 + (v011 - v010) * dz;
            const c10 = v100 + (v101 - v100) * dz;
            const c11 = v110 + (v111 - v110) * dz;
            
            // Interpolate along y-axis
            const c0 = c00 + (c01 - c00) * dy;
            const c1 = c10 + (c11 - c10) * dy;
            
            // Interpolate along x-axis
            result[comp] = c0 + (c1 - c0) * dx;
        }
        
        return result;
    }

    collectArrowData() {
        if (!this.data || !this.arrowPositions || this.arrowPositions.length === 0) {
            console.warn("No data or arrow positions available");
            return [];
        }
        
        // Pre-allocate array space to reduce memory allocation
        const arrowData = [];
        
        const [nx, ny, nz] = this.gridSize;    
        const centerIdxX = (nx - 1) / 2;
        const centerIdxY = (ny - 1) / 2;
        const centerIdxZ = (nz - 1) / 2;
        
        for (let i = 0; i < this.arrowPositions.length; i++) {
            const pos = this.arrowPositions[i];
            const gridX = pos[0] / this.cellSize[0] + centerIdxX;
            const gridY = pos[1] / this.cellSize[1] + centerIdxY;
            const gridZ = pos[2] / this.cellSize[2] + centerIdxZ;
            
            const vector = this.trilinearInterpolation(gridX, gridY, gridZ);
            
            // Calculate vector length
            const length = Math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
            if (length < 0.5) {
                continue;
            }
            
            arrowData.push({
                position: pos,
                direction: vector
            });
        }
        return arrowData;
    }

    createArrowInstances(arrowData, arrowScale) {
        const coneGeometry = new THREE.ConeGeometry(0.05, 0.2, 32);
        coneGeometry.translate(0, -0.2, 0);
        
        const cylinderGeometry = new THREE.CylinderGeometry(0.01, 0.01, 0.2, 32);
        cylinderGeometry.translate(0, -0.2, 0);
        
        const material = new THREE.MeshStandardMaterial({ 
            metalness: 0.3,
            roughness: 0.4
        });
        
        const coneMesh = new THREE.InstancedMesh(coneGeometry, material, arrowData.length);
        const cylinderMesh = new THREE.InstancedMesh(cylinderGeometry, material, arrowData.length);
        
        const coneMatrix = new THREE.Matrix4();
        const cylinderMatrix = new THREE.Matrix4();
        const arrowLength = 0.6 * arrowScale;
        const offset = arrowLength * 0.5;
        
        // Store initial positions for optimization
        this.initialCylinderPositions = new Array(arrowData.length);
        this.initialConePositions = new Array(arrowData.length);
        this.initialArrowScale = arrowScale;
        const componentIndex = this.getComponentIndex(this.component);
        
        const tempPos = new THREE.Vector3();
        const tempDir = new THREE.Vector3();
        const tempQuat = new THREE.Quaternion();
        const tempColor = new THREE.Color();
        const upVector = new THREE.Vector3(0, 1, 0);
        const offsetVector = new THREE.Vector3();
        
        for (let i = 0; i < arrowData.length; i++) {
            const { position, direction } = arrowData[i];
            
            // Use temporary objects instead of creating new ones
            tempPos.set(...position);
            tempDir.set(...direction).normalize();
            
            tempQuat.setFromUnitVectors(upVector, tempDir);
            
            offsetVector.copy(tempDir).multiplyScalar(offset);
            
            const cylinderPos = tempPos.clone().sub(offsetVector).add(tempDir.clone().multiplyScalar(0.2 * arrowScale));
            this.initialCylinderPositions[i] = cylinderPos.clone();
            cylinderMatrix.compose(cylinderPos, tempQuat, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            cylinderMesh.setMatrixAt(i, cylinderMatrix);
            
            const conePos = tempPos.clone().sub(offsetVector).add(tempDir.clone().multiplyScalar(0.4 * arrowScale));
            this.initialConePositions[i] = conePos.clone();
            coneMatrix.compose(conePos, tempQuat, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            coneMesh.setMatrixAt(i, coneMatrix);
            
            const normalizedValue = (direction[componentIndex] + 1) / 2;
            
            const color = getColor(normalizedValue, this.colormap);
            tempColor.setRGB(color.r, color.g, color.b);
            coneMesh.setColorAt(i, tempColor);
            cylinderMesh.setColorAt(i, tempColor);
        }
        
        this.arrowGroup.add(coneMesh);
        this.arrowGroup.add(cylinderMesh);
        this.arrows.push(coneMesh, cylinderMesh);
    }

    updateArrowInstances(arrowData, arrowScale) {
        if (this.arrows.length < 2) return;
        
        const coneMesh = this.arrows[0];
        const cylinderMesh = this.arrows[1];
        const coneMatrix = new THREE.Matrix4();
        const cylinderMatrix = new THREE.Matrix4();
        
        const hasScaleChanged = this.initialArrowScale !== arrowScale;
        const componentIndex = this.getComponentIndex(this.component);
        
        const tempPos = new THREE.Vector3();
        const tempDir = new THREE.Vector3();
        const tempQuat = new THREE.Quaternion();
        const tempColor = new THREE.Color();
        const upVector = new THREE.Vector3(0, 1, 0);
        const offsetVector = new THREE.Vector3();
        const tempScale = new THREE.Vector3();
        
        for (let i = 0; i < arrowData.length; i++) {
            const { direction } = arrowData[i];
            
            // Use temporary objects instead of creating new ones
            tempDir.set(...direction).normalize();
            
            tempQuat.setFromUnitVectors(upVector, tempDir);
            
            if (hasScaleChanged) {
                // Scale changed - need to recalculate positions
                const arrowLength = 0.6 * arrowScale;
                const offset = arrowLength * 0.5;
                
                const { position } = arrowData[i];
                tempPos.set(...position);
                
                offsetVector.copy(tempDir).multiplyScalar(offset);
                
                const cylinderPos = tempPos.clone().sub(offsetVector).add(tempDir.clone().multiplyScalar(0.2 * arrowScale));
                const conePos = tempPos.clone().sub(offsetVector).add(tempDir.clone().multiplyScalar(0.4 * arrowScale));
                
                tempScale.set(arrowScale, arrowScale, arrowScale);
                
                cylinderMatrix.compose(cylinderPos, tempQuat, tempScale);
                coneMatrix.compose(conePos, tempQuat, tempScale);
            } else {
                // No scale change - use stored positions (optimization)
                const cylinderPos = this.initialCylinderPositions[i];
                const conePos = this.initialConePositions[i];
                
                // Only need to update rotation (direction), position and scale remain the same
                cylinderMatrix.compose(cylinderPos, tempQuat, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
                coneMatrix.compose(conePos, tempQuat, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            }
            
            // Apply transformations
            cylinderMesh.setMatrixAt(i, cylinderMatrix);
            coneMesh.setMatrixAt(i, coneMatrix);
            
            const normalizedValue = (direction[componentIndex] + 1) / 2;
            const color = getColor(normalizedValue, this.colormap);
            tempColor.setRGB(color.r, color.g, color.b);
            coneMesh.setColorAt(i, tempColor);
            cylinderMesh.setColorAt(i, tempColor);
        }
        
        // Update scale if it changed
        if (hasScaleChanged) {
            this.initialArrowScale = arrowScale;
        }
        
        // Mark matrices and colors as needing update
        coneMesh.instanceMatrix.needsUpdate = true;
        coneMesh.instanceColor.needsUpdate = true;
        cylinderMesh.instanceMatrix.needsUpdate = true;
        cylinderMesh.instanceColor.needsUpdate = true;
    }

    clearArrows() {
        for (const arrow of this.arrows) {
            this.arrowGroup.remove(arrow);
            arrow.geometry.dispose();
            arrow.material.dispose();
        }
        this.arrows = [];
    }

    setVisible(visible) {
        this.arrowGroup.visible = visible;
    }
    
    getComponentIndex(component) {
        switch(component) {
            case 'mx': return 0;
            case 'my': return 1;
            case 'mz': return 2;
            default: return 0;
        }
    }

    dispose() {
        this.clearArrows();
        this.scene.remove(this.arrowGroup);
    }
}

export default ArrowVisualization;