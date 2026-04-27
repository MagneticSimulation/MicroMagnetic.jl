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
        this.needsResample = true;
        
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
        
        this.needsResample = true;
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

        if (this.needsResample || updatePosition || !this.arrowPositions) {
            this.calculateArrowPositions();
            this.needsResample = false;
            this.initialCylinderPositions = null;
            this.initialConePositions = null;
            this.initialArrowScale = null;
        }

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
        
        // Calculate Z sampling with boundary protection
        const clampedSampleNz = Math.max(1, Math.min(sampleNz, nz));
        const stepZ = nz / clampedSampleNz;
        const startZ = clampedSampleNz <= 1 ? (nz - 1) / 2 : stepZ / 2;

        if (sampling === 'cylindrical') {
            const { radius, ringNum, ringStep } = this.settings;
            for (let k = 0; k < clampedSampleNz; k++) {
                // Ensure gridZ is within [0, nz-1]
                let gridZ = startZ + k * stepZ;
                gridZ = Math.max(0, Math.min(nz - 1, gridZ));
                
                const z = (gridZ - (nz - 1) / 2) * cellSize[2];
                
                // Center point - ensure it's within bounds
                let centerX = (nx - 1) / 2;
                let centerY = (ny - 1) / 2;
                centerX = Math.max(0, Math.min(nx - 1, centerX));
                centerY = Math.max(0, Math.min(ny - 1, centerY));
                
                this.arrowPositions.push({
                    grid: [centerX, centerY, gridZ],
                    world: [0, 0, z]
                });
                
                const dR = radius / ringNum;
                for (let ring = 1; ring <= ringNum; ring++) {
                    const r = dR * ring;
                    const num = ringStep * ring;
                    const dTheta = (2 * Math.PI) / num;

                    for (let j = 0; j < num; j++) {
                        const theta = dTheta * j;
                        const wx = r * Math.cos(theta);
                        const wy = r * Math.sin(theta);
                        let gx = wx / cellSize[0] + (nx - 1) / 2;
                        let gy = wy / cellSize[1] + (ny - 1) / 2;
                        
                        // Ensure gx and gy are within bounds
                        gx = Math.max(0, Math.min(nx - 1, gx));
                        gy = Math.max(0, Math.min(ny - 1, gy));
                        
                        this.arrowPositions.push({
                            grid: [gx, gy, gridZ],
                            world: [wx, wy, z]
                        });
                    }
                }
            }
        } else {
            const { sampleNx, sampleNy } = this.settings;
            
            // Calculate X sampling with boundary protection
            const clampedSampleNx = Math.max(1, Math.min(sampleNx, nx));
            const stepX = nx / clampedSampleNx;
            const startX = clampedSampleNx <= 1 ? (nx - 1) / 2 : stepX / 2;
            
            // Calculate Y sampling with boundary protection
            const clampedSampleNy = Math.max(1, Math.min(sampleNy, ny));
            const stepY = ny / clampedSampleNy;
            const startY = clampedSampleNy <= 1 ? (ny - 1) / 2 : stepY / 2;

            for (let i = 0; i < clampedSampleNx; i++) {
                for (let j = 0; j < clampedSampleNy; j++) {
                    for (let k = 0; k < clampedSampleNz; k++) {
                        // Calculate grid coordinates
                        let gridX = startX + i * stepX;
                        let gridY = startY + j * stepY;
                        let gridZ = startZ + k * stepZ;
                        
                        // Ensure all grid coordinates are within valid bounds
                        gridX = Math.max(0, Math.min(nx - 1, gridX));
                        gridY = Math.max(0, Math.min(ny - 1, gridY));
                        gridZ = Math.max(0, Math.min(nz - 1, gridZ));
                    
                        this.arrowPositions.push({
                            grid: [gridX, gridY, gridZ],
                            world: [
                                (gridX - (nx - 1) / 2) * cellSize[0],
                                (gridY - (ny - 1) / 2) * cellSize[1],
                                (gridZ - (nz - 1) / 2) * cellSize[2]
                            ]
                        });
                    }
                }
            }
        }
    }

    trilinearInterpolation(x, y, z) {
        const data = this.data;
        const [nx, ny, nz] = this.gridSize;
    
        const ix = Math.max(0, Math.min(nx - 1, Math.floor(x)));
        const iy = Math.max(0, Math.min(ny - 1, Math.floor(y)));
        const iz = Math.max(0, Math.min(nz - 1, Math.floor(z)));
    
        const ix1 = Math.min(ix + 1, nx - 1);
        const iy1 = Math.min(iy + 1, ny - 1);
        const iz1 = Math.min(iz + 1, nz - 1);
        
        const tx = ix1 === ix ? 0 : (x - ix);
        const ty = iy1 === iy ? 0 : (y - iy);
        const tz = iz1 === iz ? 0 : (z - iz);
        
        const nxy = nx * ny;
        const stride = 3;
        const result = [0, 0, 0];
    
        for (let comp = 0; comp < 3; comp++) {
            const base = (iz * nxy + iy * nx + ix) * stride + comp;
        
            const v000 = data[base];
            const v001 = data[(iz1 * nxy + iy  * nx + ix ) * stride + comp];
            const v010 = data[(iz  * nxy + iy1 * nx + ix ) * stride + comp];
            const v011 = data[(iz1 * nxy + iy1 * nx + ix ) * stride + comp];
            const v100 = data[(iz  * nxy + iy  * nx + ix1) * stride + comp];
            const v101 = data[(iz1 * nxy + iy  * nx + ix1) * stride + comp];
            const v110 = data[(iz  * nxy + iy1 * nx + ix1) * stride + comp];
            const v111 = data[(iz1 * nxy + iy1 * nx + ix1) * stride + comp];
        
            const c00 = v000 + (v001 - v000) * tz;
            const c01 = v010 + (v011 - v010) * tz;
            const c10 = v100 + (v101 - v100) * tz;
            const c11 = v110 + (v111 - v110) * tz;
        
            const c0 = c00 + (c01 - c00) * ty;
            const c1 = c10 + (c11 - c10) * ty;
        
            result[comp] = c0 + (c1 - c0) * tx;
        }
    
        return result;
    }

    collectArrowData() {
        if (!this.data || !this.arrowPositions || this.arrowPositions.length === 0) {
            console.warn("No data or arrow positions available");
            return [];
        }
        
        const arrowData = [];
        
        for (let i = 0; i < this.arrowPositions.length; i++) {
            const { grid, world } = this.arrowPositions[i];
            const vector = this.trilinearInterpolation(grid[0], grid[1], grid[2]);
            
            const length = Math.sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
            if (length < 0.5) {
                continue;
            }
            
            arrowData.push({
                position: world,
                direction: vector
            });
        }
        return arrowData;
    }

    createArrowInstances(arrowData, arrowScale) {
        const cylinderHeight = 0.4;
        const cylinderGeometry = new THREE.CylinderGeometry(0.06, 0.06, cylinderHeight, 32);
        cylinderGeometry.translate(0, -cylinderHeight/2, 0); // Move the top to origin

        const coneHeight = 0.3;
        const coneGeometry = new THREE.ConeGeometry(0.15, coneHeight, 32);
        coneGeometry.translate(0, coneHeight/2, 0); // Move the base to origin
        const totalLength = cylinderHeight + coneHeight; // Total length is 0.8
        
        const material = new THREE.MeshStandardMaterial({ 
            metalness: 0.3,
            roughness: 0.4
        });
        
        const coneMesh = new THREE.InstancedMesh(coneGeometry, material, arrowData.length);
        const cylinderMesh = new THREE.InstancedMesh(cylinderGeometry, material, arrowData.length);
        
        const coneMatrix = new THREE.Matrix4();
        const cylinderMatrix = new THREE.Matrix4();
                
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

        const arrowLength = totalLength * arrowScale; 
        const halfLength = arrowLength / 2;
        
        for (let i = 0; i < arrowData.length; i++) {
            const { position, direction } = arrowData[i];
            
            // Use temporary objects instead of creating new ones
            tempPos.set(...position);
            tempDir.set(...direction).normalize();
            
            tempQuat.setFromUnitVectors(upVector, tempDir);
            

            const startOfArrow = tempPos.clone().addScaledVector(tempDir, -halfLength);
            
            // The cylinder's local origin is at its top, so position the cylinder top at startOfArrow + cylinderHeight * direction * arrowScale
            // This makes the cylinder extend from startOfArrow to startOfArrow + cylinderHeight * direction * arrowScale
            const cylinderPos = startOfArrow.clone().addScaledVector(tempDir, cylinderHeight * arrowScale);
            
            // The cone's local origin is at its base, so position the cone base at startOfArrow + cylinderHeight * direction * arrowScale
            // This makes the cone extend from the end of the cylinder to startOfArrow + totalLength * direction * arrowScale
            const conePos = startOfArrow.clone().addScaledVector(tempDir, cylinderHeight * arrowScale);
    
            this.initialCylinderPositions[i] = cylinderPos.clone();
            cylinderMatrix.compose(cylinderPos, tempQuat, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            cylinderMesh.setMatrixAt(i, cylinderMatrix);
    
            this.initialConePositions[i] = conePos.clone();
            coneMatrix.compose(conePos, tempQuat, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            coneMesh.setMatrixAt(i, coneMatrix);    
    
            const normalizedValue = (direction[componentIndex] + 1) / 2;
            const color = getColor(normalizedValue, this.colormap);
            tempColor.setRGB(color.r, color.g, color.b);
            coneMesh.setColorAt(i, tempColor);
            cylinderMesh.setColorAt(i, tempColor);
        }
        
        this.arrowGroup.add(cylinderMesh);
        this.arrowGroup.add(coneMesh);
        this.arrows.push(cylinderMesh, coneMesh);
    }

    updateArrowInstances(arrowData, arrowScale) {
        if (this.arrows.length < 2) return;
        
        const cylinderMesh = this.arrows[0];
        const coneMesh = this.arrows[1];
        const coneMatrix = new THREE.Matrix4();
        const cylinderMatrix = new THREE.Matrix4();
        
        const hasScaleChanged = this.initialArrowScale !== arrowScale;
        const componentIndex = this.getComponentIndex(this.component);
        
        const tempPos = new THREE.Vector3();
        const tempDir = new THREE.Vector3();
        const tempQuat = new THREE.Quaternion();
        const tempColor = new THREE.Color();
        const upVector = new THREE.Vector3(0, 1, 0);
        const tempScale = new THREE.Vector3();
        
        // Arrow dimensions (match createArrowInstances)
        const cylinderHeight = 0.5;
        const coneHeight = 0.3;
        const totalLength = cylinderHeight + coneHeight;
        
        for (let i = 0; i < arrowData.length; i++) {
            const { position, direction } = arrowData[i];
            
            // Use temporary objects instead of creating new ones
            tempDir.set(...direction).normalize();
            tempPos.set(...position);
            
            tempQuat.setFromUnitVectors(upVector, tempDir);
            
            if (hasScaleChanged) {
                // Scale changed - need to recalculate positions using the same logic as createArrowInstances
                const arrowLength = totalLength * arrowScale;
                const halfLength = arrowLength / 2;
                
                // Calculate the start position of the arrow (the entire arrow is centered at tempPos)
                const startOfArrow = tempPos.clone().addScaledVector(tempDir, -halfLength);
                
                // The cylinder's local origin is at its top, so position the cylinder top at startOfArrow + cylinderHeight * direction * arrowScale
                const cylinderPos = startOfArrow.clone().addScaledVector(tempDir, cylinderHeight * arrowScale);
                
                // The cone's local origin is at its base, so position the cone base at startOfArrow + cylinderHeight * direction * arrowScale
                const conePos = startOfArrow.clone().addScaledVector(tempDir, cylinderHeight * arrowScale);
                
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
            
            // Use the same safety check for color as in createArrowInstances
            let color;
            try {
                color = getColor(normalizedValue, this.colormap);
                // Ensure color has r, g, b properties
                if (!color || typeof color !== 'object' || 
                    color.r === undefined || color.g === undefined || color.b === undefined) {
                    // Fallback to white if color is invalid
                    color = { r: 1, g: 1, b: 1 };
                }
            } catch (error) {
                console.error('Error getting color:', error);
                color = { r: 1, g: 1, b: 1 };
            }
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