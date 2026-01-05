import * as THREE from 'three';

class ArrowVisualization {
    constructor(scene) {
        this.scene = scene;
        this.arrowGroup = new THREE.Group();
        this.scene.add(this.arrowGroup);
        
        this.arrows = [];
        this.arrowPositions = null;
        this.data = null;
        this.gridSize = null;
        this.dimensions = null;
        
        this.settings = {
            arrowSize: 1.0,
            samplingMethod: 'cartesian',
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
        this.settings.radius = Math.min(dimensions[1], dimensions[0]) / 2.0;
    }

    updateMagnetization(data, options = {}) {
        // Guard against null data (e.g., when GUI changes trigger before data is loaded)
        if (!data || data.length === 0) {
            return;
        }
        
        if (!this.gridSize || !this.dimensions) {
            console.error('ArrowVisualizer: gridSize and dimensions not set');
            return;
        }
        
        const [nx, ny, nz] = this.gridSize;
        
        if (data.length !== nx * ny * nz * 3) {
            console.error('ArrowVisualizer: Data length mismatch');
            return;
        }
        
        this.data = data;
        
        const cellSize = [
            this.dimensions[0] / this.gridSize[0],
            this.dimensions[1] / this.gridSize[1],
            this.dimensions[2] / this.gridSize[2]
        ];
        
        const baseScale = Math.min(...cellSize) * 0.8;
        const arrowScale = baseScale * (options.arrowScaleFactor || 1.0) * this.settings.arrowSize;
        
        // Ensure arrowPositions are calculated before collecting arrow data
        if (!this.arrowPositions || this.arrowPositions.length === 0) {
            this.calculateArrowPositions();
        }
        
        const arrowData = this.collectArrowData();
        
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

    calculateArrowPositions() {
        const [nx, ny, nz] = this.gridSize;
        const [dimX, dimY, dimZ] = this.dimensions;
        const cellSize = [dimX / nx, dimY / ny, dimZ / nz];

        this.arrowPositions = [];
        const { samplingMethod, sampleNz } = this.settings;
        const stepZ = (nz - 1) / Math.max(1, sampleNz - 1);

        if (samplingMethod === 'cylindrical') {
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
        const [nx, ny, nz] = this.gridSize;
        const data = this.data;
    
        const clamp = (value, max) => Math.max(0, Math.min(max - 1e-10, value));
        const xClamped = clamp(x, nx - 1);
        const yClamped = clamp(y, ny - 1);
        const zClamped = clamp(z, nz - 1);
    
 
        const x0 = Math.floor(xClamped);
        const x1 = Math.min(x0 + 1, nx - 1);
        const y0 = Math.floor(yClamped);
        const y1 = Math.min(y0 + 1, ny - 1);
        const z0 = Math.floor(zClamped);
        const z1 = Math.min(z0 + 1, nz - 1);
    
        const dx = xClamped - x0;
        const dy = yClamped - y0;
        const dz = zClamped - z0;
    
        const nxy = nx*ny;
        const idx000 = (z0 * nxy + y0 * nx + x0) * 3;
        const idx001 = (z1 * nxy + y0 * nx + x0) * 3;
        const idx010 = (z0 * nxy + y1 * nx + x0) * 3;
        const idx011 = (z1 * nxy + y1 * nx + x0) * 3;
        const idx100 = (z0 * nxy + y0 * nx + x1) * 3;
        const idx101 = (z1 * nxy + y0 * nx + x1) * 3;
        const idx110 = (z0 * nxy + y1 * nx + x1) * 3;
        const idx111 = (z1 * nxy + y1 * nx + x1) * 3;
    
        const result = [0, 0, 0];
        for (let comp = 0; comp < 3; comp++) {
            const v000 = data[idx000 + comp];
            const v001 = data[idx001 + comp];
            const v010 = data[idx010 + comp];
            const v011 = data[idx011 + comp];
            const v100 = data[idx100 + comp];
            const v101 = data[idx101 + comp];
            const v110 = data[idx110 + comp];
            const v111 = data[idx111 + comp];
            
            const c00 = v000 * (1 - dz) + v001 * dz;
            const c01 = v010 * (1 - dz) + v011 * dz;
            const c10 = v100 * (1 - dz) + v101 * dz;
            const c11 = v110 * (1 - dz) + v111 * dz;
    
            const c0 = c00 * (1 - dy) + c01 * dy;
            const c1 = c10 * (1 - dy) + c11 * dy;

            result[comp] = c0 * (1 - dx) + c1 * dx;
        }
    
        return result;
    }

    collectArrowData() {
        if (!this.data || !this.arrowPositions || this.arrowPositions.length === 0) {
            console.warn("No data or arrow positions available");
            return [];
        }
        const arrowData = [];
        
        const [nx, ny, nz] = this.gridSize;    
        const [dimX, dimY, dimZ] = this.dimensions;
        this.cellSize = [dimX / nx, dimY / ny, dimZ / nz];
        
        const centerIdxX = (nx - 1) / 2;
        const centerIdxY = (ny - 1) / 2;
        const centerIdxZ = (nz - 1) / 2;
        
        for (let pos of this.arrowPositions) {
            const gridX = pos[0] / this.cellSize[0] + centerIdxX;
            const gridY = pos[1] / this.cellSize[1] + centerIdxY;
            const gridZ = pos[2] / this.cellSize[2] + centerIdxZ;
            
            const vector = this.trilinearInterpolation(gridX, gridY, gridZ);
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
            color: 0x0077ff,
            metalness: 0.3,
            roughness: 0.4
        });
        
        const coneMesh = new THREE.InstancedMesh(coneGeometry, material, arrowData.length);
        const cylinderMesh = new THREE.InstancedMesh(cylinderGeometry, material, arrowData.length);
        
        const coneMatrix = new THREE.Matrix4();
        const cylinderMatrix = new THREE.Matrix4();
        const color = new THREE.Color();
        const up = new THREE.Vector3(0, 1, 0);
        const arrowLength = 0.6 * arrowScale;
        const offset = arrowLength * 0.5;
        
        for (let i = 0; i < arrowData.length; i++) {
            const { position, direction } = arrowData[i];
            const pos = new THREE.Vector3(...position);
            const dir = new THREE.Vector3(...direction).normalize();
            
            const quaternion = new THREE.Quaternion();
            quaternion.setFromUnitVectors(up, dir);
            
            const offsetVector = dir.clone().multiplyScalar(offset);
            
            const cylinderPos = pos.clone().sub(offsetVector).add(dir.clone().multiplyScalar(0.2 * arrowScale));
            cylinderMatrix.compose(cylinderPos, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            cylinderMesh.setMatrixAt(i, cylinderMatrix);
            
            const conePos = pos.clone().sub(offsetVector).add(dir.clone().multiplyScalar(0.4 * arrowScale));
            coneMatrix.compose(conePos, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            coneMesh.setMatrixAt(i, coneMatrix);
            
            color.setRGB(Math.abs(dir.x), Math.abs(dir.y), Math.abs(dir.z));
            coneMesh.setColorAt(i, color);
            cylinderMesh.setColorAt(i, color);
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
        const color = new THREE.Color();
        const up = new THREE.Vector3(0, 1, 0);
        
        for (let i = 0; i < arrowData.length; i++) {
            const { position, direction } = arrowData[i];
            const pos = new THREE.Vector3(...position);
            const dir = new THREE.Vector3(...direction).normalize();
            
            const quaternion = new THREE.Quaternion();
            quaternion.setFromUnitVectors(up, dir);
            
            const offset = dir.clone().multiplyScalar(0.3 * arrowScale);
            
            const cylinderPos = pos.clone().sub(offset).add(dir.clone().multiplyScalar(0.2 * arrowScale));
            cylinderMatrix.compose(cylinderPos, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            cylinderMesh.setMatrixAt(i, cylinderMatrix);
            
            const conePos = pos.clone().sub(offset).add(dir.clone().multiplyScalar(0.4 * arrowScale));
            coneMatrix.compose(conePos, quaternion, new THREE.Vector3(arrowScale, arrowScale, arrowScale));
            coneMesh.setMatrixAt(i, coneMatrix);
            
            color.setRGB(Math.abs(dir.x), Math.abs(dir.y), Math.abs(dir.z));
            coneMesh.setColorAt(i, color);
            cylinderMesh.setColorAt(i, color);
        }
        
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

    setArrowSize(size) {
        this.settings.arrowSize = size;
        if (this.data) {
            this.updateMagnetization(this.data, this.gridSize, this.dimensions);
        }
    }

    setSampling(method, nx, ny, nz) {
        this.settings.samplingMethod = method;
        this.settings.sampleNz = nz;
        if (method === 'cartesian') {
            this.settings.sampleNx = nx;
            this.settings.sampleNy = ny;
        }else{
            this.settings.ringNum = nx;
            this.settings.ringStep = ny;
        }
        this.calculateArrowPositions();
    }

    setVisible(visible) {
        this.arrowGroup.visible = visible;
    }

    dispose() {
        this.clearArrows();
        this.scene.remove(this.arrowGroup);
    }
}

export default ArrowVisualization;