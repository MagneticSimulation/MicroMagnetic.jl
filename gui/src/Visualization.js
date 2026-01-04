import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { GUI } from 'lil-gui';
import FDMeshVisualization from './FDMeshVisualization.js';
import VolumeVisualization from './VolumeVisualization.js';
import ArrowVisualization from './ArrowVisualization.js';

/**
 * Main Visualization class coordinating all visualizers
 */
class Visualization {
    constructor(containerId) {
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.controls = null;
        this.animationId = null;
        
        // Visualizers
        this.fdMeshVisualization = null;
        this.volumeVisualization = null;
        this.arrowVisualization = null;
        
        // Grid helper
        this.gridHelper = null;
        
        // GUI control properties
        this.controlsData = {
            showGrid: true,
            showMesh: true,
            magnetizationDisplay: 'arrows', // 'arrows' or 'volume'
            arrowSampling: 'cartesian',
            arrowSize: 1.0,
            sampleDensity: {
                nx: 10,
                ny: 10,
                nz: 10
            },
            surfacePosition: {
                x: 0,
                y: 0,
                z: 0
            },
            surfaceDirection: 'z'
        };
        
        this.gui = null;
        this.container = document.getElementById(containerId);
        
        if (this.container) {
            this.init(this.container);
        } else {
            console.error(`Container with ID '${containerId}' not found`);
        }
    }

    /**
     * Initialize visualization scene
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
        
        // Add grid helper
        this.gridHelper = new THREE.GridHelper(100, 100, 0x666666, 0x555555);
        this.gridHelper.rotation.x = Math.PI / 2;
        this.scene.add(this.gridHelper);
    
        // Add lighting
        this.scene.add(new THREE.AmbientLight(0xffffff, 0.8));
        const light = new THREE.HemisphereLight(0xffffff, 0x222222, 1.5);
        this.scene.add(light);

        // Initialize visualizers
        this.fdMeshVisualization = new FDMeshVisualization(this.scene);
        this.volumeVisualization = new VolumeVisualization(this.scene);
        this.arrowVisualization = new ArrowVisualization(this.scene);

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

        // Initialize GUI
        this.initGUI();

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
     * Display FD Mesh
     * @param {Object} meshData - Mesh data
     */
    displayFDMesh(meshData) {
        const mesh = this.fdMeshVisualization.displayFDMesh(meshData);
        
        // Adjust grid helper position
        const gridInfo = this.fdMeshVisualization.getGridInfo();
        this.gridHelper.position.set(0, 0, -gridInfo.dimensions[2] / 2);
        
        // Adjust camera
        this.fdMeshVisualization.adjustCameraToFitMesh(this.camera, this.controls);
        
        return mesh;
    }

    /**
     * Update magnetization data
     * @param {Object} data - Magnetization data
     * @param {Object} options - Visualization options
     */
    updateMagnetization(data, options = {}) {
        const gridInfo = this.fdMeshVisualization.getGridInfo();
        this.arrowVisualization.updateMagnetization(
            data, 
            gridInfo.gridSize, 
            gridInfo.dimensions, 
            options
        );
    }

    /**
     * Display custom surface
     * @param {Array} surfaceData - Surface data
     */
    displayCustomSurface(surfaceData) {
        const gridInfo = this.fdMeshVisualization.getGridInfo();
        this.volumeVisualization.displayVolume(   
            surfaceData,
            gridInfo.gridSize,
            gridInfo.dimensions
        );
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
     */
    getStatus() {
        return {
            initialized: !!this.scene,
            dimensions: this.fdMeshVisualization.getGridInfo().dimensions,
            cells: this.fdMeshVisualization.getGridInfo().gridSize,
            arrowCount: this.arrowVisualization ? this.arrowVisualization.arrows.length : 0,
            meshElementCount: this.fdMeshVisualization.meshGroup.children.length
        };
    }

    /**
     * Initialize GUI controls
     */
    initGUI() {
        const lilGuiContainer = document.getElementById('lilgui-panel');
        
        this.gui = new GUI({
            title: 'Visualization Controls',
            width: 300,
            container: lilGuiContainer
        });

        // Display settings
        const displayFolder = this.gui.addFolder('Display Settings');
        displayFolder.add(this.controlsData, 'showGrid').name('Show Grid').onChange((value) => {
            if (this.gridHelper) {
                this.gridHelper.visible = value;
            }
        });
        displayFolder.add(this.controlsData, 'showMesh').name('Show Mesh').onChange((value) => {
            this.fdMeshVisualization.setVisible(value);
        });
        displayFolder.open();

        // Magnetization display settings
        const magnetizationFolder = this.gui.addFolder('Magnetization Display');
        magnetizationFolder.add(this.controlsData, 'magnetizationDisplay', ['arrows', 'surface'])
            .name('Display Type')
            .onChange((value) => {
                this.arrowVisualization.setVisible(value === 'arrows');
                this.volumeVisualization.setVisible(value === 'volume');
            });

        // Arrow settings
        const arrowFolder = magnetizationFolder.addFolder('Arrow Settings');
        arrowFolder.add(this.controlsData, 'arrowSampling', ['cartesian', 'polar'])
            .name('Sampling Method');
        
        const densityFolder = arrowFolder.addFolder('Sample Density');
        densityFolder.add(this.controlsData.sampleDensity, 'nx', 1, 50, 1).name('NX');
        densityFolder.add(this.controlsData.sampleDensity, 'ny', 1, 50, 1).name('NY');
        densityFolder.add(this.controlsData.sampleDensity, 'nz', 1, 50, 1).name('NZ');
        densityFolder.open();
        
        arrowFolder.add(this.controlsData, 'arrowSize', 0.1, 5.0, 0.1)
            .name('Arrow Size')
            .onChange((value) => {
                this.arrowVisualization.setArrowSize(value);
            });
        arrowFolder.open();

        // Surface settings
        const surfaceFolder = magnetizationFolder.addFolder('Surface Settings');
        const positionFolder = surfaceFolder.addFolder('Position');
        positionFolder.add(this.controlsData.surfacePosition, 'x', -10, 10, 0.1).name('X');
        positionFolder.add(this.controlsData.surfacePosition, 'y', -10, 10, 0.1).name('Y');
        positionFolder.add(this.controlsData.surfacePosition, 'z', -10, 10, 0.1).name('Z');
        positionFolder.open();
        
        surfaceFolder.add(this.controlsData, 'surfaceDirection', ['x', 'y', 'z']).name('Direction');
        surfaceFolder.open();

        magnetizationFolder.open();
    }

    /**
     * Dispose resources
     */
    dispose() {
        cancelAnimationFrame(this.animationId);
        
        if (this.fdMeshVisualization) {
            this.fdMeshVisualization.dispose();
        }
        
        if (this.volumeVisualization) {
            this.volumeVisualization.dispose();
        }
        
        if (this.arrowVisualization) {
            this.arrowVisualization.dispose();
        }

        if (this.renderer) {
            this.renderer.dispose();
        }

        if (this.gui) {
            this.gui.destroy();
        }

        console.log('Visualization disposed');
    }
}

export default Visualization;