import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { GUI } from 'lil-gui';
import FDMeshVisualization from './FDMeshVisualization.js';
import VolumeVisualization from './VolumeVisualization.js';
import ArrowVisualization from './ArrowVisualization.js';
import SurfaceVisualization from './SurfaceVisualization.js';
import { getAvailableColormaps } from './colormaps.js';

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

        this.spin = null;
        this.fdMesh = null;
        
        // Visualizers
        this.fdMeshVisualization = null;
        this.volumeVisualization = null;
        this.arrowVisualization = null;
        this.surfaceVisualization = null;
        
        // Grid helper
        this.gridHelper = null;
        this._showGrid = true;
        
        // Arrow configuration (sampling settings)
        this.arrowConfig = {
            sampling: 'cartesian', // 'cartesian' | 'cylindrical'
            sampleNx: 10,
            sampleNy: 10,
            sampleNz: 10,
            arrowSize: 1.0,
            visible: true,
            component: 'mx',      // 'mx' | 'my' | 'mz'
            colormap: 'viridis'   // colormap name
        };
        
        // Surface configuration
        this.surfaceConfig = {
            type: 'surface',      // 'surface' | 'isosurface'
            component: 'mx',      // 'mx' | 'my' | 'mz'
            isoValue: 0.5,
            position: 0,          // surface position
            direction: 'z',       // 'x' | 'y' | 'z'
            visible: true,
            colormap: 'viridis'   // colormap name
        };
        
        // State flag
        this.hasSpinData = false;
        
        // GUI references
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
        this.fdMeshVisualization._visible = true;
        
        this.volumeVisualization = new VolumeVisualization(this.scene);
        this.volumeVisualization._visible = true;
        
        this.arrowVisualization = new ArrowVisualization(this.scene);
        this.surfaceVisualization = new SurfaceVisualization(this.scene);

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
        this.fdMesh = meshData;
        this.clearAllVisualizations()

        this.fdMeshVisualization.displayFDMesh(meshData);
        
        // Adjust grid helper position
        const gridInfo = this.fdMeshVisualization.getGridInfo();
        this.gridHelper.position.set(0, 0, -gridInfo.dimensions[2] / 2 - 0.01);
        
        // Adjust camera
        this.fdMeshVisualization.adjustCameraToFitMesh(this.camera, this.controls);
    
        // Set grid info for arrow visualization
        this.arrowVisualization.setGridInfo(gridInfo.gridSize, gridInfo.dimensions);
        this.surfaceVisualization.setGridInfo(gridInfo.gridSize, gridInfo.dimensions);

        this.updateSamplingRange();
        this.updateSurfacePositionRange();
    }
    
    /**
     * Update surface position range based on mesh size
     */
    updateSurfacePositionRange() {
        if (!this.gui || !this.fdMesh) return;
        if (this.gui.posSurface) {
            const maxValue = this.surfaceConfig.direction === 'x' ? this.fdMesh.nx - 1 :
                             this.surfaceConfig.direction === 'y' ? this.fdMesh.ny - 1 :
                             this.fdMesh.nz - 1;
            this.gui.posSurface.max(maxValue);
        }
    }
    
    /**
     * Update sampling density range based on mesh size
     */
    updateSamplingRange() {
        if (!this.gui || !this.fdMesh) return;
        
        const isCartesian = this.arrowConfig.sampling === 'cartesian';
        const sampling = (n) => Math.min(n, 2*Math.round(Math.sqrt(n)));
        
        if (this.gui.sampleNx) {
            this.gui.sampleNx.name(isCartesian ? 'Nx' : 'Nr');
        }
        if (this.gui.sampleNy) {
            this.gui.sampleNy.name(isCartesian ? 'Ny' : 'Nphi');
        }

        if (isCartesian) {
            this.gui.sampleNx?.max(this.fdMesh.nx);
            this.gui.sampleNx?.setValue(sampling(this.fdMesh.nx));
            
            this.gui.sampleNy?.max(this.fdMesh.ny);
            this.gui.sampleNy?.setValue(sampling(this.fdMesh.ny));
            
            this.gui.sampleNz?.max(this.fdMesh.nz);
            this.gui.sampleNz?.setValue(sampling(this.fdMesh.nz));
        } else {
            this.gui.sampleNx?.max(30).setValue(8);
            this.gui.sampleNy?.max(20).min(4).setValue(8);
        }
    }

    /**
     * Toggle arrow visibility - GUI callback
     * @param {boolean} visible - Whether arrows are visible
     */
    updateArrowVisibility(visible) {
        this.arrowConfig.visible = visible;
        this.arrowVisualization.setVisible(visible);
        
        if (visible) {
            this.gui.arrowFolder?.show();
        } else {
            this.gui.arrowFolder?.hide();
        }
    }

    /**
     * Update arrow configuration - GUI callback
     * @param {Object} config - Arrow configuration {sampling, nx, ny, nz, size, component, colormap}
     */
    updateArrowConfig(config) {
        const updatePosition = config.sampling !== undefined || 
                             config.sampleNx !== undefined || 
                             config.sampleNy !== undefined || 
                             config.sampleNz !== undefined;
        // Check if any configuration changed
        if (config.sampling !== undefined) {
            this.updateSamplingRange();
        }

        // Re-render arrows if any configuration changed and we have spin data
        if (this.hasSpinData) {
            this.arrowVisualization.updateMagnetization(
                this.spin, 
                updatePosition,
                this.arrowConfig
            );
        }
    }

    /**
     * Update magnetization data - core interface
     * @param {Object} data - Magnetization data (spin)
     */
    updateMagnetization(data) {
        this.spin = data;
        this.hasSpinData = true;
        
        // Update arrow visualization using unified interface
        this.arrowVisualization.updateMagnetization(
            this.spin, 
            false,
            this.arrowConfig
        );
        
        this.surfaceVisualization.updateVisualization(
            this.spin,
            this.surfaceConfig
        );
        
        // Enable GUI controls
        this.enableVisualizationControls(true);
    }
    
    /**
     * Update surface configuration - GUI callback
     * @param {Object} config - Surface configuration
     */
    updateSurfaceConfig(config) {
        if (config.type !== undefined) this.surfaceConfig.type = config.type;
        if (config.component !== undefined) this.surfaceConfig.component = config.component;
        if (config.isoValue !== undefined) this.surfaceConfig.isoValue = config.isoValue;
        if (config.direction !== undefined) {
            this.surfaceConfig.direction = config.direction;
            
            // Reset position to 0 when direction changes
            this.surfaceConfig.position = 0;
            if (this.gui?.posSurface) {
                this.gui.posSurface.setValue(0);
            }
            
            // Update position range for new direction
            this.updateSurfacePositionRange();
            
        }
        if (config.position !== undefined) {
            this.surfaceConfig.position = config.position;
        }
        if (config.colormap !== undefined) {
            this.surfaceConfig.colormap = config.colormap;
        }
        
        if (this.hasSpinData) {
            this.surfaceVisualization.updateVisualization(
                this.spin,
                this.surfaceConfig
            );
        }
    }
    
    /**
     * Toggle surface visibility - GUI callback
     * @param {boolean} visible - Whether surface is visible
     */
    updateSurfaceVisibility(visible) {
        this.surfaceConfig.visible = visible;
        this.surfaceVisualization.setVisible(visible);
        
        if (visible) {
            this.gui.surfaceFolder?.show();
        } else {
            this.gui.surfaceFolder?.hide();
        }
    }

    /**
     * Enable/disable visualization controls based on spin data availability
     * @param {boolean} enabled - Whether controls should be enabled
     */
    enableVisualizationControls(enabled) {
        if (!this.gui) return;
        
        if (this.gui.showArrows) {
            enabled ? this.gui.showArrows.enable() : this.gui.showArrows.disable();
        }
        if (this.gui.showSurface) {
            enabled ? this.gui.showSurface.enable() : this.gui.showSurface.disable();
        }
        
        // Show/hide folder controls
        if (enabled && this.arrowConfig.visible) {
            this.gui.arrowFolder?.show();
        }
        if (enabled && this.surfaceConfig.visible) {
            this.gui.surfaceFolder?.show();
        }
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
     * Clear all visualizations
     */
    clearAllVisualizations() {
        if (this.arrowVisualization) {
            this.arrowVisualization.clearArrows();
        }
        
        if (this.volumeVisualization) {
            this.volumeVisualization.clearVolume();
        }
        
        if (this.fdMeshVisualization) {
            this.fdMeshVisualization.clearMesh();
        }

        if (this.surfaceVisualization) {
            this.surfaceVisualization.clearSurface();
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
            container: lilGuiContainer
        });

        // Display settings
        const displayFolder = this.gui.addFolder('Display Settings');
        displayFolder.add(this, '_showGrid', true).name('Show Grid').onChange((value) => {
            if (this.gridHelper) {
                this.gridHelper.visible = value;
            }
        });
        displayFolder.add(this.fdMeshVisualization, '_visible', true).name('Show Mesh').onChange((value) => {
            this.fdMeshVisualization.setVisible(value);
        });
        displayFolder.add(this.volumeVisualization, '_visible', true).name('Show Volume').onChange((value) => {
            this.volumeVisualization.setVisible(value);
        });
        displayFolder.open();

        // Magnetization display settings
        const magnetizationFolder = this.gui.addFolder('Magnetization Display');
        
        // Arrow visibility toggle
        this.gui.showArrows = magnetizationFolder.add(this.arrowConfig, 'visible').name('Show Arrows')
            .onChange((value) => this.updateArrowVisibility(value));
        this.gui.showArrows.disable(); // Disabled until spin data arrives
        
        // Arrow settings folder
        const arrowFolder = magnetizationFolder.addFolder('Arrow Settings').hide();
        this.gui.arrowFolder = arrowFolder;
        
        // Sampling method
        arrowFolder.add(this.arrowConfig, 'sampling', ['cartesian', 'cylindrical'])
            .name('Sampling Method')
            .onChange((value) => this.updateArrowConfig({ sampling: value }));
        
        // Sample density
        const densityFolder = arrowFolder.addFolder('Sample Density');
        this.gui.sampleNx = densityFolder.add(this.arrowConfig, 'sampleNx', 1, 50, 1)
            .name('Nx').onFinishChange(() => this.updateArrowConfig({ sampleNx: this.arrowConfig.sampleNx }));
        this.gui.sampleNy = densityFolder.add(this.arrowConfig, 'sampleNy', 1, 50, 1)
            .name('Ny').onFinishChange(() => this.updateArrowConfig({ sampleNy: this.arrowConfig.sampleNy }));
        this.gui.sampleNz = densityFolder.add(this.arrowConfig, 'sampleNz', 1, 50, 1)
            .name('Nz').onFinishChange(() => this.updateArrowConfig({ sampleNz: this.arrowConfig.sampleNz }));
        this.updateSamplingRange();
        
        // Arrow size
        arrowFolder.add(this.arrowConfig, 'arrowSize', 0.1, 5.0, 0.1)
            .name('Arrow Size')
            .onChange((value) => this.updateArrowConfig({ arrowSize: value }));

        // Arrow component
        arrowFolder.add(this.arrowConfig, 'component', ['mx', 'my', 'mz'])
            .name('Component')
            .onChange((value) => this.updateArrowConfig({ component: value }));

        // Arrow colormap
        arrowFolder.add(this.arrowConfig, 'colormap', getAvailableColormaps())
            .name('Colormap')
            .onChange((value) => this.updateArrowConfig({ colormap: value }));

        // Surface visibility toggle
        this.gui.showSurface = magnetizationFolder.add(this.surfaceConfig, 'visible').name('Show Surface/Isosurface')
            .onChange((value) => this.updateSurfaceVisibility(value));
        this.gui.showSurface.disable(); // Disabled until spin data arrives
        
        // Surface settings folder
        const surfaceFolder = magnetizationFolder.addFolder('Surface/Isosurface').hide();
        this.gui.surfaceFolder = surfaceFolder;
        
        // Surface type
        surfaceFolder.add(this.surfaceConfig, 'type', ['surface', 'isosurface'])
            .name('Type')
            .onChange((value) => {
                this.updateSurfaceConfig({ type: value });
                if (value === 'surface') {
                    this.gui.posSurface.show();
                    this.gui.dirSurface.show();
                    this.gui.isoSurfaceValue.hide();
                } else {
                    this.gui.posSurface.hide();
                    this.gui.dirSurface.hide();
                    this.gui.isoSurfaceValue.show();
                }
            });

        // Surface component
        surfaceFolder.add(this.surfaceConfig, 'component', ['mx', 'my', 'mz'])
            .name('Component')
            .onChange((value) => this.updateSurfaceConfig({ component: value }));

        this.gui.posSurface = surfaceFolder.add(this.surfaceConfig, 'position', 0, 10, 1)
            .name('Position Index')
            .onChange(() => this.updateSurfaceConfig({ position: this.surfaceConfig.position }));
        this.gui.dirSurface = surfaceFolder.add(this.surfaceConfig, 'direction', ['x', 'y', 'z'])
            .name('Direction')
            .onChange((value) => this.updateSurfaceConfig({ direction: value }));
        // Colormap selector
        surfaceFolder.add(this.surfaceConfig, 'colormap', getAvailableColormaps())
            .name('Colormap')
            .onChange((value) => this.updateSurfaceConfig({ colormap: value }));

        // Isosurface settings - directly in surfaceFolder
        this.gui.isoSurfaceValue = surfaceFolder.add(this.surfaceConfig, 'isoValue', -1, 1, 0.01)
            .name('Isosurface').onChange((value) => this.updateSurfaceConfig({ isoValue: value }));

        // Initialize control visibility based on initial type using lil-gui methods
        const initialType = this.surfaceConfig.type;
        if (initialType === 'surface') {
            this.gui.posSurface.show();
            this.gui.dirSurface.show();
            this.gui.isoSurfaceValue.hide();
        } else {
            this.gui.posSurface.hide();
            this.gui.dirSurface.hide();
            this.gui.isoSurfaceValue.show();
        }

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

        if (this.surfaceVisualization) {
            this.surfaceVisualization.dispose();
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