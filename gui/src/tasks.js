
const Ms_demo = [
    { title: "Create Mesh", code: `mesh = FDMesh(; nx=100, ny=100, nz=10, dx=2e-9, dy=2e-9, dz=2e-9);`, cellType: "mesh" },
    { title: "Create Simulation", code: `sim = Sim(mesh; name="Ms_demo")`, cellType: "sim" },
    { title: "Set Ms (const)", code: `set_Ms(sim, 8e5)`, cellType: "ms" },
    { title: "Set Ms (function)", cellType: "ms",  code:
`function circular_Ms(x, y, z)
    if x^2 + y^2 <= (50e-9)^2
        return 8e5
    end
    return 0.0
end
set_Ms(sim, circular_Ms)`, },
    { title: "Set Ms (shape)", cellType: "ms", code:
`shape = Cylinder(center=(-100e-9, 0, 0), radius = 100e-9);
# only set the Ms inside the given shape\nset_Ms(sim, shape, 8e5)`, },
    { title: "Set Ms (shapes)", cellType: "ms",  code:
`shape = Cylinder(center=(-100e-9, 0, 0), radius = 50e-9);
set_region(mesh, shape, 1);

shape = Box(center=(100e-9, 0, 0), size=(100e-9, 100e-9, 50e-9));
set_region(mesh, shape, 2);

set_Ms(sim, region_map(1 => 8e5, 2 => 3e5))`, },
    { title: "Set Ms (voronoi)", cellType: "ms",  code:
`grain_ids, gb_mask, seeds = voronoi(mesh, min_dist=20, seed=100, threshold=0.28)
function spatial_Ms(i,j,k,dx,dy,dz)
    if gb_mask[i,j,k]
        return 0
    end
    return 8e5
end
set_Ms(sim, spatial_Ms)`},
]

const vortex = [
    { title: "Create Mesh", code: `mesh = FDMesh(nx=100, ny=100, nz=10, dx=2e-9, dy=2e-9, dz=2e-9);`, cellType: "mesh" },
    { title: "Set Region", code: `shape = Cylinder(radius = 100e-9);\nset_region(mesh, shape, 1);`, cellType: "set_region" },
    { title: "Create Simulation", code: `sim = Sim(mesh; name="vortex")`, cellType: "sim" },
    { title: "Set Saturation Magnetization", code: `set_Ms(sim, region_map(1 => 8e5))`, cellType: "ms" },
    { title: "Initialize Magnetization", code: `function init_fun(x, y, z)
    r = sqrt(x^2 + y^2)
    if r < 20e-9
        return (0, 0, 1)
    end
    return (y / r, -x / r, 0)
end
    init_m0(sim, init_fun)`, cellType: "m0" },
    { title: "Set SD driver", code: `set_driver(sim, driver="SD")`, cellType: "sd" },
    { title: "Add Exchange Interaction", code: `add_exch(sim, 1.3e-11)`, cellType: "exch", isInteraction: true },
    { title: "Add Demagnetization", code: `add_demag(sim)`, cellType: "demag", isInteraction: true },
    { title: "Relax System", code: `relax(sim; stopping_dmdt=0.01)`, cellType: "relax" },
];

const std4 = [
    { title: "Create Mesh", code: `mesh = FDMesh(; nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9);`, cellType: "mesh" },
    { title: "Create Simulation", code: `sim = Sim(mesh; name="std4")`, cellType: "sim" },
    { title: "Set Saturation Magnetization", code: `Ms = 8.6e5; # A/m\nset_Ms(sim, Ms)`, cellType: "ms" },
    { title: "Initialize Magnetization", code: `init_m0(sim, (1, 0.25, 0.1))`, cellType: "m0" },
    { title: "Set SD driver", code: `set_driver(sim, driver="SD")`, cellType: "sd" },
    { title: "Add Exchange Interaction", code: `A = 1.3e-11; # J/m\nadd_exch(sim, A)`, cellType: "exch", isInteraction: true },
    { title: "Add Demagnetization", code: `add_demag(sim)`, cellType: "demag", isInteraction: true },
    { title: "Relax System", code: `relax(sim; stopping_dmdt=0.01)`, cellType: "relax" },
    { title: "Set LLG driver", code: `set_driver(sim; driver="LLG", alpha=0.02, gamma=2.211e5)`, cellType: "llg" },
    { title: "Add Zeeman Interaction", code: `add_zeeman(sim, (-24.6mT, 4.3mT, 0))`, cellType: "zeeman", isInteraction: true },
    { title: "Run Simulation", code: `run_sim(sim; steps=100, dt=1e-11, save_m_every=1)`, cellType: "run" }
];

export const examples ={
    "vortex": {
        title: "Vortex",
        type: "fd",
        task: "relax",
        steps: vortex
    },
    "std4": {
        title: "std4",
        type: "fd",
        task: "relax_dyn",
        steps: std4
    },
    "Ms_demo": {
        title: "Ms_demo",
        type: "fd",
        task: "demo",
        steps: Ms_demo
    },
}

import { Cell } from './Cell.js';
import { CellManager } from './CellManager.js';

/**
 * Create a CellManager instance for a specific example
 * @param {string} containerSelector - CSS selector for the container element
 * @param {string} exampleName - Name of the example
 * @param {GUIManager} guiManager - GUI manager instance
 * @returns {CellManager} - CellManager instance containing example steps
 */
export function createExampleTaskManager(containerSelector, exampleName, guiManager) {
    if (!examples[exampleName]) {
        console.error(`Example '${exampleName}' not found`);
        return null;
    }
    
    const example = examples[exampleName];
    const taskManager = new CellManager(containerSelector, example.title, guiManager);
            
    // Clear any existing cells
    taskManager.cells = [];
    if (taskManager.contentArea) {
        taskManager.contentArea.innerHTML = '';
    }
    
    // Add each step as a separate cell
    example.steps.forEach((step, index) => {
        const cell = new Cell(step.code, `${index + 1}. ${step.title}`, null, null, step.cellType, step.isInteraction || false);
        taskManager.addCell(cell);
    });
    
    return taskManager;
}