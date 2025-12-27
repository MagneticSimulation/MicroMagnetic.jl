/**
 * tasks.js
 * 
 * Task template definitions and management for MicroMagnetic.jl GUI
 */

import { Cell } from './cell.js';
import { CellManager } from './CellManager.js';

/**
 * Relax task template definition
 */
class RelaxTask {
    constructor() {
        // Define relax task steps as individual cells
        this.taskSteps = [
            { title: "Create Mesh", code: `mesh = FDMesh(; nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9);` },
            { title: "Create Simulation", code: `sim = Sim(mesh; driver="SD", name="std4")` },
            { title: "Set Saturation Magnetization", code: `Ms = 8.6e5; # A/m\nset_Ms(sim, Ms)` },
            { title: "Add Exchange Interaction", code: `A = 1.3e-11; # J/m\nadd_exch(sim, A)` },
            { title: "Add Demagnetization", code: `add_demag(sim)` },
            { title: "Initialize Magnetization", code: `init_m0(sim, (1, 0.25, 0.1))` },
            { title: "Relax System", code: `relax(sim; stopping_dmdt=0.01)` }
        ];
    }

    /**
     * Create a CellManager instance for relax task
     * @param {string} containerSelector - CSS selector for the container element
     * @returns {CellManager} - CellManager instance containing relax task steps
     */
    createTaskManager(containerSelector) {
        // Create a new CellManager instance for this task
        const taskManager = new CellManager(containerSelector);
        
        // Set the title to "relax"
        taskManager.title = "relax";
        
        // Add each step as a separate cell
        this.taskSteps.forEach((step, index) => {
            // Create a Cell instance for each step
            const cell = new Cell(step.code, `${index + 1}. ${step.title}`);
            taskManager.addCell(cell);
        });
        
        return taskManager;
    }
}

/**
 * Create a relax task manager
 * @param {string} containerSelector - CSS selector for the container element
 * @returns {CellManager} - CellManager instance for relax task
 */
export function createRelaxTaskManager(containerSelector) {
    const relaxTask = new RelaxTask();
    return relaxTask.createTaskManager(containerSelector);
}

/**
 * Initialize task templates system
 */
export function initTaskTemplates() {
    console.log('Task templates system initialized');
}