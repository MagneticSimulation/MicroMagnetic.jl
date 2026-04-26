import GUIManager from './GUIManager.js';
import Visualization from './Visualization.js';
import SimStatePanel from './SimStatePanel.js';
import SnippetPanel from './SnippetPanel.js';
import CodeEditorPanel from './CodeEditorPanel.js';
import PlotPanel from './PlotPanel.js';

// Export to global scope
window.Visualization = Visualization ;
window.GUIManager = GUIManager;

// Initialize GUI when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    console.log('Initializing MicroMagneticGUI...');

    // Create PlotPanel for data visualization
    const plotPanel = new PlotPanel();

    // Create and initialize visualization
    const containerId = 'visualization-container';
    if (document.getElementById(containerId)) {
        // Create visualization instance with container ID and plotPanel
        const visualization = new Visualization(containerId, plotPanel);

        // Create GUI manager instance
        const guiManager = new GUIManager(visualization);
        guiManager.initWebSocketClient();

        // Create SimStatePanel
        const simStatePanel = new SimStatePanel('sim-state-panel');
        guiManager.simStatePanel = simStatePanel;

        // Create CodeEditorPanel
        const codeEditorPanel = new CodeEditorPanel('code-editor-panel', guiManager);
        guiManager.setCodeEditorPanel(codeEditorPanel);

        // Create SnippetPanel and pass the codeEditorPanel
        const snippetPanel = new SnippetPanel('snippet-panel', codeEditorPanel);

        // Attach plotPanel to guiManager
        guiManager.plotPanel = plotPanel;

        console.log('MicroMagneticGUI initialized successfully!');
    } else {
        console.error('Visualization container not found');
    }
});