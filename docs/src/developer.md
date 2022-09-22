# Developer's guide



### Mesh
```@raw html
<div class="mermaid">
graph LR;
    Mesh --> MeshCPU
    MeshCPU --> FEMesh
    MeshCPU --> FDMesh
    Mesh --> MeshGPU
    MeshGPU --> FDMeshGPU
    MeshGPU --> CubicMeshGPU
    MeshGPU --> TriangularMeshGPU
</div>
```

### Sim
```@raw html
<div class="mermaid">
graph LR
   AbstractSim --> MicroSim
   AbstractSim --> NEB
   AbstractSim --> AbstractSimGPU
   AbstractSimGPU --> MicroSimGPU
   AbstractSimGPU --> AtomicSimGPU
   AbstractSimGPU --> MonteCarloGPU
   AbstractSimGPU --> NEB_GPU
</div>
```

### Driver
```@raw html
<div class="mermaid">
graph LR;
    Driver --> LLG
    Driver --> LLG_STT
    Driver --> EnergyMinimization
    DriverGPU --> LLG_GPU
    DriverGPU --> LLG_STT_GPU
    DriverGPU --> LLG_STT_CPU_GPU
    DriverGPU --> EnergyMinimization_GPU
</div>
```

### Implemented Energies
```@raw html
<div class="mermaid">
graph LR;
    MicroEnergy --> Exchange
    MicroEnergy --> VectorExchange
    MicroEnergy --> ExchangeRKKY
    MicroEnergy --> BulkDMI
    MicroEnergy --> SpatialBulkDMI
    MicroEnergy --> Zeeman
    MicroEnergy --> Anisotropy
    MicroEnergy --> CubicAnisotropy

    MicroEnergyGPU --> ExchangeGPU
    MicroEnergyGPU --> VectorExchangeGPU
    MicroEnergyGPU --> ExchangeRKKYGPU
    MicroEnergyGPU --> BulkDMIGPU
    MicroEnergyGPU --> SpatialBulkDMIGPU
    MicroEnergyGPU --> ZeemanGPU
    MicroEnergyGPU --> AnisotropyGPU
    MicroEnergyGPU --> CubicAnisotropyGPU
    MicroEnergyGPU --> StochasticFieldGPU
</div>
```

### Driver
```@raw html
<div class="mermaid">
graph LR;
    AbstractSim --> MicroSim
    AbstractSim --> NEB
    AbstractSim --> AbstractSimGPU
    AbstractSimGPU --> MicroSimGPU
    AbstractSimGPU --> AtomicSimGPU
    AbstractSimGPU --> MonteCarloGPU
    AbstractSimGPU --> NEB_GPU
</div>
```

