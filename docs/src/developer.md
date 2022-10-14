# Developer's guide

To develop JuMag, simply using
```julia
(@v1.7) pkg> dev JuMag
```
in the Julia REPL. The folder `$JULIA_DEPOT_PATH/dev/JuMag` should be created and you can modifiy the codes in it. 
For example, we could open the file `src/JuMag.jl` and add a function `dev_test` 
```julia
function dev_test()
    return "This is a newly added function!"
end
```
We can check it in a new Julia REPL:
```julia
julia> using JuMag
[ Info: Precompiling JuMag [8b6b6816-cea2-582c-a99f-83810c20db0f]

julia> JuMag.dev_test()
"This is a newly added function!"
```

After the modification, we can push our codes into github using 
```bash
git commit -m "we added a dev function" -a
git push
```

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

