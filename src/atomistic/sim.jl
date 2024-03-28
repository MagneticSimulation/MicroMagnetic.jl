
"""
    set_mu_s(sim::AtomicSimGPU, Ms::NumberOrArrayOrFunction)

Set magnetic moment mu_s of the studied system. For example,

```julia
   set_mu_s(sim, 2*mu_B)
```
or
```julia
function circular_shape(i,j,k,dx,dy,dz)
    if (i-50.5)^2 + (j-50.5)^2 <= 50^2
        return 2*mu_B
    end
    return 0.0
end
set_mu_s(sim, circular_shape)
```
"""
function set_mu_s(sim::AtomicSimGPU, init::NumberOrArrayOrFunction)
    Float = _cuda_using_double.x ? Float64 : Float32
    Ms = zeros(Float, sim.n_nodes)
    init_scalar!(Ms, sim.mesh, init)
    copyto!(sim.mu_s, Ms)
    return true
end

function set_mu_s_kagome(sim::AtomicSimGPU, Ms::Number)
    mesh = sim.mesh
    Float = _cuda_using_double.x ? Float64 : Float32
    mu_s = zeros(Float, sim.n_nodes)
    for k = 1:mesh.nz, j = 1:mesh.ny, i = 1:mesh.nx
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        mu_s[id] = Ms
        if i%2==0 && j%2==0
            mu_s[id] = 0.0
        end
    end
    copyto!(sim.mu_s, mu_s)
    return true
end


"""
    add_exch(sim::AtomicSimGPU, J::Array; name="exch")

Add exchange energy to the system. The length of J should be equal to the length of neigbours.
"""
function add_exch(sim::AtomicSimGPU, J::Array; name="exch")
    n_nodes = sim.n_nodes
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*n_nodes)
    energy = zeros(Float, n_nodes)
    n_ngbs = sim.mesh.n_ngbs
    Js = CUDA.zeros(Float, n_ngbs)

    if length(J) != n_ngbs
        @error("The length of given Js is $(length(Js)) but we need an array with $a.")
    else
        copyto!(Js, [Float(i) for i in J])
    end

    exch = HeisenbergExchange(Js, field, energy, Float(0.0), name)
    push!(sim.interactions, exch)
    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_",name), "J", o::AbstractSim->o.interactions[id].total_energy)
        push!(sim.saver.items, item)
    end
    return exch
end


"""
    add_exch(sim::AtomicSimGPU, J::Number; name="exch")

Add exchange energy to the system.
"""
function add_exch(sim::AtomicSimGPU, J::Number; name="exch")
    n_ngbs = sim.mesh.n_ngbs
    Js = zeros(n_ngbs)
    Js .= J
    add_exch(sim, Js, name=name)
end


"""
    add_next_exch()

Add next-nearest neigbours exchange energy to the system
"""

function add_next_exch(sim::AtomicSimGPU, J::Number; name="next_exch")
    nn_ngbs = sim.mesh.nn_ngbs
    Js = zeros(nn_ngbs)
    Js .= J
    add_next_exch(sim, Js, name=name)
end

function add_next_exch(sim::AtomicSimGPU, J::Array; name="next_exch")
    n_nodes = sim.n_nodes
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*n_nodes)
    energy = zeros(Float, n_nodes)
    nn_ngbs = sim.mesh.nn_ngbs
    Js = CUDA.zeros(Float, nn_ngbs)

    if length(J) != nn_ngbs
        @error("The length of given Js is $(length(Js)) but we need an array with $a.")
    else
        copyto!(Js, [Float(i) for i in J])
    end

    next_exch = NextHeisenbergExchange(Js, field, energy, Float(0.0), name)
    push!(sim.interactions, next_exch)
    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_",name), "J", o::AbstractSim->o.interactions[id].total_energy)
        push!(sim.saver.items, item)
    end
    return next_exch
end

"""
    add_next_next_exch()

Add next-next-nearest neigbours exchange energy to the system
"""
function add_next_next_exch(sim::AtomicSimGPU, J::Number; name="next_next_exch")
    nnn_ngbs = sim.mesh.nnn_ngbs
    Js = zeros(nnn_ngbs)
    Js .= J
    add_next_next_exch(sim, Js, name=name)
end

function add_next_next_exch(sim::AtomicSimGPU, J::Array; name="next_next_exch")
    n_nodes = sim.n_nodes
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*n_nodes)
    energy = zeros(Float, n_nodes)
    nnn_ngbs = sim.mesh.nnn_ngbs
    Js = CUDA.zeros(Float, nnn_ngbs)

    if length(J) != nnn_ngbs
        @error("The length of given Js is $(length(Js)) but we need an array with $a.")
    else
        copyto!(Js, [Float(i) for i in J])
    end

    next_next_exch = NextNextHeisenbergExchange(Js, field, energy, Float(0.0), name)
    push!(sim.interactions, next_next_exch)
    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_",name), "J", o::AbstractSim->o.interactions[id].total_energy)
        push!(sim.saver.items, item)
    end
    return next_next_exch
end


"""
    add_next_next_exch()

Add 4th next_next-next-nearest neigbours exchange energy to the system
"""
function add_next_next_next_exch(sim::AtomicSimGPU, J::Number; name="next_next_next_exch")
    nnnn_ngbs = sim.mesh.nnnn_ngbs
    Js = zeros(nnnn_ngbs)
    Js .= J
    add_next_next_next_exch(sim, Js, name=name)
end

function add_next_next_next_exch(sim::AtomicSimGPU, J::Array; name="next_next_next_exch")
    n_nodes = sim.n_nodes
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3 * n_nodes)
    energy = zeros(Float, n_nodes)
    nnnn_ngbs = sim.mesh.nnnn_ngbs
    Js = CUDA.zeros(Float, nnnn_ngbs)

    if length(J) != nnnn_ngbs
        @error("The length of given Js is $(length(Js)) but we need an array with $a.")
    else
        copyto!(Js, [Float(i) for i in J])
    end

    next_next_next_exch = NextNextNextHeisenbergExchange(Js, field, energy, Float(0.0), name)
    push!(sim.interactions, next_next_next_exch)
    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_", name), "J", o::AbstractSim -> o.interactions[id].total_energy)
        push!(sim.saver.items, item)
    end
    return next_next_next_exch
end





@doc raw"""
    add_dmi(sim::AtomicSimGPU, D::Real; name="dmi")

Add bulk dmi energy to the system. The bulk dmi is defined as
```math
\mathcal{H}_\mathrm{dmi} =  \sum_{\langle i, j\rangle}  \mathbf{D}_{i j} \cdot\left(\mathbf{m}_{i} \times \mathbf{m}_{j}\right)
```
where $\mathbf{D}_{i j}$ is the DM vector. For the bulk dmi as implemented in this function, 
$\mathbf{D}_{i j} = D \hat{r}_{ij}$.
"""
function add_dmi(sim::AtomicSimGPU, D::Real; name="dmi")
    n_nodes = sim.n_nodes
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*n_nodes)
    energy = zeros(Float, n_nodes)
    n_ngbs = sim.mesh.n_ngbs

    # TODO: implement the effective field for both TriangularMeshGPU and CubicMeshGPU
    # TODO: check the sign of D
    if isa(sim.mesh, TriangularMeshGPU)
        error("bulk dmi for triangular mesh has not implemented yet")
    elseif isa(sim.mesh, CubicMeshGPU)
        dmi = HeisenbergBulkDMI(Float(D), field, energy, Float(0.0), name)
    elseif isa(sim.mesh, CylindricalTubeMeshGPU)
        nr = sim.mesh.nr
        coords = sim.mesh.coordinates
        ngbs = Array(sim.mesh.ngbs)
        #Dij stores D*r_{ij}, i.e., r_{n_r, 1}, r12, r23, ..., r_{(n_r-1)n_r},  r_{n_r, 1}
        Dij = zeros(Float, 3, nr+1)
        for i = 1:nr
            j = ngbs[1, i]  # the left one
            rx = coords[1, i] - coords[1, j] 
            ry = coords[2, i] - coords[2, j]
            rz = coords[3, i] - coords[3, j]
            r = sqrt(rx*rx+ry*ry+rz*rz)
            Dij[1, i] = D*rx/r
            Dij[2, i] = D*ry/r
            Dij[3, i] = D*rz/r
        end
        Dij[:, nr+1] .= Dij[:, 1]
        dmi = HeisenbergTubeBulkDMI(Float(D), CuArray(Dij), field, energy, Float(0.0), name)
    end
        
    push!(sim.interactions, dmi)
    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_",name), "J", o::AbstractSim->o.interactions[id].total_energy)
        push!(sim.saver.items, item)
    end
    return dmi
end


"""
    add_exch_kagome(sim::AtomicSimGPU, Jxy::Number, Jz::Number; name="exch")

Add exchange energy to the system.
"""
function add_exch_kagome(sim::AtomicSimGPU, Jxy::Number, Jz::Number; name="exch")
    n_ngbs = sim.mesh.n_ngbs
    Js = zeros(n_ngbs)

    if n_ngbs!=8
        error("The number of neigbours is not 8.")
    end

    Js[1:6] .= Jxy
    Js[7:8] .= Jz
    add_exch(sim, Js, name=name)
end

"""
    add_anis_kagome(sim::AtomicSimGPU, Ku::Float64; ax1=(-0.5,-sqrt(3)/2,0), ax2=(1,0,0), ax3=(-0.5,sqrt(3)/2,0), name="anis")

Add Anisotropy for kagome system, where the energy density is given by

```math
E_\\mathrm{anis} = - K_{u} (\\vec{m} \\cdot \\hat{u})^2
```
"""
function add_anis_kagome(sim::AtomicSimGPU, Ku::Float64; ax1=(-0.5,-sqrt(3)/2,0), ax2=(1,0,0), ax3=(-0.5,sqrt(3)/2,0), name="anis")
  n_nodes = sim.n_nodes
  T = _cuda_using_double.x ? Float64 : Float32
  field = zeros(T, 3*n_nodes)
  energy = zeros(T, n_nodes)
  anis =  KagomeAnisotropy(T(Ku), ax1, ax2, ax3, field, energy, T(0.0), name)
  push!(sim.interactions, anis)

  if sim.save_data
    id = length(sim.interactions)
    item = SaverItem(string("E_",name), "J",  o::AbstractSim->sum(o.interactions[id].energy))
    push!(sim.saver.items, item)
  end
  return anis
end

@doc raw"""
    add_anis_tube(sim::AtomicSimGPU, Ku::Float64; name="anis")

add anisotropy to the system when the tube mesh is used. The anisotropy axis $u$ 
is along with the radial direction.

```math
E_\mathrm{anis} = - K_{u} (\vec{m} \cdot \hat{u})^2
```
"""
function add_anis_tube(sim::AtomicSimGPU, Ku::Float64; name="anis")
    n_nodes = sim.n_nodes
    T = _cuda_using_double.x ? Float64 : Float32
    field = zeros(T, 3*n_nodes)
    energy = zeros(T, n_nodes)
    axes = zeros(T, 3, n_nodes)

    nr = sim.mesh.nr
    for i = 1:n_nodes
        theta = 2*pi*(i-1)/nr
        axes[1, i] = cos(theta)
        axes[2, i] = sin(theta)
    end

    anis =  TubeAnisotropy(T(Ku), CuArray(axes), field, energy, T(0.0), name)
    push!(sim.interactions, anis)
  
    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_",name), "J",  o::AbstractSim->sum(o.interactions[id].energy))
        push!(sim.saver.items, item)
    end
    return anis
  end


function add_thermal_noise(sim::AtomicSimGPU, T::NumberOrArrayOrFunction; name="thermal", k_B=k_B)
    n_nodes = sim.n_nodes
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*n_nodes)
    energy = zeros(Float, n_nodes)
    Spatial_T = CUDA.zeros(Float, n_nodes)
    eta = CUDA.zeros(Float, 3*n_nodes)
    init_scalar!(Spatial_T , sim.mesh, T)
    thermal = StochasticFieldGPU(Spatial_T, eta, field, energy, Float(0.0), -1, name, k_B)
  
    push!(sim.interactions, thermal)
  
    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_",name), "J", o::AbstractSim->o.interactions[id].exi)
        push!(sim.saver.items, item)
    end
    return thermal
end

"""
    add_magnetoelectric_laser(sim::AtomicSimGPU, lambda::Float64, E::Float64, B::Float64, omega::Float64; delta=0, direction=001, name="lasers")

Add the interaction of high-frequency lasers to multiferroic insulator Cu2OSeO3. 
The Hamiltonian is given by 

```math
\\mathcal{H}_\\mathrm{laser} =  -\\sum_{i} \\mu_s \\mathbf{m}_i \\cdot \\mathbf{B}(t) - \\sum_{i} \\mathbf{P}_i \\cdot \\mathbf{E}(t)
```

The high-frequency laser is described as 

```math
\\mathbf{E}(t) =  E ( \\sin (\\omega t + \\delta), \\cos \\omega t, 0) \\qquad
\\mathbf{B}(t) =  B ( \\cos \\omega t, -\\sin(\\omega t + \\delta), 0)
```
where δ determines the laser polarization, i.e., δ = 0 for right-circularly polarized (RCP), 
δ=π/2 for linearly polarized and δ=π for left-circularly polarized (LCP).
"""
function add_magnetoelectric_laser(sim::AtomicSimGPU, lambda::Float64, E::Float64, B::Float64, omega::Float64; delta=0, direction=001, name="lasers")
    n_nodes = sim.n_nodes
    F = _cuda_using_double.x ? Float64 : Float32
    field = zeros(F, 3*n_nodes)
    energy = zeros(F, n_nodes)

    laser = MagnetoelectricLaser(F(lambda), F(E),  F(B), F(omega), F(delta), 
                                direction, field, energy, F(0.0), name)
  
    push!(sim.interactions, laser)
  
    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_",name), "J", o::AbstractSim->o.interactions[id].total_energy)
        push!(sim.saver.items, item)
    end
    return laser
end

@doc raw"""
    add_demag(sim::AtomicSimGPU; name="demag", Nx=0, Ny=0, Nz=0 )

add dipolar interaction into the system.

```math
\mathcal{H}_{\mathrm{d}}=-\frac{\mu_0 \mu_s^2}{4 \pi} \sum_{i<j} \frac{3\left(\mathbf{m}_i \cdot \hat{\mathbf{r}}_{i j}\right)\left(\mathbf{m}_j \cdot \hat{\mathbf{r}}_{i j}\right)-\mathbf{m}_i \cdot \mathbf{m}_j}{r_{i j}^3}
```
"""
function add_demag(sim::AtomicSimGPU; name="demag", Nx=0, Ny=0, Nz=0)
    demag = init_demag_gpu(sim, Nx, Ny, Nz)
    demag.name = name
    push!(sim.interactions, demag)

    if sim.save_data
        id = length(sim.interactions)
        item = SaverItem(string("E_",name), "J", o::MicroSimGPU->o.interactions[id].total_energy)
        push!(sim.saver.items, item)
    end
    return demag
end
