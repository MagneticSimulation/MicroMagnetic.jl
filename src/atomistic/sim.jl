
"""
    set_mu_s(sim::MicroSim, Ms::NumberOrArrayOrFunction)

Set magnetic moment mu_s of the studied system. For example,

```julia
   set_mu_s(sim, 8.6e5)
```
or
```julia
function circular_shape(i,j,k,dx,dy,dz)
    if (i-50.5)^2 + (j-50.5)^2 <= 50^2
        return 8.6e5
    end
    return 0.0
end
set_mu_s(sim, circular_shape)
```
"""
function set_mu_s(sim::AtomicSimGPU, init::NumberOrArrayOrFunction)
    Float = _cuda_using_double.x ? Float64 : Float32
    Ms = zeros(Float, sim.nxyz)
    init_scalar!(Ms, sim.mesh, init)
    copyto!(sim.mu_s, Ms)
    return true
end

function set_mu_s_kagome(sim::AtomicSimGPU, Ms::Number)
    mesh = sim.mesh
    Float = _cuda_using_double.x ? Float64 : Float32
    mu_s = zeros(Float, sim.nxyz)
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
    nxyz = sim.nxyz
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*nxyz)
    energy = zeros(Float, nxyz)
    n_ngbs = sim.mesh.n_ngbs
    Js = CuArrays.zeros(Float, n_ngbs)

    if length(J) != n_ngbs
        @error("The length of given Js is $(length(Js)) but we need an array with $a.")
    else
        Js .= [Float(i) for i in J]
    end

    exch = ExchangeGPU(Js, field, energy, Float(0.0), name)
    push!(sim.interactions, exch)
    if sim.save_data
        push!(sim.saver.headers, string("E_",name))
        push!(sim.saver.units, "J")
        id = length(sim.interactions)
        push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
    end
    return exch
end


"""
    add_exch(sim::AtomicSimGPU, J::Number; name="exch")

Add exchange energy to the system.
"""
function add_exch(sim::AtomicSimGPU, J::Number; name="exch")
    n_ngbs = sim.mesh.n_ngbs
    Js = zeros(a)
    Js .= J
    add_exch(sim, Js, name=name)
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
  nxyz = sim.nxyz
  T = _cuda_using_double.x ? Float64 : Float32
  field = zeros(T, 3*nxyz)
  energy = zeros(T, nxyz)
  anis =  KagomeAnisotropy(Float(Ku), ax1, ax2, ax3, field, energy, T(0.0), name)
  push!(sim.interactions, anis)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return anis
end
