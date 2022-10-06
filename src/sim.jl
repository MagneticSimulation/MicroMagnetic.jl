"""
    Sim(mesh::Mesh; driver="LLG", name="dyn", integrator="DormandPrince")

Create a simulation instance for given mesh.

"""
function Sim(mesh::FDMesh; driver="LLG", name="dyn", integrator="DormandPrince", save_data=true)

    sim = MicroSim()

    sim.name = name
    sim.mesh = mesh
    nxyz = mesh.nx*mesh.ny*mesh.nz
    sim.nxyz = nxyz
    sim.spin = zeros(Float64,3*nxyz)
    sim.prespin = zeros(Float64,3*nxyz)
    sim.field = zeros(Float64,3*nxyz)
    sim.energy = zeros(Float64,nxyz)

    sim.Ms = zeros(Float64, nxyz)
    sim.pins = zeros(Bool, nxyz)
    sim.save_data = save_data

    if save_data
        headers = ["step", "E_total", ("m_x", "m_y", "m_z")]
        units = ["<>", "<J>",("<>", "<>", "<>")]
        results = [o::AbstractSim -> o.saver.nsteps,
                o::AbstractSim -> sum(o.energy),
                average_m]
        sim.saver = DataSaver(string(name, "_", lowercase(driver),".txt"), 0.0, 0, false, headers, units, results)
    end

    if driver in ("LLG", "LLG_STT", "LLG_STT_CPP") && save_data
        saver = sim.saver
        insert!(saver.headers, 2, "time")
        insert!(saver.units, 2, "<s>")
        insert!(saver.results, 2, o::AbstractSim -> o.saver.t)
    end
    sim.driver_name = driver
    sim.driver = create_driver(driver, integrator, nxyz)
    
   sim.interactions = []
   #println(sim)
   return sim
end

"""
    set_Ms(sim::MicroSim, Ms::NumberOrArrayOrFunction)

Set the saturation magnetization Ms of the studied system. For example,

```julia
   set_Ms(sim, 8.6e5)
```
or
```julia
function circular_Ms(i,j,k,dx,dy,dz)
    if (i-50.5)^2 + (j-50.5)^2 <= 50^2
        return 8.6e5
    end
    return 0.0
end
set_Ms(sim, circular_Ms)
```

"""
function set_Ms(sim::MicroSim, Ms::NumberOrArrayOrFunction)
    init_scalar!(sim.Ms, sim.mesh, Ms)

    if any(isnan, sim.Ms)
        error("NaN is given by the input Ms!")
    end

    return true
end

"""
    set_Ms_cylindrical(sim::MicroSim, Ms::Number; axis=ez, r1=0, r2=0)

Set the saturation magnetization Ms of the studied system to a cylindrical shape,
axis could be JuMag.ex, JuMag.ey or JuMag.ez.
"""
function set_Ms_cylindrical(sim::MicroSim, Ms::Number; axis=ez, r1=0, r2=0)
    geo = create_cylinder(sim.mesh, axis, r1=r1, r2=r2)
    for i in 1:sim.nxyz
        if geo.shape[i]
            sim.Ms[i] = Ms
        end
    end
    return true
end

"""
    set_Ms(sim::AbstractSim, geo::Geometry, Ms::Number)

Set the saturation magnetization Ms within the Geometry.
"""
function set_Ms(sim::AbstractSim, geo::Geometry, Ms::Number)
    update_scalar_geometry(sim.Ms, geo, Ms)
    return true
end


"""
    set_pinning(sim::MicroSim, ids::ArrayOrFunction)

Pinning the spins at the given ids.

```julia
function pinning_boundary(i,j,k,dx,dy,dz)
    if i == 1 || i == 100
        return true
    end
    return false
end
set_pinning(sim, pinning_boundary)
```
"""
function set_pinning(sim::MicroSim, ids::ArrayOrFunction)
    init_scalar!(sim.pins, sim.mesh, ids)
    return true
end

function set_ux(sim::AbstractSim, init_ux)
    init_scalar!(sim.driver.ux, sim.mesh, init_ux)
end

function set_ux_bounary(sim::AbstractSim, ux)
    return set_ux_bounary_implement(sim, ux)
end

function set_uy(sim::AbstractSim, init_uy)
    init_scalar!(sim.driver.uy, sim.mesh, init_uy)
end

function set_uz(sim::AbstractSim, init_uz)
    init_scalar!(sim.driver.uz, sim.mesh, init_uz)
end

function set_aj(sim::AbstractSim, init_aj)
	init_scalar!(sim.driver.aj, sim.mesh, init_aj)
end

function average_m(sim::AbstractSim)
  b = reshape(sim.spin, 3, sim.nxyz)
  mx,my,mz = 0.0,0.0,0.0
  n = 0
  ms = isa(sim, MicroSim) ? sim.Ms : sim.mu_s
  for i = 1:sim.nxyz
    if ms[i]>0
      n += 1
      mx += b[1,i]
      my += b[2,i]
      mz += b[3,i]
    end
  end
  if n == 0
    error("n should not be zero!")
  end
  return (mx/n, my/n, mz/n)
end

"""
    init_m0(sim::MicroSim, m0::TupleOrArrayOrFunction; norm=true)

Set the initial magnetization of the system. If `norm=false` the magnetization array will be not normalised.
Examples:

```julia
   init_m0(sim, (1,1,1))
```
or
```julia
   init_m0(sim, (1,1,1), norm=false)
```
or
```julia
   function uniform_m0(i,j,k,dx,dy,dz)
       return (0,0,1)
   end
   init_m0(sim, uniform_m0)
```
"""
function init_m0(sim::MicroSim, m0::TupleOrArrayOrFunction; norm=true)
  init_vector!(sim.prespin, sim.mesh, m0)
  if norm
    normalise(sim.prespin, sim.nxyz)
  end
  for i = 1:sim.nxyz
      if sim.Ms[i] == 0.0
          sim.prespin[3*i-2] = 0
          sim.prespin[3*i-1] = 0
          sim.prespin[3*i] = 0
      end
  end

  if any(isnan, sim.prespin)
    error("NaN is given by the input m0!")
  end

  sim.spin[:] .= sim.prespin[:]
end

"""
    set_driver(sim::AbstractSim; driver="LLG", integrator="DormandPrince")

Set the driver of the simulation, can be used to switch the driver.
"""
function set_driver(sim::AbstractSim; driver="LLG", integrator="DormandPrince")

    if sim.driver_name == driver
      return nothing
    end
  
    # if the driver is update, we create a new saver
    if sim.save_data
      headers = ["step", "E_total", ("m_x", "m_y", "m_z")]
      units = ["<>", "<J>",("<>", "<>", "<>")]
      results = [o::AbstractSim -> o.saver.nsteps,
                 o::AbstractSim -> o.total_energy, average_m]
      sim.saver = DataSaver(string(sim.name, "_", lowercase(driver), ".txt"), 0.0, 0, false, headers, units, results)
    end
    
    if isa(sim, AbstractSimGPU)
        sim.driver = create_driver_gpu(driver, integrator, sim.nxyz)
    else
        sim.driver = create_driver(driver, integrator, sim.nxyz)
    end
    sim.driver_name = driver
  
    if startswith(driver,"LLG") && sim.save_data
      saver = sim.saver
      insert!(saver.headers, 2, "time")
      insert!(saver.units, 2, "<s>")
      insert!(saver.results, 2, o::AbstractSim -> o.saver.t)
    end
  
end

"""
    add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")

Add a static Zeeman energy to the simulation.
"""
function add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")
    # FIXME: We need to unify the two variable names (nxyz and n_nodes)
    if isa(sim, MicroSimFEM)
        nxyz = sim.n_nodes
    else
        nxyz = sim.nxyz
    end

  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  init_vector!(field, sim.mesh, H0)

  zeeman =  Zeeman(field, energy, name)
  push!(sim.interactions, zeeman)

  if sim.save_data
      if isa(H0, Tuple)
          push!(sim.saver.headers, (string(name, "_Hx"), string(name, "_Hy"), string(name, "_Hz")))
          push!(sim.saver.units, ("<A/m>", "<A/m>", "<A/m>"))
          id = length(sim.interactions)
          fun = o::AbstractSim ->  (o.interactions[id].field[1], o.interactions[id].field[2], o.interactions[id].field[3])
          push!(sim.saver.results, fun)
      end

      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return zeeman
end

"""
    update_zeeman(sim::AbstractSim, H0::Tuple; name="zeeman")

Set the Zeeman field to H0 where H0 is TupleOrArrayOrFunction according to its name. For example,

```julia
   add_zeeman(sim, (0,0,0), name="my_H")  #create a zeeman energy with field (0,0,0) A/m
   update_zeeman(sim, (0,0,1e5), name="my_H")  #change the field to (0,0,1e5) A/m
```

"""
function update_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")
    N_spins = isa(sim, MicroSimFEM) ? sim.n_nodes : sim.nxyz
    field = zeros(Float64, 3*N_spins)
    init_vector!(field, sim.mesh, H0)

    for i in sim.interactions
        if i.name == name
            i.field[:] = field[:]
            return nothing
        end
    end
    return nothing
end

"""
    add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction, ft::Function; name="timezeeman")

Add a time varying zeeman to system.

The input `ft` is a function of time `t` and its return value should be a tuple with length 3.

Example:

```julia
  function time_fun(t)
    w = 2*pi*2.0e9
    return (sin(w*t), cos(w*t), 0)
  end

  function spatial_H(i, j, k, dx, dy, dz)
    H = 1e3
    if i<=2
        return (H, H, 0)
    end
    return (0, 0, 0)
  end

  add_zeeman(sim, spatial_H, time_fun)
```
"""
function add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction, ft::Function; name="timezeeman")
  nxyz = sim.nxyz
  init_field = zeros(Float64, 3*nxyz)
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  init_vector!(init_field, sim.mesh, H0)

  zeeman =  TimeZeeman(ft, init_field, field, energy, name)
  push!(sim.interactions, zeeman)

  if sim.save_data
    id = length(sim.interactions)
      if isa(H0, Tuple)
          push!(sim.saver.headers, (string(name, "_Hx"), string(name, "_Hy"), string(name, "_Hz")))
          push!(sim.saver.units, ("<A/m>", "<A/m>", "<A/m>"))
          fun = o::AbstractSim ->  (o.interactions[id].field[1], o.interactions[id].field[2], o.interactions[id].field[3])
          push!(sim.saver.results, fun)
      end

      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return zeeman
end

"""
    add_exch_vector(sim::AbstractSim, A::TupleOrArrayOrFunction; name="exch")

Add a vector form exchange energy to the system. The exchange constant of 3 directions can be different.
For example:
```julia
add_exc_vector(sim, (2e-12,5e-12,0))
```
"""
function add_exch_vector(sim::AbstractSim, A::TupleOrArrayOrFunction; name="exch_vector")
  nxyz = sim.nxyz
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  Spatial_A = zeros(Float64, 3*nxyz)
  init_vector!(Spatial_A , sim.mesh, A)
  exch = Vector_Exchange(Spatial_A , field, energy, name)
  push!(sim.interactions, exch)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return exch
end

"""
    add_exch(sim::AbstractSim, A::NumberOrArrayOrFunction; name="exch")

Add exchange energy to the system.
"""
function add_exch(sim::AbstractSim, A::NumberOrArrayOrFunction; name="exch")
    if isa(sim, MicroSimFEM)
        nxyz = sim.n_nodes
        Spatial_A = zeros(Float64, sim.n_cells)
        K_mat = spzeros(3*nxyz, 3*nxyz)
    else
        nxyz = sim.nxyz
        Spatial_A = zeros(Float64, sim.nxyz)
    end
  
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)

  init_scalar!(Spatial_A , sim.mesh, A)
  if isa(sim, MicroSim) 
    exch = Exchange(Spatial_A , field, energy, name)
  elseif isa(sim, MicroSimFEM)
    exch = ExchangeFEM(Spatial_A , field, energy, K_mat, false, name)
  else
	exch = HeisenbergExchange(A, field, energy, name)
  end
  push!(sim.interactions, exch)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return exch
end

"""
    add_exch(sim::AbstractSim, geo::Geometry, A::Number; name="exch")

Add exchange interaction within the Geometry, or update corresponding A if other exch is added.
"""
function add_exch(sim::AbstractSim, geo::Geometry, A::Number; name="exch")
  for interaction in sim.interactions
      if interaction.name == name
         update_scalar_geometry(interaction.A, geo, A)
         return nothing
      end
  end
  nxyz = sim.nxyz
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  Spatial_A = zeros(Float64, nxyz)
  update_scalar_geometry(Spatial_A, geo, A)
  if isa(sim, MicroSim)
    exch = Exchange(Spatial_A , field, energy, name)
  else
  exch = HeisenbergExchange(A, field, energy, name)
  end
  push!(sim.interactions, exch)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return exch
end

function add_exch_rkky(sim::AbstractSim, sigma::Float64, Delta::Float64; name="rkky")
  nxyz = sim.nxyz
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  exch = ExchangeRKKY(sigma, Delta, field, energy, name)

  push!(sim.interactions, exch)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return exch
end

"""
    add_dmi(sim::AbstractSim, D::Tuple{Real, Real, Real}; name="dmi")

Add DMI to the system. Example:

```julia
   add_dmi(sim, (1e-3, 1e-3, 0))
```
"""
function add_dmi(sim::AbstractSim, D::Tuple{Real, Real, Real}; name="dmi")
  nxyz = sim.nxyz
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  dmi =  BulkDMI(Float64(D[1]), Float64(D[2]), Float64(D[3]), field, energy, name)
  push!(sim.interactions, dmi)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return dmi
end

"""
    add_dmi(sim::AbstractSim, D::Real; name="dmi", type="bulk")

Add DMI to the system. `type` could be "bulk" or "interfacial"
Examples:

```julia
   add_dmi(sim, 1e-3, type="interfacial")
```
or
```julia
   add_dmi(sim, 1e-3, type="bulk")
```
"""
function add_dmi(sim::AbstractSim, D::Real; name="dmi", type="bulk")
    if type == "interfacial"
        return add_dmi_interfacial(sim, D, name=name)
    elseif type == "bulk"
        return add_dmi(sim, (D,D,D), name=name)
    else
        error("Supported DMI type:", "interfacial", "bulk")
    end

end

function add_dmi_interfacial(sim::AbstractSim, D::Real; name="dmi")
    nxyz = sim.nxyz
    field = zeros(Float64, 3*nxyz)
    energy = zeros(Float64, nxyz)
    dmi =  InterfacialDMI(Float64(D), field, energy, name)
    push!(sim.interactions, dmi)

    if sim.save_data
        push!(sim.saver.headers, string("E_",name))
        push!(sim.saver.units, "J")
        id = length(sim.interactions)
        push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
    end
    return dmi
end

"""
    add_dmi(sim::AbstractSim;  name="dmi")

Add DMI to the system. Example:

```julia
   add_dmi(sim, (1e-3, 1e-3, 0))
```
"""
function add_dmi(sim::AbstractSim, name="dmi")
  nxyz = sim.nxyz
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  z = 0.0
  dmi =  DMI(z,z,z,z,z,z,z,z,z,z,z,z,z,z, field, energy, name)
  push!(sim.interactions, dmi)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return dmi
end

"""
    add_demag(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0)

Add Demag to the system. `Nx`, `Ny` and `Nz` can be used to describe the macro boundary conditions which means that
the given mesh is repeated `2Nx+1`, `2Ny+1 and `2Nz+1` times in `x`, `y` and `z` direction, respectively.
"""
function add_demag(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0)
  demag = init_demag(sim, Nx, Ny, Nz)
  demag.name = name
  push!(sim.interactions, demag)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::MicroSim->sum(o.interactions[id].energy))
  end
  return demag
end

"""
    add_anis(sim::AbstractSim, Ku::NumberOrArrayOrFunction; axis=(0,0,1), name="anis")

Add Anisotropy to the system, where the energy density is given by

```math
E_\\mathrm{anis} = - K_{u} (\\vec{m} \\cdot \\hat{u})^2
```
"""
function add_anis(sim::AbstractSim, Ku::NumberOrArrayOrFunction; axis=(0,0,1), name="anis")
  nxyz = sim.nxyz
  Kus =  zeros(Float64, nxyz)
  init_scalar!(Kus, sim.mesh, Ku)
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  lt = sqrt(axis[1]^2+axis[2]^2+axis[3]^2)
  naxis = (axis[1]/lt, axis[2]/lt, axis[3]/lt)
  anis =  Anisotropy(Kus, naxis, field, energy, name)
  push!(sim.interactions, anis)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return anis
end

"""
    add_anis(sim::AbstractSim, geo::Geometry, Ku::Number; axis=(0,0,1), name="anis")

Add Anisotropy within the Geometry, or update corresponding Ku if other anis is added.
"""
function add_anis(sim::AbstractSim, geo::Geometry, Ku::Number; axis=(0,0,1), name="anis")
  lt = sqrt(axis[1]^2+axis[2]^2+axis[3]^2)
  naxis = (axis[1]/lt, axis[2]/lt, axis[3]/lt)
  for interaction in sim.interactions
      if interaction.name == name
         update_scalar_geometry(interaction.Ku, geo, Ku)
         return nothing
      end
  end
  nxyz = sim.nxyz
  Kus =  zeros(Float64, nxyz)
  update_scalar_geometry(Kus, geo, Ku)
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  anis =  Anisotropy(Kus, naxis, field, energy, name)
  push!(sim.interactions, anis)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return anis
end
"""
    update_anis(sim::MicroSim, Ku::NumberOrArrayOrFunction; name = "anis")

update anisotropy constant Ku according to its name.

Example:
```julia
    mesh = FDMesh(nx=200, ny=200, nz=12, dx=5e-9, dy=5e-9, dz=5e-9)
    sim = Sim(mesh)
    add_anis(sim, 3e4, axis = (0,0,1), name="K1")
    add_anis(sim, 1e5, axis = (1,0,0), name="K2")
    update_anis(sim, 5e4, name="K2")  #update anisotropy K2
```
"""
function update_anis(sim::MicroSim, Ku::NumberOrArrayOrFunction; name = "anis")
  nxyz = sim.nxyz
  Kus =  zeros(Float64, nxyz)
  init_scalar!(Kus, sim.mesh, Ku)
  for i in sim.interactions
    if i.name == name
      i.Ku[:] = Kus[:]
      return nothing
    end
  end
  return nothing
end

"""
    add_cubic_anis(sim::AbstractSim, Kc::Float64; axis1::Any=nothing, axis2::Any=nothing, name="cubic")

add cubic anisotropy with default axis (1,0,0) , (0,1,0), and (0,0,1)
use axis1=(1,1,0), axis2=(1,-1,0) to set a pair of normal axis, and the third axis will be calculated automatically
"""
function add_cubic_anis(sim::AbstractSim, Kc::Float64; axis1::Any=nothing, axis2::Any=nothing, name="cubic")
  axis1 = axis1 == nothing ? (1,0,0) : axis1
  axis2 = axis2 == nothing ? (0,1,0) : axis2
  norm1 = sqrt(axis1[1]^2+axis1[2]^2+axis1[3]^2)
  norm2 = sqrt(axis2[1]^2+axis2[2]^2+axis2[3]^2)
  naxis1,naxis2 = axis1./norm1,axis2./norm2
  if abs.(sum(naxis1.*naxis2)) > 1e-10
    @error("cubic axis not normal!")
    return nothing
  end
  naxis3= cross_product(axis1,axis2)
  axis = zeros(Float64,9)
  for i = 1:3
    axis[i] = naxis1[i]
    axis[i+3] = naxis2[i]
    axis[i+6] = naxis3[i]
  end

  nxyz = sim.nxyz
  field = zeros(Float64, 3*nxyz)
  energy = zeros(Float64, nxyz)
  anis =  CubicAnisotropy(axis, Kc, field, energy, name)
  push!(sim.interactions, anis)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return anis
end

"""
    relax(sim::AbstractSim; maxsteps=10000, stopping_dmdt=0.01, save_m_every = 10, save_ovf_every=-1, ovf_format = "binary", ovf_folder="ovfs", save_vtk_every=-1, vtk_folder="vtks",fields::Array{String, 1} = String[])

Relax the system using `LLG` or `SD` driver. This function works for both micromagnetic (FD and FE) and atomistic simulations, in both CPU and GPU. 

`maxsteps` is the maximum steps allowed to run. 

`stopping_dmdt` is the main stop condition, both for both for `LLG` and `SD` drivers. For standard micromagnetic simulaition, 
the typical value of `stopping_dmdt` is in the range of [0.01, 1].  In the `SD` driver, the time is not strictly defined. 
To make it comparable for the `LLG` driver, we multiply a factor of `gamma`. However, for the atomistic model 
with dimensionless unit, this factor should not be used. In this situation, `using_time_factor` should be set to `false`.

The magnetization (spins) can be stored in ovfs or vtks. ovf format can be chosen in "binary"(float64),"binary8"(float64), "binary4"(float32), "text"

Fields can be stored in vtks as well

```julia
relax(sim, save_vtk_every = 10, fields = ["demag", "exch", "anis"])
```
"""
function relax(sim::AbstractSim; maxsteps=10000, stopping_dmdt=0.01, using_time_factor=true, save_m_every = -1, 
    save_ovf_every=-1, ovf_format = "binary", ovf_folder="ovfs", save_vtk_every=-1, vtk_folder="vtks", fields::Array{String, 1} = String[])

    # to dertermine which driver is used.
    llg_driver = false
    if isa(sim.driver, LLG) || (_cuda_available.x && isa(sim.driver, LLG_GPU))
        llg_driver = true
    end

    time_factor =  using_time_factor ? 2.21e5/2 : 1.0

    N_spins = isa(sim, MicroSimFEM) ? sim.n_nodes : sim.nxyz

    if _cuda_available.x && (isa(sim, MicroSimGPU) || isa(sim, AtomicSimGPU))
        T = _cuda_using_double.x ? Float64 : Float32
        dm = CUDA.zeros(T, 3*sim.nxyz)
    else
        dm = zeros(Float64,3*N_spins)
    end


    dmdt_factor = (2 * pi / 360) * 1e9
    if _cuda_available.x && isa(sim, AtomicSimGPU)
        dmdt_factor = 1.0
    end

    step = 0
    driver = sim.driver
    @info @sprintf("Running Driver : %s.", typeof(driver))
    for i=1:maxsteps

      run_step(sim, driver)

      step_size = llg_driver ? driver.ode.step : driver.tau/time_factor

      compute_dm!(dm, sim.prespin, sim.spin, N_spins)
      max_dmdt = maximum(dm)/step_size

      t = llg_driver ? sim.driver.ode.t : 0.0
      if llg_driver
          @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e",
                        i, step_size, t, max_dmdt/dmdt_factor)
      else
          @info @sprintf("step =%5d  step_size=%10.6e    max_dmdt=%10.6e",
                        i, step_size, max_dmdt/dmdt_factor)
      end

      if save_m_every>0 && i%save_m_every == 0
          compute_system_energy(sim, sim.spin, t)
          write_data(sim)
      end
    if save_vtk_every > 0 && i%save_vtk_every == 0
        isdir(vtk_folder) || mkdir(vtk_folder)
        save_vtk_points(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, i)), fields = fields)
    end
    if save_ovf_every > 0 && i%save_ovf_every == 0
        isdir(ovf_folder) || mkdir(ovf_folder)
        save_ovf(sim, joinpath(ovf_folder, @sprintf("%s_%d", sim.name, i)), dataformat = ovf_format)
    end

    if sim.save_data
        if llg_driver 
            sim.saver.t = driver.ode.t
        end

        sim.saver.nsteps += 1
    end

    if max_dmdt < stopping_dmdt*dmdt_factor
      @info @sprintf("max_dmdt is less than stopping_dmdt=%g, Done!", stopping_dmdt)
      if save_m_every>0
          compute_system_energy(sim, sim.spin, t)
          write_data(sim)
      end
      if save_vtk_every > 0
        isdir(vtk_folder) || mkdir(vtk_folder)
        save_vtk_points(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, i)), fields = fields)
      end
      if save_ovf_every > 0
        isdir(ovf_folder) || mkdir(ovf_folder)
        save_ovf(sim, joinpath(ovf_folder, @sprintf("%s_%d", sim.name, i)), dataformat = ovf_format)
      end
      step = i
      break
    end
  end

  if step == maxsteps
      if save_vtk_every > 0
        isdir(vtk_folder) || mkdir(vtk_folder)
        save_vtk(sim, joinpath(vtk_folder, @sprintf("%s_%d", sim.name, step)), fields = fields)
      end
      if save_ovf_every > 0
        isdir(ovf_folder) || mkdir(ovf_folder)
        save_ovf(sim, joinpath(ovf_folder, @sprintf("%s_%d", sim.name, step)), dataformat = ovf_format)
      end
  end
  return nothing
end


function run_until(sim::AbstractSim, t_end::Float64, integrator::IntegratorCayley, save_data::Bool)
      if t_end < integrator.t - integrator.step
          println("Run_until: t_end >= integrator.t - integrator.step")
          return
      elseif t_end == integrator.t
          integrator.omega_t[:] = integrator.omega[:]
          omega_to_spin(integrator.omega_t, sim.prespin, sim.spin, sim.nxyz)
          sim.saver.t = t_end
          sim.saver.nsteps += 1
          if save_data
              compute_system_energy(sim, sim.spin, t_end)
              write_data(sim)
        end
          return
      elseif t_end > integrator.t - integrator.step && integrator.step > 0 && t_end < integrator.t
          interpolation_dopri5(integrator, t_end)
          omega_to_spin(integrator.omega_t, sim.prespin, sim.spin, sim.nxyz)
          sim.saver.t = t_end
          sim.saver.nsteps += 1
          if save_data
              compute_system_energy(sim, sim.spin, t_end)
              write_data(sim)
          end
          return
      end

      # so we have t_end > self.t
      if integrator.step_next<=0
          integrator.step_next = compute_init_step(sim, t_end - integrator.t)
      end

      while integrator.t < t_end
          ratio = (t_end - integrator.t)/integrator.step_next
          if ratio<1.2 && ratio>0.8
              integrator.step_next = t_end - integrator.t
          end

          advance_step(sim, integrator)
      end

      interpolation_dopri5(integrator, t_end)
      omega_to_spin(integrator.omega_t, sim.prespin, sim.spin, sim.nxyz)
      sim.saver.t = t_end
      sim.saver.nsteps += 1
      if save_data
          compute_system_energy(sim, sim.spin, t_end)
          write_data(sim)
      end
      return nothing
end


function run_until(sim::AbstractSim, t_end::Float64, integrator::Integrator, save_data::Bool)
    if t_end < integrator.t - integrator.step
        @info("Run_until: t_end >= integrator.t - integrator.step")
        return
    elseif t_end == integrator.t
        sim.saver.t = t_end
        sim.saver.nsteps += 1
        if save_data
            compute_system_energy(sim, sim.spin, t_end)
            write_data(sim)
        end
        return
    end

    # so we have t_end > self.t
    if integrator.step_next<=0
        integrator.step_next = compute_init_step_DP(sim, t_end - integrator.t)
    end

    while integrator.t < t_end
        if integrator.step_next + integrator.t> t_end
            integrator.step_next = t_end - integrator.t
        end
        advance_step(sim, integrator)
    end

    sim.saver.t = t_end
    sim.saver.nsteps += 1
    if save_data
        compute_system_energy(sim, sim.spin, t_end)
        write_data(sim)
    end
end

function run_until(sim::AbstractSim, t_end::Float64; save_data=true)
    run_until(sim, t_end, sim.driver.ode, save_data)
end

"""
    create_sim(mesh; args...)

Create a micromagnetic simulation instance with given arguments. 

- `mesh`: a mesh has to be provided to start the simulation. The mesh could be [`FDMesh`](@ref), [`FEMesh`](@ref),
[`FDMeshGPU`](@ref), [`CubicMeshGPU`](@ref), or [`TriangularMeshGPU`](@ref).

# Arguments
- `name` : the simulation name, should be a string.
- `driver` : the driver name, should be a string. By default, the driver is "SD".
- `alpha` : the Gilbert damping in the LLG equation, should be a number.
- `beta` : the nonadiabatic strength in the LLG equation with spin transfer torques (zhang-li model), should be a number.
- `gamma` : the gyromagnetic ratio, default value = 2.21e5.
- `ux`, `uy` or `uz`: the strengths of the spin transfer torque.
- `Ms`: the saturation magnetization, should be [`NumberOrArrayOrFunction`](@ref). By default, Ms=8e5
- `mu_s`: the magnetic moment, should be [`NumberOrArrayOrFunction`](@ref). By default, mu_s=2*mu_B
- `A` or `J`: the exchange constant, should be [`NumberOrArrayOrFunction`](@ref).
- `D` : the DMI constant, should be [`NumberOrArrayOrFunction`](@ref).
- `dmi_type` : the type of DMI, could be "bulk" or "interfacial".
- `Ku`: the anisotropy constant, should be [`NumberOrArrayOrFunction`](@ref).
- `axis`: the anisotropy axis, should be a tuple, such as (0,0, 1)
- `demag` : include demagnetization or not, should be a boolean, i.e., true or false.
- `H`: the external field, should be a tuple or function, i.e., [`TupleOrArrayOrFunction`](@ref). 
- `m0` : the initial magnetization, should be a tuple or function, i.e., [`TupleOrArrayOrFunction`](@ref). 
"""
function create_sim(mesh; args...)
    #convert args to a dict
    args = Dict(args)
    #Ms=8e5, A=1.3e-11, D=0, Ku=0, axis=(0,0,1), H=(0,0,0), m0=(0,0,1), demag=false, driver="SD", name="relax"
    driver = haskey(args, :driver) ? args[:driver] : "SD"
    name = haskey(args, :name) ? args[:name] : "unnamed"

    #Create the mesh using given driver and name
    sim = Sim(mesh, driver=driver, name = name)

    #If the simulation is the standard micromagnetic simulation.
    if isa(mesh, FDMesh) || isa(mesh, FDMeshGPU) || isa(mesh, FEMesh) || isa(mesh, FDMeshGPU)

        # we set the Ms anyway
        Ms = haskey(args, :Ms) ? args[:Ms] : 8e5
        set_Ms(sim, Ms)

        # add the exchange if A is given
        if haskey(args, :A)
            add_exch(sim, args[:A])
        end

        # add the DMI if D is given
        if haskey(args, :D)
            dmi_type = haskey(args, :dmi_type) ? args[:dmi_type] : "bulk"
            add_dmi(sim, args[:D], type=dmi_type)
        end

        # add the demag
        if haskey(args, :demag) && args[:demag]
            add_demag(sim)
        end

        for key in [:Ms, :A, :D, :demag]
            haskey(args, key) && delete!(args, key)
        end

    
    #If the simulation is atomistic
    elseif isa(mesh, AtomicMeshGPU)

        mu_s = haskey(args, :mu_s) ? args[:mu_s] : 2*mu_B
        set_mu_s(sim, mu_s)

        # add the exchange if A is given
        if haskey(args, :J)
            add_exch(sim, args[:J])
        end
        
        # add the DMI if D is given
        if haskey(args, :D)
            add_dmi(sim, args[:D])
        end

        for key in [:mu_s, :J, :D]
            haskey(args, key) && delete!(args, key)
        end
    else
        error("This info is for debug.")
    end

    # add the anisotropy
    if haskey(args, :Ku)
        axis = haskey(args, :axis) ? args[:axis] : (0,0,1)
        add_anis(sim, args[:Ku], axis=axis)
        haskey(args, :axis) && delete!(args, :axis)
    end

    # add the external field
    if haskey(args, :H)
        add_zeeman(sim, args[:H])
    end

    if haskey(args, :alpha) && startswith(driver, "LLG")
        sim.driver.alpha = args[:alpha]
        delete!(args, :alpha)
    end

    if haskey(args, :gamma) && startswith(driver, "LLG")
        sim.driver.gamma = args[:gamma]
        delete!(args, :gamma)
    end

    if haskey(args, :beta) && startswith(driver, "LLG_STT")
        sim.driver.beta = args[:beta]
        delete!(args, :beta)
    end

    if haskey(args, :ux) && startswith(driver, "LLG_STT")
        set_ux(sim, args[:ux])
        delete!(args, :ux)
    end

    if haskey(args, :uy) && startswith(driver, "LLG_STT")
        set_uy(sim, args[:uy])
        delete!(args, :uy)
    end

    if haskey(args, :uz) && startswith(driver, "LLG_STT")
        set_uz(sim, args[:uz])
        delete!(args, :uz)
    end

    # set m0 anyway
    m0_value = haskey(args, :m0) ? args[:m0] : (0.8, 0.6, 0)
    init_m0(sim, m0_value)
  
    for key in [:driver, :name, :Ku, :H, :m0]
        haskey(args, key) && delete!(args, key)
    end
    
    for key in args
        @warn @sprintf("Key '%s' is not used.", key)
    end
    
    return sim
    
  end

