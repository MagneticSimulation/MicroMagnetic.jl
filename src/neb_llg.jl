using Printf

mutable struct NEB_LLG <: AbstractNEB
  N::Int64  #number of images
  nxyz::Int64 # nxyz = neb.sim.nxyz*(neb.N)
  sim::AbstractSim
  init_m::Any
  intervals::Any
  driver::NEBDriver
  images::Array{Float64, 2}  #size(images) = (3*nxyz, N)
  pre_images::Array{Float64, 2}
  spin::Array{Float64, 1} #A pointer to images
  prespin::Array{Float64, 1} #A pointer to pre_images
  name::String
  field::Array{Float64, 2}
  energy::Array{Float64, 2}
  saver::NEBSaver
  k::Float64
end

mutable struct NEB_LLG_Driver <: NEBDriver
  ode::Integrator
end

function NEB_LLG(sim::AbstractSim, init_m::Any, intervals::Any; name="NEB", k=1.0e5, driver="NEB_LLG")
  nxyz = sim.mesh.nx*sim.mesh.ny*sim.mesh.nz
  N = sum(intervals)+length(intervals)+1
  images = zeros(Float64,3*sim.nxyz,N)
  pre_images = zeros(Float64,3*sim.nxyz,N)
  energy = zeros(Float64,nxyz,N)
  Ms = zeros(Float64,nxyz)
  driver = create_neb_driver(driver, nxyz, N)
  field = zeros(Float64,3*nxyz, N)
  energy = zeros(nxyz,N)

  headers = ["steps"]
  units = ["<>"]
  results = Any[o::NEB_LLG -> o.driver.ode.nsteps]
  saver = NEBSaver(string(name, ".txt"), 0, false, headers, units, results)
  #interactions = []
  spins = reshape(images, 3*sim.nxyz*N)
  prespins = reshape(pre_images, 3*sim.nxyz*N)
  neb = NEB_LLG(N, sim.nxyz*N, sim, init_m, intervals, driver, images, pre_images, spins, prespins,
             name, field,energy,saver,k)
  init_images(neb)
  return neb
end

function init_images(neb::NEB_LLG)
  sim=neb.sim
  nxyz=sim.nxyz
  N=neb.N
  m=zeros(3*nxyz,N)
  pics= neb.init_m
  intervals=neb.intervals
  if length(pics)==length(intervals)+1
    n=1
    for i=1:(length(pics)-1)
      m1=zeros(3*nxyz)
      m2=zeros(3*nxyz)
      init_vector!(m1, sim.mesh, pics[i])
      init_vector!(m2, sim.mesh, pics[i+1])
      normalise(m1,nxyz)
      normalise(m2,nxyz)
      M=interpolate_m(m1,m2,Int(intervals[i]))
      for j=1:intervals[i]+1
          m[:,n]=M[:,j]
          n+=1
      end
    end
    m3=zeros(3*nxyz)
    init_vector!(m3, sim.mesh, pics[length(pics)])
    normalise(m3,nxyz)
    m[:,N]=m3[:]
  else
    println("Input error!")
  end
  neb.pre_images .= m
  neb.images .= neb.pre_images
  for n = 1:N
    name = @sprintf("E_total_%g",n)
    push!(neb.saver.headers,name)
    push!(neb.saver.units, "J")
    fun =  o::NEB_LLG -> sum(o.energy[:,n])
    push!(neb.saver.results, fun)
  end
end


function effective_field_NEB(neb::NEB_LLG, spin::Array{Float64, 1}, t::Float64)

  sim = neb.sim
  nxyz = sim.nxyz
  #neb.spin[:] = spin[:], we already copy spin to neb.images in dopri5
  images = neb.images
  N = neb.N
  fill!(neb.field, 0.0)
  fill!(neb.energy, 0.0)
  for n = 1:neb.N
      for interaction in sim.interactions
          effective_field(interaction, sim, images[:,n], 0.0)
          neb.field[:,n] .+= interaction.field[:]
          if (n==1)||(n==neb.N)
            neb.field[:,n] .= 0
          end
          neb.energy[:,n] .+= interaction.energy[:]
  end
  end
  t = zeros(3*nxyz,N)
  for n = 1:neb.N
    if (n==1)||(n==neb.N)
      continue
    end
    E1 = sum(neb.energy[:,n-1])
    E2 = sum(neb.energy[:,n])
    E3 = sum(neb.energy[:,n+1])
    dEmax = max(abs(E3-E2),abs(E2-E1))
    dEmin = min(abs(E3-E2),abs(E2-E1))
    tip = images[:,n+1]-images[:,n]
    tim = images[:,n]-images[:,n-1]
    if (E1>E2)&&(E2>E3)
      t[:,n] = tim
    elseif (E3>E2)&&(E2>E1)
      t[:,n] = tip
    elseif E3>E1
      t[:,n] =dEmax*tip+dEmin*tim
    else
      t[:,n] =dEmin*tip+dEmax*tim
    end
    for i = 1:nxyz
      j = 3*i - 2
      fx,fy,fz = cross_product(images[j,n],images[j+1,n],images[j+2,n], t[j,n],t[j+1,n],t[j+2,n])
      t[j,n],t[j+1,n],t[j+2,n] = cross_product(fx,fy,fz, images[j,n],images[j+1,n],images[j+2,n])
    end
    norm_t = model(t[:,n])
    t[:,n] = t[:,n]/norm_t
    f = neb.field[:,n]'*t[:,n]
    neb.field[:,n] .-= f*t[:,n]
    neb.field[:,n] .+= neb.k*(model(tip)-model(tim))*t[:,n]
  end
  return 0
end


function neb_llg_call_back(neb::AbstractNEB, dm_dt::Array{Float64, 1}, spin::Array{Float64, 1}, t::Float64)
  effective_field_NEB(neb, spin, t)
  neb.field[:, 1] .= 0
  neb.field[:, neb.N] .= 0
  field = reshape(neb.field, 3*neb.nxyz)
  llg_rhs(dm_dt, spin, field, 1.0, 2.21e5, false, neb.nxyz)
  return nothing
end


function relax_NEB_LLG(neb::NEB_LLG, maxsteps::Int64, stopping_dmdt::Float64, save_m_every::Int64, save_vtk_every::Int64, vtk_folder::String,save_ovf_every::Int64,ovf_folder::String)
    N = neb.N
    sim = neb.sim

    dmdt_factor = (2 * pi / 360) * 1e9

    driver = neb.driver

    if save_m_every>0
        compute_system_energy(neb,  0.0)
        write_data(neb)
    end
    if save_vtk_every > 0
        save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, 0)))
    end
    if save_ovf_every > 0
        save_ovf(neb, joinpath(ovf_folder, @sprintf("%s_%d", neb.name, 0)))
    end

    rk_data = neb.driver.ode
    for i=1:maxsteps

        advance_step(neb, rk_data)

        step_size = rk_data.step
        max_dmdt = compute_dmdt(neb.prespin, neb.spin, neb.nxyz, step_size)
        @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e",
                     i, rk_data.step, rk_data.t, max_dmdt/dmdt_factor)

        if save_m_every>0 && i%save_m_every == 0
            compute_system_energy(neb,  0.0)
            write_data(neb)
        end
        if save_vtk_every > 0 && i%save_vtk_every == 0
            save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, i)))
        end
        if save_ovf_every > 0 && i%save_ovf_every == 0
            save_ovf(neb, joinpath(ovf_folder, @sprintf("%s_%d", neb.name, i)))
        end

        if max_dmdt < stopping_dmdt*dmdt_factor
            @info @sprintf("max_dmdt is less than stopping_dmdt=%g, Done!", stopping_dmdt)

            if save_m_every>0
                compute_system_energy(neb, 0.0)
                write_data(neb)
            end
            if save_vtk_every > 0
                save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, i)))
            end
            if save_ovf_every > 0
                save_ovf(neb, joinpath(ovf_folder, @sprintf("%s_%d", neb.name, i)))
            end
            break
        end
    end
    return nothing
end
