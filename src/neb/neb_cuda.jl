using LinearAlgebra
using CUDA

mutable struct NEB_GPU{T<:AbstractFloat} <: AbstractSim
    N::Int64 #number of images in this process
    nxyz::Int64 # nxyz = neb.sim.nxyz*N
    clib_image:: Int64
    sim::AbstractSim
    driver::NEBDriver
    image_l::CuArray{T, 1}
    spin::CuArray{T, 1}
    prespin::CuArray{T, 1}
    image_r::CuArray{T, 1}
    field::CuArray{T, 1}
    tangent::CuArray{T, 1}
    energy::CuArray{T, 1}
    energy_cpu::Array{T, 1}
    distance::Array{T, 1}
    saver_energy::DataSaver
    saver_distance::DataSaver
    spring_constant::Float64
    name::String
    NEB_GPU{T}() where {T<:AbstractFloat} = new()
end

function NEB_GPU(sim::AbstractSim, init_m::Any, intervals::Any; name="NEB", spring_constant=1.0e5, driver="LLG", clib_image=-1)
    Float = _cuda_using_double.x ? Float64 : Float32
    neb = NEB_GPU{Float}()

    mesh = sim.mesh
    nxyz = sim.mesh.nx*sim.mesh.ny*sim.mesh.nz
    N = sum(intervals) + length(intervals) - 1
    if N < 1
        error("The number of free images should larger than 1 (typically 10~20).")
    end
    neb.N = N
    neb.nxyz = sim.nxyz*neb.N
    neb.sim = sim
    neb.name = name
    neb.clib_image = clib_image

    neb.driver = create_neb_driver_gpu(driver, nxyz,N)

    neb.image_l = CUDA.zeros(Float, 3*sim.nxyz)
    neb.spin = CUDA.zeros(Float, 3*sim.nxyz*neb.N)
    neb.prespin = CUDA.zeros(Float, 3*sim.nxyz*neb.N)
    neb.image_r = CUDA.zeros(Float, 3*sim.nxyz)
    neb.field = CUDA.zeros(Float, 3*sim.nxyz*neb.N)
    neb.tangent = CUDA.zeros(Float, 3*sim.nxyz*neb.N)
    neb.energy = CUDA.zeros(Float, neb.N+2)
    neb.energy_cpu = zeros(Float, neb.N+2)
    neb.distance = zeros(Float, neb.N+1)

    step_item = SaverItem("step", "", o::NEB -> o.driver.nsteps)
    init_saver(neb)
    init_images(neb, init_m, intervals)
    compute_distance(neb)
    compute_system_energy(neb)
    compute_system_energy_two_ends(neb)

    return neb
end


function create_neb_driver_gpu(driver::String, nxyz::Int64, N::Int64)
    if driver == "LLG"
        tol = 1e-5
        dopri5 = DormandPrinceGPU(nxyz*N, neb_llg_call_back_gpu, tol)
        return NEB_LLG_Driver(0, dopri5)
    else
      error("Only Driver 'LLG' is supported!")
    end
  end

function init_saver(neb::NEB_GPU)
    step_item = SaverItem("step", "", o::NEB_GPU -> o.saver_energy.nsteps)
    neb.saver_energy = DataSaver(string(neb.name, ".txt"), false, 0.0, 0, [step_item])
    neb.saver_distance = DataSaver(string(@sprintf("%s_distance",neb.name), ".txt"), false, 0.0, 0, [step_item])
    N=neb.N
    for n = 1:N+2
        name = @sprintf("E_total_%g",n)
        item = SaverItem(name, "J", o::NEB_GPU -> o.energy_cpu[n])
        push!(neb.saver_energy.items, item)
    end

    for n = 1:N+1
        name =  @sprintf("distance_%d", n)
        item = SaverItem(name, "", o::NEB_GPU -> o.distance[n])
        push!(neb.saver_distance.items, item)
    end
end

function init_m0_Ms(sim::AbstractSim, m0::TupleOrArrayOrFunction)
  Float = _cuda_using_double.x ? Float64 : Float32
  spin = zeros(Float, 3*sim.nxyz)
  init_vector!(spin, sim.mesh, m0)
  normalise(spin, sim.nxyz)
  Ms = Array(sim.Ms)
  for i = 1:sim.nxyz
      if abs(Ms[i]) < eps(Float)
          spin[3*i-2] = 0
          spin[3*i-1] = 0
          spin[3*i] = 0
      end
  end
  return spin
end

function init_images(neb::NEB_GPU, init_m::Any, intervals::Any)
    sim = neb.sim
    nxyz = sim.nxyz

    if length(init_m) != length(intervals)+1
        error("Imput init_m or intervals wrong!")
    end

    dof = 3*nxyz
    image_id = 1
    local_image_id=1
    N = neb.N

    for id = 1:length(intervals)
        m1 = init_m0_Ms(sim, init_m[id])
        m2 = init_m0_Ms(sim, init_m[id+1])
        M = interpolate_m(m1, m2, Int(intervals[id]))
        #M = interpolate_m_spherical(m1,m2,Int(intervals[i]))
        if  image_id == 1
            CUDA.@sync copyto!(neb.image_l, M[:, 1])
        end
        for j=1:intervals[id]+1
            #images[:,n]=M[:,j]
            if  2 <= image_id <= N+1
                k = local_image_id-1
                b = view(neb.spin, k*dof+1:(k+1)*dof)
                CUDA.@sync copyto!(b, M[:, j])
                local_image_id += 1
            end
            image_id += 1
        end
    end
    m = init_m0_Ms(sim, init_m[end])
    CUDA.@sync copyto!(neb.image_r, m)

    CUDA.@sync copyto!(neb.prespin, neb.spin)
    return nothing
end

function effective_field_NEB(neb::NEB_GPU, spin::CuArray{T, 1})  where {T<:AbstractFloat}

    #neb.spin[:] = spin[:], we already copy spin to neb.spin in dopri5
    compute_micromagnetic_field(neb)

    compute_field_related_to_tangent(neb)
    return 0
  end




function compute_micromagnetic_field(neb::NEB_GPU)
    CUDA.@sync fill!(neb.field, 0.0)
    sim = neb.sim
    dof = 3*sim.nxyz
    for n = 1:neb.N
        f = view(neb.field, (n-1)*dof+1:n*dof)
        m = view(neb.spin, (n-1)*dof+1:n*dof)
        total_energy = 0
        for interaction in sim.interactions
            effective_field(interaction, sim, m, 0.0)
            f .+= sim.field
            total_energy += sum(sim.energy)
        end
        neb.energy_cpu[n+1] = total_energy
    end

    compute_system_energy_two_ends(neb)
    CUDA.@sync copyto!(neb.energy, neb.energy_cpu)
   return nothing
end

function compute_field_related_to_tangent(neb::NEB_GPU)
   compute_tangents(neb)
   compute_distance(neb)

   dof = 3*neb.sim.nxyz
   for n = 1:neb.N
       f = view(neb.field, (n-1)*dof+1:n*dof)
       t = view(neb.tangent, (n-1)*dof+1:n*dof)
      
       if n == neb.clib_image
        ft =  - 2*LinearAlgebra.dot(f, t)
       else
        ft = neb.spring_constant*(neb.distance[n+1]-neb.distance[n]) - LinearAlgebra.dot(f, t)
       end
       f .+=  ft.*t

   end

   return nothing
end

#let compute the total energy for the two ends directly for now.
#for larger systems it should be better to use MPI for sending/recieving
#energys.
function compute_system_energy_two_ends(neb::NEB_GPU)
  sim = neb.sim
  #left image
  total_energy = 0
  for interaction in sim.interactions
      effective_field(interaction, sim, neb.image_l, 0.0)
      total_energy +=  sum(sim.energy)
  end
  neb.energy_cpu[1] = total_energy

  #right image
  total_energy = 0
  for interaction in sim.interactions
      effective_field(interaction, sim, neb.image_r, 0.0)
      total_energy +=  sum(sim.energy)
  end
  neb.energy_cpu[end] = total_energy

  return 0
end

function compute_system_energy(neb::NEB_GPU)
  sim = neb.sim

  dof = 3*sim.nxyz
  for n = 1:neb.N
      sim.total_energy = 0
      b = view(neb.spin, (n-1)*dof+1:n*dof)
      for interaction in sim.interactions
          effective_field(interaction, sim, b, 0.0)
          interaction.total_energy = sum(sim.energy)
          sim.total_energy += interaction.total_energy
      end
      neb.energy_cpu[n+1] = sim.total_energy
  end
  return 0
end

function neb_llg_call_back_gpu(neb::NEB_GPU, dm_dt::CuArray{T, 1}, spin::CuArray{T, 1}, t::Float64)  where {T<:AbstractFloat}
    effective_field_NEB(neb, spin)
    neb_llg_rhs_gpu(dm_dt, spin, neb.field, neb.sim.pins, 2.21e5, neb.nxyz, neb.sim.nxyz)
    return nothing
  end


function relax(neb::NEB_GPU; maxsteps=10000, stopping_dmdt=0.05, save_m_every = 10, save_ovf_every=-1, ovf_format = "binary", ovf_folder="ovfs", save_vtk_every=-1, vtk_folder="vtks")
    if save_vtk_every>0 && !isdir(vtk_folder)
       mkdir(vtk_folder)
    end
    if save_ovf_every>0 && !isdir(ovf_folder)
       mkdir(ovf_folder)
    end
    N = neb.N
    sim = neb.sim

    dmdt_factor = (2 * pi / 360) * 1e9

    driver = neb.driver

    if save_m_every>0
        compute_system_energy(neb)
        compute_distance(neb)
        write_data(neb, neb.saver_energy)
        write_data(neb, neb.saver_distance)
    end
    if save_vtk_every > 0
        save_vtk(neb, joinpath(vtk_folder, neb.name))
    end
    if save_ovf_every > 0
        save_ovf(neb, joinpath(ovf_folder, neb.name))
    end

    rk_data = neb.driver.ode
    for i=1:maxsteps
        if !advance_step(neb, rk_data)
            @info("Advance step failed, end.")
            break
        end
        #advance_step(neb, rk_data)

        neb.saver_energy.nsteps += 1

        step_size = rk_data.step

        #we borrow neb.field to store dm
        compute_dm!(neb.field, neb.prespin, neb.spin, neb.nxyz)
        max_dmdt = maximum(view(neb.field, 1:neb.nxyz))/step_size

        @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e  time=%s",
                            i, rk_data.step, rk_data.t, max_dmdt/dmdt_factor, Dates.now())

        if save_m_every>0 && i%save_m_every == 0
            compute_system_energy(neb)
            compute_distance(neb)
            write_data(neb, neb.saver_energy)
            write_data(neb, neb.saver_distance)
        end
        if save_vtk_every > 0 && i%save_vtk_every == 0
            save_vtk(neb, joinpath(vtk_folder, neb.name))
        end
        if save_ovf_every > 0 && i%save_ovf_every == 0
            save_ovf(neb, joinpath(ovf_folder,  neb.name))
        end

        if max_dmdt < stopping_dmdt*dmdt_factor
            @info @sprintf("max_dmdt is less than stopping_dmdt=%g, Done!", stopping_dmdt)

            if save_m_every>0
                compute_system_energy(neb)
                write_data(neb, neb.saver_energy)
                write_data(neb, neb.saver_distance)
            end
            if save_vtk_every > 0
                save_vtk(neb, joinpath(vtk_folder, neb.name))
            end
            if save_ovf_every > 0
                save_ovf(neb, joinpath(vtk_folder, neb.name))
            end
            break
        end
    end
    return nothing
end

#TODO: merge save_ovf and save_vtk???
function save_ovf(neb::NEB_GPU, fname::String; 
    type::DataType = _cuda_using_double.x ? Float64 : Float32)

  sim = neb.sim
  dof = 3*sim.nxyz

  id_start = 0
  sim.spin .= neb.image_l
  save_ovf(sim, @sprintf("%s_0",fname); type = type)

  for n = 1:neb.N
      b = view(neb.spin, (n-1)*dof+1:n*dof)
      sim.spin .= b
      save_ovf(sim, @sprintf("%s_%d",fname,n); type = type)
  end
  
    sim.spin .= neb.image_r
    save_ovf(sim, @sprintf("%s_%d",fname,neb.N+1); type = type)
end

function save_vtk(neb::NEB_GPU, fname::String)
  sim = neb.sim
  dof = 3*sim.nxyz

    sim.spin .= neb.image_l
    save_vtk(sim,  @sprintf("%s_0",fname))
  for n = 1:neb.N
      b = view(neb.spin, (n-1)*dof+1:n*dof)
      sim.spin .= b
      save_vtk(sim, @sprintf("%s_%d",fname,n))
  end
    sim.spin .= neb.image_r
    save_vtk(sim, @sprintf("%s_%d",fname,neb.N+1))
end
