using MPI
using LinearAlgebra
using CuArrays

mutable struct NEB_GPU_MPI{T<:AbstractFloat} <: AbstractSim
    comm_rank::Int64
    comm_size::Int64
    N_total::Int64  #number of total images
    N::Int64 #number of images in this process
    nxyz::Int64 # nxyz = neb.sim.nxyz*N
    sim::AbstractSim
    driver::NEBDriver
    image_l::CuArray{T, 1}
    image_lc::Array{T, 1} #cuda-awared MPI seems do not work well for now.
    spin::CuArray{T, 1}
    prespin::CuArray{T, 1}
    image_r::CuArray{T, 1}
    image_rc::Array{T, 1}
    field::CuArray{T, 1}
    tangent::CuArray{T, 1}
    energy::CuArray{T, 1}
    energy_cpu::Array{T, 1}
    distance::Array{T, 1}
    saver_energy::DataSaver
    saver_distance::DataSaver
    spring_constant::Float64
    name::String
    NEB_GPU_MPI{T}() where {T<:AbstractFloat} = new()
end

function NEB_MPI(sim::AbstractSim, init_m::Any, intervals::Any; name="NEB", spring_constant=1.0e5, driver="LLG")
    MPI.Init()
    comm = MPI.COMM_WORLD

    Float = _cuda_using_double.x ? Float64 : Float32
    neb = NEB_GPU_MPI{Float}()
    comm_size = MPI.Comm_size(comm)
    comm_rank = MPI.Comm_rank(comm)
    neb.comm_rank = comm_rank
    neb.comm_size = comm_size

    mesh = sim.mesh
    nxyz = sim.mesh.nx*sim.mesh.ny*sim.mesh.nz
    N = sum(intervals) + length(intervals) - 1
    if N < 1
        error("The number of free images should larger than 1 (typically 10~20).")
    end
    neb.N_total = N
    neb.N = ceil(Int64, N/comm_size)
    if comm_rank == comm_size - 1
        neb.N = N - comm_rank*ceil(Int64, N/comm_size)
    end
    neb.nxyz = sim.nxyz*neb.N
    neb.sim = sim
    neb.name = name

    neb.driver = create_neb_driver_mpi(driver, nxyz, N)

    neb.image_l = CuArrays.zeros(Float, 3*sim.nxyz)
    neb.image_lc = zeros(Float, 3*sim.nxyz)
    neb.spin = CuArrays.zeros(Float, 3*sim.nxyz*neb.N)
    neb.prespin = CuArrays.zeros(Float, 3*sim.nxyz*neb.N)
    neb.image_r = CuArrays.zeros(Float, 3*sim.nxyz)
    neb.image_rc = zeros(Float, 3*sim.nxyz)
    neb.field = CuArrays.zeros(Float, 3*sim.nxyz*neb.N)
    neb.tangent = CuArrays.zeros(Float, 3*sim.nxyz*neb.N)
    neb.energy = CuArrays.zeros(Float, neb.N+2)
    neb.energy_cpu = zeros(Float, neb.N+2)
    neb.distance = zeros(Float, neb.N+1)

    init_saver(neb)
    init_images(neb, init_m, intervals)

    update_images_between_processes(neb)
    compute_system_energy_two_ends(neb)

    return neb
end

function init_saver(neb::NEB_GPU_MPI)

    headers = ["steps"]
    units = ["<>"]
    results = Any[o::NEB_GPU_MPI -> o.driver.nsteps]
    neb.saver_energy = DataSaver(@sprintf("%s_energy_%d.txt", neb.name, neb.comm_rank),
                                 0.0, 0, false, headers, units, results)

    distance_headers= ["steps"]
    distance_units = ["<>"]
    distance_results = Any[o::NEB_GPU_MPI -> o.driver.nsteps]

    neb.saver_distance = DataSaver(@sprintf("%s_distance_%d.txt", neb.name, neb.comm_rank),
                                   0.0, 0, false, distance_headers, distance_units, distance_results)
    for n = 1:neb.N
        name = @sprintf("E_total_%g",n)
        push!(neb.saver_energy.headers,name)
        push!(neb.saver_energy.units, "<J>")
        fun =  o::NEB_GPU_MPI -> o.energy_cpu[n]
        push!(neb.saver_energy.results, fun)
    end

    for n = 1:neb.N-1
        name =  @sprintf("distance_%d",n)
        push!(neb.saver_distance.headers,name)
        push!(neb.saver_distance.units, "<>")
        fun =  o::NEB_GPU_MPI -> o.distance[n]
        push!(neb.saver_distance.results, fun)
    end

end

function create_neb_driver_mpi(driver::String, nxyz::Int64, N::Int64)
  if driver == "LLG"
      tol = 1e-5
      dopri5 = DormandPrinceGPU(nxyz*N, neb_llg_call_back_gpu_mpi, tol)
      return NEB_LLG_Driver(0, dopri5)
  else
    error("Only Driver 'LLG' is supported!")
  end
end

function init_images(neb::NEB_GPU_MPI, init_m::Any, intervals::Any)
    sim = neb.sim
    nxyz = sim.nxyz

    if length(init_m) != length(intervals)+1
        error("Imput init_m or intervals wrong!")
    end

    dof = 3*nxyz
    rank = neb.comm_rank
    size = neb.comm_size

    image_id = 1
    local_image_id = 1

    start, stop = 0, 0
    if rank == size - 1
        stop = neb.N_total
        start = stop - neb.N + 1
    else
        start = rank*neb.N + 1
        stop = (rank+1)*neb.N
    end

    for id = 1:length(intervals)
          m1 = zeros(3*nxyz)
          m2 = zeros(3*nxyz)
          init_vector!(m1, sim.mesh, init_m[id])
          init_vector!(m2, sim.mesh, init_m[id+1])
          normalise(m1, nxyz)
          normalise(m2, nxyz)
          M = interpolate_m(m1, m2, Int(intervals[id]))
          #M = interpolate_m_spherical(m1,m2,Int(intervals[i]))
          if  rank == 0 && image_id == 1
              CuArrays.@sync copyto!(neb.image_l, M[:, 1])
          end
          for j=1:intervals[id]+1
              #images[:,n]=M[:,j]
               if  start+1 <= image_id <= stop+1
                   k = local_image_id-1
                   b = view(neb.spin, k*dof+1:(k+1)*dof)
                   CuArrays.@sync copyto!(b, M[:, j])
                   local_image_id += 1
               end
              image_id += 1
          end
    end

    if rank == size - 1
        m = zeros(3*nxyz)
        init_vector!(m, sim.mesh, init_m[end])
        normalise(m, nxyz)
        CuArrays.@sync copyto!(neb.image_r, m)
    end

    CuArrays.@sync copyto!(neb.prespin, neb.spin)

    println("We are here!!!")
    return nothing
end


function effective_field_NEB(neb::NEB_GPU_MPI, spin::CuArray{T, 1})  where {T<:AbstractFloat}

  #neb.spin[:] = spin[:], we already copy spin to neb.spin in dopri5

  update_images_between_processes(neb)

  compute_micromagnetic_field(neb)

  compute_field_related_to_tangent(neb)
  MPI.Barrier(MPI.COMM_WORLD)
  return 0
end

function update_images_between_processes(neb::NEB_GPU_MPI)
    nxyz = neb.sim.nxyz
    dof = 3*nxyz
    N = neb.N
    rank = neb.comm_rank
    size = neb.comm_size
    left_id = rank - 1
    right_id = rank < size - 1 ? rank+1 : -1

    MPI.Barrier(MPI.COMM_WORLD)

    if rank > 0
        b = view(neb.spin, 1:dof)
        CuArrays.@sync copyto!(neb.image_lc, b)
        sreq = MPI.Isend(neb.image_lc, left_id, rank+64, MPI.COMM_WORLD)
        #MPI.Wait!(sreq)
    end

    if rank < size - 1
        rreq = MPI.Irecv!(neb.image_rc, right_id, right_id+64, MPI.COMM_WORLD)
        MPI.Wait!(rreq)
        CuArrays.@sync copyto!(neb.image_r, neb.image_rc)
    end

    MPI.Barrier(MPI.COMM_WORLD)

    if rank < size - 1
        b = view(neb.spin, (N-1)*dof+1:N*dof)
        CuArrays.@sync copyto!(neb.image_rc, b)
        sreq = MPI.Isend(neb.image_rc, right_id, rank+128, MPI.COMM_WORLD)
        #MPI.Wait!(sreq)
    end

    if rank > 0
        rreq = MPI.Irecv!(neb.image_lc, left_id, left_id+128, MPI.COMM_WORLD)
        MPI.Wait!(rreq)
        CuArrays.@sync copyto!(neb.image_l, neb.image_lc)
    end

    MPI.Barrier(MPI.COMM_WORLD)
    return nothing
end

function compute_micromagnetic_field(neb::NEB_GPU_MPI)
    CuArrays.@sync fill!(neb.field, 0.0)
    sim = neb.sim
    dof = 3*sim.nxyz
    for n = 1:neb.N
        f = view(neb.field, (n-1)*dof+1:n*dof)
        m = view(neb.spin, (n-1)*dof+1:n*dof)
        sim.total_energy = 0
        for interaction in sim.interactions
            effective_field(interaction, sim, m, 0.0)
            f .+= sim.field
            interaction.total_energy = sum(sim.energy)
            sim.total_energy += interaction.total_energy
        end
        neb.energy_cpu[n+1] = sim.total_energy
    end

    compute_system_energy_two_ends(neb)
    CuArrays.@sync copyto!(neb.energy, neb.energy_cpu)
   return nothing
end

function compute_field_related_to_tangent(neb::NEB_GPU_MPI)
   compute_tangents(neb)
   compute_distance(neb)

   dof = 3*neb.sim.nxyz
   for n = 1:neb.N
       f = view(neb.field, (n-1)*dof+1:n*dof)
       t = view(neb.tangent, (n-1)*dof+1:n*dof)
       ft = neb.spring_constant*(neb.distance[n+1]-neb.distance[n]) - LinearAlgebra.dot(f, t)
       f .+=  ft.*t
   end

   return nothing
end

#let compute the total energy for the two ends directly for now.
#for larger systems it should be better to use MPI for sending/recieving
#energys.
function compute_system_energy_two_ends(neb::NEB_GPU_MPI)
  sim = neb.sim
  #left image
  sim.total_energy = 0
  for interaction in sim.interactions
      effective_field(interaction, sim, neb.image_l, 0.0)
      interaction.total_energy = sum(sim.energy)
      sim.total_energy += interaction.total_energy
  end
  neb.energy_cpu[1] = sim.total_energy

  #right image
  sim.total_energy = 0
  for interaction in sim.interactions
      effective_field(interaction, sim, neb.image_r, 0.0)
      interaction.total_energy = sum(sim.energy)
      sim.total_energy += interaction.total_energy
  end
  neb.energy_cpu[end] = sim.total_energy

  return 0
end

function compute_system_energy(neb::NEB_GPU_MPI)
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


function neb_llg_call_back_gpu_mpi(neb::NEB_GPU_MPI, dm_dt::CuArray{T, 1}, spin::CuArray{T, 1}, t::Float64)  where {T<:AbstractFloat}
  effective_field_NEB(neb, spin)
  #println("spin: ", neb.spin)
  #println("fields: ", neb.field)
  neb_llg_rhs_gpu(dm_dt, spin, neb.field, 2.21e5, neb.nxyz)
  return nothing
end


function relax(neb::NEB_GPU_MPI; maxsteps=10000, stopping_dmdt=0.05, save_m_every = 10, save_ovf_every=-1, ovf_format = "binary", ovf_folder="ovfs", save_vtk_every=-1, vtk_folder="vtks")
    if save_vtk_every>0 && !isdir(vtk_folder) && neb.comm_rank == 0
        mkdir(vtk_folder)
    end
    if save_ovf_every>0 && !isdir(ovf_folder) && neb.comm_rank == 0
      mkdir(ovf_folder)
    end
    MPI.Barrier(MPI.COMM_WORLD) # make sure that the two folders are generated.

    N = neb.N
    sim = neb.sim

    dmdt_factor = (2 * pi / 360) * 1e9

    driver = neb.driver

    if save_m_every>0
        compute_system_energy(neb)
        write_data(neb, neb.saver_energy)
        write_data(neb, neb.saver_distance)
    end
    if save_vtk_every > 0
        save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, 0)))
    end
    if save_ovf_every > 0
        save_ovf(neb, joinpath(ovf_folder, @sprintf("%s_%d", neb.name, 0)))
    end

    MPI.Barrier(MPI.COMM_WORLD)
    rk_data = neb.driver.ode
    for i=1:maxsteps

        if !advance_step(neb, rk_data)
            @info("Advance step failed, end.")
            break
        end

        neb.driver.nsteps = rk_data.nsteps

        step_size = rk_data.step

        #we borrow neb.field to store dm
        compute_dm!(neb.field, neb.prespin, neb.spin, neb.nxyz)
        max_dmdt = maximum(view(neb.field, 1:neb.nxyz))/step_size

        all_max_dmdt = [0.0]

        if neb.comm_rank == 0
            all_max_dmdt[1] = MPI.Reduce(max_dmdt, max, 0, MPI.COMM_WORLD)
        else
            MPI.Reduce(max_dmdt, max, 0, MPI.COMM_WORLD)
        end
        MPI.Bcast!(all_max_dmdt, 0, MPI.COMM_WORLD)

        if neb.comm_rank == 0
            @info @sprintf("step =%5d  step_size=%10.6e  sim.t=%10.6e  max_dmdt=%10.6e  time=%s",
                            i, rk_data.step, rk_data.t, all_max_dmdt[1]/dmdt_factor, Dates.now())
        end

        if save_m_every>0 && i%save_m_every == 0
            compute_system_energy(neb)
            write_data(neb, neb.saver_energy)
            write_data(neb, neb.saver_distance)
        end
        if save_vtk_every > 0 && i%save_vtk_every == 0
            save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, i)))
        end
        if save_ovf_every > 0 && i%save_ovf_every == 0
            save_ovf(neb, joinpath(ovf_folder, @sprintf("%s_%d", neb.name, i)))
        end

        if all_max_dmdt[1] < stopping_dmdt*dmdt_factor
            @info @sprintf("max_dmdt is less than stopping_dmdt=%g, Done!", stopping_dmdt)

            if save_m_every>0
                compute_system_energy(neb)
                write_data(neb, neb.saver_energy)
                write_data(neb, neb.saver_distance)
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


function save_ovf(neb::NEB_GPU_MPI, name::String)
  sim = neb.sim
  dof = 3*sim.nxyz
  for n = 1:neb.N
      b = view(neb.spin, (n-1)*dof+1:n*dof)
      sim.spin .= b
      fname=@sprintf("%s_n_%d.ovf", name, n+neb.N*neb.comm_rank)
      save_ovf(sim, fname)
  end
end

function save_vtk(neb::NEB_GPU_MPI, name::String)
  sim = neb.sim
  dof = 3*sim.nxyz
  for n=1:neb.N
      b = view(neb.spin, (n-1)*dof+1:n*dof)
      sim.spin .= b
      fname=@sprintf("%s_n_%d", name, n+neb.N*neb.comm_rank)
      save_vtk(sim, fname, fields=fields)
  end
end
