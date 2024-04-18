using LinearAlgebra
using Printf
abstract type NEBDriver end

mutable struct NEB <: AbstractSim
    num_images::Int64  #number of images
    n_total::Int64 # n_total = neb.sim.n_total*num_images
    sim::AbstractSim
    driver::NEBDriver
    images  #size(images) = (3*n_total, N)
    pre_images
    spin #A pointer to images
    prespin #A pointer to pre_images
    field
    energy
    distance
    total_distance
    tangents
    name::String
    saver_energy::DataSaver
    saver_distance::DataSaver
    spring_constant::Float64
    NEB() = new()
end

"""
    NEB(sim::AbstractSim, init_images::Tuple, frames_between_images::Tuple; 
name="NEB", spring_constant=1.0e5, driver="LLG")

Create a NEB instance.
args:
    sim: Sim instance that include the interactions.
    init_images::Tuple:

"""
function NEB(sim::AbstractSim, given_images::Tuple, frames_between_images::Tuple; 
    name="NEB", spring_constant=1.0e5, driver="LLG")

    n_total = sim.mesh.nx*sim.mesh.ny*sim.mesh.nz
    num_images = sum(frames_between_images)+length(given_images)

    neb = NEB()
    neb.num_images = num_images
    neb.n_total = n_total * num_images
    neb.sim = sim
    neb.name = name
    neb.pre_images = create_zeros(3*sim.n_total, num_images)
    neb.field = create_zeros(3*sim.n_total, num_images)
    neb.distance = create_zeros(sim.n_total, num_images-1)
    neb.total_distance = create_zeros(num_images-1)
    neb.energy = create_zeros(num_images)
    neb.spring_constant = spring_constant
    neb.tangents = create_zeros(3*sim.n_total, num_images)

    neb.images = init_images(given_images, frames_between_images, sim.n_total)

    neb.spin = reshape(neb.images, 3*sim.n_total*num_images)
    neb.prespin = reshape(neb.pre_images, 3*sim.n_total*num_images)

    effective_field_NEB!(neb)
    #init_saver(neb)
    #neb.driver = create_neb_driver(driver, n_total, num_images)
    return neb
end

function init_images(given_images::Tuple, frames_between_images::Tuple, n_total)
    num_images = length(given_images) + sum(frames_between_images)
    images = create_zeros(3*n_total, num_images)
    images[:, 1] .= given_images[1]
    image_count = 2
    for i = 1:(length(given_images)-1)
        if frames_between_images[i] > 0
            for j = 1:frames_between_images[i]
                images[:, image_count] .= slerp(given_images[i], given_images[i+1], j/(frames_between_images[i]+1), n_total)
                image_count += 1
            end
            images[:, image_count] .= given_images[i+1]
            image_count += 1
        end
    end
    return images
end

function create_neb_driver(driver::String, n_total::Int64, N::Int64) #TODO: FIX ME
  if driver=="SD"
      gk = zeros(Float64,3*n_total,N)
      return NEB_SD(gk, 0.0, 1e-5, 1e-14, 0)
  elseif driver == "LLG"
      tol = 1e-5
      dopri5 = DormandPrince(n_total*N, neb_llg_call_back, tol)
      return NEB_LLG_Driver(0, dopri5)
  else
    error("Only Driver 'SD' or 'LLG' is supported!")
  end
end

function init_saver(neb::NEB)
    step_item = SaverItem("step", "", o::NEB -> o.driver.nsteps)
    neb.saver_energy = DataSaver(string(neb.name, ".txt"), false, 0.0, 0, [step_item])
    neb.saver_distance = DataSaver(string(@sprintf("%s_distance",neb.name), ".txt"), false, 0.0, 0, [step_item])
    N = neb.N
    for n = 1:N
        name = @sprintf("E_total_%g",n)
        item = SaverItem(name, "J", o::NEB -> o.energy[n])
        push!(neb.saver_energy.items, item)
    end
    
    for n = 1:N-1  #TODO: using a single function?
        name =  @sprintf("distance_%d",n)
        item = SaverItem(name, "", o::NEB -> o.distance[n])
        push!(neb.saver_distance.items, item)
    end
end

function effective_field_NEB!(neb::NEB)
    sim = neb.sim
    for n = 2:neb.num_images-1
        effective_field(sim, neb.images[:,n], 0.0)
        neb.field[:,n] .= sim.field
        neb.energy[n] = sum(sim.energy)
    end
    compute_tangents!(neb.images, neb.energy,neb.tangents, neb.num_images)
    compute_distance!(neb.images, neb.distance, neb.total_distance, neb.num_images, sim.n_total)

    dev = default_backend[]
    for i = 2:neb.num_images-1
        reduce_tangent_kernel!(dev, groupsize[])(neb.field[:,i], neb.tangents[:,i], ndrange=sim.n_total)
    end
    KernelAbstractions.synchronize(dev)
    for i = 2:neb.num_images-1
        neb.field[:,i] .+= neb.spring_constant*(neb.total_distance[i]-neb.total_distance[i-1])*neb.tangents[:,i]
    end
end

function compute_tangents!(images, energy, tangents, num_images)
    Threads.@threads for i=2:num_images-1
        if energy[i+1] > energy[i] > energy[i-1]
            tangents[:, i] .= images[:, i+1] - images[:, i]
        elseif energy[i+1] <= energy[i] <= energy[i-1]
            tangents[:, i] .= images[:, i] - images[:, i-1]
        elseif (energy[i+1] > energy[i]) && (energy[i-1] > energy[i])
                v1 = images[:, i+1] - images[:, i]
                v2 = images[:, i] - images[:, i-1]
                w1 = max(abs(E[i+1]-E[i]), abs(E[i-1]-E[i]))
                w2 = min(abs(E[i+1]-E[i]), abs(E[i-1]-E[i]))
            if energy[i+1] > energy[i-1]
                tangents[:, i] .= w1*v1+w2*v2
            else
                tangents[:, i] .= w1*v2+w2*v1
            end
        end
    end
end

function compute_distance!(images, distance, total_distance, num_images::Int, n_total::Int)
    dev = default_backend[]
    kernel! = compute_distance_kernel!(dev, groupsize[])
    for i=1:num_images-1
        kernel!(images[:, i], images[:, i+1], distance[:, i], ndrange=n_total)
    end
    KernelAbstractions.synchronize(dev)
    Threads.@threads for i=1:num_images-1
        total_distance[i] = LinearAlgebra.norm(distance[:, i])
    end
    KernelAbstractions.synchronize(dev)
end

function relax(neb::NEB; maxsteps=10000, stopping_dmdt=0.01, stopping_torque=0.1, 
  save_m_every = 10, save_vtk_every=-1, vtk_folder="vtks", save_ovf_every=-1, ovf_folder="ovfs", type=Float64)

  is_relax_NEB = false
  if isa(neb.driver, NEB_SD)
      is_relax_NEB = true
  end

  if save_vtk_every>0 && !isdir(vtk_folder)
      mkdir(vtk_folder)
  end
  if save_ovf_every>0 && !isdir(ovf_folder)
    mkdir(ovf_folder)
  end

  if is_relax_NEB
      relax_NEB(neb, maxsteps, Float64(stopping_torque), save_m_every, save_vtk_every, vtk_folder,save_ovf_every, ovf_folder,type)
  else
      relax_NEB_LLG(neb, maxsteps, Float64(stopping_dmdt), save_m_every, save_vtk_every, vtk_folder,save_ovf_every, ovf_folder,type)
  end
  return nothing
end

function save_vtk(neb::NEB, name::String, fields::Array{String, 1} = String[])
  sim = neb.sim
  for n=1:neb.N
    sim.spin[:]=neb.images[:,n]
    fname=@sprintf("%s_n_%d",name,n)
    save_vtk(sim, fname, fields=fields)
  end
end

function save_ovf(neb::NEB, name::String; type::DataType=Float64)
  sim = neb.sim
  for n=1:neb.N
    sim.spin[:]=neb.images[:,n]
    fname=  @sprintf("%s_n_%d",name,n)
    save_ovf(sim, fname, type=type)
  end
end

# function compute_system_energy(neb::NEB)
#   sim = neb.sim
#   images = neb.images
#   #sim.total_energy = 0
#   fill!(neb.energy, 0.0)

#   for n = 1:neb.N
#       effective_field(sim, images[:,n], 0.0)
#       neb.energy[n] = sum(sim.energy)
#   end

#   return 0
# end


# function compute_tangents(t::Array{Float64, 2},images::Array{Float64, 2},energy::Array{Float64, 1},N::Int,n_total::Int)
#     Threads.@threads  for n = 1:N
#       if (n==1)||(n==N)
#         continue
#       end
#       E1 = energy[n-1]
#       E2 = energy[n]
#       E3 = energy[n+1]
#       dEmax = max(abs(E3-E2),abs(E2-E1))
#       dEmin = min(abs(E3-E2),abs(E2-E1))
#       tip = images[:,n+1]-images[:,n]
#       tim = images[:,n]-images[:,n-1]
  
#       if (E1>E2)&&(E2>E3)
#         t[:,n] = tim
#       elseif (E3>E2)&&(E2>E1)
#         t[:,n] = tip
#       elseif E3>E1
#         t[:,n] = dEmax*tip+dEmin*tim
#       elseif E3<E1
#         t[:,n] = dEmin*tip+dEmax*tim
#       else
#         t[:,n] = tim + tip
#       end
  
#       for i = 1:n_total
#         j = 3*i - 2
#         fx,fy,fz = cross_product(images[j,n],images[j+1,n],images[j+2,n], t[j,n],t[j+1,n],t[j+2,n])
#         t[j,n],t[j+1,n],t[j+2,n] = cross_product(fx,fy,fz, images[j,n],images[j+1,n],images[j+2,n])
#       end
  
#       norm_t = LinearAlgebra.norm(t[:,n])
#       t[:,n] = t[:,n]/norm_t
  
#     end
  
#      return nothing
#   end
