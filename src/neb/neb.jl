using LinearAlgebra
using Printf
abstract type NEBDriver end

mutable struct NEB <: AbstractSim
  N::Int64  #number of images
  nxyz::Int64 # nxyz = neb.sim.nxyz*N
  sim::AbstractSim
  init_m::Any
  driver::NEBDriver
  images::Array{Float64, 2}  #size(images) = (3*nxyz, N)
  pre_images::Array{Float64, 2}
  spin::Array{Float64, 1} #A pointer to images
  prespin::Array{Float64, 1} #A pointer to pre_images
  name::String
  field::Array{Float64, 2}
  energy::Array{Float64, 1}
  distance::Array{Float64, 1}
  saver_energy::DataSaver
  saver_distance::DataSaver
  spring_constant::Float64
  tangents::Array{Float64, 2}
end

function NEB(sim::AbstractSim, init_m::Any, intervals::Any; name="NEB", spring_constant=1.0e5, driver="LLG")
  nxyz = sim.mesh.nx*sim.mesh.ny*sim.mesh.nz
  N = sum(intervals)+length(intervals)+1
  images = zeros(Float64,3*sim.nxyz,N)
  tangents = zeros(Float64,3*sim.nxyz,N)
  pre_images = zeros(Float64,3*sim.nxyz,N)
  energy = zeros(Float64, N)
  Ms = zeros(Float64,nxyz)
  driver = create_neb_driver(driver, nxyz, N)
  field = zeros(Float64,3*nxyz, N)
  distance=zeros(Float64,N-1)

  headers = ["steps"]
  units = ["<>"]
  results = Any[o::NEB -> o.driver.nsteps]
  distance_headers= ["steps"]
  distance_units = ["<>"]
  distance_results = Any[o::NEB -> o.driver.nsteps]
  saver = DataSaver(string(name, ".txt"), 0.0, 0, false, headers, units, results)
  saver_distance = DataSaver(string(@sprintf("%s_distance",name), ".txt"),
                                    0.0, 0, false, distance_headers, distance_units, distance_results)


  spins = reshape(images, 3*sim.nxyz*N)
  prespins = reshape(pre_images, 3*sim.nxyz*N)

  neb = NEB(N, sim.nxyz*N, sim, init_m, driver, images, pre_images, spins, prespins,
             name, field, energy, distance, saver, saver_distance, spring_constant, tangents)
  init_images(neb, intervals)
  compute_distance(neb.distance,neb.images,N,nxyz)
  compute_system_energy(neb)
  return neb
end

function create_neb_driver(driver::String, nxyz::Int64, N::Int64) #TODO: FIX ME
  if driver=="SD"
      gk = zeros(Float64,3*nxyz,N)
      return NEB_SD(gk, 0.0, 1e-5, 1e-14, 0)
  elseif driver == "LLG"
      tol = 1e-5
      dopri5 = DormandPrince(nxyz*N, neb_llg_call_back, tol)
      return NEB_LLG_Driver(0, dopri5)
  else
    error("Only Driver 'SD' or 'LLG' is supported!")
  end
end

function init_images(neb::NEB, intervals::Any)
  sim=neb.sim
  nxyz=sim.nxyz
  N=neb.N
  images = neb.images
  pics = neb.init_m
  if length(pics)==length(intervals)+1
    n=1
    for i=1:length(intervals)
        init_m0(sim, pics[i])
        m1 = Array(sim.spin)
        init_m0(sim, pics[i+1])
        m2 = Array(sim.spin)
        M = interpolate_m(m1, m2, Int(intervals[i]))
        #M = interpolate_m_spherical(m1,m2,Int(intervals[i]))
        for j=1:intervals[i]+1
          images[:,n]=M[:,j]
          n += 1
        end
    end
    init_m0(sim, pics[end])
    images[:,N] .= Array(sim.spin)
  else
    println("Input error!")
  end

  #normalise(neb.spin, neb.nxyz)
  neb.pre_images[:]=neb.images[:]

  for n = 1:N
    name = @sprintf("E_total_%g",n)
    push!(neb.saver_energy.headers,name)
    push!(neb.saver_energy.units, "<J>")
    fun =  o::NEB -> o.energy[n]
    push!(neb.saver_energy.results, fun)
  end

  for n = 1:N-1  #TODO: using a single function?
    name =  @sprintf("distance_%d",n)
    push!(neb.saver_distance.headers,name)
    push!(neb.saver_distance.units, "<>")
    fun =  o::NEB -> o.distance[n]
    push!(neb.saver_distance.results, fun)
  end
end

function effective_field_NEB(neb::NEB, spin::Array{Float64, 1})
    sim = neb.sim
    nxyz = sim.nxyz
    #neb.spin[:] = spin[:], we already copy spin to neb.images in dopri5
    images = neb.images
    N = neb.N
    fill!(neb.field, 0.0)
    fill!(neb.energy, 0.0)


    for n = 2:neb.N-1
        effective_field(sim, images[:,n], 0.0)
        neb.field[:,n] .= sim.field
        neb.energy[n] = sum(sim.energy)
    end


    compute_tangents(neb.tangents,neb.images,neb.energy,N,nxyz)

    for n=1:N-1
      f = neb.field[:,n]'*neb.tangents[:,n]
      neb.field[:,n] .-= f*neb.tangents[:,n]
    end

    compute_distance(neb.distance,neb.images,N,nxyz)

    for n=2:N-1
      neb.field[:,n] .+= neb.spring_constant*(neb.distance[n]-neb.distance[n-1])*neb.tangents[:,n]
    end
end


function compute_tangents(t::Array{Float64, 2},images::Array{Float64, 2},energy::Array{Float64, 1},N::Int,nxyz::Int)
  Threads.@threads  for n = 1:N
    if (n==1)||(n==N)
      continue
    end
    E1 = energy[n-1]
    E2 = energy[n]
    E3 = energy[n+1]
    dEmax = max(abs(E3-E2),abs(E2-E1))
    dEmin = min(abs(E3-E2),abs(E2-E1))
    tip = images[:,n+1]-images[:,n]
    tim = images[:,n]-images[:,n-1]

    if (E1>E2)&&(E2>E3)
      t[:,n] = tim
    elseif (E3>E2)&&(E2>E1)
      t[:,n] = tip
    elseif E3>E1
      t[:,n] = dEmax*tip+dEmin*tim
    elseif E3<E1
      t[:,n] = dEmin*tip+dEmax*tim
    else
      t[:,n] = tim + tip
    end

    for i = 1:nxyz
      j = 3*i - 2
      fx,fy,fz = cross_product(images[j,n],images[j+1,n],images[j+2,n], t[j,n],t[j+1,n],t[j+2,n])
      t[j,n],t[j+1,n],t[j+2,n] = cross_product(fx,fy,fz, images[j,n],images[j+1,n],images[j+2,n])
    end

    norm_t = LinearAlgebra.norm(t[:,n])
    t[:,n] = t[:,n]/norm_t

  end

   return nothing
end

function compute_distance(distance::Array{Float64,1},images::Array{Float64, 2},N::Int,nxyz::Int)
  Threads.@threads for n=1:N-1
    m1 = images[:,n]
    m2 = images[:,n+1]
    l =zeros(nxyz)
    for i=1:nxyz
        j=3*i-2
        m1xm2=cross_product(m1[j],m1[j+1],m1[j+2],m2[j],m2[j+1],m2[j+2])
        l[i]=atan(norm(m1xm2),m1[j]*m2[j]+m1[j+1]*m2[j+1]+m1[j+2]*m2[j+2])
    end
    distance[n] = LinearAlgebra.norm(l)
  end
end

function compute_system_energy(neb::NEB)
  sim = neb.sim
  images = neb.images
  #sim.total_energy = 0
  fill!(neb.energy, 0.0)

  for n = 1:neb.N
      effective_field(sim, images[:,n], 0.0)
      neb.energy[n] = sum(sim.energy)
  end

  return 0
end

function relax(neb::NEB; maxsteps=10000, stopping_dmdt=0.01, stopping_torque=0.1, save_m_every = 10, save_vtk_every=-1, vtk_folder="vtks",save_ovf_every=-1, ovf_folder="ovfs",ovf_format="binary")

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
      relax_NEB(neb, maxsteps, Float64(stopping_torque), save_m_every, save_vtk_every, vtk_folder,save_ovf_every, ovf_folder,ovf_format)
  else
      relax_NEB_LLG(neb, maxsteps, Float64(stopping_dmdt), save_m_every, save_vtk_every, vtk_folder,save_ovf_every, ovf_folder,ovf_format)
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

function save_ovf(neb::NEB, name::String;ovf_format::String="binary")
  sim = neb.sim
  for n=1:neb.N
    sim.spin[:]=neb.images[:,n]
    fname=@sprintf("%s_n_%d",name,n)
    save_ovf(sim,fname,dataformat=ovf_format)
  end
end

#rotate m1 in the m1_m2 plane by theta. If m1 is parallel to m2 (m1 x m2 = 0), the plane
#is determined by m1 and ez (unless m1 itself is ez).
function rotate_to(m1::Array{T,1}, m2::Array{T,1}, theta::T) where T<:AbstractFloat
    normalise(m1, 1)
    normalise(m2, 1)
    if m1 == m2
      return m1
    end
    m = get_rotation_axis(m1,m2)
    return rotation_operator(m1,m,theta)
end

function get_rotation_axis(m1::Array{T,1}, m2::Array{T,1})where T<:AbstractFloat
    m = LinearAlgebra.cross(m1, m2)
    if m == [0.0,0,0]
        if m1 == [0,0,1.0]
            m = [1.0,0,0]
        else
            m = LinearAlgebra.cross(m1,[0,0,1.0])
        end
    end
    normalise(m, 1)
    return m
end
##rotate m1 with normalised axis m, angle theta
function rotation_operator(m1::Array{T,1}, m::Array{T,1}, theta::T) where T<:AbstractFloat
  st, ct = sin(theta), cos(theta)
  a11 = ct+(1-ct)*m[1]^2
  a12 = (1-ct)*m[1]*m[2] - st*m[3]
  a13 = (1-ct)*m[1]*m[3] + st*m[2]
  a21 = (1-ct)*m[2]*m[1] + st*m[3]
  a22 = ct+(1-ct)*m[2]*m[2]
  a23 = (1-ct)*m[2]*m[3] - st*m[1]
  a31 = (1-ct)*m[1]*m[3] - st*m[2]
  a32 = (1-ct)*m[2]*m[3] + st*m[1]
  a33 = ct+(1-ct)*m[3]*m[3]
  op=[a11 a12 a13; a21 a22 a23; a31 a32 a33]

  return op*m1
end

#Interpolate magnetization between image m1 and image m2, where the total images
#number is N+1
function interpolate_m(m1::Array{T,1}, m2::Array{T,1}, N::Int) where T<:AbstractFloat
    nxyz = Int(length(m1)/3)
    m = zeros(3*nxyz,N+1)
    b1 = reshape(m1,3,nxyz)
    b2 = reshape(m2,3,nxyz)

    for i=1:N+1
        for j=1:nxyz
            k=3*j-2
            amp = b1[:,j]'*b2[:,j]
            if amp > 1
                amp = 1.0
            end
            if amp < -1
                amp = -1.0
            end
            theta = acos(amp)
            dtheta = theta/(N+1)
            angle = (i-1)*dtheta
            m[k,i],m[k+1,i],m[k+2,i]=rotate_to(b1[:,j],b2[:,j],angle)
        end
    end
    return m
end

function cartesian2spherical(m, nxyz)
    R, theta, phi  = zeros(nxyz), zeros(nxyz), zeros(nxyz)
    for j=1:nxyz
        k=3*j-2
        r = sqrt(m[k]^2+m[k+1]^2)
        theta[j] = atan(r, m[k+2])
        phi[j] = atan(m[k+1], m[k])
        R[j] = sqrt(r^2+m[k+2]^2)
    end
    return R, theta, phi
end

#Interpolate magnetization between image m1 and image m2, where the total images
#number is N+1
function interpolate_m_spherical(m1::Array{T,1}, m2::Array{T,1}, N::Int) where T<:AbstractFloat
    nxyz = Int(length(m1)/3)

    R1, theta1, phi1 = cartesian2spherical(m1, nxyz)
    R2, theta2, phi2  = cartesian2spherical(m2, nxyz)

    m = zeros(3*nxyz,N+1)

    for i=1:N+1
        v = view(m, :, i)
        for j = 1:nxyz
            k = 3*j-2
            theta = theta1[j] + (i-1)/(N+1)*(theta2[j] - theta1[j])
            phi = phi1[j] + (i-1)/(N+1)*(phi2[j] - phi1[j])
            v[k] = R1[j]*sin(theta)*cos(phi)
            v[k+1] = R1[j]*sin(theta)*sin(phi)
            v[k+2] = R1[j]*cos(theta)
        end
    end
    return m
end
