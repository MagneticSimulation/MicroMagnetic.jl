using Random

mutable struct MonteCarlo{TF<:AbstractFloat} <:AbstractSimGPU
  mesh::TriangularMesh
  mu_s::Float64
  J::Float64
  D::Float64
  bulk_dmi::Bool
  Ku::Float64
  T::Float64
  Hx::Float64
  Hy::Float64
  Hz::Float64
  saver::DataSaver
  spin::CuArray{TF, 1}
  nextspin::CuArray{TF, 1}
  rnd::CuArray{TF, 1}
  energy::CuArray{TF, 1}
  total_energy::TF
  nxyz::Int64
  steps::Int64
  name::String
end

"""
Create a MonteCarlo simulation.

J, D and Ku in units of Joule
H in units of Tesla.

"""
function MonteCarlo(mesh::TriangularMesh; mu_s=1.0*mu_B, J=1*meV, D=0, Ku=0, H=(0,0,0), T=10, bulk_dmi=false, name="dyn")
  nxy = mesh.nx*mesh.ny
  Float = _cuda_using_double.x ? Float64 : Float32
  spin = cuzeros(Float, 3*nxy)
  nextspin = cuzeros(Float,3*nxy)
  rnd = cuzeros(Float,3*nxy)
  energy = cuzeros(Float,nxy)

  headers = ["step", "time", "E_total", ("m_x", "m_y", "m_z")]
  units = ["<>", "<s>", "<J>",("<>", "<>", "<>")]
  results = [o::AbstractSim -> o.saver.nsteps,
             o::AbstractSim -> o.saver.t,
             o::AbstractSim -> o.total_energy, average_m]
  saver = DataSaver(string(name, ".txt"), 0.0, 0, false, headers, units, results)

  return MonteCarlo(mesh, mu_s, J/k_B, D/k_B, bulk_dmi, Ku/k_B, Float64(T),
                    H[1]*mu_s/k_B, H[2]*mu_s/k_B, H[3]*mu_s/k_B,
                    saver, spin, nextspin, rnd, energy, Float(0.0), nxy, 0, name)
end

#function init_vector!(v::Array{T, 1}, mesh::TriangularMesh, init::Tuple{Real,Real,Real}) where {T<:AbstractFloat}
#  nxyz = mesh.nx*mesh.ny
#  b = reshape(v, 3, nxyz)
#  b[1, :] .= init[1]
#  b[2, :] .= init[2]
#  b[3, :] .= init[3]
#  return nothing
#end

function init_vector!(v::Array{T, 1}, mesh::TriangularMesh, init::Function) where {T<:AbstractFloat}
  nxyz = mesh.nx*mesh.ny
  dx,dy = mesh.dx, mesh.dy
  b = reshape(v, 3, nxyz)
  for i = 1:mesh.nx
    for j = 1:mesh.ny
        id = index(i, j, mesh.nx, mesh.ny)
        b[:, id] .=  init(i,j,dx,dy)
      end
  end
  return nothing
end

function init_m0(sim::MonteCarlo, m0::Any; norm=true)
  Float = _cuda_using_double.x ? Float64 : Float32
  spin = zeros(Float, 3*sim.nxyz)
  init_vector!(spin, sim.mesh, m0)
  if norm
    normalise(spin, sim.nxyz)
  end
  copyto!(sim.spin, spin)
  return true
end

function uniform_random_sphere_kernel!(m::CuDeviceArray{T, 1}, rnd::CuDeviceArray{T, 1}, nxyz::Int64) where {T<:AbstractFloat}
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= nxyz
        j = 3*index - 2
        @inbounds phi = rnd[j]*2*pi;
        @inbounds ct = 2*rnd[j+1] - 1;
        st = CUDAnative.sqrt(1-ct*ct);
        @inbounds m[j] = st*CUDAnative.cos(phi);
        @inbounds m[j+1] = st*CUDAnative.sin(phi);
        @inbounds m[j+2] = ct;
    end
   return nothing
end

function uniform_random_sphere(spin::CuArray{T, 1}, rnd::CuArray{T, 1}, N::Int64) where{T<:AbstractFloat}

    rand!(rnd)
    blk, thr = CuArrays.cudims(N)
    @cuda blocks=blk threads=thr uniform_random_sphere_kernel!(spin, rnd, N)

    return  nothing
end

function run_step_kernel!(m::CuDeviceArray{T, 1}, next_m::CuDeviceArray{T, 1}, rnd::CuDeviceArray{T, 1},
                              energy::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              J::T, D::T, bulk_dmi::Bool, Ku::T, Hx::T, Hy::T, Hz::T, temp::T,
                              nx::Int64, ny::Int64, bias::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    #for bulk dmi
    Rx = (T(1), T(1/2), T(-1/2), T(-1), T(-1/2), T(1/2))
    Ry = (T(0), T(sqrt(3)/2), T(sqrt(3)/2), T(0), T(-sqrt(3)/2), T(-sqrt(3)/2))

    delta_E = T(0.0)
    update = false
    #bias should be 0, 1 and 2.

    if index <= nx*ny
        a, b = Tuple(CuArrays.CartesianIndices((nx,ny))[index])
        if (a+b)%3 != bias
            return nothing
        end
        i = 3*index - 2
        @inbounds dmx = next_m[i] - m[i]
        @inbounds dmy = next_m[i+1] - m[i+1]
        @inbounds dmz = next_m[i+2] - m[i+2]

        delta_E = -(dmx*Hx+dmy*Hy+dmz*Hz) #zeeman
        @inbounds delta_E += Ku*(m[i+2]*m[i+2] - next_m[i+2]*next_m[i+2]) #Anisotropy

        for j = 1:6
            id = ngbs[j, index]
            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                delta_E -= J*(dmx*sx + dmy*sy + dmz*sz) #exchange

                Dx = bulk_dmi ? D*Rx[j] : -D*Ry[j]
                Dy = bulk_dmi ? D*Ry[j] : D*Rx[j]
                delta_E += volume(dmx, dmy, dmz, sx, sy, sz, Dx, Dy, T(0)) #DMI
            end
        end
        if delta_E < 0
            update = true
        else
            if rnd[i+2] < CUDAnative.exp(-delta_E/temp)
                update = true
            end
        end

        if update == true
          @inbounds m[i] = next_m[i];
          @inbounds m[i+1] = next_m[i+1];
          @inbounds m[i+2] = next_m[i+2]
        end

    end
   return nothing
end

function run_single_step(sim::MonteCarlo, bias::Int64)
  F = _cuda_using_double.x ? Float64 : Float32
  blk, thr = CuArrays.cudims(sim.nxyz)
  @cuda blocks=blk threads=thr run_step_kernel!(sim.spin, sim.nextspin,
                               sim.rnd, sim.energy, sim.mesh.ngbs,
                               F(sim.J), F(sim.D), sim.bulk_dmi, F(sim.Ku),
                               F(sim.Hx), F(sim.Hy), F(sim.Hz), F(sim.T),
                               sim.mesh.nx, sim.mesh.ny, bias)
  return nothing
end


function run_step(sim::MonteCarlo)

  uniform_random_sphere(sim.nextspin, sim.rnd, sim.nxyz)
  run_single_step(sim, 0)
  run_single_step(sim, 1)
  run_single_step(sim, 2)
  sim.steps += 1

  return  nothing
end

function total_energy_kernel!(m::CuDeviceArray{T, 1}, energy::CuDeviceArray{T, 1}, ngbs::CuDeviceArray{Int32, 2},
                              J::T, D::T, bulk_dmi::Bool, Ku::T, Hx::T, Hy::T, Hz::T,
                              nxy::Int64) where {T<:AbstractFloat}

    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    #for bulk dmi
    Rx = (T(1), T(1/2), T(-1/2), T(-1), T(-1/2), T(1/2))
    Ry = (T(0), T(sqrt(3)/2), T(sqrt(3)/2), T(0), T(-sqrt(3)/2), T(-sqrt(3)/2))

    E = T(0.0)

    if index <= nxy
        i = 3*index - 2
        @inbounds mx = m[i]
        @inbounds my = m[i+1]
        @inbounds mz = m[i+2]

        E = -(mx*Hx+my*Hy+mz*Hz) #zeeman
        E -= Ku*(mz*mz) #Anisotropy

        for j = 1:6
            id = ngbs[j, index]
            if id>0
                k = 3*id-2
                @inbounds sx = m[k]
                @inbounds sy = m[k+1]
                @inbounds sz = m[k+2]
                E -= 0.5*J*(mx*sx + my*sy + mz*sz) #exchange

                Dx = bulk_dmi ? D*Rx[j] : -D*Ry[j]
                Dy = bulk_dmi ? D*Ry[j] : D*Rx[j]
                E += 0.5*volume(mx, my, mz, sx, sy, sz, Dx, Dy, T(0)) #DMI
            end
        end

        @inbounds energy[index] = E

    end
   return nothing
end

function compute_system_energy(sim::MonteCarlo)
    F = _cuda_using_double.x ? Float64 : Float32
    blk, thr = CuArrays.cudims(sim.nxyz)
    @cuda blocks=blk threads=thr total_energy_kernel!(sim.spin, sim.energy,
                                 sim.mesh.ngbs,
                                 F(sim.J), F(sim.D), sim.bulk_dmi, F(sim.Ku),
                                 F(sim.Hx), F(sim.Hy), F(sim.Hz),
                                 sim.mesh.nxyz)
    return sum(sim.energy)*k_B
end

function run_sim(sim::MonteCarlo; maxsteps=10000, save_m_every = 10, save_vtk_every=-1)
    for i=1:maxsteps

        if save_m_every>0

            if sim.steps%save_m_every == 0
                @info @sprintf("step=%5d  total_energy=%g", sim.steps, compute_system_energy(sim))
            #compute_system_energy(sim, sim.spin, 0.0)
            #write_data(sim)
            end
        end

        if save_vtk_every > 0
            if sim.steps%save_vtk_every == 0
              save_vtk(sim, @sprintf("%s_%d", sim.name, sim.steps))
            end
        end

        run_step(sim)

    end
end
