using Random
using WriteVTK

function MonteCarlo(mesh::Mesh; name="mc")
    Float = _cuda_using_double.x ? Float64 : Float32
    sim = MonteCarlo{Float}()
    sim.mesh = mesh
    nxyz = mesh.nxyz

    sim.nxyz     = nxyz
    sim.spin     = CUDA.zeros(Float,3*nxyz)
    sim.nextspin = CUDA.zeros(Float,3*nxyz)
    sim.rnd      = CUDA.zeros(Float,3*nxyz)
    sim.energy   = CUDA.zeros(Float,nxyz)
    sim.energy2  = CUDA.zeros(Float,nxyz)
    sim.delta_E  = CUDA.zeros(Float,nxyz)
    sim.Qdens    = CUDA.zeros(Float,nxyz)
    sim.spin2    = CUDA.zeros(Float, 3*nxyz)
    sim.EnMean   = 0.0
    sim.En2Mean  = 0.0
    sim.Q        = 0.0
    sim.mxMean   = 0.0
    sim.myMean   = 0.0
    sim.mzMean   = 0.0
    sim.mx2Mean  = 0.0
    sim.my2Mean  = 0.0
    sim.mz2Mean  = 0.0
    sim.steps    = 0
    sim.name     = name
    sim.T        = 300

    sim.exch = ExchangeMC(CUDA.zeros(Float, mesh.TotRdN))
    sim.zee = ZeemanMC(CUDA.zeros(Float, 3))
    sim.anis = Anis6Fold2DMC(CUDA.zeros(Float, 2))

    # headers = ["step", "time", "E_total", ("m_x", "m_y", "m_z")]
    # units = ["<>", "<s>", "<J>",("<>", "<>", "<>")]
    # results = [o::AbstractSim -> o.saver.nsteps,
    #          o::AbstractSim -> o.saver.t,
    #          o::AbstractSim -> o.total_energy, average_m]
    # saver = DataSaver(string(name, ".txt"), 0.0, 0, false, headers, units, results)
    # sim.saver = saver

    return sim
end


function set_exch(sim::MonteCarlo; Jxy=0, Jz=0)
    T = _cuda_using_double.x ? Float64 : Float32
    J1 = Jxy/k_B
    J2 = Jz/k_B
    # J3 = J3/k_B
    J = Array([T(J1), T(J2)])
    copyto!(sim.exch.J, J)
    return nothing
end

#Hx, Hy, Hz in energy unit， just as J and D
function set_zeeman(sim::MonteCarlo; Hx=0, Hy=0, Hz=0)
    T = _cuda_using_double.x ? Float64 : Float32
    H = Array([T(Hx/k_B), T( Hy/k_B), T(Hz/k_B)])
    copyto!(sim.zee.H, H)
    return nothing
end

function set_anis_6fold2d(sim::MonteCarlo; K1=0,K2=0)
  T = _cuda_using_double.x ? Float64 : Float32
  copyto!(sim.anis.K, Array([T(K1/k_B),T(K2/k_B)]))
  return nothing
end

function init_m0(sim::MonteCarlo, m0::Union{Function, Tuple, Array}; norm=true)
  Float = _cuda_using_double.x ? Float64 : Float32
  spin = zeros(Float, 3*sim.nxyz)
  init_vector!(spin, sim.mesh, m0)
  copyto!(sim.spin, spin)
  if norm
    normalise(sim.spin, sim.nxyz)
  end
  # shape = Array(sim.shape)
  # for i in 1:sim.nxyz
  #     if !shape[i]
  #         spin[3*i-2] = NaN32
  #         spin[3*i-1] = NaN32
  #         spin[3*i] = NaN32
  #     end
  # end
  return true
end

function run_step_triangular(sim::MonteCarlo)

  uniform_random_sphere(sim.nextspin, sim.rnd, sim.nxyz)

  blk, thr = cudims(sim.nxyz)
  for bias = 0:4
    @cuda blocks=blk threads=thr  run_monte_carlo_kernel!(sim.spin,
                                        sim.nextspin, sim.rnd,
                                        sim.delta_E, sim.T,
                                        sim.zee.H, sim.exch.J, sim.anis.K,
                                        sim.mesh.TotRdN,sim.mesh.ngbs,sim.mesh.n_ngbs,
                                        sim.mesh.nx, sim.mesh.ny, sim.mesh.nz,bias)
  end
  sim.steps += 1

  return  nothing
end

function calcEnAndStd(sim::MonteCarlo;Volu=false)
  blk, thr = cudims(sim.nxyz)
  @cuda blocks=blk threads=thr   getEnergyAndStd!(sim.spin,sim.spin2,
                              sim.energy, sim.energy2, sim.Qdens,
                              sim.zee.H,sim.exch.J, sim.anis.K,
                              sim.mesh.TotRdN,sim.mesh.ngbs,sim.mesh.n_ngbs,
                              sim.mesh.nx, sim.mesh.ny, sim.mesh.nz,Volu)
  sim.EnMean  =Float64(sum(sim.energy))/sim.nxyz
  sim.En2Mean =Float64(sum(sim.energy2))/sim.nxyz
  b1 = reshape(sim.spin, 3, sim.nxyz)
  b2 = reshape(sim.spin2,3, sim.nxyz)
  sim.mxMean, sim.myMean, sim.mzMean = Tuple(sum(b1, dims=2)./sim.nxyz) 
  sim.mx2Mean,sim.my2Mean,sim.mz2Mean= Tuple(sum(b2, dims=2)./sim.nxyz) 
  sim.Q=Float64(sum(sim.Qdens))
  return  nothing
end

function run_sim(sim::MonteCarlo; maxsteps=10000, save_m_every = 10, save_vtk_every=-1, save_ovf_every=-1, ovf_format="binary8")
    # cubic = isa(sim.mesh, CubicMeshGPU) ? true : false
    T = _cuda_using_double.x ? Float64 : Float32
    for i=1:maxsteps

        if save_m_every>0
            if sim.steps%save_m_every == 0
                calcEnAndStd(sim)
                @info @sprintf("step=%10d T= %6g total_energy= %6g  mx= %6g  my= %6g  mz= %6g  Q= %6g", sim.steps, sim.T,sim.energy,sim.mxMean,sim.myMean,sim.mzMean,sim.Q)
            end
        end

        if save_ovf_every > 0
            if sim.steps%save_ovf_every == 0
                save_ovf(sim, @sprintf("%s_%d", sim.name, sim.steps), dataformat = ovf_format)
            end
        end

        if save_vtk_every > 0
            if sim.steps%save_vtk_every == 0
              save_vtk(sim, @sprintf("%s_%d", sim.name, sim.steps))
            end
        end

        # if cubic
        #    run_step_cubic(sim)
        # else
           run_step_triangular(sim)
        # end

    end
end

function uniform_random_sphere(spin::CuArray{T, 1}, rnd::CuArray{T, 1}, N::Int64) where{T<:AbstractFloat}

    CUDA.rand!(rnd)
    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr uniform_random_sphere_kernel!(spin, rnd, N)

    return  nothing
end

function uniform_random_circle_xy(spin::CuArray{T, 1}, rnd::CuArray{T, 1}, N::Int64) where{T<:AbstractFloat}

    CUDA.rand!(rnd)
    blk, thr = cudims(N)
    @cuda blocks=blk threads=thr uniform_random_circle_xy_kernel!(spin, rnd, N)

    return  nothing
end
