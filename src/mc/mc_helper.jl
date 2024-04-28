abstract type MonteCarloEnergy end

mutable struct ExchangeDMI{T<:AbstractFloat} <: MonteCarloEnergy
    Jx::T #nearest exchange constant
    Jy::T
    Jz::T
    D::AbstractArray{T,2} #nearest DMI vector
end

mutable struct ZeemanMC{T<:AbstractFloat} <: MonteCarloEnergy
    Hx::T
    Hy::T
    Hz::T
end

mutable struct AnisotropyMC{T<:AbstractFloat} <: MonteCarloEnergy
    Ku::T
    axis::Tuple{T,T,T}
    Kc::T
end

mutable struct KagomeAnisotropyMC{T<:AbstractFloat} <: MonteCarloEnergy
    Ku::T
end

mutable struct MonteCarlo{TF<:AbstractFloat} <: AbstractSim
    mesh::Mesh
    shape::AbstractArray{Bool,1}
    exch::MonteCarloEnergy
    zeeman::MonteCarloEnergy
    anis::MonteCarloEnergy
    saver::DataSaver
    spin::AbstractArray{TF,1}
    nextspin::AbstractArray{TF,1}
    rnd::AbstractArray{TF,1}
    energy::AbstractArray{TF,1}
    delta_E::AbstractArray{TF,1}
    n_total::Int64
    steps::Int64
    name::String
    T::TF
    mc_2d::Bool
    MonteCarlo{T}() where {T<:AbstractFloat} = new()
end

function uniform_random_sphere(spin::AbstractArray{T,1}, rnd::AbstractArray{T,1},
                               N::Int64) where {T<:AbstractFloat}
    rand!(rnd)
    kernel! = uniform_random_sphere_kernel!(default_backend[], groupsize[])
    kernel!(spin, rnd; ndrange=N)
    KernelAbstractions.synchronize(default_backend[])

    return nothing
end

function uniform_random_circle_xy(spin::AbstractArray{T,1}, rnd::AbstractArray{T,1},
                                  N::Int64) where {T<:AbstractFloat}
    rand!(rnd)
    kernel! = uniform_random_circle_xy_kernel!(default_backend[], groupsize[])
    kernel!(spin, rnd; ndrange=N)
    KernelAbstractions.synchronize(default_backend[])

    return nothing
end

function compute_dE_zeeman_anisotropy_energy(sim::MonteCarlo, za::AnisotropyMC, bias::Int64)
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    cubic = isa(mesh, CubicMesh)

    ze = sim.zeeman
    kernel! = dE_zeeman_anisotropy_energy_kernel!(default_backend[], groupsize[])
    kernel!(sim.spin, sim.nextspin, sim.shape, sim.delta_E, ze.Hx, ze.Hy, ze.Hz, za.Ku,
            za.Kc, za.axis[1], za.axis[2], za.axis[3], bias, cubic; ndrange=(nx, ny, nz))
    KernelAbstractions.synchronize(default_backend[])
    return nothing
end

function compute_zeeman_anisotropy_energy(sim::MonteCarlo, za::AnisotropyMC)
    ze = sim.zeeman
    kernel! = zeeman_anisotropy_energy_kernel!(default_backend[], groupsize[])
    kernel!(sim.spin, sim.shape, sim.energy, ze.Hx, ze.Hy, ze.Hz, za.Ku, za.Kc, za.axis[1],
            za.axis[2], za.axis[3]; ndrange=sim.n_total)
    KernelAbstractions.synchronize(default_backend[])
    return nothing
end

function add_dE_exch_dmi_energy(sim::MonteCarlo, exch::ExchangeDMI, bias::Int64)
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    cubic = isa(mesh, CubicMesh)

    ex = exch
    kernel! = add_dE_exch_dmi_energy_kernel!(default_backend[], groupsize[])
    kernel!(sim.spin, sim.nextspin, sim.shape, sim.delta_E, mesh.ngbs, mesh.n_ngbs, ex.Jx,
            ex.Jy, ex.Jz, ex.D, bias, cubic; ndrange=(nx, ny, nz))
    KernelAbstractions.synchronize(default_backend[])
    return nothing
end

function add_exch_dmi_energy(sim::MonteCarlo, exch::ExchangeDMI)
    mesh = sim.mesh
    ex = exch
    kernel! = add_exch_dmi_energy_kernel!(default_backend[], groupsize[])
    kernel!(sim.spin, sim.shape, sim.energy, mesh.ngbs, mesh.n_ngbs, ex.Jx, ex.Jy, ex.Jz,
            ex.D; ndrange=sim.n_total)
    KernelAbstractions.synchronize(default_backend[])
    return nothing
end

function compute_system_energy(sim::MonteCarlo)
    compute_zeeman_anisotropy_energy(sim, sim.anis)
    add_exch_dmi_energy(sim, sim.exch)
    return sum(sim.energy)
end

function run_step_bias(sim::MonteCarlo, bias::Int64)
    mesh = sim.mesh
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz

    compute_dE_zeeman_anisotropy_energy(sim, sim.anis, bias) #compute zeeman and anis energy
    add_dE_exch_dmi_energy(sim, sim.exch, bias) #include exch and DMI 

    cubic = isa(mesh, CubicMesh)

    kernel! = run_monte_carlo_kernel!(default_backend[], groupsize[])
    kernel!(sim.spin, sim.nextspin, sim.rnd, sim.shape, sim.delta_E, sim.T, bias, cubic;
            ndrange=(nx, ny, nz))
    KernelAbstractions.synchronize(default_backend[])

    return nothing
end
