function effective_field(zee::Zeeman, sim::AtomisticSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total

    groupsize = 512
    zeeman_kernel!(backend[], groupsize)(spin, zee.field, zee.energy, sim.mu_s, T(1);
                                         ndrange=N)
    KernelAbstractions.synchronize(backend[])
    return nothing
end

function effective_field(zee::TimeZeeman, sim::AtomisticSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    tx, ty, tz = zee.time_fun(t)
    groupsize = 512
    time_zeeman_kernel!(backend[], groupsize)(spin, zee.field, zee.init_field, zee.energy,
                                              sim.mu_s, T(1), T(tx), T(ty), T(tz);
                                              ndrange=N)
    KernelAbstractions.synchronize(backend[])
    return nothing
end

function effective_field(anis::Anisotropy, sim::AtomisticSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    axis = anis.axis

    groupsize = 512
    anisotropy_kernel!(backend[], groupsize)(spin, anis.field, anis.energy, anis.Ku,
                                             axis[1], axis[2], axis[3], sim.mu_s, T(1);
                                             ndrange=N)

    KernelAbstractions.synchronize(backend[])

    return nothing
end

function effective_field(anis::TubeAnisotropy, sim::AtomisticSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    groupsize = 512
    spatial_anisotropy_kernel!(backend[], groupsize)(spin, anis.field, anis.energy, anis.Ku,
                                                     anis.axes, sim.mu_s, T(1); ndrange=N)
    KernelAbstractions.synchronize(backend[])
    return nothing
end

function effective_field(exch::HeisenbergExchange, sim::AtomisticSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh
    groupsize = 512

    # The exchange interaction for nearest neighbours
    atomistic_exchange_kernel!(backend[], groupsize)(exch.field, exch.energy, exch.Js1,
                                                     spin, sim.mu_s, mesh.ngbs, mesh.n_ngbs,
                                                     T(0); ndrange=N)
    KernelAbstractions.synchronize(backend[])

    # The exchange interaction for next-nearest neighbours
    if hasproperty(mesh, n_ngbs2) && length(exch.Js2) == mesh.n_ngbs2
        atomistic_exchange_kernel!(backend[], groupsize)(exch.field, exch.energy, exch.Js2,
                                                         spin, sim.mu_s, mesh.ngbs2,
                                                         mesh.n_ngbs2, T(1); ndrange=N)
        KernelAbstractions.synchronize(backend[])
    end

    # The exchange interaction for next-next-nearest neighbours
    if hasproperty(mesh, n_ngbs3) && length(exch.Js3) == mesh.n_ngbs3
        atomistic_exchange_kernel!(backend[], groupsize)(exch.field, exch.energy, exch.Js3,
                                                         spin, sim.mu_s, mesh.ngbs3,
                                                         mesh.n_ngbs3, T(1); ndrange=N)
        KernelAbstractions.synchronize(backend[])
    end

    # The exchange interaction for next-next-next-nearest neighbours
    if hasproperty(mesh, n_ngbs4) && length(exch.Js4) == mesh.n_ngbs4
        atomistic_exchange_kernel!(backend[], groupsize)(exch.field, exch.energy, exch.Js4,
                                                         spin, sim.mu_s, mesh.ngbs4,
                                                         mesh.n_ngbs4, T(1); ndrange=N)
        KernelAbstractions.synchronize(backend[])
    end

    return nothing
end

function effective_field(dmi::HeisenbergBulkDMI, sim::AtomisticSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh

    groupsize = 512
    atomistic_dmi_kernel!(backend[], groupsize)(dmi.field, dmi.energy, dmi.Dij, spin,
                                                sim.mu_s, mesh.ngbs, mesh.n_ngbs; ndrange=N)

    KernelAbstractions.synchronize(backend[])

    return nothing
end

function effective_field(dmi::HeisenbergTubeBulkDMI, sim::AtomisticSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh

    groupsize = 512
    tube_bulk_dmi_kernel!(backend[], groupsize)(dmi.field, dmi.energy, dmi.D, dmi.Dij, spin,
                                                sim.mu_s, mesh.ngbs, mesh.n_ngbs, mesh.nr;
                                                ndrange=N)
    KernelAbstractions.synchronize(backend[])

    return nothing
end

function effective_field(stochastic::StochasticField, sim::AtomisticSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    integrator = sim.driver.ode

    if integrator.nsteps > stochastic.nsteps
        randn!(stochastic.eta)
        stochastic.nsteps = integrator.nsteps
    end

    mu0 = 4 * pi * 1e-7
    dt = integrator.step
    gamma = sim.driver.gamma
    alpha = sim.driver.alpha
    k_B = stochastic.k_B
    factor = 2 * alpha * k_B / (gamma * dt)

    volume = 1.0 / mu0 # we need this factor to make the energy density correctly
    groupsize = 512
    stochastic_field_kernel!(backend[], groupsize)(spin, stochastic.field,
                                                   stochastic.energy, sim.mu_s,
                                                   stochastic.eta, stochastic.temperature,
                                                   factor, volume; ndrange=N)

    KernelAbstractions.synchronize(backend[])
    return nothing
end
