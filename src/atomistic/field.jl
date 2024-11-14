function effective_field(zee::Zeeman, sim::AtomisticSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total

    kernal = zeeman_kernel!(default_backend[], groupsize[])
    kernal(spin, zee.field, zee.energy, sim.mu_s, T(1); ndrange=N)
    
    return nothing
end

function effective_field(zee::TimeZeeman, sim::AtomisticSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    tx, ty, tz = zee.time_fun(t)
    kernal = time_zeeman_kernel!(default_backend[], groupsize[])
    kernal(spin, zee.field, zee.init_field, zee.energy, sim.mu_s, T(1), T(tx), T(ty), T(tz);
           ndrange=N)
    
    return nothing
end

function effective_field(anis::Anisotropy, sim::AtomisticSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    axis = anis.axis

    kernal = anisotropy_kernel!(default_backend[], groupsize[])
    kernal(spin, anis.field, anis.energy, anis.Ku, axis[1], axis[2], axis[3], sim.mu_s,
           T(1); ndrange=N)

    return nothing
end

function effective_field(anis::TubeAnisotropy, sim::AtomisticSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    kernal = spatial_anisotropy_kernel!(default_backend[], groupsize[])
    kernal(spin, anis.field, anis.energy, anis.Ku, anis.axes, sim.mu_s, T(1); ndrange=N)

    return nothing
end

function effective_field(exch::HeisenbergExchange, sim::AtomisticSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh

    # The exchange interaction for nearest neighbours
    kernal = atomistic_exchange_kernel!(default_backend[], groupsize[])
    kernal(exch.field, exch.energy, exch.Js1, spin, sim.mu_s, mesh.ngbs, mesh.n_ngbs, T(0);
           ndrange=N)
    

    # The exchange interaction for next-nearest neighbours
    if hasproperty(mesh, :n_ngbs2) && length(exch.Js2) == mesh.n_ngbs2
        kernal(exch.field, exch.energy, exch.Js2, spin, sim.mu_s, mesh.ngbs2, mesh.n_ngbs2,
               T(1); ndrange=N)
    end

    # The exchange interaction for next-next-nearest neighbours
    if hasproperty(mesh, :n_ngbs3) && length(exch.Js3) == mesh.n_ngbs3
        kernal(exch.field, exch.energy, exch.Js3, spin, sim.mu_s, mesh.ngbs3, mesh.n_ngbs3,
               T(1); ndrange=N)
    end

    # The exchange interaction for next-next-next-nearest neighbours
    if hasproperty(mesh, :n_ngbs4) && length(exch.Js4) == mesh.n_ngbs4
        kernal(exch.field, exch.energy, exch.Js4, spin, sim.mu_s, mesh.ngbs4, mesh.n_ngbs4,
               T(1); ndrange=N)
    end

    return nothing
end

function effective_field(exch::SpatialHeisenberg, sim::AtomisticSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh

    kernal = atomistic_spatial_exchange_kernel!(default_backend[], groupsize[])
    kernal(exch.field, exch.energy, exch.Js, spin, sim.mu_s, mesh.ngbs, mesh.n_ngbs;
           ndrange=N)
    return nothing
end

function effective_field(dmi::HeisenbergDMI, sim::AtomisticSim, spin::AbstractArray{T,1},
                         t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh

    kernal = atomistic_dmi_kernel!(default_backend[], groupsize[])
    kernal(dmi.field, dmi.energy, dmi.Dij, spin, sim.mu_s, mesh.ngbs, mesh.n_ngbs;
           ndrange=N)

    return nothing
end

function effective_field(dmi::SpatialHeisenbergDMI, sim::AtomisticSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh

    kernal = atomistic_spatial_dmi_kernel!(default_backend[], groupsize[])
    kernal(dmi.field, dmi.energy, dmi.Dij, spin, sim.mu_s, mesh.ngbs, mesh.n_ngbs;
           ndrange=N)

    return nothing
end

function effective_field(dmi::HeisenbergCantedDMI, sim::AtomisticSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh

    kernal = atomistic_canted_dmi_kernel!(default_backend[], groupsize[])
    kernal(dmi.field, dmi.energy, dmi.Dij, spin, sim.mu_s, mesh.ngbs, mesh.n_ngbs;
           ndrange=(mesh.nx, mesh.ny, mesh.nz))

    return nothing
end

function effective_field(dmi::HeisenbergTubeBulkDMI, sim::AtomisticSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh

    kernal = tube_bulk_dmi_kernel!(default_backend[], groupsize[])
    kernal(dmi.field, dmi.energy, dmi.D, dmi.Dij, spin, sim.mu_s, mesh.ngbs, mesh.n_ngbs,
           mesh.nr; ndrange=N)
    
    return nothing
end

function effective_field(stochastic::StochasticField, sim::AtomisticSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    integrator = sim.driver.integrator

    if integrator.nsteps > stochastic.nsteps
        randn!(stochastic.eta)
        stochastic.nsteps = integrator.nsteps
    end

    mu0 = 4 * pi * 1e-7
    dt = integrator.step
    gamma = sim.driver.gamma
    alpha = sim.driver.alpha
    k_B = stochastic.k_B
    factor = 2 * mu0 * alpha * k_B / (gamma * dt) * stochastic.scaling_fun(t)

    stochastic_field_kernel!(default_backend[], groupsize[])(spin, stochastic.field,
                                                             stochastic.energy, sim.mu_s,
                                                             stochastic.eta,
                                                             stochastic.temperature, factor,
                                                             T(1); ndrange=N)

    return nothing
end
