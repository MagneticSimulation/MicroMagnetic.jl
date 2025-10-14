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
    if zee.is_scalar
       tx = zee.time_fun(t)
       zee.time_fx = tx
       zee.time_fy = tx
       zee.time_fz = tx
    else
       tx, ty, tz = zee.time_fun(t)
       zee.time_fx = tx
       zee.time_fy = ty
       zee.time_fz = tz
    end

    kernel = time_zeeman_kernel!(default_backend[], groupsize[])
    kernel(spin, zee.field, zee.init_field, zee.energy, sim.mu_s, T(1), T(zee.time_fx), T(zee.time_fy), T(zee.time_fz);
           ndrange=N)

    return nothing
end

function effective_field(anis::Anisotropy, sim::AtomisticSim, spin::AbstractArray{T,1},
                         t::Float64; output=nothing) where {T<:AbstractFloat}
    N = sim.n_total
    axis = anis.axis

    heff = output === nothing ? anis.field : output
    kernal = anisotropy_kernel!(default_backend[])
    kernal(spin, heff, anis.energy, anis.Ku, axis[1], axis[2], axis[3], sim.mu_s,
           T(1); ndrange=N)

    return nothing
end

function effective_field(anis::HexagonalAnisotropy, sim::AtomisticSim, spin::AbstractArray{T,1},
                        t::Float64) where {T<:AbstractFloat}
    N = sim.n_total

    hexagonal_anisotropy_kernel!(default_backend[])(spin, anis.field, anis.energy, anis.K1, 
              anis.K2, anis.K3, sim.mu_s, T(1); ndrange=N)

    return nothing
end

function effective_field(anis::CubicAnisotropy, sim::AtomisticSim, spin::AbstractArray{T,1},
                        t::Float64; output=nothing) where {T<:AbstractFloat}
    N = sim.n_total
    
    a1x, a1y, a1z = anis.axis1
    a2x, a2y, a2z = anis.axis2
    a3x, a3y, a3z = anis.axis3
    heff = output === nothing ? anis.field : output
    cubic_anisotropy_kernel!(default_backend[])(spin, heff, anis.energy, anis.Kc,
                              T(a1x), T(a1y), T(a1z), T(a2x), T(a2y),
                              T(a2z), T(a3x), T(a3y), T(a3z), sim.mu_s,
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
                         spin::AbstractArray{T,1}, t::Float64; output=nothing) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh

    heff = output == nothing ? exch.field : output

    # The exchange interaction for nearest neighbours
    kernal = atomistic_exchange_kernel!(default_backend[])
    kernal(heff, exch.energy, exch.Js1, spin, sim.mu_s, mesh.ngbs, mesh.n_ngbs, T(0);
           ndrange=N)

    # The exchange interaction for next-nearest neighbours
    if hasproperty(mesh, :n_ngbs2) && length(exch.Js2) == mesh.n_ngbs2
        kernal(heff, exch.energy, exch.Js2, spin, sim.mu_s, mesh.ngbs2, mesh.n_ngbs2,
               T(1); ndrange=N)
    end

    # The exchange interaction for next-next-nearest neighbours
    if hasproperty(mesh, :n_ngbs3) && length(exch.Js3) == mesh.n_ngbs3
        kernal(heff, exch.energy, exch.Js3, spin, sim.mu_s, mesh.ngbs3, mesh.n_ngbs3,
               T(1); ndrange=N)
    end

    # The exchange interaction for next-next-next-nearest neighbours
    if hasproperty(mesh, :n_ngbs4) && length(exch.Js4) == mesh.n_ngbs4
        kernal(heff, exch.energy, exch.Js4, spin, sim.mu_s, mesh.ngbs4, mesh.n_ngbs4,
               T(1); ndrange=N)
    end

    return nothing
end

function effective_field(exch::BiquadraticExchange, sim::AtomisticSim,
                         spin::AbstractArray{T,1}, t::Float64) where {T<:AbstractFloat}
    N = sim.n_total
    mesh = sim.mesh

    kernel! = atomistic_exchange_bq_kernel!(default_backend[])
    kernel!(exch.field, exch.energy, exch.Ks, spin, sim.mu_s, mesh.ngbs, mesh.n_ngbs;
            ndrange=N)
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
    scaling_factor = stochastic.scaling_fun(t)
    factor::T = 2 * alpha * k_B / (gamma * dt) * scaling_factor

    stochastic.scaling_factor = scaling_factor

    stochastic_field_kernel!(default_backend[], groupsize[])(spin, stochastic.field,
                                                             stochastic.energy, sim.mu_s,
                                                             stochastic.eta,
                                                             stochastic.temperature, factor,
                                                             T(1); ndrange=N)

    return nothing
end
