# implement the GPSM Scheme B from "Two improved Gauss-Seidel projection methods for Landau-Lifshitz-Gilbert equation"
# They solve the dimensionally reduced LLG equation
# dm/dt = - m x h - alpha * m x m x h
# where h = epsilon*Delta m + hs, with epsilon*Delta m the exchange field and hs include stray field, zeeman field, and anisotropy field, etc. 
# LLG equation is 
# dm/dt = - gamma_L m x H - alpha * gamma_L * m x (m x H)
# where gamma_L = gamma_G / (1 + alpha * alpha) and H is the total effective field
# Hex =2*A/(mu_0*Ms) Delta m, so we have epsilon = 2*A/(mu_0*Ms)*gamma_L and hs = gamma_L * H

mutable struct GPSM{T<:AbstractFloat} <: Integrator
    t::Float64
    step::T
    nsteps::Int64
    G::Any # I - epsilon*dt*Laplaian
    initialized::Bool
    Laplaian::AbstractSparseMatrix{T}
    gn1::AbstractArray{T,1}
    gn2::AbstractArray{T,1}
    gn3::AbstractArray{T,1}
    rhs::AbstractArray{T,1}
    exch::Union{Nothing, MicroEnergy}
    GPSM{T}() where {T<:AbstractFloat} = new()
end

function GPSM(n_total::Int64, step::Float64)
    T = Float[]
    gpsm = GPSM{T}()
    gpsm.step = step
    gpsm.nsteps = 0
    gpsm.t = 0
    gpsm.initialized = false
    gpsm.gn1 = create_zeros(n_total)
    gpsm.gn2 = create_zeros(n_total)
    gpsm.gn3 = create_zeros(n_total)
    gpsm.rhs = create_zeros(n_total)
    return gpsm
end

function initialize!(integrator::GPSM, sim::AbstractSim)
    if integrator.initialized
        return nothing
    end
    idx = findfirst(x -> isa(x, ExchangeFE) || isa(x, UniformExchange), sim.interactions)
    if idx !== nothing
        integrator.exch = splice!(sim.interactions, idx)
    else
        error("GPSM integrator requires Exchange interaction!")
    end
    integrator.Laplaian = build_exch_matrix(integrator.exch, sim)
    alpha = sim.driver.alpha
    gamma_G = sim.driver.gamma 
    dt = integrator.step * gamma_G / (1 + alpha * alpha)
    integrator.G = I - dt * integrator.Laplaian
    integrator.G = cholesky(integrator.G)
    
    integrator.initialized = true

    alpha = sim.driver.alpha
    gamma_G = sim.driver.gamma 
    dt = integrator.step * gamma_G / (1 + alpha * alpha)
    N = sim.n_total

    effective_field(sim, sim.spin, integrator.t)
    spin = reshape(sim.spin, 3, N)
    m1 = @view spin[1, :]
    m2 = @view spin[2, :]
    m3 = @view spin[3, :]

    field = reshape(sim.field, 3, N)
    f1 = @view field[1, :]
    f2 = @view field[2, :]
    f3 = @view field[3, :]
    
    rhs = integrator.rhs
    rhs .= m1 .+ dt .* f1
    x = integrator.G \ rhs
    integrator.gn1 .= x

    rhs .= m2 .+ dt .* f2
    integrator.gn2 .= integrator.G \ rhs
    
    rhs .= m3 .+ dt .* f3
    integrator.gn3 .= integrator.G \ rhs
    return nothing
end

function advance_step(sim::AbstractSim, integrator::GPSM)

    if !integrator.initialized
        initialize!(integrator, sim)
    end

    sim.prespin .= sim.spin

    alpha = sim.driver.alpha
    gamma_G = sim.driver.gamma 
    dt = integrator.step * gamma_G / (1 + alpha * alpha)

    N = sim.n_total
    spin = reshape(sim.spin, 3, N)
    m1 = @view spin[1, :]
    m2 = @view spin[2, :]
    m3 = @view spin[3, :]

    field = reshape(sim.field, 3, N)
    f1 = @view field[1, :]
    f2 = @view field[2, :]
    f3 = @view field[3, :]

    g1, g2, g3 = integrator.gn1, integrator.gn2, integrator.gn3
    rhs = integrator.rhs
    
    # compute m1star
    mg = m1.*g1 .+ m2.*g2 .+ m3.*g3
    mm = m1.*m1 .+ m2.*m2 .+ m3.*m3
    m1star = m1 .- (m2 .* g3 - m3 .* g2) - alpha .* mg .* m1 + alpha .* mm .* g1

    # update the effective field and g1
    m1 .= m1star
    effective_field(sim, sim.spin, integrator.t)
    rhs .= m1star + dt .* f1
    g1 .= integrator.G \ rhs

    # compute m2star
    mg = m1.*g1 .+ m2.*g2 .+ m3.*g3
    mm = m1.*m1 .+ m2.*m2 .+ m3.*m3
    m2star = m2 .- (m3 .* g1 - m1 .* g3) .- alpha .* mg .* m2 + alpha .* mm .* g2

    # compute the effective field and g2
    m2 .= m2star
    effective_field(sim, sim.spin, integrator.t)
    rhs .= m2star + dt .* f2
    g2 .= integrator.G \ rhs

    # compute m3star
    mg = m1.*g1 .+ m2.*g2 .+ m3.*g3
    mm = m1.*m1 .+ m2.*m2 .+ m3.*m3
    m3star = m3 .- (m1 .* g2 - m2 .* g1) .- alpha .* mg .* m3 + alpha .* mm .* g3
    
    # compute the effective field and g3
    m3 .= m3star
    effective_field(sim, sim.spin, integrator.t)
    rhs .= m3star + dt .* f3
    g3 .= integrator.G \ rhs
    
    normalise(sim.spin, sim.n_total)
    
    integrator.nsteps += 1
    integrator.t = integrator.nsteps * integrator.step
    
    return nothing
end


