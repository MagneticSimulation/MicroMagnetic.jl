# Licensed under the MicroMagnetic.jl license.

# Internal optimisers for the transition module.  The Atlas implementation
# follows Spirit 2.2.0: stereographic 2N coordinates, one LBFGS history shared
# by the complete chain, memory three, an RMS rotation cap of 0.05, and chart
# changes at s_z*a_3 < -0.6.
#
# Spirit reference: G. P. Müller et al., Phys. Rev. B 99, 224414 (2019),
# https://doi.org/10.1103/PhysRevB.99.224414.

abstract type _TransitionSolverState end

mutable struct _VPState <: _TransitionSolverState
    time_step::Float64
    mass::Float64
end

mutable struct _FIREState <: _TransitionSolverState
    time_step::Float64
    alpha::Float64
    positive_steps::Int
end

mutable struct _LBFGSAtlasState <: _TransitionSolverState
    memory::Int
    maxmove::Float64
    local_iter::Int
    coordinates::Any
    directions::Any
    residuals::Any
    residuals_last::Any
    q::Any
    updates::Vector{Any}
    gradient_updates::Vector{Any}
    rho::Vector{Float64}
    alpha::Vector{Float64}
    factors::Any
end

const _TRANSITION_SOLVERS = (:vp, :lbfgs_atlas, :fire)

function _check_solver_symbol(solver::Symbol)
    solver in _TRANSITION_SOLVERS ||
        error("solver must be one of $_TRANSITION_SOLVERS, got :$solver.")
    return nothing
end

@kernel function _atlas_gradient_kernel!(
    residuals, @Const(images), @Const(forces), @Const(coordinates),
    @Const(active), inverse_energy_scale, n_sites)

    index = @index(Global)
    site = mod(index - 1, n_sites) + 1
    image = div(index - 1, n_sites) + 1
    j = 3 * site - 2
    q = 2 * site - 1
    if active[site]
        x, y, z = images[j, image], images[j + 1, image], images[j + 2, image]
        fx = forces[j, image] * inverse_energy_scale
        fy = forces[j + 1, image] * inverse_energy_scale
        fz = forces[j + 2, image] * inverse_energy_scale
        a3 = coordinates[site, image]
        j00 = y * y + z * (z + a3)
        j10 = -x * y
        j01 = j10
        j11 = x * x + z * (z + a3)
        j02 = -x * (z + a3)
        j12 = -y * (z + a3)
        residuals[q, image] = -(j00 * fx + j01 * fy + j02 * fz)
        residuals[q + 1, image] = -(j10 * fx + j11 * fy + j12 * fz)
    else
        residuals[q, image] = zero(eltype(residuals))
        residuals[q + 1, image] = zero(eltype(residuals))
    end
end

@kernel function _atlas_rotate_kernel!(
    images, @Const(coordinates), @Const(directions), @Const(active), n_sites)

    index = @index(Global)
    site = mod(index - 1, n_sites) + 1
    image = div(index - 1, n_sites) + 1
    if active[site]
        j = 3 * site - 2
        q = 2 * site - 1
        x, y, z = images[j, image], images[j + 1, image], images[j + 2, image]
        dx, dy = directions[q, image], directions[q + 1, image]
        a3 = coordinates[site, image]
        gamma = one(z) + z * a3
        denom = (x * x + y * y) / gamma +
                2 * (dx * x + dy * y) + gamma * (dx * dx + dy * dy)
        scale = inv(gamma + denom)
        images[j, image] = 2 * (x + dx * gamma) * scale
        images[j + 1, image] = 2 * (y + dy * gamma) * scale
        images[j + 2, image] = a3 * (gamma - denom) * scale
    end
end

@kernel function _atlas_chart_factors_kernel!(
    factors, coordinates, @Const(images), n_sites)

    index = @index(Global)
    site = mod(index - 1, n_sites) + 1
    image = div(index - 1, n_sites) + 1
    z = images[3 * site, image]
    a3 = coordinates[site, image]
    if z * a3 < zero(z)
        new_a3 = z > zero(z) ? one(z) : -one(z)
        coordinates[site, image] = new_a3
        factors[site, image] = (one(z) - new_a3 * z) / (one(z) + new_a3 * z)
    else
        factors[site, image] = one(z)
    end
end

@kernel function _atlas_metric_kernel!(
    metric, @Const(coordinates), @Const(images), n_sites)
    index = @index(Global)
    site = mod(index - 1, n_sites) + 1
    image = div(index - 1, n_sites) + 1
    metric[site, image] = images[3 * site, image] * coordinates[site, image]
end

@kernel function _atlas_scale_kernel!(vectors, @Const(factors), n_sites)
    index = @index(Global)
    site = mod(index - 1, n_sites) + 1
    image = div(index - 1, n_sites) + 1
    q = 2 * site - 1
    factor = factors[site, image]
    vectors[q, image] *= factor
    vectors[q + 1, image] *= factor
end

function _atlas_state(backend, images; memory::Int=3, maxmove::Real=0.05)
    memory > 0 || error("LBFGS Atlas memory must be positive.")
    maxmove > 0 || error("LBFGS Atlas maxmove must be positive.")
    n_sites = div(size(images, 1), 3)
    n_images = size(images, 2)
    host = Float64.(Array(images))
    coordinate_host = Matrix{Float64}(undef, n_sites, n_images)
    for image in 1:n_images, site in 1:n_sites
        coordinate_host[site, image] = host[3 * site, image] > 0 ? 1.0 : -1.0
    end
    coordinates = create_zeros(n_sites, n_images)
    _copy_transition!(coordinates, coordinate_host)
    residuals = create_zeros(2 * n_sites, n_images)
    return _LBFGSAtlasState(
        memory, Float64(maxmove), 0, coordinates,
        similar(residuals), residuals, similar(residuals), similar(residuals),
        [create_zeros(size(residuals)...) for _ in 1:memory],
        [create_zeros(size(residuals)...) for _ in 1:memory],
        zeros(memory), zeros(memory), create_zeros(n_sites, n_images))
end

function _reset_solver_state!(state::_LBFGSAtlasState)
    state.local_iter = 0
    fill!(state.directions, 0)
    fill!(state.residuals, 0)
    fill!(state.residuals_last, 0)
    fill!(state.q, 0)
    foreach(x -> fill!(x, 0), state.updates)
    foreach(x -> fill!(x, 0), state.gradient_updates)
    fill!(state.rho, 0)
    fill!(state.alpha, 0)
    return state
end

_reset_solver_state!(state::_VPState) = state
function _reset_solver_state!(state::_FIREState)
    state.alpha = 0.1
    state.positive_steps = 0
    return state
end

function _atlas_residuals!(state, backend, images, forces, energy_scale)
    n_sites = div(size(images, 1), 3)
    kernel! = _atlas_gradient_kernel!(get_backend(images), groupsize[])
    kernel!(state.residuals, images, forces, state.coordinates,
            backend.active, inv(Float64(energy_scale)), n_sites;
            ndrange=n_sites * size(images, 2))
    return state.residuals
end

function _atlas_first_direction!(state)
    copyto!(state.residuals_last, state.residuals)
    @. state.directions = -state.residuals
    state.local_iter = 1
    return state.directions
end

function _atlas_search_direction!(state::_LBFGSAtlasState)
    state.local_iter == 0 && return _atlas_first_direction!(state)
    slot = mod(state.local_iter - 1, state.memory) + 1
    copyto!(state.updates[slot], state.directions)
    @. state.gradient_updates[slot] = state.residuals - state.residuals_last
    inverse_rho = Float64(dot(state.gradient_updates[slot], state.updates[slot]))
    if !(inverse_rho > 1e-300)
        _reset_solver_state!(state)
        return _atlas_first_direction!(state)
    end
    state.rho[slot] = inv(inverse_rho)
    copyto!(state.q, state.residuals)
    count = min(state.local_iter, state.memory)
    recent = [mod(slot - k - 1, state.memory) + 1 for k in 0:(count - 1)]
    for index in recent
        state.alpha[index] = state.rho[index] * Float64(dot(state.updates[index], state.q))
        @. state.q -= state.alpha[index] * state.gradient_updates[index]
    end
    yy = Float64(dot(state.gradient_updates[slot], state.gradient_updates[slot]))
    gamma = yy * state.rho[slot] > 1e-300 ? inv(yy * state.rho[slot]) : 1.0
    @. state.directions = gamma * state.q
    for index in reverse(recent)
        beta = state.rho[index] * Float64(dot(state.gradient_updates[index], state.directions))
        @. state.directions += (state.alpha[index] - beta) * state.updates[index]
    end
    copyto!(state.residuals_last, state.residuals)
    @. state.directions = -state.directions
    state.local_iter += 1
    return state.directions
end

function _atlas_transform_charts!(state, backend, images)
    n_sites = size(state.coordinates, 1)
    n_total = n_sites * size(state.coordinates, 2)
    metric! = _atlas_metric_kernel!(get_backend(images), groupsize[])
    metric!(state.factors, state.coordinates, images, n_sites; ndrange=n_total)
    needs_change = Float64(minimum(state.factors)) < -0.6
    needs_change || return false
    factors! = _atlas_chart_factors_kernel!(get_backend(images), groupsize[])
    factors!(state.factors, state.coordinates, images, n_sites; ndrange=n_total)
    scale! = _atlas_scale_kernel!(get_backend(images), groupsize[])
    for vectors in (state.directions, state.residuals_last,
                    state.updates..., state.gradient_updates...)
        scale!(vectors, state.factors, n_sites; ndrange=n_total)
    end
    for index in eachindex(state.rho)
        inverse_rho = Float64(dot(state.updates[index],
                                  state.gradient_updates[index]))
        state.rho[index] = inverse_rho > 1e-300 ? inv(inverse_rho) : 0.0
    end
    return true
end

function _atlas_step!(state::_LBFGSAtlasState, backend, images, forces;
                      energy_scale::Real=meV, stationary=falses(size(images, 2)))
    _atlas_residuals!(state, backend, images, forces, energy_scale)
    for image in eachindex(stationary)
        stationary[image] && fill!(view(state.residuals, :, image), 0)
    end
    _atlas_search_direction!(state)
    for image in eachindex(stationary)
        stationary[image] && fill!(view(state.directions, :, image), 0)
    end
    n_sites = div(size(images, 1), 3)
    rms = 0.0
    for image in axes(images, 2)
        rms = max(rms, sqrt(Float64(dot(view(state.directions, :, image),
                                        view(state.directions, :, image))) / n_sites))
    end
    scaling = rms > state.maxmove ? state.maxmove / rms : 1.0
    state.directions .*= scaling
    kernel! = _atlas_rotate_kernel!(get_backend(images), groupsize[])
    kernel!(images, state.coordinates, state.directions, backend.active, n_sites;
            ndrange=n_sites * size(images, 2))
    _atlas_transform_charts!(state, backend, images)
    return scaling * rms
end

function _init_gneb_solver_state(solver, band, time_step, mass, maximum_step)
    _check_solver_symbol(solver)
    solver === :vp && return _VPState(Float64(time_step), Float64(mass))
    solver === :fire && return _FIREState(min(Float64(time_step), 0.05), 0.1, 0)
    return _atlas_state(band.backend, band.images; memory=3, maxmove=0.05)
end

function _gneb_lbfgs_step!(band, state::_LBFGSAtlasState)
    stationary = copy(band.stationary)
    band.climbing_image > 0 && (stationary[band.climbing_image] = false)
    return _atlas_step!(state, band.backend, band.images, band.forces;
                        energy_scale=band.energy_scale, stationary=stationary)
end
