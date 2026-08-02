export SaddlePoint, find_saddle

# Geodesic minimum-mode following for spin systems follows the SPS method of
# G. P. Müller et al., Phys. Rev. Lett. 121, 197202 (2018),
# https://doi.org/10.1103/PhysRevLett.121.197202.

"""
    SaddlePoint

Result of a minimum-mode-following search. `first_order` is true only when the
final state is stationary and its checked Hessian spectrum contains exactly
one negative mode. Energy, Hessian eigenvalues, residuals, and the projected
energy gradient are stored in joules.
"""
struct SaddlePoint
    spin::Vector{Float64}
    energy::Float64
    unstable_mode::Vector{Float64}
    eigenvalues::Vector{Float64}
    residuals::Vector{Float64}
    maximum_gradient::Float64
    negative_mode_count::Int
    converged::Bool
    first_order::Bool
    steps::Int
    history::Vector{NamedTuple}
end

struct _RelaxedTransitionState
    spin::Vector{Float64}
    energy::Float64
    maximum_gradient::Float64
    maximum_dmdt::Float64
    converged::Bool
    trajectory::Matrix{Float64}
end

function _select_seed_mode(initial_mode, modes::HessianModes)
    seed = Float64.(Array(initial_mode))
    norm(seed) > 0 || error("The initial mode is zero.")
    seed ./= norm(seed)
    overlaps = [dot(seed, view(modes.eigenvectors, :, j))
                for j in axes(modes.eigenvectors, 2)]
    index = argmax(abs.(overlaps))
    mode = copy(view(modes.eigenvectors, :, index))
    overlaps[index] < 0 && (mode .*= -1)
    return mode, index, abs(overlaps[index])
end

function _track_transition_mode(previous_mode, modes, previous_spin, spin, backend)
    transported = similar(spin)
    _transport_tangent!(transported, previous_mode, previous_spin, spin, backend)
    transported ./= norm(transported)
    overlaps = Float64[]
    candidates = Any[]
    for column in axes(modes.eigenvectors, 2)
        candidate = _device_vector(view(modes.eigenvectors, :, column))
        _normalize_tangent!(candidate, spin, backend)
        push!(candidates, candidate)
        push!(overlaps, Float64(dot(transported, candidate)))
    end
    index = argmax(abs.(overlaps))
    mode = candidates[index]
    overlaps[index] < 0 && (mode .*= -1)
    return mode, index, abs(overlaps[index])
end

# Compute the MMF climbing force at a spin: the gradient with its component
# along the unstable mode inverted, ``f = -g + 2 <g, mode> mode``. This is the
# same formula used by the climbing image in GNEB and by Spirit's MMF method.
# The result is projected onto the spin tangent space and returned in a fresh
# device vector. Used by the LBFGS Atlas saddle solver below.
function _mmf_climbing_force(backend, spin, gradient, mode)
    force = similar(gradient)
    @. force = -gradient
    parallel = Float64(dot(gradient, mode))
    @. force = force + 2 * parallel * mode
    _project_tangent!(force, spin, backend)
    return force
end

# One Spirit-style LBFGS Atlas step for the MMF climbing force.
function _mmf_lbfgs_step!(state::_LBFGSAtlasState, backend, spin, gradient,
                          mode, maximum_step)
    climbing_force = _mmf_climbing_force(backend, spin, gradient, mode)
    updated = copy(spin)
    state.maxmove = Float64(maximum_step)
    step = _atlas_step!(state, backend, reshape(updated, :, 1),
                        reshape(climbing_force, :, 1);
                        energy_scale=meV, stationary=falses(1))
    return step, updated
end

function _find_saddle(backend::_TransitionBackend, minimum_spin, initial_mode;
                      direction=1, maximum_steps=4000, step_size=0.03,
                      minimum_step=1e-10, step_shrink=0.5, step_growth=1.1,
                      maximum_stagnant_steps=100,
                      gradient_tolerance=5e-27,
                      eigenvalue_tolerance=1e-28,
                      minimum_mode_overlap=0.35,
                      transverse_relaxation=0.0, mode_update_every=20,
                      log_every=20,
                      n_modes=6, finite_difference_step=2e-3,
                      krylov_dimension=max(28, n_modes + 2),
                      maximum_restarts=5, residual_tolerance=1e-28,
                      relative_residual_tolerance=5e-3,
                      breakdown_tolerance=1e-30, random_seed=314159,
                      symmetry=nothing, on_step=nothing,
                      # Solver selection for the climbing phase.
                      solver::Symbol=:lbfgs_atlas,
                      tangent_basis::Symbol=:full_2N)
    direction in (-1, 1) || error("direction must be -1 or 1.")
    maximum_steps > 0 || error("maximum_steps must be positive.")
    mode_update_every > 0 || error("mode_update_every must be positive.")
    log_every > 0 || error("log_every must be positive.")
    0 < step_shrink < 1 || error("step_shrink must be between 0 and 1.")
    step_growth >= 1 || error("step_growth must be at least 1.")
    _check_solver_symbol(solver)
    solver in (:lbfgs_atlas, :fire) ||
        error("MMF solver must be :lbfgs_atlas or :fire.")
    effective_solver = solver

    spin = _device_vector(minimum_spin)
    if !isnothing(symmetry)
        source = copy(spin)
        _project_transition_spin!(spin, source, symmetry, backend)
    end
    previous_spin = copy(spin)
    previous_mode = nothing
    previous_modes = nothing
    selected_mode = _device_vector(initial_mode)
    selected_eigenvalues = Float64[]
    selected_residuals = Float64[]
    history = NamedTuple[]
    converged = false
    maximum_gradient = Inf
    cached_modes = nothing
    cached_index = 0
    force_mode_refresh = false
    best_negative_gradient = Inf
    best_negative_spin = copy(spin)
    best_negative_mode = copy(selected_mode)
    last_negative_improvement = 0
    completed_steps = 0
    fire_velocity = create_zeros(length(spin))
    fire_time_step = step_size
    fire_alpha = 0.1
    fire_positive_steps = 0
    climbing_solver_state = Ref{_TransitionSolverState}()
    climbing_solver_initialized = false
    maximum_fire_step = max(step_size, 0.05)

    for step in 1:maximum_steps
        refresh_mode = force_mode_refresh || step == 1 ||
                       (step - 1) % mode_update_every == 0
        mode, index, overlap, modes = if refresh_mode
            starts = isnothing(previous_modes) ? initial_mode : previous_modes
            refreshed = _compute_hessian_modes(
                backend, spin; n_modes=n_modes, initial_modes=starts,
                finite_difference_step=finite_difference_step,
                krylov_dimension=krylov_dimension,
                maximum_restarts=maximum_restarts,
                residual_tolerance=residual_tolerance,
                relative_residual_tolerance=relative_residual_tolerance,
                breakdown_tolerance=breakdown_tolerance,
                random_seed=random_seed + step - 1,
                symmetry=symmetry,
                tangent_basis=tangent_basis)
            refreshed_mode, refreshed_index, refreshed_overlap =
                if isnothing(previous_mode)
                    host_mode, host_index, host_overlap =
                        _select_seed_mode(initial_mode, refreshed)
                    (_device_vector(host_mode), host_index, host_overlap)
                else
                    _track_transition_mode(
                        previous_mode, refreshed, previous_spin, spin, backend)
                end
            cached_modes = refreshed
            cached_index = refreshed_index
            force_mode_refresh = false
            (refreshed_mode, refreshed_index, refreshed_overlap, refreshed)
        else
            transported = similar(previous_mode)
            _transport_tangent!(
                transported, previous_mode, previous_spin, spin, backend)
            _normalize_tangent!(transported, spin, backend)
            (transported, cached_index, 1.0, cached_modes)
        end
        selected_mode = mode
        if refresh_mode && climbing_solver_initialized &&
           effective_solver === :lbfgs_atlas
            _reset_solver_state!(climbing_solver_state[])
        end
        selected_eigenvalues = modes.eigenvalues
        selected_residuals = modes.residuals
        eigenvalue = modes.eigenvalues[index]

        gradient = create_zeros(length(spin))
        energy = _transition_energy_gradient!(gradient, backend, spin)
        if !isnothing(symmetry)
            source = copy(gradient)
            _project_transition_vector!(
                gradient, source, spin, symmetry, backend)
        end
        maximum_gradient = _maximum_site_norm(gradient, backend)
        completed_steps = step

        if eigenvalue < -eigenvalue_tolerance
            if maximum_gradient < best_negative_gradient
                best_negative_gradient = maximum_gradient
                last_negative_improvement = step
                copyto!(best_negative_spin, spin)
                best_negative_mode = copy(mode)
            elseif step - last_negative_improvement >= maximum_stagnant_steps
                copyto!(spin, best_negative_spin)
                copyto!(previous_spin, best_negative_spin)
                previous_mode = copy(best_negative_mode)
                previous_modes = nothing
                fill!(fire_velocity, 0)
                fire_time_step = max(minimum_step,
                                     fire_time_step * step_shrink)
                fire_alpha = 0.1
                fire_positive_steps = 0
                last_negative_improvement = step
                force_mode_refresh = true
                # Reset the climbing solver as well so stale history does not
                # pull the search back toward the stagnated iterate.
                if climbing_solver_initialized
                    _reset_solver_state!(climbing_solver_state[])
                end
                maximum_fire_step = max(minimum_step,
                                         maximum_fire_step * step_shrink)
                continue
            end
        else
            best_negative_gradient = Inf
            last_negative_improvement = step
        end

        converged = eigenvalue < -eigenvalue_tolerance &&
                    maximum_gradient <= gradient_tolerance
        update_step = 0.0
        search_direction = create_zeros(length(spin))
        updated = similar(spin)
        if !converged
            force = -gradient
            if eigenvalue >= -eigenvalue_tolerance
                # Positive-eigenvalue region: push away from the minimum along
                # the selected mode. This phase is MMF-specific and always
                # uses the legacy steepest-ascent + transverse-SD strategy.
                fill!(fire_velocity, 0)
                fire_time_step = step_size
                fire_alpha = 0.1
                fire_positive_steps = 0
                search_direction .= direction .* mode
                transverse_force = force .- dot(force, mode) .* mode
                transverse_scale = _maximum_site_norm(transverse_force, backend)
                transverse_scale > 0 &&
                    (search_direction .+= transverse_relaxation .*
                                          transverse_force ./ transverse_scale)
                _project_tangent!(search_direction, spin, backend)
                direction_scale =
                    _scale_by_maximum_site!(search_direction, backend)
                direction_scale > 0 ||
                    error("Saddle search direction vanished at step $step.")
                update_step = overlap < minimum_mode_overlap ?
                              0.5 * step_size : step_size
                _geodesic_step!(updated, spin, search_direction, update_step,
                                backend)
            else
                # Negative-eigenvalue (climbing) region: dispatch to the
                # requested solver. The legacy :fire path is unchanged.
                if effective_solver === :lbfgs_atlas
                    if !climbing_solver_initialized
                        climbing_solver_state[] = _atlas_state(
                            backend, reshape(spin, :, 1);
                            memory=3, maxmove=maximum_fire_step)
                        climbing_solver_initialized = true
                    end
                    update_step, updated = _mmf_lbfgs_step!(
                        climbing_solver_state[], backend, spin, gradient,
                        mode, maximum_fire_step)
                else
                    # MicroMagnetic FIRE extension.
                    search_direction .= force .- 2 * dot(force, mode) .* mode
                    _project_tangent!(search_direction, spin, backend)
                    search_direction ./= meV
                    fire_velocity .+= fire_time_step .* search_direction
                    power = Float64(dot(fire_velocity, search_direction))
                    speed = Float64(norm(fire_velocity))
                    force_norm = Float64(norm(search_direction))
                    if power > 0
                        if speed > 0 && force_norm > 0
                            @. fire_velocity =
                                (1 - fire_alpha) * fire_velocity +
                                fire_alpha * speed / force_norm * search_direction
                        end
                        fire_positive_steps += 1
                        if fire_positive_steps > 5
                            fire_time_step =
                                min(0.5, fire_time_step * step_growth)
                            fire_alpha *= 0.99
                        end
                    else
                        fill!(fire_velocity, 0)
                        fire_time_step =
                            max(minimum_step, fire_time_step * step_shrink)
                        fire_alpha = 0.1
                        fire_positive_steps = 0
                    end
                    copyto!(search_direction, fire_velocity)
                    direction_scale =
                        _maximum_site_norm(search_direction, backend)
                    if direction_scale > 0
                        update_step = min(
                            fire_time_step * direction_scale, maximum_fire_step)
                        search_direction ./= direction_scale
                        overlap < minimum_mode_overlap &&
                            (update_step *= 0.5)
                        _geodesic_step!(updated, spin, search_direction,
                                        update_step, backend)
                    else
                        update_step = 0.0
                        copyto!(updated, spin)
                    end
                end
            end
        else
            copyto!(updated, spin)
        end

        item = (step=step, energy=energy, eigenvalue=eigenvalue,
                maximum_gradient=maximum_gradient, mode_overlap=overlap,
                eigen_residual=modes.residuals[index], step_size=update_step)
        push!(history, item)
        if step == 1 || step % log_every == 0 || converged
            @info "Saddle search progress" step energy eigenvalue maximum_gradient overlap update_step solver=effective_solver
        end
        isnothing(on_step) || on_step(step, Array(spin), item)
        converged && break

        copyto!(previous_spin, spin)
        previous_mode = copy(mode)
        previous_modes = copy(modes.eigenvectors)
        if !isnothing(symmetry)
            source = copy(updated)
            _project_transition_spin!(
                updated, source, symmetry, backend)
        end
        eigenvalue < -eigenvalue_tolerance && effective_solver === :fire &&
            _transport_tangent!(
                fire_velocity, fire_velocity, spin, updated, backend)
        copyto!(spin, updated)
    end

    validation_modes = _compute_hessian_modes(
        backend, spin; n_modes=max(2, n_modes),
        initial_modes=Array(selected_mode),
        finite_difference_step=finite_difference_step,
        krylov_dimension=max(krylov_dimension, max(2, n_modes) + 2),
        maximum_restarts=maximum_restarts,
        residual_tolerance=residual_tolerance,
        relative_residual_tolerance=relative_residual_tolerance,
        breakdown_tolerance=breakdown_tolerance,
        random_seed=random_seed + maximum_steps,
        symmetry=symmetry,
        tangent_basis=tangent_basis)
    gradient = create_zeros(length(spin))
    _transition_gradient!(gradient, backend, spin)
    if !isnothing(symmetry)
        source = copy(gradient)
        _project_transition_vector!(
            gradient, source, spin, symmetry, backend)
    end
    maximum_gradient = _maximum_site_norm(gradient, backend)
    negative_count = count(<(-eigenvalue_tolerance), validation_modes.eigenvalues)
    stationary = maximum_gradient <= gradient_tolerance
    first_order = stationary && negative_count == 1
    return SaddlePoint(
        Float64.(Array(spin)), _transition_energy(backend, spin),
        Float64.(Array(selected_mode)), validation_modes.eigenvalues,
        validation_modes.residuals, maximum_gradient, negative_count,
        converged, first_order, completed_steps, history)
end

"""
    find_saddle(create_sim, minimum_spin, mode; direction=1, ...)

Follow one Hessian mode away from a minimum and then invert the force component
along the tracked minimum mode until a first-order saddle is reached.

This is the systematic saddle-point search (SPS) method introduced by G. P.
Müller et al., *Phys. Rev. Lett.* **121**, 197202 (2018),
https://doi.org/10.1103/PhysRevLett.121.197202. The geodesic minimum-mode
following implementation also follows the scalable formulation discussed by
H. Schrautzer et al., *Phys. Rev. B* **112**, 104433 (2025),
https://doi.org/10.1103/z673-hhnp.
"""
function find_saddle(create_sim, minimum_spin, mode; kwargs...)
    return _find_saddle(_TransitionBackend(create_sim), minimum_spin, mode; kwargs...)
end

function _relax_transition_state(backend::_TransitionBackend, spin;
                                 maximum_steps=20_000,
                                 stopping_dmdt=1e-3,
                                 using_time_factor=true,
                                 maximum_trajectory_frames=400)
    sim = backend.create_sim()
    _set_transition_spin!(sim, spin)
    maximum_trajectory_frames >= 2 ||
        error("maximum_trajectory_frames must be at least 2.")
    time_factor = using_time_factor ? 2.21e5 / 2 : 1.0
    dmdt_factor = using_time_factor ? (2 * pi / 360) * 1e9 : 1.0
    dm = create_zeros(3 * sim.n_total)
    driver = sim.driver
    llg_driver = driver isa LLG
    stride = max(1, cld(maximum_steps, maximum_trajectory_frames - 1))
    frames = Vector{Vector{Float64}}()
    push!(frames, Float64.(Array(sim.spin)))
    maximum_dmdt = Inf
    converged = false

    for step in 1:maximum_steps
        run_step(sim, driver)
        if llg_driver && driver.integrator isa AdaptiveRK
            copyto!(sim.prespin, sim.spin)
            copyto!(sim.spin, view(driver.integrator.y_current,
                                   1:(3 * sim.n_total)))
        end
        step_size = llg_driver ? driver.integrator.step :
                    driver.tau / time_factor
        step_size > 0 || error("Relaxation step size must be positive.")
        compute_dm!(dm, sim.prespin, sim.spin, sim.n_total)
        maximum_dmdt = Float64(maximum(dm) / step_size / dmdt_factor)
        converged = maximum_dmdt < stopping_dmdt
        if step % stride == 0 || converged || step == maximum_steps
            push!(frames, Float64.(Array(sim.spin)))
        end
        converged && break
    end

    relaxed_spin = Float64.(Array(sim.spin))
    if frames[end] != relaxed_spin
        push!(frames, relaxed_spin)
    end
    device_spin = _device_vector(relaxed_spin)
    gradient = create_zeros(length(device_spin))
    _transition_gradient!(gradient, backend, device_spin)
    return _RelaxedTransitionState(
        relaxed_spin, _transition_energy(backend, device_spin),
        _maximum_site_norm(gradient, backend), maximum_dmdt,
        converged, hcat(frames...))
end

function _relax_saddle_branches(backend::_TransitionBackend, saddle::SaddlePoint;
                                perturbation=0.03, maximum_steps=20_000,
                                stopping_dmdt=1e-3, using_time_factor=true)
    spin = _device_vector(saddle.spin)
    mode = _device_vector(saddle.unstable_mode)
    _normalize_tangent!(mode, spin, backend)
    _scale_by_maximum_site!(mode, backend)
    minus_spin, plus_spin = similar(spin), similar(spin)
    _geodesic_step!(minus_spin, spin, mode, -perturbation, backend)
    _geodesic_step!(plus_spin, spin, mode, perturbation, backend)
    minus_state = _relax_transition_state(
        backend, minus_spin; maximum_steps=maximum_steps,
        stopping_dmdt=stopping_dmdt, using_time_factor=using_time_factor)
    plus_state = _relax_transition_state(
        backend, plus_spin; maximum_steps=maximum_steps,
        stopping_dmdt=stopping_dmdt, using_time_factor=using_time_factor)
    return minus_state, plus_state
end
