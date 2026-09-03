export TransitionPath, TransitionResult, find_transitions, plot_transition_paths

# The automatic exploration workflow implements the systematic saddle-point
# search strategy demonstrated for skyrmion escape, collapse, and duplication
# by G. P. Müller et al., Phys. Rev. Lett. 121, 197202 (2018),
# https://doi.org/10.1103/PhysRevLett.121.197202.

"""
    TransitionPath

A validated GNEB path. `energies`, `forward_barrier`, and `reverse_barrier`
are in joules. `reaction_coordinate` is the normalized geodesic path length
and is dimensionless. All three force values are the maximum absolute GNEB-force
site norm over every active spin and image, in joules.
"""
struct TransitionPath
    images::Matrix{Float64}
    energies::Vector{Float64}
    reaction_coordinate::Vector{Float64}
    saddle_index::Int
    forward_barrier::Float64
    reverse_barrier::Float64
    maximum_force::Float64
    ordinary_maximum_force::Float64
    climbing_maximum_force::Float64
    ordinary_converged::Bool
    climbing_converged::Bool
    saddle::SaddlePoint
    mode_index::Int
    direction::Int
    valid::Bool
    zero_vector_count::Int
    maximum_norm_error::Float64
end

"""
    TransitionResult

Result of `find_transitions`. `modes` contains the Hessian eigensystem used for
the search, `transitions` contains only validated paths, and `attempts` records
the outcome of every explored mode and direction.
"""
struct TransitionResult
    modes::HessianModes
    transitions::Vector{TransitionPath}
    attempts::Vector{NamedTuple}
end

"""
    plot_transition_paths(
        paths; labels=nothing, energy_unit=meV, energy_label=nothing,
        output=nothing, smooth=true)

Plot one or more transition energy profiles. `paths` may be a
`TransitionResult`, a vector of `TransitionPath` objects, a dictionary, or a
vector of `label => path` pairs. Following the current MicroMagnetic NEB
example, energies are displayed relative to the minimum energy of each path.
`energy_unit` is expressed in joules and defaults to `meV`; the y-axis label is
updated consistently for J, eV, or meV. `energy_label` overrides the complete
y-axis label. Shape-preserving interpolation affects only the line shown in the
figure, while markers and output tables retain the computed GNEB images and
energies. Pass `smooth=false` to connect raw points.
"""
function plot_transition_paths end

@kernel function _safe_slerp_kernel!(output, @Const(spin_a), @Const(spin_b),
                                     fraction, @Const(active), angle_tolerance)
    site = @index(Global)
    j = 3 * site - 2
    if !active[site]
        output[j] = spin_a[j]
        output[j + 1] = spin_a[j + 1]
        output[j + 2] = spin_a[j + 2]
    else
        an = sqrt(spin_a[j]^2 + spin_a[j + 1]^2 + spin_a[j + 2]^2)
        bn = sqrt(spin_b[j]^2 + spin_b[j + 1]^2 + spin_b[j + 2]^2)
        ax, ay, az = spin_a[j] / an, spin_a[j + 1] / an, spin_a[j + 2] / an
        bx, by, bz = spin_b[j] / bn, spin_b[j + 1] / bn, spin_b[j + 2] / bn
        cosine = clamp(ax * bx + ay * by + az * bz, -one(an), one(an))
        angle = acos(cosine)
        x, y, z = ax, ay, az
        if angle <= angle_tolerance
            x, y, z = ax, ay, az
        elseif pi - angle <= angle_tolerance
            ux, uy, uz = if abs(ax) <= abs(ay) && abs(ax) <= abs(az)
                (one(ax), zero(ax), zero(ax))
            elseif abs(ay) <= abs(az)
                (zero(ax), one(ax), zero(ax))
            else
                (zero(ax), zero(ax), one(ax))
            end
            radial = ux * ax + uy * ay + uz * az
            ux, uy, uz = ux - radial * ax, uy - radial * ay, uz - radial * az
            un = sqrt(ux^2 + uy^2 + uz^2)
            ux, uy, uz = ux / un, uy / un, uz / un
            x = cospi(fraction) * ax + sinpi(fraction) * ux
            y = cospi(fraction) * ay + sinpi(fraction) * uy
            z = cospi(fraction) * az + sinpi(fraction) * uz
        else
            denominator = sin(angle)
            wa = sin((one(fraction) - fraction) * angle) / denominator
            wb = sin(fraction * angle) / denominator
            x, y, z = wa * ax + wb * bx, wa * ay + wb * by, wa * az + wb * bz
        end
        local_norm = sqrt(x^2 + y^2 + z^2)
        output[j] = x / local_norm
        output[j + 1] = y / local_norm
        output[j + 2] = z / local_norm
    end
end

@kernel function _path_validation_kernel!(norm_errors, zero_flags, finite_flags,
                                          @Const(path), @Const(active),
                                          site_count, norm_tolerance)
    index = @index(Global)
    site = mod(index - 1, site_count) + 1
    image = div(index - 1, site_count) + 1
    if active[site]
        j = 3 * site - 2
        x, y, z = path[j, image], path[j + 1, image], path[j + 2, image]
        finite = isfinite(x) && isfinite(y) && isfinite(z)
        local_norm = finite ? sqrt(x^2 + y^2 + z^2) : zero(eltype(path))
        norm_errors[index] = finite ? abs(local_norm - one(local_norm)) :
                             typemax(eltype(norm_errors))
        zero_flags[index] = finite && local_norm <= norm_tolerance ?
                            one(eltype(zero_flags)) : zero(eltype(zero_flags))
        finite_flags[index] = finite ? one(eltype(finite_flags)) :
                              zero(eltype(finite_flags))
    else
        norm_errors[index] = zero(eltype(norm_errors))
        zero_flags[index] = zero(eltype(zero_flags))
        finite_flags[index] = one(eltype(finite_flags))
    end
end

function _path_validation(path, active; norm_tolerance=1e-8)
    site_count = length(active)
    size(path, 1) == 3 * site_count ||
        throw(DimensionMismatch("Path row count and active mask disagree."))
    image_count = size(path, 2)
    device_path = create_zeros(size(path)...)
    _copy_transition!(device_path, path)
    active_device = kernel_array(Vector{Bool}(active))
    entry_count = site_count * image_count
    norm_errors = create_zeros(entry_count)
    zero_flags = create_zeros(entry_count)
    finite_flags = create_zeros(entry_count)
    kernel! = _path_validation_kernel!(get_backend(device_path), groupsize[])
    kernel!(norm_errors, zero_flags, finite_flags, device_path, active_device,
            site_count, eltype(device_path)(norm_tolerance); ndrange=entry_count)
    finite = Int(round(sum(finite_flags))) == entry_count
    zero_count = Int(round(sum(zero_flags)))
    maximum_error = Float64(maximum(norm_errors))
    return (valid=finite && zero_count == 0 && maximum_error <= norm_tolerance,
            finite=finite, zero_vector_count=zero_count,
            maximum_norm_error=maximum_error)
end

function _safe_interpolate_device(
    backend::_TransitionBackend, a, b, fraction;
    angle_tolerance=1e-10)

    output = create_zeros(length(a))
    kernel! = _safe_slerp_kernel!(get_backend(output), groupsize[])
    kernel!(output, a, b, eltype(a)(fraction), backend.active,
            eltype(a)(angle_tolerance); ndrange=length(backend.active_cpu))
    return output
end

function _safe_interpolate(backend::_TransitionBackend, spin_a, spin_b, fraction;
                           angle_tolerance=1e-10)
    length(spin_a) == length(spin_b) ||
        throw(DimensionMismatch("Spin arrays must have equal length."))
    0 <= fraction <= 1 || error("fraction must be in [0, 1].")
    a, b = _device_vector(spin_a), _device_vector(spin_b)
    endpoints = hcat(Float64.(Array(a)), Float64.(Array(b)))
    validation = _path_validation(endpoints, backend.active_cpu)
    validation.valid ||
        error("Cannot interpolate invalid endpoint spins: $validation")
    return _safe_interpolate_device(
        backend, a, b, fraction; angle_tolerance=angle_tolerance)
end

function _safe_interpolated_path(
    backend::_TransitionBackend, spin_a, spin_b, image_count;
    angle_tolerance=1e-10)

    image_count >= 3 || error("A path needs at least three images.")
    length(spin_a) == length(spin_b) ||
        throw(DimensionMismatch("Spin arrays must have equal length."))
    a, b = _device_vector(spin_a), _device_vector(spin_b)
    endpoints = hcat(Float64.(Array(a)), Float64.(Array(b)))
    validation = _path_validation(endpoints, backend.active_cpu)
    validation.valid ||
        error("Cannot interpolate invalid endpoint spins: $validation")
    path = zeros(Float64, length(a), image_count)
    path[:, 1] .= view(endpoints, :, 1)
    path[:, end] .= view(endpoints, :, 2)
    for image in 2:(image_count - 1)
        fraction = (image - 1) / (image_count - 1)
        interpolated = _safe_interpolate_device(
            backend, a, b, fraction; angle_tolerance=angle_tolerance)
        path[:, image] .= Array(interpolated)
    end
    full_validation = _path_validation(path, backend.active_cpu)
    full_validation.valid ||
        error("Interpolated path is invalid: $full_validation")
    return path
end

function _resample_transition_trajectory(backend::_TransitionBackend,
                                         trajectory, image_count)
    image_count >= 2 || error("A resampled trajectory needs at least two images.")
    size(trajectory, 2) >= 2 ||
        error("A trajectory needs at least two recorded states.")
    cumulative = zeros(Float64, size(trajectory, 2))
    for index in 2:length(cumulative)
        cumulative[index] = cumulative[index - 1] +
            _spin_rms_distance(
                view(trajectory, :, index - 1),
                view(trajectory, :, index), backend.active_cpu)
    end
    cumulative[end] > 0 || error("Cannot resample a zero-length trajectory.")
    output = zeros(Float64, size(trajectory, 1), image_count)
    output[:, 1] .= view(trajectory, :, 1)
    output[:, end] .= view(trajectory, :, size(trajectory, 2))
    segment = 1
    for image in 2:(image_count - 1)
        target = cumulative[end] * (image - 1) / (image_count - 1)
        while segment < length(cumulative) - 1 &&
              cumulative[segment + 1] < target
            segment += 1
        end
        span = cumulative[segment + 1] - cumulative[segment]
        fraction = span > 0 ? (target - cumulative[segment]) / span : 0.0
        interpolated = _safe_interpolate(
            backend, view(trajectory, :, segment),
            view(trajectory, :, segment + 1), fraction)
        output[:, image] .= Array(interpolated)
    end
    validation = _path_validation(output, backend.active_cpu)
    validation.valid ||
        error("Resampled relaxation trajectory is invalid: $validation")
    return output
end

function _collect_transition_path(
    band::GNEB, saddle, mode_index, direction, active,
    ordinary_maximum_force, climbing_maximum_force,
    ordinary_converged, climbing_converged)
    _update_gneb!(band)
    images = gneb_images(band)
    energies = copy(band.energies)
    reaction_coordinate = copy(band.reaction_coordinate)
    reaction_coordinate[end] > 0 &&
        (reaction_coordinate ./= reaction_coordinate[end])
    saddle_index = argmax(energies)
    validation = _path_validation(images, active)
    maximum_energy = energies[saddle_index]
    return TransitionPath(
        images, energies, reaction_coordinate, saddle_index,
        maximum_energy - energies[1], maximum_energy - energies[end],
        band.maximum_force, ordinary_maximum_force,
        climbing_maximum_force, ordinary_converged, climbing_converged,
        saddle, mode_index, direction,
        validation.valid, validation.zero_vector_count,
        validation.maximum_norm_error)
end

function _reverse_transition_path(path::TransitionPath)
    image_count = size(path.images, 2)
    reaction_coordinate =
        1 .- reverse(path.reaction_coordinate)
    return TransitionPath(
        reverse(path.images; dims=2), reverse(path.energies),
        reaction_coordinate, image_count + 1 - path.saddle_index,
        path.reverse_barrier, path.forward_barrier,
        path.maximum_force, path.ordinary_maximum_force,
        path.climbing_maximum_force, path.ordinary_converged,
        path.climbing_converged, path.saddle, path.mode_index,
        path.direction, path.valid, path.zero_vector_count,
        path.maximum_norm_error)
end

function _orient_transition_path(path::TransitionPath, reference, active)
    start_distance = _spin_rms_distance(
        view(path.images, :, 1), reference, active)
    end_distance = _spin_rms_distance(
        view(path.images, :, size(path.images, 2)), reference, active)
    return end_distance < start_distance ?
           _reverse_transition_path(path) : path
end

function _duplicate_transition_path(
    candidate::TransitionPath, transitions, active,
    endpoint_tolerance, saddle_tolerance, barrier_tolerance)

    for path in transitions
        same_orientation =
            _spin_rms_distance(
                view(candidate.images, :, 1), view(path.images, :, 1),
                active) < endpoint_tolerance &&
            _spin_rms_distance(
                view(candidate.images, :, size(candidate.images, 2)),
                view(path.images, :, size(path.images, 2)),
                active) < endpoint_tolerance
        reverse_orientation =
            _spin_rms_distance(
                view(candidate.images, :, 1),
                view(path.images, :, size(path.images, 2)),
                active) < endpoint_tolerance &&
            _spin_rms_distance(
                view(candidate.images, :, size(candidate.images, 2)),
                view(path.images, :, 1), active) < endpoint_tolerance
        (same_orientation || reverse_orientation) || continue
        same_saddle = _spin_rms_distance(
            candidate.saddle.spin, path.saddle.spin,
            active) < saddle_tolerance
        equivalent_barrier =
            barrier_tolerance > 0 &&
            abs(candidate.saddle.energy - path.saddle.energy) <
            barrier_tolerance
        (same_saddle || equivalent_barrier) && return true
    end
    return false
end

function _explore_transition_endpoints!(
    create_sim, root_minimum, root_modes, transitions, attempts,
    backend, search_keywords;
    exploration_depth, output_folder, on_attempt, n_transitions,
    distinct_minimum_tolerance, duplicate_tolerance,
    equivalent_barrier_tolerance)

    exploration_depth > 0 || return false
    n_transitions > 0 && length(transitions) >= n_transitions &&
        return true
    endpoints = Vector{Vector{Float64}}()
    for path in transitions
        endpoint = copy(view(path.images, :, size(path.images, 2)))
        _spin_rms_distance(
            endpoint, root_minimum, backend.active_cpu) >=
            distinct_minimum_tolerance || continue
        any(existing -> _spin_rms_distance(
                endpoint, existing, backend.active_cpu) <
                duplicate_tolerance, endpoints) && continue
        push!(endpoints, endpoint)
    end

    for (endpoint_offset, endpoint) in enumerate(endpoints)
        minimum_index = endpoint_offset + 1
        child_folder = isnothing(output_folder) ? nothing :
                       joinpath(
                           output_folder,
                           @sprintf("minimum_%02d_search", minimum_index))
        endpoint_modes = _compute_hessian_modes(
            backend, endpoint;
            n_modes=search_keywords.n_modes,
            finite_difference_step=
                search_keywords.finite_difference_step,
            krylov_dimension=search_keywords.krylov_dimension,
            maximum_restarts=search_keywords.maximum_restarts,
            residual_tolerance=search_keywords.residual_tolerance,
            relative_residual_tolerance=
                search_keywords.relative_residual_tolerance,
            breakdown_tolerance=search_keywords.breakdown_tolerance,
            random_seed=search_keywords.random_seed + 100_000 *
                        minimum_index,
            tangent_basis=search_keywords.tangent_basis)
        endpoint_symmetries = _detect_transition_symmetries(
            backend, (endpoint,))
        projected_blocks = Matrix{Float64}[]
        for symmetry in endpoint_symmetries
            block = _symmetry_projected_modes(
                endpoint_modes, endpoint, symmetry, backend)
            size(block, 2) > 0 && push!(projected_blocks, block)
        end
        child_initial_modes = isempty(projected_blocks) ?
                              nothing : reduce(hcat, projected_blocks)
        child_symmetry = isempty(endpoint_symmetries) ?
                         nothing : first(endpoint_symmetries)
        child = find_transitions(
            create_sim, endpoint;
            search_keywords...,
            initial_modes=child_initial_modes,
            symmetry=child_symmetry,
            exploration_depth=exploration_depth - 1,
            resume_result=nothing, on_attempt=nothing,
            n_transitions=0, output_folder=child_folder)
        for attempt in child.attempts
            push!(attempts, merge(
                attempt,
                (minimum_index=minimum_index,
                 detail="minimum $minimum_index: $(attempt.detail)")))
        end
        for child_path in child.transitions
            path = _orient_transition_path(
                child_path, root_minimum, backend.active_cpu)
            _duplicate_transition_path(
                path, transitions, backend.active_cpu,
                duplicate_tolerance, duplicate_tolerance,
                equivalent_barrier_tolerance) && continue
            push!(transitions, path)
            if !isnothing(output_folder)
                transition_folder = joinpath(
                    output_folder,
                    @sprintf("transition_%02d", length(transitions)))
                _write_transition(path, transition_folder, create_sim)
                _write_combined_energies(transitions, output_folder)
                _write_attempts(attempts, output_folder)
            end
            if !isnothing(on_attempt)
                on_attempt(TransitionResult(
                    root_modes, copy(transitions), copy(attempts)))
            end
            n_transitions > 0 &&
                length(transitions) >= n_transitions && return true
        end
    end
    return false
end

function _spin_rms_distance(spin_a, spin_b, active)
    length(spin_a) == length(spin_b) ||
        throw(DimensionMismatch("Spin arrays must have equal length."))
    distance_squared, count = 0.0, 0
    for site in eachindex(active)
        active[site] || continue
        j = 3 * site - 2
        ax, ay, az = spin_a[j], spin_a[j + 1], spin_a[j + 2]
        bx, by, bz = spin_b[j], spin_b[j + 1], spin_b[j + 2]
        sine = sqrt((ay * bz - az * by)^2 + (az * bx - ax * bz)^2 +
                    (ax * by - ay * bx)^2)
        cosine = clamp(ax * bx + ay * by + az * bz, -1.0, 1.0)
        distance_squared += atan(sine, cosine)^2
        count += 1
    end
    return sqrt(distance_squared / max(count, 1))
end

function _reaction_spacing_ratio(coordinate)
    spacing = diff(coordinate)
    isempty(spacing) && return Inf
    minimum_spacing = minimum(spacing)
    minimum_spacing > 0 || return Inf
    return maximum(spacing) / minimum_spacing
end

function _save_transition_spin(create_sim, spin, basename)
    mkpath(dirname(basename))
    sim = create_sim()
    _set_transition_spin!(sim, spin)
    save_ovf(sim, basename)
    return basename * ".ovf"
end

function _write_hessian_modes(modes, folder)
    mkpath(folder)
    open(joinpath(folder, "hessian_modes.txt"), "w") do io
        println(io, "# mode eigenvalue_J residual_J relative_residual converged")
        for index in eachindex(modes.eigenvalues)
            scale = abs(modes.eigenvalues[index])
            relative = scale > 0 ? modes.residuals[index] / scale : Inf
            @printf(io, "%d %.16e %.16e %.16e %s\n", index,
                    modes.eigenvalues[index], modes.residuals[index],
                    relative,
                    string(modes.converged[index]))
        end
    end
end

function _write_attempts(attempts, folder)
    mkpath(folder)
    open(joinpath(folder, "attempts.txt"), "w") do io
        println(io, "# minimum mode direction status detail")
        for attempt in attempts
            detail = replace(attempt.detail, '\n' => ' ')
            minimum_index = get(attempt, :minimum_index, 1)
            println(io, minimum_index, ' ', attempt.mode_index, ' ',
                    attempt.direction, ' ', attempt.status, ' ', detail)
        end
    end
end

function _write_mmf_history(saddle, folder, mode_index, direction)
    mkpath(folder)
    suffix = direction > 0 ? "plus" : "minus"
    filename = joinpath(
        folder, @sprintf("mode_%02d_%s_mmf_history.txt", mode_index, suffix))
    open(filename, "w") do io
        @printf(io, "# converged=%s first_order=%s negative_modes=%d\n",
                string(saddle.converged), string(saddle.first_order),
                saddle.negative_mode_count)
        println(io, "# step energy_J eigenvalue_J maximum_gradient_J mode_overlap eigen_residual_J step_size")
        for item in saddle.history
            @printf(io, "%d %.16e %.16e %.16e %.16e %.16e %.16e\n",
                    item.step, item.energy, item.eigenvalue,
                    item.maximum_gradient, item.mode_overlap,
                    item.eigen_residual, item.step_size)
        end
    end
    return filename
end

function _write_transition(path::TransitionPath, folder, create_sim)
    mkpath(folder)
    open(joinpath(folder, "energy.txt"), "w") do io
        print(io, "# steps")
        for image in eachindex(path.energies)
            print(io, " E_total_", image)
        end
        println(io)
        print(io, "# <>")
        for _ in eachindex(path.energies)
            print(io, " <J>")
        end
        println(io)
        @printf(io, "0")
        for energy in path.energies
            @printf(io, " %.16e", energy)
        end
        println(io)
    end
    open(joinpath(folder, "distance.txt"), "w") do io
        print(io, "# steps")
        for image in 1:(length(path.reaction_coordinate) - 1)
            print(io, " distance_", image)
        end
        println(io)
        print(io, "# <>")
        for _ in 1:(length(path.reaction_coordinate) - 1)
            print(io, " <>")
        end
        println(io)
        @printf(io, "0")
        for distance in diff(path.reaction_coordinate)
            @printf(io, " %.16e", distance)
        end
        println(io)
    end
    _save_transition_spin(create_sim, path.saddle.spin,
                          joinpath(folder, "saddle"))
    image_folder = joinpath(folder, "images")
    mkpath(image_folder)
    for image in axes(path.images, 2)
        _save_transition_spin(
            create_sim, view(path.images, :, image),
            joinpath(image_folder, @sprintf("image_%03d", image - 1)))
    end
end

function _write_combined_energies(transitions, folder)
    open(joinpath(folder, "transitions_energy.txt"), "w") do io
        println(io, "# transition mode direction image reaction_coordinate energy_J energy_above_minimum_J energy_above_minimum_meV")
        for (index, path) in enumerate(transitions)
            reference = minimum(path.energies)
            for image in eachindex(path.energies)
                relative = path.energies[image] - reference
                @printf(io, "%d %d %d %d %.16e %.16e %.16e %.16e\n",
                        index, path.mode_index, path.direction, image - 1,
                        path.reaction_coordinate[image], path.energies[image],
                        relative, relative / meV)
            end
        end
    end
end

function _validated_gneb_saddle(
    backend, band, image, original_saddle;
    n_modes, finite_difference_step, krylov_dimension,
    maximum_restarts, residual_tolerance, relative_residual_tolerance,
    breakdown_tolerance, random_seed, force_tolerance,
    eigenvalue_tolerance, tangent_basis)

    spin = Float64.(Array(view(band.images, :, image)))
    modes = _compute_hessian_modes(
        backend, spin; n_modes=max(2, n_modes),
        initial_modes=original_saddle.unstable_mode,
        finite_difference_step=finite_difference_step,
        krylov_dimension=max(krylov_dimension, max(2, n_modes) + 2),
        maximum_restarts=maximum_restarts,
        residual_tolerance=residual_tolerance,
        relative_residual_tolerance=relative_residual_tolerance,
        breakdown_tolerance=breakdown_tolerance,
        random_seed=random_seed,
        tangent_basis=tangent_basis)
    device_spin = _device_vector(spin)
    gradient = create_zeros(length(device_spin))
    energy = _transition_energy_gradient!(gradient, backend, device_spin)
    maximum_gradient = _maximum_site_norm(gradient, backend)
    negative_count = count(<(-eigenvalue_tolerance), modes.eigenvalues)
    first_order = maximum_gradient <= force_tolerance &&
                  negative_count == 1
    unstable_index = argmin(modes.eigenvalues)
    return SaddlePoint(
        spin, energy,
        copy(modes.eigenvectors[:, unstable_index]),
        modes.eigenvalues, modes.residuals, maximum_gradient,
        negative_count, first_order, first_order,
        original_saddle.steps + band.steps, original_saddle.history)
end

"""
    find_transitions(create_sim, minimum_spin; ...)

Search Hessian modes in both directions, validate first-order saddles, relax
the two branches, and refine each accepted connection with ordinary and
climbing-image GNEB. Failed candidates are recorded without stopping the
remaining search. The dedicated GNEB solver uses joules for both the
Riemannian energy gradient and `spring_constant`. Set `exploration_depth > 0`
to repeat the same mode search from newly discovered endpoint minima; any
compatible discrete symmetry is detected there and used only as a numerical
subspace constraint, never as a preassigned transition label.

The exploration strategy follows the systematic saddle-point search (SPS)
method of G. P. Müller et al., *Phys. Rev. Lett.* **121**, 197202 (2018),
https://doi.org/10.1103/PhysRevLett.121.197202. Accepted saddles are refined
with the GNEB method of P. F. Bessarab, V. M. Uzdin, and H. Jónsson,
*Comput. Phys. Commun.* **196**, 335-347 (2015),
https://doi.org/10.1016/j.cpc.2015.07.001.

When `output_folder` is set, the default output contains only the relaxed
minimum, Hessian modes, the attempts table, accepted transition folders, and a
combined energy table. Set `save_diagnostics=true` to additionally write MMF
and GNEB iteration histories under `diagnostics/`.
"""
function find_transitions(
    create_sim, minimum_spin;
    n_modes::Int=6, initial_modes=nothing, directions=(-1, 1),
    tracking_n_modes::Int=max(2, n_modes),
    tracking_mode_buffer::Int=1,
    resume_result=nothing, on_attempt=nothing,
    exploration_depth::Int=0,
    symmetry=nothing,
    n_transitions::Int=0, images::Int=21, output_folder=nothing,
    finite_difference_step::Real=2e-3,
    krylov_dimension::Int=max(28, n_modes + 2),
    tracking_krylov_dimension::Int=max(28, tracking_n_modes + 2),
    maximum_restarts::Int=5, residual_tolerance::Real=1e-28,
    relative_residual_tolerance::Real=5e-3,
    breakdown_tolerance::Real=1e-30, random_seed::Int=314159,
    saddle_maximum_steps::Int=4000, saddle_step_size::Real=0.03,
    saddle_mode_update_every::Int=20,
    gradient_tolerance::Real=5e-27,
    saddle_entry_gradient_tolerance::Real=2e-3meV,
    eigenvalue_tolerance::Real=1e-28,
    transverse_relaxation::Real=0.0,
    branch_perturbation::Real=0.03,
    branch_maximum_steps::Int=20_000,
    branch_stopping_dmdt::Real=1e-5,
    minimum_match_tolerance::Real=0.15,
    distinct_minimum_tolerance::Real=0.05,
    duplicate_tolerance::Real=0.05,
    equivalent_barrier_tolerance::Real=0.1meV,
    maximum_spacing_ratio::Real=3.0,
    saddle_energy_tolerance::Real=0.1meV,
    spring_constant::Real=meV,
    ordinary_maximum_steps::Int=20_000,
    ordinary_force_tolerance::Real=1e-4meV,
    climbing_maximum_steps::Int=20_000,
    climbing_force_tolerance=nothing,
    gneb_time_step::Real=0.1,
    gneb_mass::Real=1.0,
    gneb_maximum_step::Real=0.03,
    require_gneb_convergence::Bool=true,
    save_data_every::Int=20,
    save_diagnostics::Bool=false,
    using_time_factor::Bool=true,
    # Spirit-aligned defaults; FIRE remains a MicroMagnetic extension.
    gneb_solver::Symbol=:vp,
    saddle_solver::Symbol=:lbfgs_atlas,
    spring_force_ratio::Real=0.0,
    path_shortening_constant::Real=0.0,
    tangent_basis::Symbol=:full_2N)

    images >= 5 || error("images must be at least 5.")
    all(direction -> direction in (-1, 1), directions) ||
        error("directions may contain only -1 and 1.")
    tracking_mode_buffer >= 0 ||
        error("tracking_mode_buffer must be nonnegative.")
    exploration_depth >= 0 ||
        error("exploration_depth must be nonnegative.")
    backend = _TransitionBackend(create_sim)
    active_spins = max(count(backend.active_cpu), 1)
    resolved_climbing_force_tolerance = isnothing(climbing_force_tolerance) ?
        min(5e-5meV,
            0.5 * Float64(saddle_energy_tolerance) / active_spins) :
        Float64(climbing_force_tolerance)
    resolved_climbing_force_tolerance > 0 ||
        error("climbing_force_tolerance must be positive.")
    isnothing(climbing_force_tolerance) && @info(
        "Automatic climbing convergence tolerance",
        maximum_site_force_J=resolved_climbing_force_tolerance,
        maximum_site_force_meV=resolved_climbing_force_tolerance / meV,
        active_spins=active_spins)
    spring_constant > 0 ||
        error("spring_constant must be positive.")
    minimum = Float64.(Array(minimum_spin))
    !isnothing(resume_result) && !(resume_result isa TransitionResult) &&
        error("resume_result must be a TransitionResult.")
    modes = if !isnothing(resume_result)
        resume_result.modes
    elseif isnothing(initial_modes)
        _compute_hessian_modes(
            backend, minimum; n_modes=n_modes,
            finite_difference_step=finite_difference_step,
            krylov_dimension=krylov_dimension,
            maximum_restarts=maximum_restarts,
            residual_tolerance=residual_tolerance,
            relative_residual_tolerance=relative_residual_tolerance,
            breakdown_tolerance=breakdown_tolerance,
            random_seed=random_seed, symmetry=symmetry,
            tangent_basis=tangent_basis)
    else
        _evaluate_hessian_seeds(
            backend, minimum, initial_modes;
            finite_difference_step=finite_difference_step,
            residual_tolerance=residual_tolerance,
            relative_residual_tolerance=relative_residual_tolerance,
            symmetry=symmetry)
    end
    candidate_count = size(modes.eigenvectors, 2)
    transitions = isnothing(resume_result) ?
                  TransitionPath[] : copy(resume_result.transitions)
    attempts = isnothing(resume_result) ?
               NamedTuple[] : copy(resume_result.attempts)
    completed_attempts = Set(
        (attempt.mode_index, attempt.direction) for attempt in attempts
        if get(attempt, :minimum_index, 1) == 1 &&
           !(attempt.status in (:failed, :incomplete)))

    if !isnothing(output_folder)
        mkpath(output_folder)
        _write_hessian_modes(modes, output_folder)
        _save_transition_spin(create_sim, minimum,
                              joinpath(output_folder, "minimum"))
    end

    stop_search = false
    for mode_index in 1:candidate_count
        for direction in directions
        (mode_index, direction) in completed_attempts && continue
        status, detail = :failed, ""
        try
            tracked_mode_count = min(
                tracking_n_modes,
                max(2, mode_index + tracking_mode_buffer))
            saddle = _find_saddle(
                backend, minimum, view(modes.eigenvectors, :, mode_index);
                direction=direction, maximum_steps=saddle_maximum_steps,
                step_size=saddle_step_size,
                mode_update_every=saddle_mode_update_every,
                gradient_tolerance=max(
                    gradient_tolerance, saddle_entry_gradient_tolerance),
                eigenvalue_tolerance=eigenvalue_tolerance,
                transverse_relaxation=transverse_relaxation,
                n_modes=tracked_mode_count,
                finite_difference_step=finite_difference_step,
                krylov_dimension=max(
                    tracking_krylov_dimension, tracked_mode_count + 2),
                maximum_restarts=maximum_restarts,
                residual_tolerance=residual_tolerance,
                relative_residual_tolerance=relative_residual_tolerance,
                breakdown_tolerance=breakdown_tolerance,
                random_seed=random_seed + 1000 * mode_index +
                            (direction > 0 ? 1 : 2),
                symmetry=symmetry,
                solver=saddle_solver,
                tangent_basis=tangent_basis)
            !isnothing(output_folder) && save_diagnostics &&
                _write_mmf_history(
                    saddle, joinpath(output_folder, "diagnostics"),
                    mode_index, direction)
            saddle.negative_mode_count == 1 &&
                saddle.maximum_gradient <= saddle_entry_gradient_tolerance ||
                error("candidate did not reach the index-one GNEB entry " *
                      "region (negative_modes=$(saddle.negative_mode_count), " *
                      "maximum_gradient=$(saddle.maximum_gradient) J)")

            minus_state, plus_state = _relax_saddle_branches(
                backend, saddle; perturbation=branch_perturbation,
                maximum_steps=branch_maximum_steps,
                stopping_dmdt=branch_stopping_dmdt,
                using_time_factor=using_time_factor)
            minus_distance = _spin_rms_distance(
                minimum, minus_state.spin, backend.active_cpu)
            plus_distance = _spin_rms_distance(
                minimum, plus_state.spin, backend.active_cpu)
            near_state, far_state = minus_distance <= plus_distance ?
                                    (minus_state, plus_state) :
                                    (plus_state, minus_state)
            near_distance = min(minus_distance, plus_distance)
            endpoint_distance = _spin_rms_distance(
                minus_state.spin, plus_state.spin, backend.active_cpu)
            minus_state.converged ||
                error("negative saddle branch relaxation did not converge " *
                      "(maximum_dmdt=$(minus_state.maximum_dmdt))")
            plus_state.converged ||
                error("positive saddle branch relaxation did not converge " *
                      "(maximum_dmdt=$(plus_state.maximum_dmdt))")
            near_distance <= minimum_match_tolerance ||
                error("neither saddle branch returned to the supplied minimum " *
                      "(distances=$minus_distance,$plus_distance)")
            endpoint_distance >= distinct_minimum_tolerance ||
                error("both saddle branches relaxed to the same minimum " *
                      "(distance=$endpoint_distance)")

            duplicate = any(transitions) do path
                same_endpoint = _spin_rms_distance(
                    path.images[:, end], far_state.spin,
                    backend.active_cpu) < duplicate_tolerance
                same_saddle = _spin_rms_distance(
                    path.saddle.spin, saddle.spin,
                    backend.active_cpu) < duplicate_tolerance
                equivalent_barrier =
                    equivalent_barrier_tolerance > 0 &&
                    abs(path.saddle.energy - saddle.energy) <
                    equivalent_barrier_tolerance
                return same_endpoint && (same_saddle || equivalent_barrier)
            end
            duplicate && error("transition duplicates an accepted path")

            left_distance = _spin_rms_distance(
                minimum, saddle.spin, backend.active_cpu)
            right_distance = _spin_rms_distance(
                saddle.spin, far_state.spin, backend.active_cpu)
            interior = images - 3
            left_frames = clamp(
                round(Int, interior * left_distance /
                           max(left_distance + right_distance, eps())),
                1, interior - 1)
            right_frames = interior - left_frames
            near_trajectory = hcat(saddle.spin, near_state.trajectory)
            far_trajectory = hcat(saddle.spin, far_state.trajectory)
            left_trajectory = reverse(near_trajectory; dims=2)
            left_trajectory[:, 1] .= minimum
            far_trajectory[:, end] .= far_state.spin
            left_path = _resample_transition_trajectory(
                backend, left_trajectory, left_frames + 2)
            right_path = _resample_transition_trajectory(
                backend, far_trajectory, right_frames + 2)
            initial_path = hcat(
                left_path, view(right_path, :, 2:size(right_path, 2)))

            prefix = !isnothing(output_folder) && save_diagnostics ?
                     joinpath(output_folder, "diagnostics",
                              @sprintf("mode_%02d_%s_gneb", mode_index,
                                       direction > 0 ? "plus" : "minus")) : ""
            saddle_column = left_frames + 2
            band = GNEB(
                create_sim, initial_path; name=prefix * "_ordinary",
                spring_constant=spring_constant,
                stationary_images=[saddle_column])
            relax_gneb!(
                band; climbing=false,
                maximum_steps=ordinary_maximum_steps,
                force_tolerance=ordinary_force_tolerance,
                time_step=gneb_time_step,
                mass=gneb_mass,
                maximum_step=gneb_maximum_step,
                save_data_every=save_data_every,
                symmetry=symmetry,
                solver=gneb_solver,
                spring_force_ratio=spring_force_ratio,
                path_shortening_constant=path_shortening_constant)
            ordinary_maximum_force = band.maximum_force
            ordinary_converged = band.converged
            require_gneb_convergence && !ordinary_converged &&
                error("anchored ordinary GNEB did not converge " *
                      "(maximum_force=$ordinary_maximum_force J)")

            # Spirit keeps the ordinary images mobile during climbing so the
            # path tangent follows the refined saddle instead of being frozen
            # at its ordinary-band value.
            set_climbing_image!(band, saddle_column)
            climbing_entry_force = band.maximum_force
            climbing_entry_energy = maximum(band.energies)
            climbing_unstable = Ref(false)
            stop_unstable_climbing = function (current_band, item)
                force_limit = max(1e-3meV, 10 * climbing_entry_force)
                energy_limit = climbing_entry_energy + meV
                if item.maximum_force > force_limit ||
                   item.maximum_energy > energy_limit
                    climbing_unstable[] = true
                    return :stop
                end
                return nothing
            end
            climbing_solver = gneb_solver === :lbfgs_atlas ?
                              gneb_solver : :vp
            band.name = prefix * "_climbing_" * String(climbing_solver)
            relax_gneb!(
                band; climbing=true,
                maximum_steps=climbing_maximum_steps,
                force_tolerance=resolved_climbing_force_tolerance,
                time_step=gneb_time_step,
                mass=gneb_mass,
                maximum_step=gneb_maximum_step,
                save_data_every=save_data_every,
                solver=climbing_solver,
                symmetry=symmetry,
                on_step=stop_unstable_climbing,
                spring_force_ratio=spring_force_ratio,
                path_shortening_constant=path_shortening_constant)
            climbing_unstable[] &&
                error("climbing-image GNEB became unstable relative to " *
                      "its converged ordinary band")
            climbing_maximum_force = band.maximum_force
            climbing_converged = band.converged
            require_gneb_convergence && !climbing_converged &&
                error("climbing-image GNEB did not converge " *
                      "(maximum_force=$climbing_maximum_force J)")
            refined_saddle = _validated_gneb_saddle(
                backend, band, saddle_column, saddle;
                n_modes=tracked_mode_count,
                finite_difference_step=finite_difference_step,
                krylov_dimension=max(
                    tracking_krylov_dimension, tracked_mode_count + 2),
                maximum_restarts=maximum_restarts,
                residual_tolerance=residual_tolerance,
                relative_residual_tolerance=relative_residual_tolerance,
                breakdown_tolerance=breakdown_tolerance,
                random_seed=random_seed + 10_000 * mode_index +
                            (direction > 0 ? 1 : 2),
                force_tolerance=resolved_climbing_force_tolerance,
                eigenvalue_tolerance=eigenvalue_tolerance,
                tangent_basis=tangent_basis)
            !require_gneb_convergence || refined_saddle.first_order ||
                error("climbing GNEB image is not a validated first-order " *
                      "saddle (negative_modes=" *
                      "$(refined_saddle.negative_mode_count))")
            path = _collect_transition_path(
                band, refined_saddle, mode_index, direction,
                backend.active_cpu,
                ordinary_maximum_force, climbing_maximum_force,
                ordinary_converged, climbing_converged)
            path.valid || error("refined GNEB path contains invalid spins")
            spacing_ratio = _reaction_spacing_ratio(path.reaction_coordinate)
            spacing_ratio <= maximum_spacing_ratio ||
                error("refined GNEB image spacing ratio is $spacing_ratio")
            saddle_barrier = refined_saddle.energy - path.energies[1]
            barrier_error = abs(path.forward_barrier - saddle_barrier)
            barrier_error <= saddle_energy_tolerance ||
                error("climbing GNEB barrier disagrees with the MMF saddle " *
                      "by $(barrier_error / meV) meV")
            push!(transitions, path)
            status, detail = :accepted, "first-order saddle and GNEB path accepted"
            if !isnothing(output_folder)
                transition_folder = joinpath(
                    output_folder, @sprintf("transition_%02d", length(transitions)))
                _write_transition(path, transition_folder, create_sim)
                _write_combined_energies(transitions, output_folder)
            end
        catch error_value
            detail = sprint(showerror, error_value)
            incomplete = occursin("did not converge", detail) ||
                         occursin("not a validated first-order saddle", detail) ||
                         occursin("did not reach the index-one", detail) ||
                         occursin("climbing GNEB image is not", detail)
            status = occursin("duplicates", detail) ? :duplicate :
                     incomplete ? :incomplete : :rejected
        end
        push!(attempts, (minimum_index=1, mode_index=mode_index,
                          direction=direction,
                          status=status, detail=detail))
        !isnothing(output_folder) && _write_attempts(attempts, output_folder)
        if !isnothing(on_attempt)
            on_attempt(TransitionResult(modes, copy(transitions),
                                        copy(attempts)))
        end
        if n_transitions > 0 && length(transitions) >= n_transitions
            stop_search = true
            break
        end
        end
        stop_search && break
    end
    if !stop_search && exploration_depth > 0
        search_keywords = (
            n_modes=n_modes,
            directions=directions,
            tracking_n_modes=tracking_n_modes,
            tracking_mode_buffer=tracking_mode_buffer,
            images=images,
            finite_difference_step=finite_difference_step,
            krylov_dimension=krylov_dimension,
            tracking_krylov_dimension=tracking_krylov_dimension,
            maximum_restarts=maximum_restarts,
            residual_tolerance=residual_tolerance,
            relative_residual_tolerance=relative_residual_tolerance,
            breakdown_tolerance=breakdown_tolerance,
            random_seed=random_seed,
            saddle_maximum_steps=saddle_maximum_steps,
            saddle_step_size=saddle_step_size,
            saddle_mode_update_every=saddle_mode_update_every,
            gradient_tolerance=gradient_tolerance,
            saddle_entry_gradient_tolerance=
                saddle_entry_gradient_tolerance,
            eigenvalue_tolerance=eigenvalue_tolerance,
            transverse_relaxation=transverse_relaxation,
            branch_perturbation=branch_perturbation,
            branch_maximum_steps=branch_maximum_steps,
            branch_stopping_dmdt=branch_stopping_dmdt,
            minimum_match_tolerance=minimum_match_tolerance,
            distinct_minimum_tolerance=distinct_minimum_tolerance,
            duplicate_tolerance=duplicate_tolerance,
            equivalent_barrier_tolerance=
                equivalent_barrier_tolerance,
            maximum_spacing_ratio=maximum_spacing_ratio,
            saddle_energy_tolerance=saddle_energy_tolerance,
            spring_constant=spring_constant,
            ordinary_maximum_steps=ordinary_maximum_steps,
            ordinary_force_tolerance=ordinary_force_tolerance,
            climbing_maximum_steps=climbing_maximum_steps,
            climbing_force_tolerance=resolved_climbing_force_tolerance,
            gneb_time_step=gneb_time_step,
            gneb_mass=gneb_mass,
            gneb_maximum_step=gneb_maximum_step,
            require_gneb_convergence=require_gneb_convergence,
            save_data_every=save_data_every,
            using_time_factor=using_time_factor,
            gneb_solver=gneb_solver,
            saddle_solver=saddle_solver,
            spring_force_ratio=spring_force_ratio,
            path_shortening_constant=path_shortening_constant,
            tangent_basis=tangent_basis)
        stop_search = _explore_transition_endpoints!(
            create_sim, minimum, modes, transitions, attempts,
            backend, search_keywords;
            exploration_depth=exploration_depth,
            output_folder=output_folder, on_attempt=on_attempt,
            n_transitions=n_transitions,
            distinct_minimum_tolerance=distinct_minimum_tolerance,
            duplicate_tolerance=duplicate_tolerance,
            equivalent_barrier_tolerance=
                equivalent_barrier_tolerance)
    end
    return TransitionResult(modes, transitions, attempts)
end
