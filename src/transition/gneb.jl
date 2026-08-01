export GNEB, relax_gneb!, set_climbing_image!, gneb_images, set_image_type!

# GNEB reference: P. F. Bessarab, V. M. Uzdin, and H. Jónsson,
# Comput. Phys. Commun. 196, 335-347 (2015),
# https://doi.org/10.1016/j.cpc.2015.07.001.

# Image type symbols used by the enhanced GNEB force routine. They mirror
# Spirit's GNEB_Image_Type enum:
#   :normal    – standard NEB image with spring force + perpendicular force
#   :climbing  – climbing image, force inverted along the tangent (no spring)
#   :falling   – transition image with gradient force but no spring force,
#                used to smooth the hand-off between normal and climbing
#   :stationary– image held fixed (endpoints or anchored saddle columns)
const _GNEB_IMAGE_TYPES = (:normal, :climbing, :falling, :stationary)

"""
    GNEB(create_sim, path; name="GNEB", spring_constant=meV)

Geodesic elastic band for spin systems. Unlike the legacy `NEB` simulation
type, `GNEB` works with the Riemannian energy gradient in
joules, geodesic image distances, and a dedicated manifold
velocity-projection optimizer.
The endpoints are fixed. Additional images may be held fixed with
`stationary_images`.

The force model extends the original Henkelman–Jónsson NEB force with three
optional Spirit-compatible enhancements, all enabled through `relax_gneb!`
keyword arguments and inactive by default:

  * **Path shortening** (`path_shortening_constant>0`): adds a finite-difference
    secant curvature force projected onto the plane orthogonal to both the
    tangent and the gradient. This compresses redundant path length when the
    number of images is small, matching Spirit's `GNEB::Calculate_Force`
    shortening term.
  * **Energy-weighted spring** (`spring_force_ratio ∈ [0,1]`): blends the
    pure reaction-coordinate spring with an energy-difference-weighted
    spring using cubic Hermite interpolation, matching Spirit's
    `spring_force_ratio`. ``ratio=0`` (default) recovers the pure Rx spring.
  * **Falling image** (`image_type=:falling`): an image type between normal
    and climbing that carries the gradient force but no spring force,
    smoothing the switch to climbing mode.

The geodesic force and distance follow P. F. Bessarab, V. M. Uzdin, and H.
Jónsson, *Comput. Phys. Commun.* **196**, 335-347 (2015),
https://doi.org/10.1016/j.cpc.2015.07.001. The energy-weighted tangent follows
G. Henkelman and H. Jónsson, *J. Chem. Phys.* **113**, 9978-9985 (2000),
https://doi.org/10.1063/1.1323224. Spirit-aligned extensions refer to G. P.
Müller et al., *Phys. Rev. B* **99**, 224414 (2019),
https://doi.org/10.1103/PhysRevB.99.224414.
"""
mutable struct GNEB{B,I}
    backend::B
    images::I
    gradients::I
    forces::I
    previous_forces::I
    segment_logs::I
    transported_logs::I
    tangents::I
    velocity::I
    energies::Vector{Float64}
    reaction_coordinate::Vector{Float64}
    stationary::BitVector
    image_types::Vector{Symbol}
    climbing_image::Int
    spring_constant::Float64
    energy_scale::Float64
    name::String
    converged::Bool
    steps::Int
    maximum_force::Float64
    maximum_component::Float64
    history::Vector{NamedTuple}
end

@kernel function _gneb_log_map_kernel!(
    output, @Const(spin_from), @Const(spin_to), @Const(active))

    site = @index(Global)
    j = 3 * site - 2
    if active[site]
        x1, y1, z1 = spin_from[j], spin_from[j + 1], spin_from[j + 2]
        x2, y2, z2 = spin_to[j], spin_to[j + 1], spin_to[j + 2]
        cosine = clamp(x1 * x2 + y1 * y2 + z1 * z2, -one(x1), one(x1))
        angle = acos(cosine)
        sine = sqrt(max(zero(x1), one(x1) - cosine * cosine))
        if angle <= eltype(output)(1e-12)
            output[j] = zero(x1)
            output[j + 1] = zero(x1)
            output[j + 2] = zero(x1)
        elseif sine > eltype(output)(1e-10)
            scale = angle / sine
            output[j] = scale * (x2 - cosine * x1)
            output[j + 1] = scale * (y2 - cosine * y1)
            output[j + 2] = scale * (z2 - cosine * z1)
        else
            ax, ay, az = if abs(x1) <= abs(y1) && abs(x1) <= abs(z1)
                (one(x1), zero(x1), zero(x1))
            elseif abs(y1) <= abs(z1)
                (zero(x1), one(x1), zero(x1))
            else
                (zero(x1), zero(x1), one(x1))
            end
            radial = ax * x1 + ay * y1 + az * z1
            tx, ty, tz = ax - radial * x1, ay - radial * y1,
                         az - radial * z1
            tangent_norm = sqrt(tx * tx + ty * ty + tz * tz)
            output[j] = angle * tx / tangent_norm
            output[j + 1] = angle * ty / tangent_norm
            output[j + 2] = angle * tz / tangent_norm
        end
    else
        output[j] = zero(eltype(output))
        output[j + 1] = zero(eltype(output))
        output[j + 2] = zero(eltype(output))
    end
end

@kernel function _gneb_retraction_step_kernel!(
    spin, @Const(displacement), scale, @Const(active))

    site = @index(Global)
    j = 3 * site - 2
    if active[site]
        x = spin[j] + scale * displacement[j]
        y = spin[j + 1] + scale * displacement[j + 1]
        z = spin[j + 2] + scale * displacement[j + 2]
        local_norm = sqrt(x * x + y * y + z * z)
        spin[j] = x / local_norm
        spin[j + 1] = y / local_norm
        spin[j + 2] = z / local_norm
    end
end

function _gneb_log_map!(output, spin_from, spin_to, backend)
    kernel! = _gneb_log_map_kernel!(default_backend[], groupsize[])
    kernel!(output, spin_from, spin_to, backend.active;
            ndrange=length(backend.active_cpu))
    return output
end

function GNEB(create_sim, path::AbstractMatrix; name="GNEB",
              spring_constant::Real=meV, energy_scale::Real=meV,
              stationary_images=Int[])
    size(path, 2) >= 3 || error("A GNEB path needs at least three images.")
    spring_constant > 0 || error("spring_constant must be positive.")
    energy_scale > 0 || error("energy_scale must be positive.")
    backend = _TransitionBackend(create_sim)
    validation = _path_validation(path, backend.active_cpu)
    validation.valid ||
        error("Refusing to initialize an invalid GNEB path: $validation")
    images = create_zeros(size(path)...)
    _copy_transition!(images, path)
    gradients = similar(images)
    forces = similar(images)
    previous_forces = similar(images)
    segment_logs = similar(images)
    transported_logs = similar(images)
    tangents = similar(images)
    velocity = similar(images)
    fill!(gradients, 0)
    fill!(forces, 0)
    fill!(previous_forces, 0)
    fill!(segment_logs, 0)
    fill!(transported_logs, 0)
    fill!(tangents, 0)
    fill!(velocity, 0)
    stationary = falses(size(path, 2))
    stationary[1] = true
    stationary[end] = true
    for image in stationary_images
        1 <= image <= size(path, 2) ||
            error("stationary image $image is outside the path.")
        stationary[image] = true
    end
    image_types = fill(:normal, size(path, 2))
    image_types[1] = :stationary
    image_types[end] = :stationary
    for image in 1:length(image_types)
        stationary[image] && (image_types[image] = :stationary)
    end
    band = GNEB(
        backend, images, gradients, forces, previous_forces,
        segment_logs, transported_logs, tangents, velocity,
        zeros(Float64, size(path, 2)), zeros(Float64, size(path, 2)),
        stationary, image_types, 0, Float64(spring_constant),
        Float64(energy_scale), String(name), false, 0, Inf, Inf,
        NamedTuple[])
    _update_gneb!(band; refresh_all=true)
    return band
end

function _update_gneb_energies_gradients!(band::GNEB; refresh_all=false)
    for image in axes(band.images, 2)
        !refresh_all && band.stationary[image] && continue
        spin = view(band.images, :, image)
        gradient = view(band.gradients, :, image)
        band.energies[image] =
            _transition_energy_gradient!(gradient, band.backend, spin)
    end
    return band
end

function _update_gneb_reaction_coordinate!(band::GNEB)
    fill!(band.segment_logs, 0)
    fill!(band.transported_logs, 0)
    band.reaction_coordinate[1] = 0.0
    for image in 2:length(band.reaction_coordinate)
        segment = view(band.segment_logs, :, image - 1)
        _gneb_log_map!(
            segment, view(band.images, :, image - 1),
            view(band.images, :, image), band.backend)
        band.reaction_coordinate[image] =
            band.reaction_coordinate[image - 1] +
            Float64(norm(segment))
    end
    for image in 2:(size(band.images, 2) - 1)
        _transport_tangent!(
            view(band.transported_logs, :, image),
            view(band.segment_logs, :, image - 1),
            view(band.images, :, image - 1),
            view(band.images, :, image), band.backend)
    end
    return band
end

function _gneb_tangent!(tangent, band::GNEB, image)
    current = view(band.images, :, image)
    forward = view(band.segment_logs, :, image)
    backward = view(band.transported_logs, :, image)
    forward .= view(band.images, :, image + 1) .- current
    backward .= current .- view(band.images, :, image - 1)

    energy = band.energies[image]
    energy_plus = band.energies[image + 1]
    energy_minus = band.energies[image - 1]
    if energy_plus > energy > energy_minus
        copyto!(tangent, forward)
    elseif energy_plus < energy < energy_minus
        copyto!(tangent, backward)
    elseif (energy_plus < energy && energy > energy_minus) ||
           (energy_plus > energy && energy < energy_minus)
        maximum_difference = max(
            abs(energy_plus - energy), abs(energy_minus - energy))
        minimum_difference = min(
            abs(energy_plus - energy), abs(energy_minus - energy))
        if energy_plus > energy_minus
            @. tangent = maximum_difference * forward +
                         minimum_difference * backward
        else
            @. tangent = minimum_difference * forward +
                         maximum_difference * backward
        end
    else
        @. tangent = forward + backward
    end
    _project_tangent!(tangent, current, band.backend)
    tangent_norm = norm(tangent)
    tangent_norm > 0 ||
        error("GNEB tangent vanished at image $image.")
    tangent ./= tangent_norm
    return tangent
end

function _gneb_endpoint_tangent(band, image)
    tangent = create_zeros(size(band.images, 1))
    if image == 1
        tangent .= view(band.images, :, 2) .- view(band.images, :, 1)
    else
        tangent .= view(band.images, :, image) .- view(band.images, :, image - 1)
    end
    _project_tangent!(tangent, view(band.images, :, image), band.backend)
    tangent_norm = Float64(norm(tangent))
    tangent_norm > 0 && (tangent ./= tangent_norm)
    return tangent
end

function _gneb_hermite_lengths(band::GNEB, ratio::Float64;
                               n_interpolations::Int=20)
    n_images = size(band.images, 2)
    derivatives = zeros(n_images)
    for image in 1:n_images
        tangent = image in (1, n_images) ?
                  _gneb_endpoint_tangent(band, image) :
                  view(band.tangents, :, image)
        derivatives[image] =
            -Float64(dot(view(band.gradients, :, image), tangent))
    end
    rx_samples = Float64[]
    energy_samples = Float64[]
    for image in 1:(n_images - 1)
        x0, x1 = band.reaction_coordinate[image:image + 1]
        e0, e1 = band.energies[image:image + 1]
        h = x1 - x0
        for k in 0:n_interpolations
            image > 1 && k == 0 && continue
            t = k / n_interpolations
            h00, h10 = 2t^3 - 3t^2 + 1, t^3 - 2t^2 + t
            h01, h11 = -2t^3 + 3t^2, t^3 - t^2
            push!(rx_samples, x0 + t * h)
            push!(energy_samples, h00 * e0 + h10 * h * derivatives[image] +
                                  h01 * e1 + h11 * h * derivatives[image + 1])
        end
    end
    range_rx = maximum(rx_samples) - minimum(rx_samples)
    range_energy = maximum(energy_samples) - minimum(energy_samples)
    scale_energy = max(floatmin(Float64), maximum(abs, energy_samples))
    if !(range_rx > 0) || !(range_energy > eps(Float64) * scale_energy)
        return [0.0; diff(band.reaction_coordinate)]
    end
    lengths = zeros(n_images)
    offset = 1
    for image in 2:n_images
        for _ in 1:n_interpolations
            drx = (1 - ratio) * (rx_samples[offset + 1] - rx_samples[offset]) / range_rx
            de = ratio * (energy_samples[offset + 1] - energy_samples[offset]) / range_energy
            lengths[image] += hypot(drx, de)
            offset += 1
        end
        lengths[image] *= range_rx
    end
    return lengths
end

function _update_gneb_forces!(band::GNEB; spring_force_ratio::Real=0.0,
                              path_shortening_constant::Real=0.0)
    fill!(band.forces, 0)
    fill!(band.tangents, 0)
    ratio = clamp(Float64(spring_force_ratio), 0.0, 1.0)
    for image in 2:(size(band.images, 2) - 1)
        _gneb_tangent!(view(band.tangents, :, image), band, image)
    end
    lengths = ratio > 0 ? _gneb_hermite_lengths(band, ratio) : Float64[]
    for image in 2:(size(band.images, 2) - 1)
        tangent = view(band.tangents, :, image)
        force = view(band.forces, :, image)
        force .= -view(band.gradients, :, image)
        parallel_force = Float64(dot(force, tangent))
        image_type = band.image_types[image]
        if image_type === :climbing || image == band.climbing_image
            force .-= 2 * parallel_force .* tangent
        elseif image_type === :stationary || band.stationary[image]
            fill!(force, 0)
        elseif image_type === :falling
            nothing
        else
            force .-= parallel_force .* tangent
            left_distance = band.reaction_coordinate[image] -
                            band.reaction_coordinate[image - 1]
            right_distance = band.reaction_coordinate[image + 1] -
                             band.reaction_coordinate[image]
            spring_force = ratio > 0 ?
                band.spring_constant * (lengths[image + 1] - lengths[image]) :
                band.spring_constant * (right_distance - left_distance)
            force .+= spring_force .* tangent
            if path_shortening_constant > 0
                force .+= _path_shortening_force!(
                    band, image, tangent, force,
                    Float64(path_shortening_constant))
            end
        end
        _project_tangent!(force, view(band.images, :, image), band.backend)
    end
    maximum_force = 0.0
    maximum_component = 0.0
    for image in 2:(size(band.images, 2) - 1)
        band.stationary[image] && image != band.climbing_image && continue
        image_force = view(band.forces, :, image)
        maximum_force = max(maximum_force,
                            _maximum_site_norm(image_force, band.backend))
        maximum_component = max(maximum_component,
                                Float64(maximum(abs, image_force)))
    end
    band.maximum_force = maximum_force
    band.maximum_component = maximum_component
    return band
end

function _path_shortening_force!(band::GNEB, image, tangent, gradient_force,
                                  constant)
    forward = similar(tangent)
    backward = similar(tangent)
    forward .= view(band.images, :, image + 1) .- view(band.images, :, image)
    backward .= view(band.images, :, image) .- view(band.images, :, image - 1)
    forward_norm = Float64(norm(forward))
    backward_norm = Float64(norm(backward))
    forward_norm > 0 && (forward ./= forward_norm)
    backward_norm > 0 && (backward ./= backward_norm)
    shrink = similar(tangent)
    @. shrink = forward - backward
    gradient_norm = Float64(norm(gradient_force))
    if gradient_norm > 0
        gradient_direction = gradient_force ./ gradient_norm
        shrink .-= Float64(dot(shrink, gradient_direction)) .* gradient_direction
    end
    shrink .-= Float64(dot(shrink, tangent)) .* tangent
    shrink_norm = Float64(norm(shrink))
    if shrink_norm > 0
        shrink ./= shrink_norm
        shrink .*= max(gradient_norm,
                       count(band.backend.active_cpu) * constant)
    end
    return shrink
end

function _update_gneb!(band::GNEB; refresh_all=false,
                       spring_force_ratio::Real=0.0,
                       path_shortening_constant::Real=0.0)
    _update_gneb_energies_gradients!(band; refresh_all=refresh_all)
    _update_gneb_reaction_coordinate!(band)
    _update_gneb_forces!(band; spring_force_ratio=spring_force_ratio,
                         path_shortening_constant=path_shortening_constant)
    return band
end

"""
    set_climbing_image!(band, image=:automatic)

Set the climbing image. `:automatic` selects the highest-energy interior
image. Passing `0` disables climbing. A climbing image is never stationary.
"""
function set_climbing_image!(band::GNEB, image=:automatic)
    selected = image === :automatic ?
               argmax(view(band.energies, 2:(length(band.energies) - 1))) + 1 :
               Int(image)
    0 <= selected <= size(band.images, 2) ||
        error("climbing image is outside the path.")
    selected in (1, size(band.images, 2)) &&
        error("An endpoint cannot be a climbing image.")
    # Clear any previous climbing marker from image_types.
    for k in 1:length(band.image_types)
        if band.image_types[k] === :climbing
            band.image_types[k] = band.stationary[k] ? :stationary : :normal
        end
    end
    band.climbing_image = selected
    if selected > 0
        band.stationary[selected] = false
        band.image_types[selected] = :climbing
    end
    fill!(band.velocity, 0)
    _update_gneb!(band)
    return selected
end

"""
    set_image_type!(band, image, type)

Set the type of an individual image. ``type`` must be one of
``(:normal, :climbing, :falling, :stationary)``. Endpoints can only be set
to ``:stationary``. This is the user-facing entry point for the Spirit-style
falling-image workflow.
"""
function set_image_type!(band::GNEB, image::Int, type::Symbol)
    1 <= image <= size(band.images, 2) ||
        error("image $image is outside the path.")
    type in _GNEB_IMAGE_TYPES ||
        error("image type must be one of $_GNEB_IMAGE_TYPES, got :$type.")
    if image in (1, size(band.images, 2))
        type === :stationary ||
            error("Endpoints can only be :stationary, got :$type.")
    end
    if type === :stationary
        band.stationary[image] = true
    else
        band.stationary[image] = false
    end
    if type === :climbing
        for k in eachindex(band.image_types)
            band.image_types[k] === :climbing &&
                (band.image_types[k] = band.stationary[k] ? :stationary : :normal)
        end
        band.climbing_image = image
    elseif band.climbing_image == image
        band.climbing_image = 0
    end
    band.image_types[image] = type
    fill!(band.velocity, 0)
    _update_gneb!(band)
    return type
end

function _gneb_velocity_projection_step!(
    band::GNEB, time_step, mass, maximum_step)

    # Spirit's VP solver uses the previous and current forces to update the
    # velocity, then projects the complete chain velocity onto the current
    # force. Dividing by the spring energy maps joules to Spirit's
    # dimensionless energy scale without allocating scaled band copies.
    inverse_energy_scale = inv(band.energy_scale)
    band.velocity .+= (0.5 / mass * inverse_energy_scale) .*
                     (band.previous_forces .+ band.forces)
    power = Float64(dot(band.velocity, band.forces))
    force_norm_squared = Float64(dot(band.forces, band.forces))
    if power > 0 && force_norm_squared > 0
        band.velocity .= (power / force_norm_squared) .* band.forces
    else
        fill!(band.velocity, 0)
    end

    copyto!(band.previous_forces, band.forces)
    for image in 2:(size(band.images, 2) - 1)
        if band.stationary[image] && image != band.climbing_image
            fill!(view(band.velocity, :, image), 0)
            fill!(view(band.previous_forces, :, image), 0)
            continue
        end
        velocity = view(band.velocity, :, image)
        force = view(band.forces, :, image)
        displacement = view(band.tangents, :, image)
        @. displacement = time_step * velocity +
                          (0.5 / mass) * time_step *
                          inverse_energy_scale * force
        local_scale = _maximum_site_norm(displacement, band.backend)
        scale = 1.0
        local_scale == 0 && continue
        kernel! = _gneb_retraction_step_kernel!(
            default_backend[], groupsize[])
        kernel!(
            view(band.images, :, image), displacement, Float[](scale),
            band.backend.active; ndrange=length(band.backend.active_cpu))
    end
    return band
end

function _gneb_fire_step!(
    band::GNEB, time_step, alpha, positive_steps, maximum_step)

    inverse_energy_scale = inv(band.energy_scale)
    band.velocity .+=
        (time_step * inverse_energy_scale) .* band.forces
    power = Float64(dot(band.velocity, band.forces))
    speed = Float64(norm(band.velocity))
    force_norm = Float64(norm(band.forces))
    if power > 0
        if speed > 0 && force_norm > 0
            @. band.velocity =
                (1 - alpha) * band.velocity +
                alpha * speed / force_norm * band.forces
        end
        positive_steps += 1
        if positive_steps > 5
            time_step = min(0.5, 1.1 * time_step)
            alpha *= 0.99
        end
    else
        fill!(band.velocity, 0)
        time_step *= 0.5
        alpha = 0.1
        positive_steps = 0
    end

    for image in 2:(size(band.images, 2) - 1)
        if band.stationary[image] && image != band.climbing_image
            fill!(view(band.velocity, :, image), 0)
            continue
        end
        velocity = view(band.velocity, :, image)
        displacement = view(band.tangents, :, image)
        displacement .= time_step .* velocity
        displacement_norm =
            _maximum_site_norm(displacement, band.backend)
        displacement_norm == 0 && continue
        scale = displacement_norm > maximum_step ?
                maximum_step / displacement_norm : 1.0
        previous = copy(view(band.images, :, image))
        updated = similar(displacement)
        _geodesic_step!(
            updated, previous, displacement, scale, band.backend)
        _transport_tangent!(
            velocity, velocity, previous, updated, band.backend)
        copyto!(view(band.images, :, image), updated)
    end
    return time_step, alpha, positive_steps
end

function _write_gneb_history(band::GNEB)
    isempty(band.name) && return nothing
    open(band.name * "_energy.txt", "w") do io
        print(io, "# steps")
        for image in eachindex(band.energies)
            print(io, " E_total_", image)
        end
        println(io)
        print(io, "# <>")
        for _ in eachindex(band.energies)
            print(io, " <J>")
        end
        println(io)
        for item in band.history
            @printf(io, "%d", item.step)
            for energy in item.energies
                @printf(io, " %.16e", energy)
            end
            println(io)
        end
    end
    open(band.name * "_distance.txt", "w") do io
        print(io, "# steps")
        for image in 1:(length(band.reaction_coordinate) - 1)
            print(io, " distance_", image)
        end
        println(io)
        print(io, "# <>")
        for _ in 1:(length(band.reaction_coordinate) - 1)
            print(io, " <>")
        end
        println(io)
        for item in band.history
            @printf(io, "%d", item.step)
            for distance in diff(item.reaction_coordinate)
                @printf(io, " %.16e", distance)
            end
            println(io)
        end
    end
    open(band.name * "_gneb_history.txt", "w") do io
        println(io,
                "# step maximum_site_force_J maximum_component_J climbing_image")
        println(io, "# <> <J> <J> <>")
        for item in band.history
            @printf(io, "%d %.16e %.16e %d\n", item.step,
                    item.maximum_force, item.maximum_component,
                    item.climbing_image)
        end
    end
    return nothing
end

"""
    relax_gneb!(band; climbing=false, maximum_steps=20_000,
                force_tolerance=1e-7meV, solver=:vp, symmetry=nothing, ...)

Relax a `GNEB` chain with one of three solvers:

  * `solver=:vp`         – Spirit's velocity projection (default, legacy).
  * `solver=:lbfgs_atlas`– Spirit's stereographic, chain-wide LBFGS Atlas.
  * `solver=:fire`       – MicroMagnetic manifold FIRE extension.

Convergence is the maximum site norm of the tangential GNEB force in joules,
matching Spirit. The maximum component is retained as a diagnostic. An
internally detected transition symmetry
may be supplied to keep the path in its invariant subspace.

# Spirit-compatible force enhancements (optional, off by default)

  * `spring_force_ratio=η ∈ [0,1]`   – blend Rx and energy-weighted spring.
  * `path_shortening_constant=c`     – enable Spirit path shortening for c>0 J/spin.
  * `set_image_type!(band, i, :falling)` – use a falling image (per-image).

"""
function relax_gneb!(
    band::GNEB; climbing::Bool=false, maximum_steps::Int=20_000,
    force_tolerance::Real=1e-7meV, time_step::Real=0.1,
    mass::Real=1.0, maximum_step::Real=0.03, log_every::Int=100,
    save_data_every::Int=20, solver::Symbol=:vp, symmetry=nothing,
    on_step=nothing,
    # Spirit-compatible force enhancements (off by default).
    spring_force_ratio::Real=0.0,
    path_shortening_constant::Real=0.0)

    maximum_steps > 0 || error("maximum_steps must be positive.")
    force_tolerance > 0 || error("force_tolerance must be positive.")
    time_step > 0 || error("time_step must be positive.")
    mass > 0 || error("mass must be positive.")
    maximum_step > 0 || error("maximum_step must be positive.")
    _check_solver_symbol(solver)
    climbing && band.climbing_image == 0 &&
        set_climbing_image!(band, :automatic)
    !climbing && band.climbing_image != 0 &&
        set_climbing_image!(band, 0)

    band.converged = false
    empty!(band.history)
    fill!(band.velocity, 0)
    !isnothing(symmetry) && _project_gneb_symmetry!(band, symmetry)
    _update_gneb!(band; spring_force_ratio=spring_force_ratio,
                  path_shortening_constant=path_shortening_constant)
    copyto!(band.previous_forces, band.forces)

    fire_time_step = min(Float64(time_step), 0.05)
    fire_alpha = 0.1
    fire_positive_steps = 0
    vp_time_step = Float64(time_step)
    vp_improving_steps = 0
    atlas_improving_steps = 0
    solver_state = _init_gneb_solver_state(
        solver, band, time_step, mass, maximum_step)
    vp_accepted_images = solver === :vp ? similar(band.images) : nothing
    atlas_accepted_images = solver === :lbfgs_atlas ? similar(band.images) : nothing
    if solver === :lbfgs_atlas && climbing
        # A climbing image starts close to a saddle, where Spirit's fixed
        # 0.05-radian trust radius can be unnecessarily large.  Choose the
        # initial radius from the dimensionless force and let the existing
        # accept/reject adaptation enlarge it again when that is safe.
        dimensionless_force = band.maximum_force / band.energy_scale
        solver_state.maxmove = clamp(dimensionless_force, 1e-10, 0.05)
    end

    for step in 1:maximum_steps
        band.steps = step
        current_maximum_energy = maximum(band.energies)
        band.converged = band.maximum_force < force_tolerance
        effective_time_step = solver === :fire ? fire_time_step :
                              solver === :vp ? vp_time_step : solver_state.maxmove

        item = (step=step, maximum_force=band.maximum_force,
                maximum_component=band.maximum_component,
                maximum_energy=current_maximum_energy,
                climbing_image=band.climbing_image,
                time_step=effective_time_step,
                solver=solver, energies=copy(band.energies),
                reaction_coordinate=copy(band.reaction_coordinate))
        if step == 1 ||
           (save_data_every > 0 && step % save_data_every == 0) ||
           band.converged
            push!(band.history, item)
        end
        if step == 1 || step % log_every == 0 || band.converged
            @info "GNEB progress" step maximum_force=band.maximum_force maximum_energy=current_maximum_energy climbing_image=band.climbing_image solver=solver
        end
        step_result = isnothing(on_step) ? nothing : on_step(band, item)
        step_result === :stop && break
        band.converged && break

        if solver === :vp
            copyto!(vp_accepted_images, band.images)
            _gneb_velocity_projection_step!(
                band, effective_time_step, Float64(mass),
                Float64(maximum_step))
        elseif solver === :fire
            fire_time_step, fire_alpha, fire_positive_steps =
                _gneb_fire_step!(
                    band, effective_time_step, fire_alpha,
                    fire_positive_steps, Float64(maximum_step))
        elseif solver === :lbfgs_atlas
            copyto!(atlas_accepted_images, band.images)
            _gneb_lbfgs_step!(band, solver_state)
        end
        !isnothing(symmetry) && _project_gneb_symmetry!(band, symmetry)
        _update_gneb!(band; spring_force_ratio=spring_force_ratio,
                      path_shortening_constant=path_shortening_constant)
        if solver === :vp
            new_force = band.maximum_force
            old_force = item.maximum_force
            if isfinite(old_force) && old_force > force_tolerance &&
               new_force > 1.25 * old_force
                copyto!(band.images, vp_accepted_images)
                vp_time_step = max(1e-12, 0.5 * vp_time_step)
                vp_improving_steps = 0
                fill!(band.velocity, 0)
                _update_gneb!(
                    band; spring_force_ratio=spring_force_ratio,
                    path_shortening_constant=path_shortening_constant)
                copyto!(band.previous_forces, band.forces)
            elseif new_force < 0.9 * old_force
                vp_improving_steps += 1
                if vp_improving_steps >= 5
                    vp_time_step = min(max(0.2, Float64(time_step)),
                                       1.1 * vp_time_step)
                    vp_improving_steps = 0
                end
            else
                vp_improving_steps = 0
            end
        elseif solver === :lbfgs_atlas
            new_force = band.maximum_force
            old_force = item.maximum_force
            if isfinite(old_force) && old_force > force_tolerance &&
               new_force > 1.25 * old_force
                # Reject a step outside the local Atlas trust region.  Merely
                # shrinking the next step leaves the chain in the overshot
                # configuration and can trigger false instability rejection.
                copyto!(band.images, atlas_accepted_images)
                solver_state.maxmove = max(1e-10, 0.5 * solver_state.maxmove)
                atlas_improving_steps = 0
                _reset_solver_state!(solver_state)
                _update_gneb!(
                    band; spring_force_ratio=spring_force_ratio,
                    path_shortening_constant=path_shortening_constant)
            elseif new_force < 0.8 * old_force
                atlas_improving_steps += 1
                if atlas_improving_steps >= 5
                    solver_state.maxmove = min(0.05,
                                                1.1 * solver_state.maxmove)
                    atlas_improving_steps = 0
                end
            else
                atlas_improving_steps = 0
            end
        end
    end
    save_data_every > 0 && _write_gneb_history(band)
    return band
end

"""
    gneb_images(band)

Return the complete GNEB chain on the CPU.
"""
gneb_images(band::GNEB) = Float64.(Array(band.images))
