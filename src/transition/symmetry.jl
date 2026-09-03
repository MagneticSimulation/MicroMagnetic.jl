struct _TransitionSymmetry{P}
    name::Symbol
    partner::P
    partner_cpu::Vector{Int}
    signs::NTuple{3,Float64}
end

@kernel function _project_transition_spin_symmetry_kernel!(
    output, @Const(input), @Const(partner), sign_x, sign_y, sign_z,
    @Const(active))

    site = @index(Global)
    j = 3 * site - 2
    if active[site]
        partner_site = partner[site]
        p = 3 * partner_site - 2
        x = 0.5 * (input[j] + sign_x * input[p])
        y = 0.5 * (input[j + 1] + sign_y * input[p + 1])
        z = 0.5 * (input[j + 2] + sign_z * input[p + 2])
        local_norm = sqrt(x * x + y * y + z * z)
        if local_norm > eps(eltype(output))
            output[j] = x / local_norm
            output[j + 1] = y / local_norm
            output[j + 2] = z / local_norm
        else
            output[j] = input[j]
            output[j + 1] = input[j + 1]
            output[j + 2] = input[j + 2]
        end
    else
        output[j] = input[j]
        output[j + 1] = input[j + 1]
        output[j + 2] = input[j + 2]
    end
end

@kernel function _project_transition_vector_symmetry_kernel!(
    output, @Const(input), @Const(partner), sign_x, sign_y, sign_z,
    @Const(active))

    site = @index(Global)
    j = 3 * site - 2
    if active[site]
        partner_site = partner[site]
        p = 3 * partner_site - 2
        output[j] = 0.5 * (input[j] + sign_x * input[p])
        output[j + 1] = 0.5 * (input[j + 1] + sign_y * input[p + 1])
        output[j + 2] = 0.5 * (input[j + 2] + sign_z * input[p + 2])
    else
        output[j] = zero(eltype(output))
        output[j + 1] = zero(eltype(output))
        output[j + 2] = zero(eltype(output))
    end
end

function _half_turn_transition_symmetry(backend::_TransitionBackend)
    sim = backend.evaluator
    mesh = sim.mesh
    all(property -> hasproperty(mesh, property), (:nx, :ny, :nz)) ||
        return nothing
    nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
    nx * ny * nz == length(backend.active_cpu) || return nothing
    partner_cpu = Vector{Int}(undef, nx * ny * nz)
    for k in 1:nz, j in 1:ny, i in 1:nx
        site = (k - 1) * nx * ny + (j - 1) * nx + i
        partner_cpu[site] =
            (k - 1) * nx * ny + (ny - j) * nx + (nx + 1 - i)
    end
    all(backend.active_cpu .== backend.active_cpu[partner_cpu]) ||
        return nothing
    partner = kernel_array(partner_cpu)
    return _TransitionSymmetry(
        :half_turn_z, partner, partner_cpu, (-1.0, -1.0, 1.0))
end

function _transition_symmetry_error(spin, symmetry::_TransitionSymmetry, active)
    values = Float64.(Array(spin))
    error_squared = 0.0
    count = 0
    sx, sy, sz = symmetry.signs
    for site in eachindex(active)
        active[site] || continue
        partner = symmetry.partner_cpu[site]
        j, p = 3 * site - 2, 3 * partner - 2
        dx = values[j] - sx * values[p]
        dy = values[j + 1] - sy * values[p + 1]
        dz = values[j + 2] - sz * values[p + 2]
        error_squared += dx * dx + dy * dy + dz * dz
        count += 1
    end
    return sqrt(error_squared / max(count, 1))
end

function _detect_transition_symmetries(
    backend::_TransitionBackend, states;
    tolerance=5e-3)

    symmetries = _TransitionSymmetry[]
    half_turn = _half_turn_transition_symmetry(backend)
    isnothing(half_turn) && return symmetries
    errors = [_transition_symmetry_error(
                  state, half_turn, backend.active_cpu) for state in states]
    maximum(errors) <= tolerance && push!(symmetries, half_turn)
    return symmetries
end

function _project_transition_spin!(
    output, input, symmetry::_TransitionSymmetry,
    backend::_TransitionBackend)

    sx, sy, sz = symmetry.signs
    kernel! = _project_transition_spin_symmetry_kernel!(
        default_backend[], groupsize[])
    kernel!(output, input, symmetry.partner, eltype(output)(sx), eltype(output)(sy),
            eltype(output)(sz), backend.active;
            ndrange=length(backend.active_cpu))
    return output
end

function _project_transition_vector!(
    output, input, spin, symmetry::_TransitionSymmetry,
    backend::_TransitionBackend)

    sx, sy, sz = symmetry.signs
    kernel! = _project_transition_vector_symmetry_kernel!(
        default_backend[], groupsize[])
    kernel!(output, input, symmetry.partner, eltype(output)(sx), eltype(output)(sy),
            eltype(output)(sz), backend.active;
            ndrange=length(backend.active_cpu))
    _project_tangent!(output, spin, backend)
    return output
end

function _project_gneb_symmetry!(
    band, symmetry::_TransitionSymmetry)

    for image in 2:(size(band.images, 2) - 1)
        spin = view(band.images, :, image)
        scratch = view(band.transported_logs, :, image)
        copyto!(scratch, spin)
        _project_transition_spin!(
            spin, scratch, symmetry, band.backend)

        velocity = view(band.velocity, :, image)
        copyto!(scratch, velocity)
        _project_transition_vector!(
            velocity, scratch, spin, symmetry, band.backend)
    end
    return band
end

function _symmetry_projected_modes(
    modes::HessianModes, spin, symmetry::_TransitionSymmetry,
    backend::_TransitionBackend; minimum_norm=1e-8)

    device_spin = _device_vector(spin)
    projected = Any[]
    for column in axes(modes.eigenvectors, 2)
        source = _device_vector(view(modes.eigenvectors, :, column))
        candidate = create_zeros(length(source))
        _project_transition_vector!(
            candidate, source, device_spin, symmetry, backend)
        _orthogonalize!(candidate, projected)
        candidate_norm = Float64(norm(candidate))
        candidate_norm > minimum_norm || continue
        candidate ./= candidate_norm
        push!(projected, candidate)
    end
    isempty(projected) &&
        return zeros(Float64, length(spin), 0)
    output = zeros(Float64, length(spin), length(projected))
    for column in eachindex(projected)
        output[:, column] .= Array(projected[column])
    end
    return output
end
